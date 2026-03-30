#include <CLI/CLI.hpp>
#include <rhbm_gem/core/command/CommandApi.hpp>
#include <rhbm_gem/utils/domain/ScopeTimer.hpp>

#include "internal/command/CommandOptionSupport.hpp"
#include "internal/command/CommandRegistry.hpp"
#include <rhbm_gem/utils/domain/StringHelper.hpp>

#include <functional>
#include <memory>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>

namespace {

template <typename ValueType>
std::optional<ValueType> BuildScalarDefault(const ValueType & value)
{
    return value;
}

template <>
std::optional<std::string> BuildScalarDefault(const std::string & value)
{
    if (value.empty())
    {
        return std::nullopt;
    }
    return value;
}

std::optional<std::filesystem::path> BuildPathDefault(const std::filesystem::path & value)
{
    if (value.empty())
    {
        return std::nullopt;
    }
    return value;
}

template <typename ValueType>
std::string JoinCsvValues(const std::vector<ValueType> & values, char delimiter)
{
    std::ostringstream stream;
    for (std::size_t index = 0; index < values.size(); ++index)
    {
        if (index != 0)
        {
            stream << delimiter;
        }
        if constexpr (std::is_same_v<ValueType, std::string>)
        {
            stream << values[index];
        }
        else
        {
            stream << values[index];
        }
    }
    return stream.str();
}

template <typename ValueType>
std::optional<std::string> BuildCsvDefault(const std::vector<ValueType> & values, char delimiter)
{
    if (values.empty())
    {
        return std::nullopt;
    }
    return JoinCsvValues(values, delimiter);
}

std::pair<std::string, std::vector<std::string>> ParseReferenceGroupArgument(
    const std::string & text,
    char assignment_delimiter,
    char item_delimiter)
{
    const auto delimiter_position{ text.find(assignment_delimiter) };
    if (delimiter_position == std::string::npos)
    {
        throw CLI::ValidationError(
            "--ref-group",
            "Expected value in the form <group>=key1,key2,...");
    }

    const std::string group_name{ text.substr(0, delimiter_position) };
    const std::string member_text{ text.substr(delimiter_position + 1) };
    if (group_name.empty())
    {
        throw CLI::ValidationError("--ref-group", "Reference group name cannot be empty.");
    }

    const auto members{ StringHelper::ParseListOption<std::string>(member_text, item_delimiter) };
    if (members.empty())
    {
        throw CLI::ValidationError("--ref-group", "Reference group members cannot be empty.");
    }
    return { group_name, members };
}

std::optional<rhbm_gem::CommonOption> ResolveCommonOptionByField(std::string_view python_name)
{
    using rhbm_gem::CommonOption;
    if (python_name == "thread_size") return CommonOption::Threading;
    if (python_name == "verbose_level") return CommonOption::Verbose;
    if (python_name == "database_path") return CommonOption::Database;
    if (python_name == "folder_path") return CommonOption::OutputFolder;
    return std::nullopt;
}

template <typename Request, typename FieldType>
void BindCliField(
    CLI::App * /*command*/,
    Request * /*request*/,
    const rhbm_gem::RequestObjectFieldSpec<Request, FieldType> & /*field*/)
{
}

template <typename Request, typename FieldType>
void BindCliField(
    CLI::App * command,
    Request * request,
    const rhbm_gem::RequestScalarFieldSpec<Request, FieldType> & field)
{
    const auto & current_value{ request->*(field.member) };
    if constexpr (std::is_same_v<FieldType, std::string>)
    {
        rhbm_gem::command_cli::AddStringOption(
            command,
            field.cli_flags,
            [request, member = field.member](const std::string & value) { request->*member = value; },
            field.help,
            BuildScalarDefault(current_value),
            field.required);
    }
    else
    {
        rhbm_gem::command_cli::AddScalarOption<FieldType>(
            command,
            field.cli_flags,
            [request, member = field.member](FieldType value) { request->*member = value; },
            field.help,
            BuildScalarDefault(current_value),
            field.required);
    }
}

template <typename Request>
void BindCliField(
    CLI::App * command,
    Request * request,
    const rhbm_gem::RequestPathFieldSpec<Request> & field)
{
    rhbm_gem::command_cli::AddPathOption(
        command,
        field.cli_flags,
        [request, member = field.member](const std::filesystem::path & value) { request->*member = value; },
        field.help,
        BuildPathDefault(request->*(field.member)),
        field.required);
}

template <typename Request, typename EnumType>
void BindCliField(
    CLI::App * command,
    Request * request,
    const rhbm_gem::RequestEnumFieldSpec<Request, EnumType> & field)
{
    rhbm_gem::command_cli::AddEnumOption<EnumType>(
        command,
        field.cli_flags,
        [request, member = field.member](EnumType value) { request->*member = value; },
        field.help,
        BuildScalarDefault(request->*(field.member)),
        field.required);
}

template <typename Request, typename ElementType>
void BindCliField(
    CLI::App * command,
    Request * request,
    const rhbm_gem::RequestCsvListFieldSpec<Request, ElementType> & field)
{
    rhbm_gem::command_cli::AddStringOption(
        command,
        field.cli_flags,
        [request, member = field.member, delimiter = field.delimiter](const std::string & value)
        {
            request->*member = StringHelper::ParseListOption<ElementType>(value, delimiter);
        },
        field.help,
        BuildCsvDefault(request->*(field.member), field.delimiter),
        field.required);
}

template <typename Request>
void BindCliField(
    CLI::App * command,
    Request * request,
    const rhbm_gem::RequestRefGroupFieldSpec<Request> & field)
{
    rhbm_gem::command_cli::AddRepeatableStringOption(
        command,
        field.cli_flags,
        [request,
         member = field.member,
         assignment_delimiter = field.assignment_delimiter,
         item_delimiter = field.item_delimiter](const std::string & value)
        {
            auto & group_map{ request->*member };
            const auto [group_name, members]{
                ParseReferenceGroupArgument(value, assignment_delimiter, item_delimiter)
            };
            auto & stored_members{ group_map[group_name] };
            stored_members.insert(stored_members.end(), members.begin(), members.end());
        },
        field.help,
        field.required);
}

template <typename Request>
void BindRequestFields(CLI::App * command, Request * request)
{
    Request::VisitFields([&](const auto & field)
    {
        BindCliField(command, request, field);
    });
}

void BindCommonOptions(
    CLI::App * command,
    rhbm_gem::CommonCommandRequest * common,
    rhbm_gem::CommonOptionProfile profile)
{
    const auto common_options{ rhbm_gem::CommonOptionMaskForProfile(profile) };
    rhbm_gem::CommonCommandRequest::VisitFields([&](const auto & field)
    {
        const auto common_option{ ResolveCommonOptionByField(field.python_name) };
        if (!common_option.has_value()
            || !rhbm_gem::HasCommonOption(common_options, common_option.value()))
        {
            return;
        }
        BindCliField(command, common, field);
    });
}

template <typename CommandType, typename RequestType, auto RunCommandFn>
void RegisterCommand(
    CLI::App & app,
    std::string_view name,
    std::string_view description,
    rhbm_gem::CommonOptionProfile profile)
{
    CLI::App * command{
        app.add_subcommand(
            std::string(name),
            std::string(description))
    };
    auto request{ std::make_shared<RequestType>() };
    BindCommonOptions(command, &request->common, profile);
    BindRequestFields(command, request.get());
    command->callback([request]()
    {
        ScopeTimer timer("Command CLI callback");
        const auto report{ RunCommandFn(*request) };
        if (!report.executed)
        {
            throw CLI::RuntimeError(1);
        }
    });
}

} // namespace

namespace rhbm_gem {

void ConfigureCommandCli(CLI::App & app)
{
    app.require_subcommand(1);

#define RHBM_GEM_COMMAND(COMMAND_ID, CLI_NAME, DESCRIPTION, PROFILE)                           \
    RegisterCommand<COMMAND_ID##Command, COMMAND_ID##Request, &Run##COMMAND_ID>(               \
        app,                                                                                   \
        CLI_NAME,                                                                              \
        DESCRIPTION,                                                                           \
        CommonOptionProfile::PROFILE);
#include <rhbm_gem/core/command/CommandList.def>
#undef RHBM_GEM_COMMAND
}

} // namespace rhbm_gem
