#include "CommandCli.hpp"

#include <rhbm_gem/utils/domain/ScopeTimer.hpp>
#include <rhbm_gem/utils/domain/StringHelper.hpp>

#include "CommandEnumMetadata.hpp"
#include "CommandRequestSchema.hpp"

#include <CLI/CLI.hpp>
#include <filesystem>
#include <memory>
#include <optional>
#include <sstream>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace rhbm_gem {
namespace {

template <typename ValueType>
std::optional<ValueType> MakeCliDefaultValue(const ValueType & value)
{
    return value;
}

template <>
std::optional<std::string> MakeCliDefaultValue(const std::string & value)
{
    if (value.empty())
    {
        return std::nullopt;
    }
    return value;
}

std::optional<std::filesystem::path> MakeCliDefaultValue(const std::filesystem::path & value)
{
    if (value.empty())
    {
        return std::nullopt;
    }
    return value;
}

template <typename ValueType>
std::optional<std::string> MakeCliDefaultValue(const std::vector<ValueType> & values, char delimiter)
{
    if (values.empty())
    {
        return std::nullopt;
    }

    std::ostringstream stream;
    for (std::size_t index = 0; index < values.size(); ++index)
    {
        if (index != 0)
        {
            stream << delimiter;
        }
        stream << values[index];
    }
    return stream.str();
}

void ApplyCliOptionSettings(CLI::Option & option, bool required)
{
    if (required)
    {
        option.required();
    }
}

template <typename DefaultType>
void ApplyCliOptionSettings(
    CLI::Option & option,
    bool required,
    const std::optional<DefaultType> & default_value)
{
    ApplyCliOptionSettings(option, required);
    if (default_value.has_value())
    {
        option.default_val(*default_value);
    }
}

std::optional<std::string> MakeCliPathDefaultText(const std::filesystem::path & value)
{
    const auto default_value{ MakeCliDefaultValue(value) };
    if (!default_value.has_value())
    {
        return std::nullopt;
    }
    return default_value->string();
}

std::pair<std::string, std::vector<std::string>> ParseReferenceGroupArgument(
    const std::string & text,
    const std::string & option_name,
    char assignment_delimiter,
    char item_delimiter)
{
    const auto delimiter_position{ text.find(assignment_delimiter) };
    if (delimiter_position == std::string::npos)
    {
        throw CLI::ValidationError(
            option_name,
            "Expected value in the form <group>=key1,key2,...");
    }

    const std::string group_name{ text.substr(0, delimiter_position) };
    const std::string member_text{ text.substr(delimiter_position + 1) };
    if (group_name.empty())
    {
        throw CLI::ValidationError(option_name, "Reference group name cannot be empty.");
    }

    const auto members{ string_helper::ParseListOption<std::string>(member_text, item_delimiter) };
    if (members.empty())
    {
        throw CLI::ValidationError(option_name, "Reference group members cannot be empty.");
    }
    return { group_name, members };
}

template <typename Request, typename FieldType>
void BindCliField(
    CLI::App & command,
    Request * request,
    const internal::RequestScalarFieldSpec<Request, FieldType> & field)
{
    const auto & current_value{ request->*(field.member) };
    auto & option{
        *command.add_option_function<FieldType>(
            field.cli_flags,
            [request, member = field.member](const FieldType & value)
            {
                request->*member = value;
            },
            field.help)
    };
    ApplyCliOptionSettings(option, field.required, MakeCliDefaultValue(current_value));
}

template <typename Request>
void BindCliField(
    CLI::App & command,
    Request * request,
    const internal::RequestPathFieldSpec<Request> & field)
{
    auto & option{
        *command.add_option_function<std::string>(
            field.cli_flags,
            [request, member = field.member](const std::string & value)
            {
                request->*member = std::filesystem::path{ value };
            },
            field.help)
    };
    ApplyCliOptionSettings(
        option,
        field.required,
        MakeCliPathDefaultText(request->*(field.member)));
}

template <typename Request, typename EnumType>
void BindCliField(
    CLI::App & command,
    Request * request,
    const internal::RequestEnumFieldSpec<Request, EnumType> & field)
{
    auto & option{
        *command.add_option_function<EnumType>(
            field.cli_flags,
            [request, member = field.member](const EnumType & value)
            {
                request->*member = value;
            },
            field.help)
    };
    ApplyCliOptionSettings(option, field.required, MakeCliDefaultValue(request->*(field.member)));
    option.transform(CLI::CheckedTransformer(
        internal::BuildCommandEnumCliMap<EnumType>(), CLI::ignore_case));
}

template <typename Request, typename ElementType>
void BindCliField(
    CLI::App & command,
    Request * request,
    const internal::RequestCsvListFieldSpec<Request, ElementType> & field)
{
    auto & option{
        *command.add_option_function<std::string>(
            field.cli_flags,
            [request, member = field.member, delimiter = field.delimiter](const std::string & value)
            {
                request->*member = string_helper::ParseListOption<ElementType>(value, delimiter);
            },
            field.help)
    };
    ApplyCliOptionSettings(
        option,
        field.required,
        MakeCliDefaultValue(request->*(field.member), field.delimiter));
}

template <typename Request>
void BindCliField(
    CLI::App & command,
    Request * request,
    const internal::RequestRefGroupFieldSpec<Request> & field)
{
    auto & option{
        *command.add_option_function<std::vector<std::string>>(
            field.cli_flags,
            [request,
             member = field.member,
             option_name = field.cli_flags,
             assignment_delimiter = field.assignment_delimiter,
             item_delimiter = field.item_delimiter](const std::vector<std::string> & values)
            {
                auto & group_map{ request->*member };
                for (const auto & value : values)
                {
                    const auto [group_name, members]{
                        ParseReferenceGroupArgument(
                            value,
                            option_name,
                            assignment_delimiter,
                            item_delimiter)
                    };
                    auto & stored_members{ group_map[group_name] };
                    stored_members.insert(stored_members.end(), members.begin(), members.end());
                }
            },
            field.help)
    };
    // Reference groups accumulate repeated CLI occurrences into one logical map.
    option.multi_option_policy(CLI::MultiOptionPolicy::TakeAll);
    ApplyCliOptionSettings(option, field.required);
}

template <typename RequestType>
void BindRequestFields(CLI::App & command, RequestType & request)
{
    auto & base_request{ static_cast<CommandRequestBase &>(request) };
    internal::VisitBaseRequestFields([&](const auto & field)
    {
        BindCliField(command, &base_request, field);
    });
    internal::VisitRequestFields<RequestType>([&](const auto & field)
    {
        BindCliField(command, &request, field);
    });
}

template <typename RequestType, auto RunCommandFn>
void RegisterCommand(
    CLI::App & app,
    std::string_view name,
    std::string_view description)
{
    CLI::App & command{
        *app.add_subcommand(
            std::string(name),
            std::string(description))
    };
    auto request{ std::make_shared<RequestType>() };

    BindRequestFields(command, *request);

    command.callback([request]()
    {
        ScopeTimer timer("Command CLI callback");
        const auto result{ RunCommandFn(*request) };
        if (!result.succeeded)
        {
            throw CLI::RuntimeError(1);
        }
    });
}

} // namespace

void ConfigureCommandCli(CLI::App & app)
{
    app.require_subcommand(1);

#define RHBM_GEM_COMMAND(COMMAND_ID, CLI_NAME, DESCRIPTION)                                    \
    RegisterCommand<COMMAND_ID##Request, &Run##COMMAND_ID>(                                    \
        app,                                                                                   \
        CLI_NAME,                                                                              \
        DESCRIPTION);
#include <rhbm_gem/core/command/CommandManifest.def>
#undef RHBM_GEM_COMMAND
}

} // namespace rhbm_gem
