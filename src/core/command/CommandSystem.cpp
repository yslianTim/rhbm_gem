#include <rhbm_gem/core/command/CommandSystem.hpp>

#include <rhbm_gem/utils/domain/ScopeTimer.hpp>
#include <rhbm_gem/utils/domain/StringHelper.hpp>

#include "detail/CommandRequestSchema.hpp"
#include "RHBMTestCommand.hpp"
#include "MapSimulationCommand.hpp"
#include "PotentialAnalysisCommand.hpp"
#include "PotentialDisplayCommand.hpp"
#include "ResultDumpCommand.hpp"
#ifdef RHBM_GEM_ENABLE_EXPERIMENTAL_FEATURE
#include "MapVisualizationCommand.hpp"
#include "PositionEstimationCommand.hpp"
#endif

#include <CLI/CLI.hpp>
#include <filesystem>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <string_view>
#include <type_traits>
#include <vector>

namespace rhbm_gem {
namespace {

template <typename CommandType, typename RequestType>
CommandResult RunCommand(const RequestType & request)
{
    CommandType command{};
    command.ApplyRequest(request);

    const auto & internal_issues{ command.GetValidationIssues() };
    std::vector<ValidationIssue> public_issues;
    public_issues.reserve(internal_issues.size());
    for (const auto & issue : internal_issues)
    {
        public_issues.push_back(ValidationIssue{ issue.option_name, issue.message });
    }

    CommandResult result;
    result.succeeded = command.Run();
    result.issues = public_issues;
    return result;
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
    if (field.required)
    {
        option.required();
    }
    if constexpr (std::is_same_v<FieldType, std::string>)
    {
        if (!current_value.empty())
        {
            option.default_val(current_value);
        }
    }
    else
    {
        option.default_val(current_value);
    }
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
    if (field.required)
    {
        option.required();
    }
    const auto & current_value{ request->*(field.member) };
    if (!current_value.empty())
    {
        option.default_val(current_value.string());
    }
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
    if (field.required)
    {
        option.required();
    }
    option.default_val(request->*(field.member));

    std::map<std::string, EnumType> option_map;
    for (const auto & option_definition : internal::CommandEnumTraits<EnumType>::kOptions)
    {
        for (const auto alias : option_definition.cli_aliases)
        {
            option_map.emplace(std::string(alias), option_definition.value);
        }
    }
    option.transform(CLI::CheckedTransformer(option_map, CLI::ignore_case));
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
    if (field.required)
    {
        option.required();
    }
    const auto & current_value{ request->*(field.member) };
    if (!current_value.empty())
    {
        std::ostringstream default_value;
        for (std::size_t index = 0; index < current_value.size(); ++index)
        {
            if (index != 0)
            {
                default_value << field.delimiter;
            }
            default_value << current_value[index];
        }
        option.default_val(default_value.str());
    }
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
                for (const auto & value : values)
                {
                    const auto delimiter_position{ value.find(assignment_delimiter) };
                    if (delimiter_position == std::string::npos)
                    {
                        throw CLI::ValidationError(
                            option_name,
                            "Expected value in the form <group>=key1,key2,...");
                    }

                    const std::string group_name{ value.substr(0, delimiter_position) };
                    const std::string member_text{ value.substr(delimiter_position + 1) };
                    if (group_name.empty())
                    {
                        throw CLI::ValidationError(
                            option_name,
                            "Reference group name cannot be empty.");
                    }

                    const auto members{
                        string_helper::ParseListOption<std::string>(member_text, item_delimiter)
                    };
                    if (members.empty())
                    {
                        throw CLI::ValidationError(
                            option_name,
                            "Reference group members cannot be empty.");
                    }

                    auto & group_map{ request->*member };
                    auto & stored_members{ group_map[group_name] };
                    stored_members.insert(stored_members.end(), members.begin(), members.end());
                }
            },
            field.help)
    };
    // Reference groups accumulate repeated CLI occurrences into one logical map.
    option.multi_option_policy(CLI::MultiOptionPolicy::TakeAll);
    if (field.required)
    {
        option.required();
    }
}

template <typename RequestType, auto RunCommandFn>
void RegisterCommand(CLI::App & app, std::string_view name, std::string_view description)
{
    CLI::App & command{ *app.add_subcommand(std::string(name), std::string(description)) };
    auto request{ std::make_shared<RequestType>() };
    auto & base_request{ static_cast<CommandRequestBase &>(*request) };
    internal::CommandRequestSchema<CommandRequestBase>::Visit([&](const auto & field)
    {
        BindCliField(command, &base_request, field);
    });
    internal::CommandRequestSchema<RequestType>::Visit([&](const auto & field)
    {
        BindCliField(command, request.get(), field);
    });

    command.callback([request]()
    {
        ScopeTimer timer("Command CLI callback");
        const auto result{ RunCommandFn(*request) };
        if (!result.succeeded) throw CLI::RuntimeError(1);
    });
}

} // namespace

const std::vector<CommandInfo> & ListCommands()
{
    static const std::vector<CommandInfo> commands{
#define RHBM_GEM_COMMAND(COMMAND_ID, CLI_NAME, DESCRIPTION)                                    \
        CommandInfo{ CLI_NAME, DESCRIPTION },
#include <rhbm_gem/core/command/CommandManifest.def>
#undef RHBM_GEM_COMMAND
    };
    return commands;
}

#define RHBM_GEM_COMMAND(COMMAND_ID, CLI_NAME, DESCRIPTION)                                    \
CommandResult Run##COMMAND_ID(const COMMAND_ID##Request & request)                             \
{                                                                                              \
    return RunCommand<COMMAND_ID##Command>(request);                                           \
}
#include <rhbm_gem/core/command/CommandManifest.def>
#undef RHBM_GEM_COMMAND

int RunCommandCLI(int argc, char * argv[])
{
    CLI::App app{"RHBM-GEM"};
    app.require_subcommand(1);

#define RHBM_GEM_COMMAND(COMMAND_ID, CLI_NAME, DESCRIPTION)                                    \
    RegisterCommand<COMMAND_ID##Request, &Run##COMMAND_ID>(app, CLI_NAME, DESCRIPTION);
#include <rhbm_gem/core/command/CommandManifest.def>
#undef RHBM_GEM_COMMAND

    try
    {
        app.parse(argc, argv);
    }
    catch (const CLI::ParseError & error)
    {
        return app.exit(error);
    }
    return 0;
}

} // namespace rhbm_gem
