#include <rhbm_gem/core/command/CommandSystem.hpp>

#include <rhbm_gem/utils/domain/ScopeTimer.hpp>
#include <rhbm_gem/utils/domain/StringHelper.hpp>

#include "detail/CommandCatalog.hpp"
#include "detail/CommandRequestSchema.hpp"

#include <CLI/CLI.hpp>
#include <filesystem>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <string_view>
#include <type_traits>
#include <typeinfo>
#include <unordered_map>
#include <vector>

namespace rhbm_gem {
namespace {

template <typename FieldType>
struct VectorTraits
{
    static constexpr bool kIsVector{ false };
};

template <typename ElementType, typename Allocator>
struct VectorTraits<std::vector<ElementType, Allocator>>
{
    static constexpr bool kIsVector{ true };
    using Element = ElementType;
};

template <typename FieldType>
struct IsReferenceGroupField : std::false_type
{
};

template <>
struct IsReferenceGroupField<std::unordered_map<std::string, std::vector<std::string>>> : std::true_type
{
};

template <typename EnumType, typename = void>
struct HasCommandEnumTraits : std::false_type
{
};

template <typename EnumType>
struct HasCommandEnumTraits<EnumType, std::void_t<decltype(internal::CommandEnumTraits<EnumType>::kOptions)>>
    : std::true_type
{
};

template <typename Request>
void BindPathCliField(
    CLI::App & command,
    Request * request,
    const command_internal::FieldSpec<Request, std::filesystem::path> & field)
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
    const auto & current_value{ request->*(field.member) };
    if (!current_value.empty())
    {
        option.default_val(current_value.string());
    }
}

template <typename Request, typename FieldType>
void BindReferenceGroupCliField(
    CLI::App & command,
    Request * request,
    const command_internal::FieldSpec<Request, FieldType> & field)
{
    static constexpr char kRefGroupAssignmentDelimiter{ '=' };
    static constexpr char kRefGroupItemDelimiter{ ',' };
    auto & option{
        *command.add_option_function<std::vector<std::string>>(
            field.cli_flags,
            [request,
             member = field.member,
             option_name = field.cli_flags](const std::vector<std::string> & values)
            {
                for (const auto & value : values)
                {
                    const auto delimiter_position{ value.find(kRefGroupAssignmentDelimiter) };
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
                        string_helper::ParseListOption<std::string>(
                            member_text,
                            kRefGroupItemDelimiter)
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
}

template <typename Request, typename FieldType>
void BindVectorCliField(
    CLI::App & command,
    Request * request,
    const command_internal::FieldSpec<Request, FieldType> & field)
{
    using ElementType = typename VectorTraits<FieldType>::Element;
    static constexpr char kCsvListDelimiter{ ',' };
    auto & option{
        *command.add_option_function<std::string>(
            field.cli_flags,
            [request, member = field.member](const std::string & value)
            {
                request->*member =
                    string_helper::ParseListOption<ElementType>(value, kCsvListDelimiter);
            },
            field.help)
    };
    const auto & current_value{ request->*(field.member) };
    if (!current_value.empty())
    {
        std::ostringstream default_value;
        for (std::size_t index = 0; index < current_value.size(); ++index)
        {
            if (index != 0)
            {
                default_value << kCsvListDelimiter;
            }
            default_value << current_value[index];
        }
        option.default_val(default_value.str());
    }
}

template <typename Request, typename FieldType>
void BindEnumCliField(
    CLI::App & command,
    Request * request,
    const command_internal::FieldSpec<Request, FieldType> & field)
{
    auto & option{
        *command.add_option_function<FieldType>(
            field.cli_flags,
            [request, member = field.member](const FieldType & value)
            {
                request->*member = value;
            },
            field.help)
    };
    option.default_val(request->*(field.member));

    std::map<std::string, FieldType> option_map;
    for (const auto & option_definition : internal::CommandEnumTraits<FieldType>::kOptions)
    {
        for (const auto alias : option_definition.cli_aliases)
        {
            option_map.emplace(std::string(alias), option_definition.value);
        }
    }
    option.transform(CLI::CheckedTransformer(option_map, CLI::ignore_case));
}

template <typename Request, typename FieldType>
void BindScalarCliField(
    CLI::App & command,
    Request * request,
    const command_internal::FieldSpec<Request, FieldType> & field)
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

template <typename Request, typename FieldType>
void BindCliField(
    CLI::App & command,
    Request * request,
    const command_internal::FieldSpec<Request, FieldType> & field)
{
    if constexpr (std::is_same_v<FieldType, std::filesystem::path>)
    {
        BindPathCliField(command, request, field);
    }
    else if constexpr (IsReferenceGroupField<FieldType>::value)
    {
        BindReferenceGroupCliField(command, request, field);
    }
    else if constexpr (VectorTraits<FieldType>::kIsVector)
    {
        BindVectorCliField(command, request, field);
    }
    else if constexpr (HasCommandEnumTraits<FieldType>::value)
    {
        BindEnumCliField(command, request, field);
    }
    else
    {
        BindScalarCliField(command, request, field);
    }
}

template <typename RequestType>
void RegisterCommand(CLI::App & app, std::string_view name, std::string_view description)
{
    CLI::App & command{ *app.add_subcommand(std::string(name), std::string(description)) };
    auto request{ std::make_shared<RequestType>() };
    auto & base_request{ static_cast<CommandRequestBase &>(*request) };
    command_internal::CommandRequestSchema<CommandRequestBase>::Visit([&](const auto & field)
    {
        BindCliField(command, &base_request, field);
    });
    command_internal::CommandRequestSchema<RequestType>::Visit([&](const auto & field)
    {
        BindCliField(command, request.get(), field);
    });

    command.callback([request]()
    {
        ScopeTimer timer("Command CLI callback");
        const auto result{ rhbm_gem::RunCommand(*request) };
        if (!result.succeeded) throw CLI::RuntimeError(1);
    });
}

} // namespace

const std::vector<CommandInfo> & ListCommands()
{
    static const std::vector<CommandInfo> commands = []
    {
        std::vector<CommandInfo> command_infos;
        command_internal::VisitCommandCatalog([&](const auto & entry)
        {
            command_infos.push_back(CommandInfo{ entry.cli_name, entry.description });
        });
        return command_infos;
    }();
    return commands;
}

namespace command_internal {

CommandResult RunCommandByRequestType(
    const std::type_info & request_type,
    const CommandRequestBase & request)
{
    CommandResult result;
    result.succeeded = false;
    result.issues.push_back(ValidationIssue{ "request_type", "Unsupported command request type." });
    bool matched{ false };
    VisitCommandCatalog([&](const auto & entry)
    {
        if (matched) return;
        using CommandType = typename std::decay_t<decltype(entry)>::Command;
        using RequestType = typename std::decay_t<decltype(entry)>::Request;
        if (request_type != typeid(RequestType)) return;

        result = ExecuteCommand<CommandType>(static_cast<const RequestType &>(request));
        matched = true;
    });
    return result;
}

} // namespace command_internal

int RunCommandCLI(int argc, char * argv[])
{
    CLI::App app{"RHBM-GEM"};
    app.require_subcommand(1);

    command_internal::VisitCommandCatalog([&](const auto & entry)
    {
        using RequestType = typename std::decay_t<decltype(entry)>::Request;
        RegisterCommand<RequestType>(app, entry.cli_name, entry.description);
    });

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
