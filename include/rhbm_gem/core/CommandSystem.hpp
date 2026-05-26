#pragma once

#include <typeinfo>
#include <type_traits>

#include <rhbm_gem/core/CommandTypes.hpp>

namespace rhbm_gem::core {

int RunCommandCLI(int argc, char * argv[]);

namespace command_internal {

CommandResult RunCommandByRequestType(
    const std::type_info & request_type,
    const CommandRequestBase & request);

} // namespace command_internal

template <typename RequestType>
CommandResult RunCommand(const RequestType & request)
{
    static_assert(
        std::is_base_of_v<CommandRequestBase, RequestType>,
        "RunCommand<RequestType> requires RequestType to derive from CommandRequestBase.");
    return command_internal::RunCommandByRequestType(typeid(RequestType), request);
}

} // namespace rhbm_gem::core
