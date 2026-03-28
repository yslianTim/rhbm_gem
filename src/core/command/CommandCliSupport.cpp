#include <CLI/CLI.hpp>
#include <rhbm_gem/core/command/CommandApi.hpp>
#include <rhbm_gem/utils/domain/ScopeTimer.hpp>

#include "internal/command/CommandCliSupport.hpp"
#include "internal/command/HRLModelTestCommand.hpp"
#include "internal/command/MapSimulationCommand.hpp"
#include "internal/command/MapVisualizationCommand.hpp"
#include "internal/command/PositionEstimationCommand.hpp"
#include "internal/command/PotentialAnalysisCommand.hpp"
#include "internal/command/PotentialDisplayCommand.hpp"
#include "internal/command/ResultDumpCommand.hpp"

#include <functional>
#include <memory>
#include <string>
#include <string_view>

namespace {

template <typename RequestType, auto RunCommandFn, typename BindOptionsFn>
void RegisterCommand(
    CLI::App & app,
    std::string_view name,
    std::string_view description,
    rhbm_gem::CommonOptionProfile profile,
    BindOptionsFn bind_options)
{
    CLI::App * command{
        app.add_subcommand(
            std::string(name),
            std::string(description))
    };
    auto request{ std::make_shared<RequestType>() };
    rhbm_gem::command_cli::BindCommonOptions(command, request->common, profile);
    std::invoke(bind_options, command, *request);
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
    RegisterCommand<COMMAND_ID##Request, &Run##COMMAND_ID>(                                    \
        app,                                                                                   \
        CLI_NAME,                                                                              \
        DESCRIPTION,                                                                           \
        CommonOptionProfile::PROFILE,                                                          \
        Bind##COMMAND_ID##RequestOptions);
#include <rhbm_gem/core/command/CommandList.def>
#undef RHBM_GEM_COMMAND
}

} // namespace rhbm_gem
