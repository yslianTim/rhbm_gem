#include <CLI/CLI.hpp>
#include <rhbm_gem/core/command/CommandApi.hpp>
#include <rhbm_gem/utils/domain/ScopeTimer.hpp>

#include "internal/command/CommandOptionSupport.hpp"
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

void BindCommonOptions(
    CLI::App * command,
    rhbm_gem::CommonCommandRequest & common,
    rhbm_gem::CommonOptionProfile profile)
{
    const auto common_options{ rhbm_gem::CommonOptionMaskForProfile(profile) };
    if (rhbm_gem::HasCommonOption(common_options, rhbm_gem::CommonOption::Threading))
    {
        rhbm_gem::command_cli::AddScalarOption<int>(
            command,
            "-j,--jobs",
            [&](int value) { common.thread_size = value; },
            "Number of threads",
            common.thread_size);
    }
    if (rhbm_gem::HasCommonOption(common_options, rhbm_gem::CommonOption::Verbose))
    {
        rhbm_gem::command_cli::AddScalarOption<int>(
            command,
            "-v,--verbose",
            [&](int value) { common.verbose_level = value; },
            "Verbose level",
            common.verbose_level);
    }
    if (rhbm_gem::HasCommonOption(common_options, rhbm_gem::CommonOption::Database))
    {
        rhbm_gem::command_cli::AddPathOption(
            command,
            "-d,--database",
            [&](const std::filesystem::path & value) { common.database_path = value; },
            "Database file path",
            common.database_path);
    }
    if (rhbm_gem::HasCommonOption(common_options, rhbm_gem::CommonOption::OutputFolder))
    {
        rhbm_gem::command_cli::AddPathOption(
            command,
            "-o,--folder",
            [&](const std::filesystem::path & value) { common.folder_path = value; },
            "folder path for output files",
            common.folder_path);
    }
}

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
    BindCommonOptions(command, request->common, profile);
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
