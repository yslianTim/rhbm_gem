#pragma once

#include "command/MapSimulationCommand.hpp"
#include "command/PotentialAnalysisCommand.hpp"
#include "command/PotentialDisplayCommand.hpp"
#include "command/RHBMTestCommand.hpp"
#include "command/ResultDumpCommand.hpp"
#ifdef RHBM_GEM_ENABLE_EXPERIMENTAL_FEATURE
#include "command/MapVisualizationCommand.hpp"
#include "command/PositionEstimationCommand.hpp"
#endif

#include <string_view>
#include <tuple>
#include <vector>

namespace rhbm_gem::command_internal {

template <typename CommandType>
struct CommandEntry
{
    using Command = CommandType;
    using Request = typename CommandType::RequestType;

    std::string_view cli_name;
    std::string_view description;
    std::string_view request_type_name;
};

inline constexpr auto kStableCommands = std::tuple{
    CommandEntry<PotentialAnalysisCommand>{
        "potential_analysis",
        "Run potential analysis",
        "PotentialAnalysisRequest"
    },
    CommandEntry<PotentialDisplayCommand>{
        "potential_display",
        "Run potential display",
        "PotentialDisplayRequest"
    },
    CommandEntry<ResultDumpCommand>{
        "result_dump",
        "Run result dump",
        "ResultDumpRequest"
    },
    CommandEntry<MapSimulationCommand>{
        "map_simulation",
        "Run map simulation command",
        "MapSimulationRequest"
    },
    CommandEntry<RHBMTestCommand>{
        "rhbm_test",
        "Run RHBM simulation test",
        "RHBMTestRequest"
    },
};

#ifdef RHBM_GEM_ENABLE_EXPERIMENTAL_FEATURE
inline constexpr auto kExperimentalCommands = std::tuple{
    CommandEntry<MapVisualizationCommand>{
        "map_visualization",
        "Run map visualization",
        "MapVisualizationRequest"
    },
    CommandEntry<PositionEstimationCommand>{
        "position_estimation",
        "Run atom position estimation",
        "PositionEstimationRequest"
    },
};
#endif

template <typename Visitor>
void VisitCommandCatalog(Visitor && visitor)
{
    std::apply([&visitor](const auto &... entries)
    {
        (static_cast<void>(visitor(entries)), ...);
    }, kStableCommands);
#ifdef RHBM_GEM_ENABLE_EXPERIMENTAL_FEATURE
    std::apply([&visitor](const auto &... entries)
    {
        (static_cast<void>(visitor(entries)), ...);
    }, kExperimentalCommands);
#endif
}

template <typename CommandType>
CommandResult ExecuteCommand(const typename CommandType::RequestType & request)
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

} // namespace rhbm_gem::command_internal
