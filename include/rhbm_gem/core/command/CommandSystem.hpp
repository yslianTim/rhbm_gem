#pragma once

#include <cstddef>
#include <string_view>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include <rhbm_gem/core/command/CommandTypes.hpp>

namespace rhbm_gem {

CommandResult RunPotentialAnalysis(const PotentialAnalysisRequest & request);
CommandResult RunPotentialDisplay(const PotentialDisplayRequest & request);
CommandResult RunResultDump(const ResultDumpRequest & request);
CommandResult RunMapSimulation(const MapSimulationRequest & request);
CommandResult RunRHBMTest(const RHBMTestRequest & request);

#ifdef RHBM_GEM_ENABLE_EXPERIMENTAL_FEATURE
CommandResult RunMapVisualization(const MapVisualizationRequest & request);
CommandResult RunPositionEstimation(const PositionEstimationRequest & request);
#endif

const std::vector<CommandInfo> & ListCommands();
int RunCommandCLI(int argc, char * argv[]);

namespace command {

template <typename RequestType>
struct CommandEntry
{
    using Request = RequestType;

    std::string_view cli_name;
    std::string_view description;
    std::string_view request_type_name;
    std::string_view run_function_name;
    CommandResult (*run)(const RequestType & request);
};

inline constexpr auto kStableCommands = std::tuple{
    CommandEntry<PotentialAnalysisRequest>{
        "potential_analysis",
        "Run potential analysis",
        "PotentialAnalysisRequest",
        "RunPotentialAnalysis",
        &RunPotentialAnalysis,
    },
    CommandEntry<PotentialDisplayRequest>{
        "potential_display",
        "Run potential display",
        "PotentialDisplayRequest",
        "RunPotentialDisplay",
        &RunPotentialDisplay,
    },
    CommandEntry<ResultDumpRequest>{
        "result_dump",
        "Run result dump",
        "ResultDumpRequest",
        "RunResultDump",
        &RunResultDump,
    },
    CommandEntry<MapSimulationRequest>{
        "map_simulation",
        "Run map simulation command",
        "MapSimulationRequest",
        "RunMapSimulation",
        &RunMapSimulation,
    },
    CommandEntry<RHBMTestRequest>{
        "rhbm_test",
        "Run RHBM simulation test",
        "RHBMTestRequest",
        "RunRHBMTest",
        &RunRHBMTest,
    },
};

#ifdef RHBM_GEM_ENABLE_EXPERIMENTAL_FEATURE
inline constexpr auto kExperimentalCommands = std::tuple{
    CommandEntry<MapVisualizationRequest>{
        "map_visualization",
        "Run map visualization",
        "MapVisualizationRequest",
        "RunMapVisualization",
        &RunMapVisualization,
    },
    CommandEntry<PositionEstimationRequest>{
        "position_estimation",
        "Run atom position estimation",
        "PositionEstimationRequest",
        "RunPositionEstimation",
        &RunPositionEstimation,
    },
};
#endif

namespace detail {

template <typename Tuple, typename Visitor, std::size_t... Indices>
void VisitTuple(Tuple && commands, Visitor & visitor, std::index_sequence<Indices...>)
{
    (static_cast<void>(visitor(std::get<Indices>(std::forward<Tuple>(commands)))), ...);
}

template <typename Tuple, typename Visitor>
void VisitTuple(Tuple && commands, Visitor & visitor)
{
    using TupleType = std::remove_reference_t<Tuple>;
    VisitTuple(
        std::forward<Tuple>(commands),
        visitor,
        std::make_index_sequence<std::tuple_size_v<TupleType>>{});
}

} // namespace detail

template <typename Visitor>
void VisitCommands(Visitor && visitor)
{
    detail::VisitTuple(kStableCommands, visitor);
#ifdef RHBM_GEM_ENABLE_EXPERIMENTAL_FEATURE
    detail::VisitTuple(kExperimentalCommands, visitor);
#endif
}

} // namespace command

} // namespace rhbm_gem
