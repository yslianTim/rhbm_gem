#include "internal/CommandCatalog.hpp"
#include "internal/CommandRuntimeRegistry.hpp"

namespace rhbm_gem {

const std::vector<CommandDescriptor> & CommandCatalog()
{
    static const std::vector<CommandDescriptor> catalog{
    // BEGIN GENERATED: command-catalog-entries
    CommandDescriptor{CommandId::PotentialAnalysis, "potential_analysis", "Run potential analysis", CommonOptionProfile::DatabaseWorkflow, &BindPotentialAnalysisRuntime},
    CommandDescriptor{CommandId::PotentialDisplay, "potential_display", "Run potential display", CommonOptionProfile::DatabaseWorkflow, &BindPotentialDisplayRuntime},
    CommandDescriptor{CommandId::ResultDump, "result_dump", "Run result dump", CommonOptionProfile::DatabaseWorkflow, &BindResultDumpRuntime},
    CommandDescriptor{CommandId::MapSimulation, "map_simulation", "Run map simulation command", CommonOptionProfile::FileWorkflow, &BindMapSimulationRuntime},
    CommandDescriptor{CommandId::MapVisualization, "map_visualization", "Run map visualization", CommonOptionProfile::FileWorkflow, &BindMapVisualizationRuntime},
    CommandDescriptor{CommandId::PositionEstimation, "position_estimation", "Run atom position estimation", CommonOptionProfile::FileWorkflow, &BindPositionEstimationRuntime},
    CommandDescriptor{CommandId::ModelTest, "model_test", "Run HRL model simulation test", CommonOptionProfile::FileWorkflow, &BindModelTestRuntime},
// END GENERATED: command-catalog-entries
    };

    return catalog;
}

} // namespace rhbm_gem
