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

#include <filesystem>
#include <string_view>
#include <tuple>
#include <unordered_map>
#include <vector>

namespace rhbm_gem::command_internal {

template <typename Owner, typename ValueType>
struct RequestField
{
    const char * field_name;
    const char * cli_flags;
    const char * help;
    ValueType Owner::* member;

    RequestField(
        const char * field_name_value,
        const char * cli_flags_value,
        const char * help_value,
        ValueType Owner::* member_value) :
        field_name{ field_name_value },
        cli_flags{ cli_flags_value },
        help{ help_value },
        member{ member_value }
    {
    }
};

template <typename Visitor, typename... Fields>
void VisitFieldList(Visitor && visitor, const Fields &... fields)
{
    (static_cast<void>(visitor(fields)), ...);
}

template <typename Request>
struct RequestFieldCatalog;

template <>
struct RequestFieldCatalog<CommandRequestBase>
{
    template <typename Visitor>
    static void Visit(Visitor && visitor)
    {
        using Self = CommandRequestBase;
        VisitFieldList(visitor,
            RequestField{"job_count",  "-j,--jobs",    "Number of threads",  &Self::job_count },
            RequestField{"verbosity",  "-v,--verbose", "Verbose level",      &Self::verbosity },
            RequestField{"output_dir", "-o,--folder",  "folder for outputs", &Self::output_dir});
    }
};

template <>
struct RequestFieldCatalog<PotentialAnalysisRequest>
{
    template <typename Visitor>
    static void Visit(Visitor && visitor)
    {
        using Self = PotentialAnalysisRequest;
        VisitFieldList(visitor,
            RequestField{ "database_path", "-d,--database",
                "Database file path",
                &Self::database_path },
            RequestField{ "model_file_path", "-a,--model",
                "Model file path",
                &Self::model_file_path },
            RequestField{ "map_file_path", "-m,--map",
                "Map file path",
                &Self::map_file_path },
            RequestField{ "training_report_dir", "--training-report-dir",
                "Optional output directory for alpha training reports",
                &Self::training_report_dir },
            RequestField{ "simulation_flag", "--simulation",
                "Simulation flag",
                &Self::simulation_flag },
            RequestField{ "simulated_map_resolution", "-r,--sim-resolution",
                "Set simulated map's resolution (blurring width)",
                &Self::simulated_map_resolution },
            RequestField{ "saved_key_tag", "-k,--save-key",
                "New key tag for saving ModelObject results into database",
                &Self::saved_key_tag },
            RequestField{ "training_alpha_flag", "--training-alpha",
                "Turn On/Off alpha training flag",
                &Self::training_alpha_flag },
            RequestField{ "training_alpha_min", "--training-alpha-min",
                "Minimum alpha value for training search",
                &Self::training_alpha_min },
            RequestField{ "training_alpha_max", "--training-alpha-max",
                "Maximum alpha value for training search",
                &Self::training_alpha_max },
            RequestField{ "training_alpha_step", "--training-alpha-step",
                "Step size for alpha training search",
                &Self::training_alpha_step },
            RequestField{ "asymmetry_flag", "--asymmetry",
                "Turn On/Off asymmetry flag",
                &Self::asymmetry_flag },
            RequestField{ "sampling_size", "-s,--sampling",
                "Number of sampling points per atom",
                &Self::sampling_size },
            RequestField{ "sampling_range_min", "--sampling-min",
                "Minimum sampling range",
                &Self::sampling_range_min },
            RequestField{ "sampling_range_max", "--sampling-max",
                "Maximum sampling range",
                &Self::sampling_range_max },
            RequestField{ "sampling_height", "--sampling-height",
                "Maximum sampling height",
                &Self::sampling_height },
            RequestField{ "fit_range_min", "--fit-min",
                "Minimum fitting range",
                &Self::fit_range_min },
            RequestField{ "fit_range_max", "--fit-max",
                "Maximum fitting range",
                &Self::fit_range_max },
            RequestField{ "alpha_r", "--alpha-r",
                "Alpha value for R",
                &Self::alpha_r },
            RequestField{ "alpha_g", "--alpha-g",
                "Alpha value for G",
                &Self::alpha_g });
    }
};

template <>
struct RequestFieldCatalog<PotentialDisplayRequest>
{
    template <typename Visitor>
    static void Visit(Visitor && visitor)
    {
        using Self = PotentialDisplayRequest;
        VisitFieldList(visitor,
            RequestField{ "database_path", "-d,--database",
                "Database file path",
                &Self::database_path },
            RequestField{ "painter_choice", "-p,--painter",
                "Painter choice",
                &Self::painter_choice },
            RequestField{ "model_key_tag_list", "-k,--model-keylist",
                "List of model key tag to be display",
                &Self::model_key_tag_list },
            RequestField{ "reference_model_groups", "--ref-group",
                "Reference group in the form group=key1,key2,...",
                &Self::reference_model_groups },
            RequestField{ "pick_chain_id", "--pick-chain",
                "Pick chain ID",
                &Self::pick_chain_id },
            RequestField{ "pick_residue", "--pick-residue",
                "Pick residue type",
                &Self::pick_residue },
            RequestField{ "pick_element", "--pick-element",
                "Pick element type",
                &Self::pick_element },
            RequestField{ "veto_chain_id", "--veto-chain",
                "Veto chain ID",
                &Self::veto_chain_id },
            RequestField{ "veto_residue", "--veto-residue",
                "Veto residue type",
                &Self::veto_residue },
            RequestField{ "veto_element", "--veto-element",
                "Veto element type",
                &Self::veto_element });
    }
};

template <>
struct RequestFieldCatalog<ResultDumpRequest>
{
    template <typename Visitor>
    static void Visit(Visitor && visitor)
    {
        using Self = ResultDumpRequest;
        VisitFieldList(visitor,
            RequestField{ "database_path", "-d,--database",
                "Database file path",
                &Self::database_path },
            RequestField{ "printer_choice", "-p,--printer",
                "Printer choice",
                &Self::printer_choice },
            RequestField{ "model_key_tag_list", "-k,--model-keylist",
                "List of model key tag to be display",
                &Self::model_key_tag_list },
            RequestField{ "map_file_path", "-m,--map",
                "Map file path",
                &Self::map_file_path });
    }
};

template <>
struct RequestFieldCatalog<MapSimulationRequest>
{
    template <typename Visitor>
    static void Visit(Visitor && visitor)
    {
        using Self = MapSimulationRequest;
        VisitFieldList(visitor,
            RequestField{ "model_file_path", "-a,--model",
                "Model file path",
                &Self::model_file_path },
            RequestField{ "map_file_name", "-n,--name",
                "File name for output map files",
                &Self::map_file_name },
            RequestField{ "potential_model_choice", "--potential-model",
                "Atomic potential model option",
                &Self::potential_model_choice },
            RequestField{ "partial_charge_choice", "--charge",
                "Partial charge table option",
                &Self::partial_charge_choice },
            RequestField{ "cutoff_distance", "-c,--cut-off",
                "Cutoff distance",
                &Self::cutoff_distance },
            RequestField{ "grid_spacing", "-g,--grid-spacing",
                "Grid spacing",
                &Self::grid_spacing },
            RequestField{ "blurring_width_list", "--blurring-width",
                "Blurring width (list) setting",
                &Self::blurring_width_list });
    }
};

template <>
struct RequestFieldCatalog<RHBMTestRequest>
{
    template <typename Visitor>
    static void Visit(Visitor && visitor)
    {
        using Self = RHBMTestRequest;
        VisitFieldList(visitor,
            RequestField{ "tester_choice", "-t,--tester",
                "Tester option",
                &Self::tester_choice },
            RequestField{ "fit_range_min", "--fit-min",
                "Minimum fitting range",
                &Self::fit_range_min },
            RequestField{ "fit_range_max", "--fit-max",
                "Maximum fitting range",
                &Self::fit_range_max },
            RequestField{ "alpha_r", "--alpha-r",
                "Alpha value for R",
                &Self::alpha_r },
            RequestField{ "alpha_g", "--alpha-g",
                "Alpha value for G",
                &Self::alpha_g });
    }
};

#ifdef RHBM_GEM_ENABLE_EXPERIMENTAL_FEATURE
template <>
struct RequestFieldCatalog<MapVisualizationRequest>
{
    template <typename Visitor>
    static void Visit(Visitor && visitor)
    {
        using Self = MapVisualizationRequest;
        VisitFieldList(visitor,
            RequestField{ "model_file_path", "-a,--model",
                "Model file path",
                &Self::model_file_path },
            RequestField{ "map_file_path", "-m,--map",
                "Map file path",
                &Self::map_file_path },
            RequestField{ "atom_serial_id", "-i,--atom-id",
                "Atom serial ID for visualization",
                &Self::atom_serial_id },
            RequestField{ "sampling_size", "-s,--sampling",
                "Number of sampling points per atom",
                &Self::sampling_size },
            RequestField{ "window_size", "--window-size",
                "Window size for sampling",
                &Self::window_size });
    }
};

template <>
struct RequestFieldCatalog<PositionEstimationRequest>
{
    template <typename Visitor>
    static void Visit(Visitor && visitor)
    {
        using Self = PositionEstimationRequest;
        VisitFieldList(visitor,
            RequestField{ "map_file_path", "-m,--map",
                "Map file path",
                &Self::map_file_path },
            RequestField{ "iteration_count", "--iter",
                "Iteration count for estimation",
                &Self::iteration_count },
            RequestField{ "knn_size", "--knn",
                "KNN size for estimation",
                &Self::knn_size },
            RequestField{ "alpha", "--alpha",
                "Alpha value for robust regression",
                &Self::alpha },
            RequestField{ "threshold_ratio", "--threshold",
                "Ratio of threshold of map values",
                &Self::threshold_ratio },
            RequestField{ "dedup_tolerance", "--dedup-tolerance",
                "Tolerance for deduplicating points",
                &Self::dedup_tolerance });
    }
};
#endif

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
    std::vector<CommandDiagnostic> public_issues;
    public_issues.reserve(internal_issues.size());
    for (const auto & issue : internal_issues)
    {
        public_issues.push_back(CommandDiagnostic{ issue.option_name, issue.message });
    }

    CommandResult result;
    result.succeeded = command.Run();
    result.issues = public_issues;
    return result;
}

} // namespace rhbm_gem::command_internal
