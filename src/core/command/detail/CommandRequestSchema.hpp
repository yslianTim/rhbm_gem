#pragma once

#include <filesystem>
#include <string>
#include <unordered_map>
#include <vector>

#include <rhbm_gem/core/command/CommandTypes.hpp>

namespace rhbm_gem::command_internal {

template <typename Owner, typename ValueType>
struct FieldSpec
{
    const char * field_name;
    const char * cli_flags;
    const char * help;
    ValueType Owner::* member;

    FieldSpec(
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
void VisitFields(Visitor && visitor, const Fields &... fields)
{
    (static_cast<void>(visitor(fields)), ...);
}

template <typename Request>
struct CommandRequestSchema;

template <>
struct CommandRequestSchema<CommandRequestBase>
{
    template <typename Visitor>
    static void Visit(Visitor && visitor)
    {
        using Self = CommandRequestBase;
        VisitFields(visitor,
            FieldSpec{ "job_count", "-j,--jobs",
                "Number of threads",
                &Self::job_count },
            FieldSpec{ "verbosity", "-v,--verbose",
                "Verbose level",
                &Self::verbosity },
            FieldSpec{ "output_dir", "-o,--folder",
                "folder path for output files",
                &Self::output_dir });
    }
};

template <>
struct CommandRequestSchema<PotentialAnalysisRequest>
{
    template <typename Visitor>
    static void Visit(Visitor && visitor)
    {
        using Self = PotentialAnalysisRequest;
        VisitFields(visitor,
            FieldSpec{ "database_path", "-d,--database",
                "Database file path",
                &Self::database_path },
            FieldSpec{ "model_file_path", "-a,--model",
                "Model file path",
                &Self::model_file_path },
            FieldSpec{ "map_file_path", "-m,--map",
                "Map file path",
                &Self::map_file_path },
            FieldSpec{ "training_report_dir", "--training-report-dir",
                "Optional output directory for alpha training reports",
                &Self::training_report_dir },
            FieldSpec{ "simulation_flag", "--simulation",
                "Simulation flag",
                &Self::simulation_flag },
            FieldSpec{ "simulated_map_resolution", "-r,--sim-resolution",
                "Set simulated map's resolution (blurring width)",
                &Self::simulated_map_resolution },
            FieldSpec{ "saved_key_tag", "-k,--save-key",
                "New key tag for saving ModelObject results into database",
                &Self::saved_key_tag },
            FieldSpec{ "training_alpha_flag", "--training-alpha",
                "Turn On/Off alpha training flag",
                &Self::training_alpha_flag },
            FieldSpec{ "training_alpha_min", "--training-alpha-min",
                "Minimum alpha value for training search",
                &Self::training_alpha_min },
            FieldSpec{ "training_alpha_max", "--training-alpha-max",
                "Maximum alpha value for training search",
                &Self::training_alpha_max },
            FieldSpec{ "training_alpha_step", "--training-alpha-step",
                "Step size for alpha training search",
                &Self::training_alpha_step },
            FieldSpec{ "asymmetry_flag", "--asymmetry",
                "Turn On/Off asymmetry flag",
                &Self::asymmetry_flag },
            FieldSpec{ "sampling_size", "-s,--sampling",
                "Number of sampling points per atom",
                &Self::sampling_size },
            FieldSpec{ "sampling_range_min", "--sampling-min",
                "Minimum sampling range",
                &Self::sampling_range_min },
            FieldSpec{ "sampling_range_max", "--sampling-max",
                "Maximum sampling range",
                &Self::sampling_range_max },
            FieldSpec{ "sampling_height", "--sampling-height",
                "Maximum sampling height",
                &Self::sampling_height },
            FieldSpec{ "fit_range_min", "--fit-min",
                "Minimum fitting range",
                &Self::fit_range_min },
            FieldSpec{ "fit_range_max", "--fit-max",
                "Maximum fitting range",
                &Self::fit_range_max },
            FieldSpec{ "alpha_r", "--alpha-r",
                "Alpha value for R",
                &Self::alpha_r },
            FieldSpec{ "alpha_g", "--alpha-g",
                "Alpha value for G",
                &Self::alpha_g });
    }
};

template <>
struct CommandRequestSchema<PotentialDisplayRequest>
{
    template <typename Visitor>
    static void Visit(Visitor && visitor)
    {
        using Self = PotentialDisplayRequest;
        VisitFields(visitor,
            FieldSpec{ "database_path", "-d,--database",
                "Database file path",
                &Self::database_path },
            FieldSpec{ "painter_choice", "-p,--painter",
                "Painter choice",
                &Self::painter_choice },
            FieldSpec{ "model_key_tag_list", "-k,--model-keylist",
                "List of model key tag to be display",
                &Self::model_key_tag_list },
            FieldSpec{ "reference_model_groups", "--ref-group",
                "Reference group in the form group=key1,key2,...",
                &Self::reference_model_groups },
            FieldSpec{ "pick_chain_id", "--pick-chain",
                "Pick chain ID",
                &Self::pick_chain_id },
            FieldSpec{ "pick_residue", "--pick-residue",
                "Pick residue type",
                &Self::pick_residue },
            FieldSpec{ "pick_element", "--pick-element",
                "Pick element type",
                &Self::pick_element },
            FieldSpec{ "veto_chain_id", "--veto-chain",
                "Veto chain ID",
                &Self::veto_chain_id },
            FieldSpec{ "veto_residue", "--veto-residue",
                "Veto residue type",
                &Self::veto_residue },
            FieldSpec{ "veto_element", "--veto-element",
                "Veto element type",
                &Self::veto_element });
    }
};

template <>
struct CommandRequestSchema<ResultDumpRequest>
{
    template <typename Visitor>
    static void Visit(Visitor && visitor)
    {
        using Self = ResultDumpRequest;
        VisitFields(visitor,
            FieldSpec{ "database_path", "-d,--database",
                "Database file path",
                &Self::database_path },
            FieldSpec{ "printer_choice", "-p,--printer",
                "Printer choice",
                &Self::printer_choice },
            FieldSpec{ "model_key_tag_list", "-k,--model-keylist",
                "List of model key tag to be display",
                &Self::model_key_tag_list },
            FieldSpec{ "map_file_path", "-m,--map",
                "Map file path",
                &Self::map_file_path });
    }
};

template <>
struct CommandRequestSchema<MapSimulationRequest>
{
    template <typename Visitor>
    static void Visit(Visitor && visitor)
    {
        using Self = MapSimulationRequest;
        VisitFields(visitor,
            FieldSpec{ "model_file_path", "-a,--model",
                "Model file path",
                &Self::model_file_path },
            FieldSpec{ "map_file_name", "-n,--name",
                "File name for output map files",
                &Self::map_file_name },
            FieldSpec{ "potential_model_choice", "--potential-model",
                "Atomic potential model option",
                &Self::potential_model_choice },
            FieldSpec{ "partial_charge_choice", "--charge",
                "Partial charge table option",
                &Self::partial_charge_choice },
            FieldSpec{ "cutoff_distance", "-c,--cut-off",
                "Cutoff distance",
                &Self::cutoff_distance },
            FieldSpec{ "grid_spacing", "-g,--grid-spacing",
                "Grid spacing",
                &Self::grid_spacing },
            FieldSpec{ "blurring_width_list", "--blurring-width",
                "Blurring width (list) setting",
                &Self::blurring_width_list });
    }
};

template <>
struct CommandRequestSchema<RHBMTestRequest>
{
    template <typename Visitor>
    static void Visit(Visitor && visitor)
    {
        using Self = RHBMTestRequest;
        VisitFields(visitor,
            FieldSpec{ "tester_choice", "-t,--tester",
                "Tester option",
                &Self::tester_choice },
            FieldSpec{ "fit_range_min", "--fit-min",
                "Minimum fitting range",
                &Self::fit_range_min },
            FieldSpec{ "fit_range_max", "--fit-max",
                "Maximum fitting range",
                &Self::fit_range_max },
            FieldSpec{ "alpha_r", "--alpha-r",
                "Alpha value for R",
                &Self::alpha_r },
            FieldSpec{ "alpha_g", "--alpha-g",
                "Alpha value for G",
                &Self::alpha_g });
    }
};

#ifdef RHBM_GEM_ENABLE_EXPERIMENTAL_FEATURE
template <>
struct CommandRequestSchema<MapVisualizationRequest>
{
    template <typename Visitor>
    static void Visit(Visitor && visitor)
    {
        using Self = MapVisualizationRequest;
        VisitFields(visitor,
            FieldSpec{ "model_file_path", "-a,--model",
                "Model file path",
                &Self::model_file_path },
            FieldSpec{ "map_file_path", "-m,--map",
                "Map file path",
                &Self::map_file_path },
            FieldSpec{ "atom_serial_id", "-i,--atom-id",
                "Atom serial ID for visualization",
                &Self::atom_serial_id },
            FieldSpec{ "sampling_size", "-s,--sampling",
                "Number of sampling points per atom",
                &Self::sampling_size },
            FieldSpec{ "window_size", "--window-size",
                "Window size for sampling",
                &Self::window_size });
    }
};

template <>
struct CommandRequestSchema<PositionEstimationRequest>
{
    template <typename Visitor>
    static void Visit(Visitor && visitor)
    {
        using Self = PositionEstimationRequest;
        VisitFields(visitor,
            FieldSpec{ "map_file_path", "-m,--map",
                "Map file path",
                &Self::map_file_path },
            FieldSpec{ "iteration_count", "--iter",
                "Iteration count for estimation",
                &Self::iteration_count },
            FieldSpec{ "knn_size", "--knn",
                "KNN size for estimation",
                &Self::knn_size },
            FieldSpec{ "alpha", "--alpha",
                "Alpha value for robust regression",
                &Self::alpha },
            FieldSpec{ "threshold_ratio", "--threshold",
                "Ratio of threshold of map values",
                &Self::threshold_ratio },
            FieldSpec{ "dedup_tolerance", "--dedup-tolerance",
                "Tolerance for deduplicating points",
                &Self::dedup_tolerance });
    }
};
#endif

} // namespace rhbm_gem::command_internal
