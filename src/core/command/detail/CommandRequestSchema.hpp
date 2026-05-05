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
    const char * python_name;
    const char * cli_flags;
    const char * help;
    ValueType Owner::* member;
    bool required{ false };
    char delimiter{ ',' };
    char assignment_delimiter{ '=' };
    char item_delimiter{ ',' };

    FieldSpec(
        const char * python_name_value,
        const char * cli_flags_value,
        const char * help_value,
        ValueType Owner::* member_value,
        bool required_value = false) :
        python_name{ python_name_value },
        cli_flags{ cli_flags_value },
        help{ help_value },
        member{ member_value },
        required{ required_value }
    {
    }

    FieldSpec(
        const char * python_name_value,
        const char * cli_flags_value,
        const char * help_value,
        ValueType Owner::* member_value,
        bool required_value,
        char delimiter_value) :
        FieldSpec{ python_name_value, cli_flags_value, help_value, member_value, required_value }
    {
        delimiter = delimiter_value;
        item_delimiter = delimiter_value;
    }

    FieldSpec(
        const char * python_name_value,
        const char * cli_flags_value,
        const char * help_value,
        ValueType Owner::* member_value,
        bool required_value,
        char assignment_delimiter_value,
        char item_delimiter_value) :
        FieldSpec{ python_name_value, cli_flags_value, help_value, member_value, required_value }
    {
        assignment_delimiter = assignment_delimiter_value;
        item_delimiter = item_delimiter_value;
    }
};

template <typename Request>
struct CommandRequestSchema;

template <>
struct CommandRequestSchema<CommandRequestBase>
{
    template <typename Visitor>
    static void Visit(Visitor && visitor)
    {
        using Self = CommandRequestBase;
        visitor(FieldSpec{
            "job_count",
            "-j,--jobs",
            "Number of threads",
            &Self::job_count });
        visitor(FieldSpec{
            "verbosity",
            "-v,--verbose",
            "Verbose level",
            &Self::verbosity });
        visitor(FieldSpec{
            "output_dir",
            "-o,--folder",
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
        visitor(FieldSpec{
            "database_path",
            "-d,--database",
            "Database file path",
            &Self::database_path });
        visitor(FieldSpec{
            "model_file_path",
            "-a,--model",
            "Model file path",
            &Self::model_file_path,
            true });
        visitor(FieldSpec{
            "map_file_path",
            "-m,--map",
            "Map file path",
            &Self::map_file_path,
            true });
        visitor(FieldSpec{
            "training_report_dir",
            "--training-report-dir",
            "Optional output directory for alpha training reports",
            &Self::training_report_dir });
        visitor(FieldSpec{
            "simulation_flag",
            "--simulation",
            "Simulation flag",
            &Self::simulation_flag });
        visitor(FieldSpec{
            "simulated_map_resolution",
            "-r,--sim-resolution",
            "Set simulated map's resolution (blurring width)",
            &Self::simulated_map_resolution });
        visitor(FieldSpec{
            "saved_key_tag",
            "-k,--save-key",
            "New key tag for saving ModelObject results into database",
            &Self::saved_key_tag });
        visitor(FieldSpec{
            "training_alpha_flag",
            "--training-alpha",
            "Turn On/Off alpha training flag",
            &Self::training_alpha_flag });
        visitor(FieldSpec{
            "training_alpha_min",
            "--training-alpha-min",
            "Minimum alpha value for training search",
            &Self::training_alpha_min });
        visitor(FieldSpec{
            "training_alpha_max",
            "--training-alpha-max",
            "Maximum alpha value for training search",
            &Self::training_alpha_max });
        visitor(FieldSpec{
            "training_alpha_step",
            "--training-alpha-step",
            "Step size for alpha training search",
            &Self::training_alpha_step });
        visitor(FieldSpec{
            "asymmetry_flag",
            "--asymmetry",
            "Turn On/Off asymmetry flag",
            &Self::asymmetry_flag });
        visitor(FieldSpec{
            "sampling_size",
            "-s,--sampling",
            "Number of sampling points per atom",
            &Self::sampling_size });
        visitor(FieldSpec{
            "sampling_range_min",
            "--sampling-min",
            "Minimum sampling range",
            &Self::sampling_range_min });
        visitor(FieldSpec{
            "sampling_range_max",
            "--sampling-max",
            "Maximum sampling range",
            &Self::sampling_range_max });
        visitor(FieldSpec{
            "sampling_height",
            "--sampling-height",
            "Maximum sampling height",
            &Self::sampling_height });
        visitor(FieldSpec{
            "fit_range_min",
            "--fit-min",
            "Minimum fitting range",
            &Self::fit_range_min });
        visitor(FieldSpec{
            "fit_range_max",
            "--fit-max",
            "Maximum fitting range",
            &Self::fit_range_max });
        visitor(FieldSpec{
            "alpha_r",
            "--alpha-r",
            "Alpha value for R",
            &Self::alpha_r });
        visitor(FieldSpec{
            "alpha_g",
            "--alpha-g",
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
        visitor(FieldSpec{
            "database_path",
            "-d,--database",
            "Database file path",
            &Self::database_path });
        visitor(FieldSpec{
            "painter_choice",
            "-p,--painter",
            "Painter choice",
            &Self::painter_choice,
            true });
        visitor(FieldSpec{
            "model_key_tag_list",
            "-k,--model-keylist",
            "List of model key tag to be display",
            &Self::model_key_tag_list,
            true });
        visitor(FieldSpec{
            "reference_model_groups",
            "--ref-group",
            "Reference group in the form group=key1,key2,...",
            &Self::reference_model_groups });
        visitor(FieldSpec{
            "pick_chain_id",
            "--pick-chain",
            "Pick chain ID",
            &Self::pick_chain_id });
        visitor(FieldSpec{
            "pick_residue",
            "--pick-residue",
            "Pick residue type",
            &Self::pick_residue });
        visitor(FieldSpec{
            "pick_element",
            "--pick-element",
            "Pick element type",
            &Self::pick_element });
        visitor(FieldSpec{
            "veto_chain_id",
            "--veto-chain",
            "Veto chain ID",
            &Self::veto_chain_id });
        visitor(FieldSpec{
            "veto_residue",
            "--veto-residue",
            "Veto residue type",
            &Self::veto_residue });
        visitor(FieldSpec{
            "veto_element",
            "--veto-element",
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
        visitor(FieldSpec{
            "database_path",
            "-d,--database",
            "Database file path",
            &Self::database_path });
        visitor(FieldSpec{
            "printer_choice",
            "-p,--printer",
            "Printer choice",
            &Self::printer_choice,
            true });
        visitor(FieldSpec{
            "model_key_tag_list",
            "-k,--model-keylist",
            "List of model key tag to be display",
            &Self::model_key_tag_list,
            true });
        visitor(FieldSpec{
            "map_file_path",
            "-m,--map",
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
        visitor(FieldSpec{
            "model_file_path",
            "-a,--model",
            "Model file path",
            &Self::model_file_path,
            true });
        visitor(FieldSpec{
            "map_file_name",
            "-n,--name",
            "File name for output map files",
            &Self::map_file_name });
        visitor(FieldSpec{
            "potential_model_choice",
            "--potential-model",
            "Atomic potential model option",
            &Self::potential_model_choice });
        visitor(FieldSpec{
            "partial_charge_choice",
            "--charge",
            "Partial charge table option",
            &Self::partial_charge_choice });
        visitor(FieldSpec{
            "cutoff_distance",
            "-c,--cut-off",
            "Cutoff distance",
            &Self::cutoff_distance });
        visitor(FieldSpec{
            "grid_spacing",
            "-g,--grid-spacing",
            "Grid spacing",
            &Self::grid_spacing });
        visitor(FieldSpec{
            "blurring_width_list",
            "--blurring-width",
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
        visitor(FieldSpec{
            "tester_choice",
            "-t,--tester",
            "Tester option",
            &Self::tester_choice });
        visitor(FieldSpec{
            "fit_range_min",
            "--fit-min",
            "Minimum fitting range",
            &Self::fit_range_min });
        visitor(FieldSpec{
            "fit_range_max",
            "--fit-max",
            "Maximum fitting range",
            &Self::fit_range_max });
        visitor(FieldSpec{
            "alpha_r",
            "--alpha-r",
            "Alpha value for R",
            &Self::alpha_r });
        visitor(FieldSpec{
            "alpha_g",
            "--alpha-g",
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
        visitor(FieldSpec{
            "model_file_path",
            "-a,--model",
            "Model file path",
            &Self::model_file_path,
            true });
        visitor(FieldSpec{
            "map_file_path",
            "-m,--map",
            "Map file path",
            &Self::map_file_path,
            true });
        visitor(FieldSpec{
            "atom_serial_id",
            "-i,--atom-id",
            "Atom serial ID for visualization",
            &Self::atom_serial_id });
        visitor(FieldSpec{
            "sampling_size",
            "-s,--sampling",
            "Number of sampling points per atom",
            &Self::sampling_size });
        visitor(FieldSpec{
            "window_size",
            "--window-size",
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
        visitor(FieldSpec{
            "map_file_path",
            "-m,--map",
            "Map file path",
            &Self::map_file_path,
            true });
        visitor(FieldSpec{
            "iteration_count",
            "--iter",
            "Iteration count for estimation",
            &Self::iteration_count });
        visitor(FieldSpec{
            "knn_size",
            "--knn",
            "KNN size for estimation",
            &Self::knn_size });
        visitor(FieldSpec{
            "alpha",
            "--alpha",
            "Alpha value for robust regression",
            &Self::alpha });
        visitor(FieldSpec{
            "threshold_ratio",
            "--threshold",
            "Ratio of threshold of map values",
            &Self::threshold_ratio });
        visitor(FieldSpec{
            "dedup_tolerance",
            "--dedup-tolerance",
            "Tolerance for deduplicating points",
            &Self::dedup_tolerance });
    }
};
#endif

} // namespace rhbm_gem::command_internal
