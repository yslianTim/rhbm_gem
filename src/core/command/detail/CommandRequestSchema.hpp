#pragma once

#include <filesystem>
#include <string>
#include <unordered_map>
#include <vector>

#include <rhbm_gem/core/command/CommandTypes.hpp>

namespace rhbm_gem::command_internal {

template <typename Owner, typename FieldType>
struct RequestScalarFieldSpec
{
    const char * python_name;
    const char * cli_flags;
    const char * help;
    FieldType Owner::* member;
    bool required{ false };
};

template <typename Owner, typename FieldType>
RequestScalarFieldSpec(const char *, const char *, const char *, FieldType Owner::*)
    -> RequestScalarFieldSpec<Owner, FieldType>;

template <typename Owner, typename FieldType>
RequestScalarFieldSpec(const char *, const char *, const char *, FieldType Owner::*, bool)
    -> RequestScalarFieldSpec<Owner, FieldType>;

template <typename Owner>
struct RequestPathFieldSpec
{
    const char * python_name;
    const char * cli_flags;
    const char * help;
    std::filesystem::path Owner::* member;
    bool required{ false };
};

template <typename Owner>
RequestPathFieldSpec(const char *, const char *, const char *, std::filesystem::path Owner::*)
    -> RequestPathFieldSpec<Owner>;

template <typename Owner>
RequestPathFieldSpec(const char *, const char *, const char *, std::filesystem::path Owner::*, bool)
    -> RequestPathFieldSpec<Owner>;

template <typename Owner, typename EnumType>
struct RequestEnumFieldSpec
{
    const char * python_name;
    const char * cli_flags;
    const char * help;
    EnumType Owner::* member;
    bool required{ false };
};

template <typename Owner, typename EnumType>
RequestEnumFieldSpec(const char *, const char *, const char *, EnumType Owner::*)
    -> RequestEnumFieldSpec<Owner, EnumType>;

template <typename Owner, typename EnumType>
RequestEnumFieldSpec(const char *, const char *, const char *, EnumType Owner::*, bool)
    -> RequestEnumFieldSpec<Owner, EnumType>;

template <typename Owner, typename ElementType>
struct RequestCsvListFieldSpec
{
    const char * python_name;
    const char * cli_flags;
    const char * help;
    std::vector<ElementType> Owner::* member;
    bool required{ false };
    char delimiter{ ',' };
};

template <typename Owner, typename ElementType>
RequestCsvListFieldSpec(const char *, const char *, const char *, std::vector<ElementType> Owner::*)
    -> RequestCsvListFieldSpec<Owner, ElementType>;

template <typename Owner, typename ElementType>
RequestCsvListFieldSpec(const char *, const char *, const char *, std::vector<ElementType> Owner::*, bool)
    -> RequestCsvListFieldSpec<Owner, ElementType>;

template <typename Owner, typename ElementType>
RequestCsvListFieldSpec(
    const char *,
    const char *,
    const char *,
    std::vector<ElementType> Owner::*,
    bool,
    char)
    -> RequestCsvListFieldSpec<Owner, ElementType>;

template <typename Owner>
struct RequestRefGroupFieldSpec
{
    const char * python_name;
    const char * cli_flags;
    const char * help;
    std::unordered_map<std::string, std::vector<std::string>> Owner::* member;
    bool required{ false };
    char assignment_delimiter{ '=' };
    char item_delimiter{ ',' };
};

template <typename Owner>
RequestRefGroupFieldSpec(
    const char *,
    const char *,
    const char *,
    std::unordered_map<std::string, std::vector<std::string>> Owner::*)
    -> RequestRefGroupFieldSpec<Owner>;

template <typename Owner>
RequestRefGroupFieldSpec(
    const char *,
    const char *,
    const char *,
    std::unordered_map<std::string, std::vector<std::string>> Owner::*,
    bool)
    -> RequestRefGroupFieldSpec<Owner>;

template <typename Owner>
RequestRefGroupFieldSpec(
    const char *,
    const char *,
    const char *,
    std::unordered_map<std::string, std::vector<std::string>> Owner::*,
    bool,
    char,
    char)
    -> RequestRefGroupFieldSpec<Owner>;

template <typename Request>
struct CommandRequestSchema;

template <>
struct CommandRequestSchema<CommandRequestBase>
{
    template <typename Visitor>
    static void Visit(Visitor && visitor)
    {
        using Self = CommandRequestBase;
        visitor(RequestScalarFieldSpec{
            "job_count",
            "-j,--jobs",
            "Number of threads",
            &Self::job_count });
        visitor(RequestScalarFieldSpec{
            "verbosity",
            "-v,--verbose",
            "Verbose level",
            &Self::verbosity });
        visitor(RequestPathFieldSpec{
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
        visitor(RequestPathFieldSpec{
            "database_path",
            "-d,--database",
            "Database file path",
            &Self::database_path });
        visitor(RequestPathFieldSpec{
            "model_file_path",
            "-a,--model",
            "Model file path",
            &Self::model_file_path,
            true });
        visitor(RequestPathFieldSpec{
            "map_file_path",
            "-m,--map",
            "Map file path",
            &Self::map_file_path,
            true });
        visitor(RequestPathFieldSpec{
            "training_report_dir",
            "--training-report-dir",
            "Optional output directory for alpha training reports",
            &Self::training_report_dir });
        visitor(RequestScalarFieldSpec{
            "simulation_flag",
            "--simulation",
            "Simulation flag",
            &Self::simulation_flag });
        visitor(RequestScalarFieldSpec{
            "simulated_map_resolution",
            "-r,--sim-resolution",
            "Set simulated map's resolution (blurring width)",
            &Self::simulated_map_resolution });
        visitor(RequestScalarFieldSpec{
            "saved_key_tag",
            "-k,--save-key",
            "New key tag for saving ModelObject results into database",
            &Self::saved_key_tag });
        visitor(RequestScalarFieldSpec{
            "training_alpha_flag",
            "--training-alpha",
            "Turn On/Off alpha training flag",
            &Self::training_alpha_flag });
        visitor(RequestScalarFieldSpec{
            "training_alpha_min",
            "--training-alpha-min",
            "Minimum alpha value for training search",
            &Self::training_alpha_min });
        visitor(RequestScalarFieldSpec{
            "training_alpha_max",
            "--training-alpha-max",
            "Maximum alpha value for training search",
            &Self::training_alpha_max });
        visitor(RequestScalarFieldSpec{
            "training_alpha_step",
            "--training-alpha-step",
            "Step size for alpha training search",
            &Self::training_alpha_step });
        visitor(RequestScalarFieldSpec{
            "asymmetry_flag",
            "--asymmetry",
            "Turn On/Off asymmetry flag",
            &Self::asymmetry_flag });
        visitor(RequestScalarFieldSpec{
            "sampling_size",
            "-s,--sampling",
            "Number of sampling points per atom",
            &Self::sampling_size });
        visitor(RequestScalarFieldSpec{
            "sampling_range_min",
            "--sampling-min",
            "Minimum sampling range",
            &Self::sampling_range_min });
        visitor(RequestScalarFieldSpec{
            "sampling_range_max",
            "--sampling-max",
            "Maximum sampling range",
            &Self::sampling_range_max });
        visitor(RequestScalarFieldSpec{
            "sampling_height",
            "--sampling-height",
            "Maximum sampling height",
            &Self::sampling_height });
        visitor(RequestScalarFieldSpec{
            "fit_range_min",
            "--fit-min",
            "Minimum fitting range",
            &Self::fit_range_min });
        visitor(RequestScalarFieldSpec{
            "fit_range_max",
            "--fit-max",
            "Maximum fitting range",
            &Self::fit_range_max });
        visitor(RequestScalarFieldSpec{
            "alpha_r",
            "--alpha-r",
            "Alpha value for R",
            &Self::alpha_r });
        visitor(RequestScalarFieldSpec{
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
        visitor(RequestPathFieldSpec{
            "database_path",
            "-d,--database",
            "Database file path",
            &Self::database_path });
        visitor(RequestEnumFieldSpec{
            "painter_choice",
            "-p,--painter",
            "Painter choice",
            &Self::painter_choice,
            true });
        visitor(RequestCsvListFieldSpec{
            "model_key_tag_list",
            "-k,--model-keylist",
            "List of model key tag to be display",
            &Self::model_key_tag_list,
            true });
        visitor(RequestRefGroupFieldSpec{
            "reference_model_groups",
            "--ref-group",
            "Reference group in the form group=key1,key2,...",
            &Self::reference_model_groups });
        visitor(RequestScalarFieldSpec{
            "pick_chain_id",
            "--pick-chain",
            "Pick chain ID",
            &Self::pick_chain_id });
        visitor(RequestScalarFieldSpec{
            "pick_residue",
            "--pick-residue",
            "Pick residue type",
            &Self::pick_residue });
        visitor(RequestScalarFieldSpec{
            "pick_element",
            "--pick-element",
            "Pick element type",
            &Self::pick_element });
        visitor(RequestScalarFieldSpec{
            "veto_chain_id",
            "--veto-chain",
            "Veto chain ID",
            &Self::veto_chain_id });
        visitor(RequestScalarFieldSpec{
            "veto_residue",
            "--veto-residue",
            "Veto residue type",
            &Self::veto_residue });
        visitor(RequestScalarFieldSpec{
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
        visitor(RequestPathFieldSpec{
            "database_path",
            "-d,--database",
            "Database file path",
            &Self::database_path });
        visitor(RequestEnumFieldSpec{
            "printer_choice",
            "-p,--printer",
            "Printer choice",
            &Self::printer_choice,
            true });
        visitor(RequestCsvListFieldSpec{
            "model_key_tag_list",
            "-k,--model-keylist",
            "List of model key tag to be display",
            &Self::model_key_tag_list,
            true });
        visitor(RequestPathFieldSpec{
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
        visitor(RequestPathFieldSpec{
            "model_file_path",
            "-a,--model",
            "Model file path",
            &Self::model_file_path,
            true });
        visitor(RequestScalarFieldSpec{
            "map_file_name",
            "-n,--name",
            "File name for output map files",
            &Self::map_file_name });
        visitor(RequestEnumFieldSpec{
            "potential_model_choice",
            "--potential-model",
            "Atomic potential model option",
            &Self::potential_model_choice });
        visitor(RequestEnumFieldSpec{
            "partial_charge_choice",
            "--charge",
            "Partial charge table option",
            &Self::partial_charge_choice });
        visitor(RequestScalarFieldSpec{
            "cutoff_distance",
            "-c,--cut-off",
            "Cutoff distance",
            &Self::cutoff_distance });
        visitor(RequestScalarFieldSpec{
            "grid_spacing",
            "-g,--grid-spacing",
            "Grid spacing",
            &Self::grid_spacing });
        visitor(RequestCsvListFieldSpec{
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
        visitor(RequestEnumFieldSpec{
            "tester_choice",
            "-t,--tester",
            "Tester option",
            &Self::tester_choice });
        visitor(RequestScalarFieldSpec{
            "fit_range_min",
            "--fit-min",
            "Minimum fitting range",
            &Self::fit_range_min });
        visitor(RequestScalarFieldSpec{
            "fit_range_max",
            "--fit-max",
            "Maximum fitting range",
            &Self::fit_range_max });
        visitor(RequestScalarFieldSpec{
            "alpha_r",
            "--alpha-r",
            "Alpha value for R",
            &Self::alpha_r });
        visitor(RequestScalarFieldSpec{
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
        visitor(RequestPathFieldSpec{
            "model_file_path",
            "-a,--model",
            "Model file path",
            &Self::model_file_path,
            true });
        visitor(RequestPathFieldSpec{
            "map_file_path",
            "-m,--map",
            "Map file path",
            &Self::map_file_path,
            true });
        visitor(RequestScalarFieldSpec{
            "atom_serial_id",
            "-i,--atom-id",
            "Atom serial ID for visualization",
            &Self::atom_serial_id });
        visitor(RequestScalarFieldSpec{
            "sampling_size",
            "-s,--sampling",
            "Number of sampling points per atom",
            &Self::sampling_size });
        visitor(RequestScalarFieldSpec{
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
        visitor(RequestPathFieldSpec{
            "map_file_path",
            "-m,--map",
            "Map file path",
            &Self::map_file_path,
            true });
        visitor(RequestScalarFieldSpec{
            "iteration_count",
            "--iter",
            "Iteration count for estimation",
            &Self::iteration_count });
        visitor(RequestScalarFieldSpec{
            "knn_size",
            "--knn",
            "KNN size for estimation",
            &Self::knn_size });
        visitor(RequestScalarFieldSpec{
            "alpha",
            "--alpha",
            "Alpha value for robust regression",
            &Self::alpha });
        visitor(RequestScalarFieldSpec{
            "threshold_ratio",
            "--threshold",
            "Ratio of threshold of map values",
            &Self::threshold_ratio });
        visitor(RequestScalarFieldSpec{
            "dedup_tolerance",
            "--dedup-tolerance",
            "Tolerance for deduplicating points",
            &Self::dedup_tolerance });
    }
};
#endif

} // namespace rhbm_gem::command_internal
