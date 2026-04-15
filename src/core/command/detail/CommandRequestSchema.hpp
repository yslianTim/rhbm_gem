#pragma once

#include <filesystem>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <rhbm_gem/core/command/CommandApi.hpp>

namespace rhbm_gem::internal {

template <typename Owner, typename FieldType>
struct RequestScalarFieldSpec
{
    const char * python_name;
    const char * cli_flags;
    const char * help;
    bool required;
    FieldType Owner::* member;
};

template <typename Owner>
struct RequestPathFieldSpec
{
    const char * python_name;
    const char * cli_flags;
    const char * help;
    bool required;
    std::filesystem::path Owner::* member;
};

template <typename Owner, typename EnumType>
struct RequestEnumFieldSpec
{
    const char * python_name;
    const char * cli_flags;
    const char * help;
    bool required;
    EnumType Owner::* member;
};

template <typename Owner, typename ElementType>
struct RequestCsvListFieldSpec
{
    const char * python_name;
    const char * cli_flags;
    const char * help;
    bool required;
    char delimiter;
    std::vector<ElementType> Owner::* member;
};

template <typename Owner>
struct RequestRefGroupFieldSpec
{
    const char * python_name;
    const char * cli_flags;
    const char * help;
    bool required;
    char assignment_delimiter;
    char item_delimiter;
    std::unordered_map<std::string, std::vector<std::string>> Owner::* member;
};

template <typename Owner, typename FieldType>
constexpr RequestScalarFieldSpec<Owner, FieldType> MakeScalarField(
    const char * python_name,
    const char * cli_flags,
    const char * help,
    FieldType Owner::* member,
    bool required = false)
{
    return RequestScalarFieldSpec<Owner, FieldType>{ python_name, cli_flags, help, required, member };
}

template <typename Owner>
constexpr RequestPathFieldSpec<Owner> MakePathField(
    const char * python_name,
    const char * cli_flags,
    const char * help,
    std::filesystem::path Owner::* member,
    bool required = false)
{
    return RequestPathFieldSpec<Owner>{ python_name, cli_flags, help, required, member };
}

template <typename Owner, typename EnumType>
constexpr RequestEnumFieldSpec<Owner, EnumType> MakeEnumField(
    const char * python_name,
    const char * cli_flags,
    const char * help,
    EnumType Owner::* member,
    bool required = false)
{
    return RequestEnumFieldSpec<Owner, EnumType>{ python_name, cli_flags, help, required, member };
}

template <typename Owner, typename ElementType>
constexpr RequestCsvListFieldSpec<Owner, ElementType> MakeCsvListField(
    const char * python_name,
    const char * cli_flags,
    const char * help,
    std::vector<ElementType> Owner::* member,
    bool required = false,
    char delimiter = ',')
{
    return RequestCsvListFieldSpec<Owner, ElementType>{
        python_name, cli_flags, help, required, delimiter, member
    };
}

template <typename Owner>
constexpr RequestRefGroupFieldSpec<Owner> MakeRefGroupField(
    const char * python_name,
    const char * cli_flags,
    const char * help,
    std::unordered_map<std::string, std::vector<std::string>> Owner::* member,
    bool required = false,
    char assignment_delimiter = '=',
    char item_delimiter = ',')
{
    return RequestRefGroupFieldSpec<Owner>{
        python_name,
        cli_flags,
        help,
        required,
        assignment_delimiter,
        item_delimiter,
        member
    };
}

template <typename Request>
struct CommandRequestSchema;

template <typename Visitor>
void VisitBaseRequestFields(Visitor && visitor)
{
    using Base = CommandRequestBase;
    visitor(MakeScalarField<Base>(
        "job_count",
        "-j,--jobs",
        "Number of threads",
        &Base::job_count));
    visitor(MakeScalarField<Base>(
        "verbosity",
        "-v,--verbose",
        "Verbose level",
        &Base::verbosity));
    visitor(MakePathField<Base>(
        "output_dir",
        "-o,--folder",
        "folder path for output files",
        &Base::output_dir));
}

template <typename Request, typename Visitor>
void VisitRequestFields(Visitor && visitor)
{
    CommandRequestSchema<Request>::Visit(std::forward<Visitor>(visitor));
}

template <>
struct CommandRequestSchema<PotentialAnalysisRequest>
{
    template <typename Visitor>
    static void Visit(Visitor && visitor)
    {
        using Self = PotentialAnalysisRequest;
        visitor(MakePathField<Self>(
            "database_path",
            "-d,--database",
            "Database file path",
            &Self::database_path));
        visitor(MakePathField<Self>(
            "model_file_path",
            "-a,--model",
            "Model file path",
            &Self::model_file_path,
            true));
        visitor(MakePathField<Self>(
            "map_file_path",
            "-m,--map",
            "Map file path",
            &Self::map_file_path,
            true));
        visitor(MakePathField<Self>(
            "training_report_dir",
            "--training-report-dir",
            "Optional output directory for alpha training reports",
            &Self::training_report_dir));
        visitor(MakeScalarField<Self>(
            "simulation_flag",
            "--simulation",
            "Simulation flag",
            &Self::simulation_flag));
        visitor(MakeScalarField<Self>(
            "simulated_map_resolution",
            "-r,--sim-resolution",
            "Set simulated map's resolution (blurring width)",
            &Self::simulated_map_resolution));
        visitor(MakeScalarField<Self>(
            "saved_key_tag",
            "-k,--save-key",
            "New key tag for saving ModelObject results into database",
            &Self::saved_key_tag));
        visitor(MakeScalarField<Self>(
            "training_alpha_flag",
            "--training-alpha",
            "Turn On/Off alpha training flag",
            &Self::training_alpha_flag));
        visitor(MakeScalarField<Self>(
            "training_alpha_min",
            "--training-alpha-min",
            "Minimum alpha value for training search",
            &Self::training_alpha_min));
        visitor(MakeScalarField<Self>(
            "training_alpha_max",
            "--training-alpha-max",
            "Maximum alpha value for training search",
            &Self::training_alpha_max));
        visitor(MakeScalarField<Self>(
            "training_alpha_step",
            "--training-alpha-step",
            "Step size for alpha training search",
            &Self::training_alpha_step));
        visitor(MakeScalarField<Self>(
            "asymmetry_flag",
            "--asymmetry",
            "Turn On/Off asymmetry flag",
            &Self::asymmetry_flag));
        visitor(MakeScalarField<Self>(
            "sampling_size",
            "-s,--sampling",
            "Number of sampling points per atom",
            &Self::sampling_size));
        visitor(MakeScalarField<Self>(
            "sampling_range_min",
            "--sampling-min",
            "Minimum sampling range",
            &Self::sampling_range_min));
        visitor(MakeScalarField<Self>(
            "sampling_range_max",
            "--sampling-max",
            "Maximum sampling range",
            &Self::sampling_range_max));
        visitor(MakeScalarField<Self>(
            "sampling_height",
            "--sampling-height",
            "Maximum sampling height",
            &Self::sampling_height));
        visitor(MakeScalarField<Self>(
            "fit_range_min",
            "--fit-min",
            "Minimum fitting range",
            &Self::fit_range_min));
        visitor(MakeScalarField<Self>(
            "fit_range_max",
            "--fit-max",
            "Maximum fitting range",
            &Self::fit_range_max));
        visitor(MakeScalarField<Self>(
            "alpha_r",
            "--alpha-r",
            "Alpha value for R",
            &Self::alpha_r));
        visitor(MakeScalarField<Self>(
            "alpha_g",
            "--alpha-g",
            "Alpha value for G",
            &Self::alpha_g));
    }
};

template <>
struct CommandRequestSchema<PotentialDisplayRequest>
{
    template <typename Visitor>
    static void Visit(Visitor && visitor)
    {
        using Self = PotentialDisplayRequest;
        visitor(MakePathField<Self>(
            "database_path",
            "-d,--database",
            "Database file path",
            &Self::database_path));
        visitor(MakeEnumField<Self>(
            "painter_choice",
            "-p,--painter",
            "Painter choice",
            &Self::painter_choice,
            true));
        visitor(MakeCsvListField<Self, std::string>(
            "model_key_tag_list",
            "-k,--model-keylist",
            "List of model key tag to be display",
            &Self::model_key_tag_list,
            true));
        visitor(MakeRefGroupField<Self>(
            "reference_model_groups",
            "--ref-group",
            "Reference group in the form group=key1,key2,...",
            &Self::reference_model_groups));
        visitor(MakeScalarField<Self>(
            "pick_chain_id",
            "--pick-chain",
            "Pick chain ID",
            &Self::pick_chain_id));
        visitor(MakeScalarField<Self>(
            "pick_residue",
            "--pick-residue",
            "Pick residue type",
            &Self::pick_residue));
        visitor(MakeScalarField<Self>(
            "pick_element",
            "--pick-element",
            "Pick element type",
            &Self::pick_element));
        visitor(MakeScalarField<Self>(
            "veto_chain_id",
            "--veto-chain",
            "Veto chain ID",
            &Self::veto_chain_id));
        visitor(MakeScalarField<Self>(
            "veto_residue",
            "--veto-residue",
            "Veto residue type",
            &Self::veto_residue));
        visitor(MakeScalarField<Self>(
            "veto_element",
            "--veto-element",
            "Veto element type",
            &Self::veto_element));
    }
};

template <>
struct CommandRequestSchema<ResultDumpRequest>
{
    template <typename Visitor>
    static void Visit(Visitor && visitor)
    {
        using Self = ResultDumpRequest;
        visitor(MakePathField<Self>(
            "database_path",
            "-d,--database",
            "Database file path",
            &Self::database_path));
        visitor(MakeEnumField<Self>(
            "printer_choice",
            "-p,--printer",
            "Printer choice",
            &Self::printer_choice,
            true));
        visitor(MakeCsvListField<Self, std::string>(
            "model_key_tag_list",
            "-k,--model-keylist",
            "List of model key tag to be display",
            &Self::model_key_tag_list,
            true));
        visitor(MakePathField<Self>(
            "map_file_path",
            "-m,--map",
            "Map file path",
            &Self::map_file_path));
    }
};

template <>
struct CommandRequestSchema<MapSimulationRequest>
{
    template <typename Visitor>
    static void Visit(Visitor && visitor)
    {
        using Self = MapSimulationRequest;
        visitor(MakePathField<Self>(
            "model_file_path",
            "-a,--model",
            "Model file path",
            &Self::model_file_path,
            true));
        visitor(MakeScalarField<Self>(
            "map_file_name",
            "-n,--name",
            "File name for output map files",
            &Self::map_file_name));
        visitor(MakeEnumField<Self>(
            "potential_model_choice",
            "--potential-model",
            "Atomic potential model option",
            &Self::potential_model_choice));
        visitor(MakeEnumField<Self>(
            "partial_charge_choice",
            "--charge",
            "Partial charge table option",
            &Self::partial_charge_choice));
        visitor(MakeScalarField<Self>(
            "cutoff_distance",
            "-c,--cut-off",
            "Cutoff distance",
            &Self::cutoff_distance));
        visitor(MakeScalarField<Self>(
            "grid_spacing",
            "-g,--grid-spacing",
            "Grid spacing",
            &Self::grid_spacing));
        visitor(MakeCsvListField<Self, double>(
            "blurring_width_list",
            "--blurring-width",
            "Blurring width (list) setting",
            &Self::blurring_width_list));
    }
};

template <>
struct CommandRequestSchema<HRLModelTestRequest>
{
    template <typename Visitor>
    static void Visit(Visitor && visitor)
    {
        using Self = HRLModelTestRequest;
        visitor(MakeEnumField<Self>(
            "tester_choice",
            "-t,--tester",
            "Tester option",
            &Self::tester_choice));
        visitor(MakeScalarField<Self>(
            "fit_range_min",
            "--fit-min",
            "Minimum fitting range",
            &Self::fit_range_min));
        visitor(MakeScalarField<Self>(
            "fit_range_max",
            "--fit-max",
            "Maximum fitting range",
            &Self::fit_range_max));
        visitor(MakeScalarField<Self>(
            "alpha_r",
            "--alpha-r",
            "Alpha value for R",
            &Self::alpha_r));
        visitor(MakeScalarField<Self>(
            "alpha_g",
            "--alpha-g",
            "Alpha value for G",
            &Self::alpha_g));
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
        visitor(MakePathField<Self>(
            "model_file_path",
            "-a,--model",
            "Model file path",
            &Self::model_file_path,
            true));
        visitor(MakePathField<Self>(
            "map_file_path",
            "-m,--map",
            "Map file path",
            &Self::map_file_path,
            true));
        visitor(MakeScalarField<Self>(
            "atom_serial_id",
            "-i,--atom-id",
            "Atom serial ID for visualization",
            &Self::atom_serial_id));
        visitor(MakeScalarField<Self>(
            "sampling_size",
            "-s,--sampling",
            "Number of sampling points per atom",
            &Self::sampling_size));
        visitor(MakeScalarField<Self>(
            "window_size",
            "--window-size",
            "Window size for sampling",
            &Self::window_size));
    }
};

template <>
struct CommandRequestSchema<PositionEstimationRequest>
{
    template <typename Visitor>
    static void Visit(Visitor && visitor)
    {
        using Self = PositionEstimationRequest;
        visitor(MakePathField<Self>(
            "map_file_path",
            "-m,--map",
            "Map file path",
            &Self::map_file_path,
            true));
        visitor(MakeScalarField<Self>(
            "iteration_count",
            "--iter",
            "Iteration count for estimation",
            &Self::iteration_count));
        visitor(MakeScalarField<Self>(
            "knn_size",
            "--knn",
            "KNN size for estimation",
            &Self::knn_size));
        visitor(MakeScalarField<Self>(
            "alpha",
            "--alpha",
            "Alpha value for robust regression",
            &Self::alpha));
        visitor(MakeScalarField<Self>(
            "threshold_ratio",
            "--threshold",
            "Ratio of threshold of map values",
            &Self::threshold_ratio));
        visitor(MakeScalarField<Self>(
            "dedup_tolerance",
            "--dedup-tolerance",
            "Tolerance for deduplicating points",
            &Self::dedup_tolerance));
    }
};
#endif

} // namespace rhbm_gem::internal
