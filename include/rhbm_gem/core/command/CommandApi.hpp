#pragma once

#include <filesystem>
#include <string>
#include <unordered_map>
#include <vector>

#include <rhbm_gem/core/command/CommandContract.hpp>
#include <rhbm_gem/core/command/CommandEnumClass.hpp>

namespace CLI
{
    class App;
}

namespace rhbm_gem {

template <typename Owner, typename FieldType>
struct RequestObjectFieldSpec
{
    const char * python_name;
    FieldType Owner::* member;
};

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
constexpr RequestObjectFieldSpec<Owner, FieldType> MakeObjectField(
    const char * python_name,
    FieldType Owner::* member)
{
    return RequestObjectFieldSpec<Owner, FieldType>{ python_name, member };
}

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

struct ExecutionReport
{
    bool prepared{ false };
    bool executed{ false };
    std::vector<ValidationIssue> validation_issues{};
};

struct CommonCommandRequest
{
    int thread_size{ 1 };
    int verbose_level{ 3 };
    std::filesystem::path folder_path{};

    template <typename Visitor>
    static void VisitFields(Visitor && visitor)
    {
        using Self = CommonCommandRequest;
        visitor(MakeScalarField<Self>(
            "thread_size",
            "-j,--jobs",
            "Number of threads",
            &Self::thread_size));
        visitor(MakeScalarField<Self>(
            "verbose_level",
            "-v,--verbose",
            "Verbose level",
            &Self::verbose_level));
        visitor(MakePathField<Self>(
            "folder_path",
            "-o,--folder",
            "folder path for output files",
            &Self::folder_path));
    }
};

struct PotentialAnalysisRequest
{
    CommonCommandRequest common{};
    std::filesystem::path database_path{ GetDefaultDatabasePath() };
    std::filesystem::path model_file_path{};
    std::filesystem::path map_file_path{};
    bool simulation_flag{ false };
    double simulated_map_resolution{ 0.0 };
    std::string saved_key_tag{ "model" };
    std::filesystem::path training_report_dir{};
    bool training_alpha_flag{ false };
    bool asymmetry_flag{ false };
    int sampling_size{ 1500 };
    double sampling_range_min{ 0.0 };
    double sampling_range_max{ 1.5 };
    double sampling_height{ 0.1 };
    double fit_range_min{ 0.0 };
    double fit_range_max{ 1.0 };
    double alpha_r{ 0.1 };
    double alpha_g{ 0.2 };

    template <typename Visitor>
    static void VisitFields(Visitor && visitor)
    {
        using Self = PotentialAnalysisRequest;
        visitor(MakeObjectField<Self>("common", &Self::common));
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

struct PotentialDisplayRequest
{
    CommonCommandRequest common{};
    std::filesystem::path database_path{ GetDefaultDatabasePath() };
    PainterType painter_choice{ PainterType::MODEL };
    std::vector<std::string> model_key_tag_list{};
    std::unordered_map<std::string, std::vector<std::string>> reference_model_groups{};
    std::string pick_chain_id{};
    std::string veto_chain_id{};
    std::string pick_residue{};
    std::string veto_residue{};
    std::string pick_element{};
    std::string veto_element{};

    template <typename Visitor>
    static void VisitFields(Visitor && visitor)
    {
        using Self = PotentialDisplayRequest;
        visitor(MakeObjectField<Self>("common", &Self::common));
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

struct ResultDumpRequest
{
    CommonCommandRequest common{};
    std::filesystem::path database_path{ GetDefaultDatabasePath() };
    PrinterType printer_choice{ PrinterType::GAUS_ESTIMATES };
    std::vector<std::string> model_key_tag_list{};
    std::filesystem::path map_file_path{};

    template <typename Visitor>
    static void VisitFields(Visitor && visitor)
    {
        using Self = ResultDumpRequest;
        visitor(MakeObjectField<Self>("common", &Self::common));
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

struct MapSimulationRequest
{
    CommonCommandRequest common{};
    std::filesystem::path model_file_path{};
    std::string map_file_name{ "sim_map" };
    PotentialModel potential_model_choice{ PotentialModel::FIVE_GAUS_CHARGE };
    PartialCharge partial_charge_choice{ PartialCharge::PARTIAL };
    double cutoff_distance{ 5.0 };
    double grid_spacing{ 0.5 };
    std::vector<double> blurring_width_list{ 1.50 };

    template <typename Visitor>
    static void VisitFields(Visitor && visitor)
    {
        using Self = MapSimulationRequest;
        visitor(MakeObjectField<Self>("common", &Self::common));
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

struct HRLModelTestRequest
{
    CommonCommandRequest common{};
    TesterType tester_choice{ TesterType::BENCHMARK };
    double fit_range_min{ 0.0 };
    double fit_range_max{ 1.0 };
    double alpha_r{ 0.1 };
    double alpha_g{ 0.2 };

    template <typename Visitor>
    static void VisitFields(Visitor && visitor)
    {
        using Self = HRLModelTestRequest;
        visitor(MakeObjectField<Self>("common", &Self::common));
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
struct MapVisualizationRequest
{
    CommonCommandRequest common{};
    std::filesystem::path model_file_path{};
    std::filesystem::path map_file_path{};
    int atom_serial_id{ 1 };
    int sampling_size{ 100 };
    double window_size{ 5.0 };

    template <typename Visitor>
    static void VisitFields(Visitor && visitor)
    {
        using Self = MapVisualizationRequest;
        visitor(MakeObjectField<Self>("common", &Self::common));
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

struct PositionEstimationRequest
{
    CommonCommandRequest common{};
    std::filesystem::path map_file_path{};
    int iteration_count{ 15 };
    std::size_t knn_size{ 20 };
    double alpha{ 2.0 };
    double threshold_ratio{ 0.01 };
    double dedup_tolerance{ 1.0e-2 };

    template <typename Visitor>
    static void VisitFields(Visitor && visitor)
    {
        using Self = PositionEstimationRequest;
        visitor(MakeObjectField<Self>("common", &Self::common));
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

void ConfigureCommandCli(::CLI::App & app);

#define RHBM_GEM_COMMAND(COMMAND_ID, CLI_NAME, DESCRIPTION)                                    \
    ExecutionReport Run##COMMAND_ID(const COMMAND_ID##Request & request);
#include <rhbm_gem/core/command/CommandList.def>
#undef RHBM_GEM_COMMAND

} // namespace rhbm_gem
