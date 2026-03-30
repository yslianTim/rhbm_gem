#pragma once

#include <filesystem>
#include <string>
#include <vector>

#include <rhbm_gem/core/command/CommandContract.hpp>
#include <rhbm_gem/core/command/CommandEnumClass.hpp>

namespace CLI
{
    class App;
}

namespace rhbm_gem {

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
    std::filesystem::path database_path{ GetDefaultDatabasePath() };
    std::filesystem::path folder_path{ "" };
};

struct PotentialAnalysisRequest
{
    CommonCommandRequest common{};
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
};

struct PotentialDisplayRequest
{
    CommonCommandRequest common{};
    PainterType painter_choice{ PainterType::MODEL };
    std::string model_key_tag_list{};
    std::string ref_model_key_tag_list{};
    std::string pick_chain_id{};
    std::string veto_chain_id{};
    std::string pick_residue{};
    std::string veto_residue{};
    std::string pick_element{};
    std::string veto_element{};
};

struct ResultDumpRequest
{
    CommonCommandRequest common{};
    PrinterType printer_choice{ PrinterType::GAUS_ESTIMATES };
    std::string model_key_tag_list{};
    std::filesystem::path map_file_path{};
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
    std::string blurring_width_list{ "1.50" };
};

struct HRLModelTestRequest
{
    CommonCommandRequest common{};
    TesterType tester_choice{ TesterType::BENCHMARK };
    double fit_range_min{ 0.0 };
    double fit_range_max{ 1.0 };
    double alpha_r{ 0.1 };
    double alpha_g{ 0.2 };
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
};

struct PositionEstimationRequest
{
    CommonCommandRequest common{};
    std::filesystem::path map_file_path{};
    int iteration_count{ 15 };
    int knn_size{ 20 };
    double alpha{ 2.0 };
    double threshold_ratio{ 0.01 };
    double dedup_tolerance{ 1.0e-2 };
};
#endif

void ConfigureCommandCli(::CLI::App & app);

#define RHBM_GEM_COMMAND(COMMAND_ID, CLI_NAME, DESCRIPTION, PROFILE)                           \
    ExecutionReport Run##COMMAND_ID(const COMMAND_ID##Request & request);
#include <rhbm_gem/core/command/CommandList.def>
#undef RHBM_GEM_COMMAND

} // namespace rhbm_gem
