#pragma once

#include <filesystem>
#include <string>
#include <vector>

#include <rhbm_gem/core/command/CommandBase.hpp>
#include <rhbm_gem/core/command/OptionEnumClass.hpp>

namespace rhbm_gem::gui {

struct CommonExecutionRequest
{
    int thread_size{ 1 };
    int verbose_level{ 3 };
    std::filesystem::path database_path{ GetDefaultDatabasePath() };
    std::filesystem::path folder_path{ "" };
};

struct MapSimulationRequest
{
    CommonExecutionRequest common{};
    std::filesystem::path model_file_path{};
    std::string map_file_name{ "sim_map" };
    PotentialModel potential_model_choice{ PotentialModel::FIVE_GAUS_CHARGE };
    PartialCharge partial_charge_choice{ PartialCharge::PARTIAL };
    double cutoff_distance{ 5.0 };
    double grid_spacing{ 0.5 };
    std::string blurring_width_list{ "1.50" };
};

struct PotentialAnalysisRequest
{
    CommonExecutionRequest common{};
    std::filesystem::path model_file_path{};
    std::filesystem::path map_file_path{};
    bool simulation_flag{ false };
    double simulated_map_resolution{ 0.0 };
    std::string saved_key_tag{ "model" };
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

struct ResultDumpRequest
{
    CommonExecutionRequest common{};
    PrinterType printer_choice{ PrinterType::GAUS_ESTIMATES };
    std::string model_key_tag_list{};
    std::filesystem::path map_file_path{};
};

struct ExecutionResult
{
    bool prepared{ false };
    bool executed{ false };
    std::vector<ValidationIssue> validation_issues{};
};

class GuiCommandExecutor
{
public:
    static ExecutionResult ExecuteMapSimulation(const MapSimulationRequest & request);
    static ExecutionResult ExecutePotentialAnalysis(const PotentialAnalysisRequest & request);
    static ExecutionResult ExecuteResultDump(const ResultDumpRequest & request);
};

} // namespace rhbm_gem::gui
