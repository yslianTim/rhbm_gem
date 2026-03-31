#pragma once

#include <cstdint>
#include <filesystem>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

#include <rhbm_gem/core/command/CommandEnums.hpp>
#include <rhbm_gem/core/command/CommandPaths.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>

namespace rhbm_gem {

enum class ValidationPhase : std::uint8_t
{
    Parse = 0,
    Prepare = 1
};

struct ValidationIssue
{
    std::string option_name;
    ValidationPhase phase;
    LogLevel level;
    std::string message;
    bool auto_corrected{ false };
};

enum class CommandOutcome : std::uint8_t
{
    ValidationFailed = 0,
    ExecutionFailed  = 1,
    Succeeded        = 2
};

struct CommandResult
{
    CommandOutcome outcome{ CommandOutcome::ValidationFailed };
    std::vector<ValidationIssue> issues{};
};

struct CommandInfo
{
    std::string_view name;
    std::string_view description;
};

const std::vector<CommandInfo> & ListCommands();

struct CommandRequestBase
{
    int job_count{ 1 };
    int verbosity{ 3 };
    std::filesystem::path output_dir{};
};

struct PotentialAnalysisRequest : public CommandRequestBase
{
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
};

struct PotentialDisplayRequest : public CommandRequestBase
{
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
};

struct ResultDumpRequest : public CommandRequestBase
{
    std::filesystem::path database_path{ GetDefaultDatabasePath() };
    PrinterType printer_choice{ PrinterType::GAUS_ESTIMATES };
    std::vector<std::string> model_key_tag_list{};
    std::filesystem::path map_file_path{};
};

struct MapSimulationRequest : public CommandRequestBase
{
    std::filesystem::path model_file_path{};
    std::string map_file_name{ "sim_map" };
    PotentialModel potential_model_choice{ PotentialModel::FIVE_GAUS_CHARGE };
    PartialCharge partial_charge_choice{ PartialCharge::PARTIAL };
    double cutoff_distance{ 5.0 };
    double grid_spacing{ 0.5 };
    std::vector<double> blurring_width_list{ 1.50 };
};

struct HRLModelTestRequest : public CommandRequestBase
{
    TesterType tester_choice{ TesterType::BENCHMARK };
    double fit_range_min{ 0.0 };
    double fit_range_max{ 1.0 };
    double alpha_r{ 0.1 };
    double alpha_g{ 0.2 };
};

#ifdef RHBM_GEM_ENABLE_EXPERIMENTAL_FEATURE
struct MapVisualizationRequest : public CommandRequestBase
{
    std::filesystem::path model_file_path{};
    std::filesystem::path map_file_path{};
    int atom_serial_id{ 1 };
    int sampling_size{ 100 };
    double window_size{ 5.0 };
};

struct PositionEstimationRequest : public CommandRequestBase
{
    std::filesystem::path map_file_path{};
    int iteration_count{ 15 };
    std::size_t knn_size{ 20 };
    double alpha{ 2.0 };
    double threshold_ratio{ 0.01 };
    double dedup_tolerance{ 1.0e-2 };
};
#endif

CommandResult RunPotentialAnalysis(const PotentialAnalysisRequest & request);
CommandResult RunPotentialDisplay(const PotentialDisplayRequest & request);
CommandResult RunResultDump(const ResultDumpRequest & request);
CommandResult RunMapSimulation(const MapSimulationRequest & request);
CommandResult RunHRLModelTest(const HRLModelTestRequest & request);

#ifdef RHBM_GEM_ENABLE_EXPERIMENTAL_FEATURE
CommandResult RunMapVisualization(const MapVisualizationRequest & request);
CommandResult RunPositionEstimation(const PositionEstimationRequest & request);
#endif

} // namespace rhbm_gem
