#pragma once

#include <array>
#include <cstddef>
#include <cstdlib>
#include <filesystem>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace rhbm_gem {

enum class PainterType : int
{
    GAUS       = 0,
    MODEL      = 1,
    COMPARISON = 2,
    DEMO       = 3,
    ATOM       = 4
};

enum class PrinterType : int
{
    ATOM_POSITION  = 0,
    MAP_VALUE      = 1,
    GAUS_ESTIMATES = 2,
    ATOM_OUTLIER   = 3
};

enum class PotentialModel : int
{
    SINGLE_GAUS      = 0,
    FIVE_GAUS_CHARGE = 1,
    SINGLE_GAUS_USER = 2
};

enum class PartialCharge : int
{
    NEUTRAL = 0,
    PARTIAL = 1,
    AMBER   = 2
};

enum class TesterType : int
{
    BENCHMARK          = 0,
    DATA_OUTLIER       = 1,
    MEMBER_OUTLIER     = 2,
    MODEL_ALPHA_DATA   = 3,
    MODEL_ALPHA_MEMBER = 4,
    NEIGHBOR_DISTANCE  = 5
};

inline std::filesystem::path GetDefaultDataRootPath()
{
    if (const char * configured_root{ std::getenv("RHBM_GEM_DATA_DIR") };
        configured_root != nullptr && configured_root[0] != '\0')
    {
        return std::filesystem::path(configured_root);
    }

    if (const char * home{ std::getenv("HOME") };
        home != nullptr && home[0] != '\0')
    {
        return std::filesystem::path(home) / ".rhbmgem" / "data";
    }

    return std::filesystem::path(".rhbmgem") / "data";
}

inline std::filesystem::path GetDefaultDatabasePath()
{
    return GetDefaultDataRootPath() / "database.sqlite";
}

struct ValidationIssue
{
    std::string option_name;
    std::string message;
};

struct CommandResult
{
    bool succeeded{ false };
    std::vector<ValidationIssue> issues{};
};

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
    double training_alpha_min{ 0.0 };
    double training_alpha_max{ 1.0 };
    double training_alpha_step{ 0.1 };
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

struct RHBMTestRequest : public CommandRequestBase
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

namespace internal {

template <typename EnumType, std::size_t AliasCount>
struct CommandEnumDefinition
{
    EnumType value;
    std::string_view binding_token;
    std::array<std::string_view, AliasCount> cli_aliases;
};

template <typename EnumType>
struct CommandEnumTraits;

template <>
struct CommandEnumTraits<PainterType>
{
    inline static constexpr std::array<CommandEnumDefinition<PainterType, 2>, 5> kOptions{{
        { PainterType::GAUS, "GAUS", { "0", "gaus" } },
        { PainterType::ATOM, "ATOM", { "4", "atom" } },
        { PainterType::MODEL, "MODEL", { "1", "model" } },
        { PainterType::COMPARISON, "COMPARISON", { "2", "comparison" } },
        { PainterType::DEMO, "DEMO", { "3", "demo" } }
    }};
};

template <>
struct CommandEnumTraits<PrinterType>
{
    inline static constexpr std::array<CommandEnumDefinition<PrinterType, 2>, 4> kOptions{{
        { PrinterType::ATOM_POSITION, "ATOM_POSITION", { "0", "atom_pos" } },
        { PrinterType::MAP_VALUE, "MAP_VALUE", { "1", "map" } },
        { PrinterType::GAUS_ESTIMATES, "GAUS_ESTIMATES", { "2", "gaus" } },
        { PrinterType::ATOM_OUTLIER, "ATOM_OUTLIER", { "3", "atom_out" } }
    }};
};

template <>
struct CommandEnumTraits<PotentialModel>
{
    inline static constexpr std::array<CommandEnumDefinition<PotentialModel, 2>, 3> kOptions{{
        { PotentialModel::SINGLE_GAUS, "SINGLE_GAUS", { "0", "single" } },
        { PotentialModel::FIVE_GAUS_CHARGE, "FIVE_GAUS_CHARGE", { "1", "five" } },
        { PotentialModel::SINGLE_GAUS_USER, "SINGLE_GAUS_USER", { "2", "user" } }
    }};
};

template <>
struct CommandEnumTraits<PartialCharge>
{
    inline static constexpr std::array<CommandEnumDefinition<PartialCharge, 2>, 3> kOptions{{
        { PartialCharge::NEUTRAL, "NEUTRAL", { "0", "neutral" } },
        { PartialCharge::PARTIAL, "PARTIAL", { "1", "partial" } },
        { PartialCharge::AMBER, "AMBER", { "2", "amber" } }
    }};
};

template <>
struct CommandEnumTraits<TesterType>
{
    inline static constexpr std::array<CommandEnumDefinition<TesterType, 2>, 6> kOptions{{
        { TesterType::BENCHMARK, "BENCHMARK", { "0", "benchmark" } },
        { TesterType::DATA_OUTLIER, "DATA_OUTLIER", { "1", "data_outlier" } },
        { TesterType::MEMBER_OUTLIER, "MEMBER_OUTLIER", { "2", "member_outlier" } },
        { TesterType::MODEL_ALPHA_DATA, "MODEL_ALPHA_DATA", { "3", "alpha_data" } },
        { TesterType::MODEL_ALPHA_MEMBER, "MODEL_ALPHA_MEMBER", { "4", "alpha_member" } },
        { TesterType::NEIGHBOR_DISTANCE, "NEIGHBOR_DISTANCE", { "5", "neighbor_distance" } }
    }};
};

} // namespace internal

} // namespace rhbm_gem
