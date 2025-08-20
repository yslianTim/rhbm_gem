#include "PositionEstimationCommand.hpp"
#include "MapObject.hpp"
#include "DataObjectManager.hpp"
#include "KDTreeAlgorithm.hpp"
#include "ScopeTimer.hpp"
#include "ArrayStats.hpp"
#include "Logger.hpp"
#include "CommandRegistry.hpp"
#include "ChimeraXHelper.hpp"
#include "FilePathHelper.hpp"

#include <random>
#include <memory>
#include <vector>
#include <array>
#include <unordered_set>
#include <cmath>
#include <cstdint>

#ifdef USE_OPENMP
#include <omp.h>
#endif

namespace {
CommandRegistrar<PositionEstimationCommand> registrar_model_test{
    "position_estimation",
    "Run atom position estimation"};
}

PositionEstimationCommand::PositionEstimationCommand(void) :
    CommandBase(), m_options{}, m_selected_voxel_list{}, m_position_list{},
    m_kd_tree_root{ nullptr }
{
    Logger::Log(LogLevel::Debug, "PositionEstimationCommand::PositionEstimationCommand() called.");
}

void PositionEstimationCommand::RegisterCLIOptionsExtend(CLI::App * cmd)
{
    Logger::Log(LogLevel::Debug, "PositionEstimationCommand::RegisterCLIOptionsExtend() called.");
    cmd->add_option("-m,--map", m_options.map_file_path,
        "Map file path")->required();
    cmd->add_option("--alpha", m_options.alpha,
        "Alpha value for robust regression")->default_val(m_options.alpha);
    cmd->add_option("--iter", m_options.iteration_count,
        "Iteration count for estimation")->default_val(m_options.iteration_count);
    cmd->add_option("--knn", m_options.knn_size,
        "KNN size for estimation")->default_val(m_options.knn_size);
    cmd->add_option("--threshold", m_options.threshold_ratio,
        "Ratio of threshold of map values")->default_val(m_options.threshold_ratio);
}

bool PositionEstimationCommand::Execute(void)
{
    Logger::Log(LogLevel::Debug, "PositionEstimationCommand::Execute() called.");
    auto data_manager{ GetDataManagerPtr() };
    data_manager->SetDatabaseManager(m_options.database_path);
    try
    {
        data_manager->ProcessFile(m_options.map_file_path, "map");
    }
    catch (const std::exception & e)
    {
        Logger::Log(LogLevel::Error,
            "PositionEstimationCommand::Execute() : " + std::string(e.what()));
        return false;
    }

    auto map_object{ data_manager->GetTypedDataObject<MapObject>("map") };
    map_object->MapValueArrayNormalization();
    
    RunMapValueConvergence(map_object.get());
    OutputPointList();
    return true;
}

bool PositionEstimationCommand::ValidateOptions(void) const
{
    Logger::Log(LogLevel::Debug, "PositionEstimationCommand::ValidateOptions() called");

    return true;
}

void PositionEstimationCommand::RunMapValueConvergence(MapObject * map_object)
{
    Logger::Log(LogLevel::Debug, "PositionEstimationCommand::RunMapValueConvergence() called");
    map_object->Display();

    BuildVoxelList(map_object);

    std::vector<VoxelNode> query_point_list(m_selected_voxel_list);
    for (int t = 0; t < m_options.iteration_count; t++)
    {
        Logger::Log(LogLevel::Info, "Iteration: " + std::to_string(t + 1));
        UpdatePointList(query_point_list);
    }
    
    BuildUniquePointList(query_point_list);
}

void PositionEstimationCommand::BuildVoxelList(MapObject * map_object)
{
    Logger::Log(LogLevel::Debug, "PositionEstimationCommand::BuildVoxelList() called");
    ScopeTimer timer("PositionEstimationCommand::BuildVoxelList");
    m_selected_voxel_list.clear();
    auto array_size{ map_object->GetMapValueArraySize() };
    m_selected_voxel_list.reserve(array_size);

    auto map_value_max{ map_object->GetMapValueMax() };
    auto threshold{ map_value_max * m_options.threshold_ratio };

    auto process_voxel = [&](size_t index, std::vector<VoxelNode> & list)
    {
        auto position{ map_object->GetGridPosition(index) };
        float value{ map_object->GetMapValue(index) };
        if (value <= threshold) return;
        list.emplace_back(position, value);
    };

#ifdef USE_OPENMP
    #pragma omp parallel num_threads(m_options.thread_size)
    {
        std::vector<VoxelNode> thread_local_list;
        thread_local_list.reserve(array_size / static_cast<size_t>(m_options.thread_size) + 1);

        #pragma omp for schedule(static)
        for (size_t i = 0; i < array_size; i++)
        {
            process_voxel(i, thread_local_list);
        }

        #pragma omp critical
        {
            m_selected_voxel_list.insert(
                m_selected_voxel_list.end(),
                thread_local_list.begin(), thread_local_list.end());
        }
    }
    #else
    for (size_t i = 0; i < array_size; i++)
    {
        process_voxel(i, m_selected_voxel_list);
    }
#endif

    if (m_selected_voxel_list.size() < m_selected_voxel_list.capacity())
    {
        m_selected_voxel_list.shrink_to_fit();
    }

    m_kd_tree_root = KDTreeAlgorithm<VoxelNode>::BuildKDTree(m_selected_voxel_list);

    Logger::Log(LogLevel::Info,
        "Number of selected voxels to be estimated from map = "
        + std::to_string(m_selected_voxel_list.size()) + " / "
        + std::to_string(array_size) + " voxels."
    );
}

void PositionEstimationCommand::UpdatePointList(std::vector<VoxelNode> & query_point_list)
{
    Logger::Log(LogLevel::Debug, "PositionEstimationCommand::UpdatePointList() called");
    std::vector<std::array<float, 3>> updated_position_list;
    updated_position_list.resize(query_point_list.size(), {0.0f, 0.0f, 0.0f});
    auto update_point_position = [&](size_t index)
    {
        auto knn_list{
            KDTreeAlgorithm<VoxelNode>::KNearestNeighbors(
                m_kd_tree_root.get(), &query_point_list[index], m_options.knn_size)
        };

        float weight_sum{ 0.0f };
        std::array<float, 3> point_position_update{ 0.0f, 0.0f, 0.0f };
        for (size_t j = 0; j < m_options.knn_size; j++)
        {
            auto w{ std::exp(m_options.alpha * std::log(knn_list[j]->GetValue())) };
            auto & query_point_position{ knn_list[j]->GetPosition() };
            point_position_update[0] += w * query_point_position[0];
            point_position_update[1] += w * query_point_position[1];
            point_position_update[2] += w * query_point_position[2];
            weight_sum += w;
        }
        point_position_update[0] /= weight_sum;
        point_position_update[1] /= weight_sum;
        point_position_update[2] /= weight_sum;
        updated_position_list[index] = point_position_update;
        query_point_list[index].SetPosition(updated_position_list[index]);
    };

#ifdef USE_OPENMP
    #pragma omp parallel for num_threads(m_options.thread_size)
#endif
    for (size_t i = 0; i < query_point_list.size(); i++)
    {
        update_point_position(i);
    }
}

void PositionEstimationCommand::BuildUniquePointList(const std::vector<VoxelNode> & point_list)
{
    Logger::Log(LogLevel::Debug, "PositionEstimationCommand::BuildUniquePointList() called");
    ScopeTimer timer("PositionEstimationCommand::BuildUniquePointList");
    if (point_list.empty())
    {
        Logger::Log(LogLevel::Warning, "Point list is empty. Nothing to degenerate.");
        return;
    }
    m_position_list.clear();
    m_position_list.reserve(point_list.size());

    auto pack_key = [](int x, int y, int z) noexcept -> uint64_t
    {
        // Pack three signed 21-bit values into a single 64-bit integer.
        // Layout: [63..42] X, [41..21] Y, [20..0] Z.
        constexpr uint64_t mask{ (1ULL << 21) - 1ULL };
        constexpr uint64_t offset{ 1ULL << 20 }; // allow negative values
        uint64_t ux{
            static_cast<uint64_t>(static_cast<int64_t>(x) + static_cast<int64_t>(offset)) & mask
        };
        uint64_t uy{
            static_cast<uint64_t>(static_cast<int64_t>(y) + static_cast<int64_t>(offset)) & mask
        };
        uint64_t uz{
            static_cast<uint64_t>(static_cast<int64_t>(z) + static_cast<int64_t>(offset)) & mask
        };
        return (ux << 42) | (uy << 21) | uz;
    };

    auto point_size_origin{ point_list.size() };
    auto tolerance{ 1.0e-2f };
    auto inv_tolerance{ 1.0f / tolerance };
    std::unordered_set<uint64_t> seen;
    seen.reserve(point_list.size());
    for (const auto & point : point_list)
    {
        const auto & position{ point.GetPosition() };
        uint64_t key{
            pack_key(
                static_cast<int>(std::floor(position[0] * inv_tolerance)),
                static_cast<int>(std::floor(position[1] * inv_tolerance)),
                static_cast<int>(std::floor(position[2] * inv_tolerance)))
        };
        if (seen.emplace(key).second)
        {
            m_position_list.emplace_back(position);
        }
    }

    auto removed_count{ point_size_origin - m_position_list.size() };
    Logger::Log(LogLevel::Info,
        "Number of removed points = "
        + std::to_string(removed_count) +" / "+ std::to_string(point_size_origin) +
        ", remaining number of unique points = " + std::to_string(m_position_list.size()));
}

void PositionEstimationCommand::OutputPointList(void) const
{
    Logger::Log(LogLevel::Debug, "PositionEstimationCommand::OutputPointList() called");
    if (m_position_list.empty())
    {
        Logger::Log(LogLevel::Warning, "No points to output.");
        return;
    }
    Logger::Log(LogLevel::Info,
        "Outputting point position list: " + std::to_string(m_position_list.size()) + " points.");

    auto map_file_name{ FilePathHelper::GetFileName(m_options.map_file_path, false) };
    auto output_file{
        FilePathHelper::EnsureTrailingSlash(m_options.folder_path)
        + "point_list_" + map_file_name + ".cmm" };
    ChimeraXHelper::WriteCMMPoints(m_position_list, output_file, 0.05f);
    Logger::Log(LogLevel::Info, "Output file: " + output_file);
}
