#include "PositionEstimationCommand.hpp"
#include "MapObject.hpp"
#include "DataObjectManager.hpp"
#include "KDTreeAlgorithm.hpp"
#include "ScopeTimer.hpp"
#include "ArrayStats.hpp"
#include "Logger.hpp"
#include "CommandRegistry.hpp"

#include <random>
#include <memory>
#include <vector>

#ifdef USE_OPENMP
#include <omp.h>
#endif

namespace {
CommandRegistrar<PositionEstimationCommand> registrar_model_test{
    "position_estimation",
    "Run atom position estimation"};
}

PositionEstimationCommand::PositionEstimationCommand(void) :
    CommandBase(), m_options{}, m_selected_voxel_list{}
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
    cmd->add_option("-k,--knn", m_options.knn_size,
        "KNN size for estimation")->default_val(m_options.knn_size);
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

    //data_manager->SaveDataObject("map", m_options.saved_key_tag);
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

    std::vector<VoxelNode> voxel_list{ m_selected_voxel_list };
    for (int t = 0; t < m_options.iteration_count; t++)
    {
        Logger::Log(LogLevel::Info, "Iteration: " + std::to_string(t + 1));
        UpdateVoxelPosition(voxel_list);
    }
    DisplayVoxelPositionChange(voxel_list);
    //DegenerateVoxelList(voxel_list);
    //DisplayVoxelList(voxel_list);
}

void PositionEstimationCommand::BuildVoxelList(MapObject * map_object)
{
    Logger::Log(LogLevel::Debug, "PositionEstimationCommand::BuildVoxelList() called");
    ScopeTimer timer("PositionEstimationCommand::BuildVoxelList");
    m_selected_voxel_list.clear();
    auto array_size{ map_object->GetMapValueArraySize() };
    m_selected_voxel_list.reserve(array_size);

    auto map_value_max{ map_object->GetMapValueMax() };
    auto threahold{ map_value_max * 0.01f };

#ifdef USE_OPENMP
#pragma omp parallel
    {
        std::vector<VoxelNode> thread_local_list;
        thread_local_list.reserve(array_size / static_cast<size_t>(omp_get_num_threads()) + 1);

#pragma omp for schedule(static)
        for (size_t i = 0; i < array_size; i++)
        {
            auto position{ map_object->GetGridPosition(i) };
            float value{ map_object->GetMapValue(i) };
            if (value <= threahold) continue; // Skip negative values & keep first 99% of map values
            thread_local_list.emplace_back(position, value);
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
        auto position{ map_object->GetGridPosition(i) };
        float value{ map_object->GetMapValue(i) };
        if (value <= threahold) continue; // Skip negative values & keep first 99% of map values
        m_selected_voxel_list.emplace_back(position, value);
    }
#endif

    if (m_selected_voxel_list.size() < m_selected_voxel_list.capacity())
    {
        m_selected_voxel_list.shrink_to_fit();
    }

    Logger::Log(LogLevel::Info,
        "Number of selected voxels to be estimated from map = "
        + std::to_string(m_selected_voxel_list.size()) + " / "
        + std::to_string(array_size) + " voxels."
    );
}

void PositionEstimationCommand::UpdateVoxelPosition(std::vector<VoxelNode> & voxel_list)
{
    Logger::Log(LogLevel::Debug, "PositionEstimationCommand::UpdateVoxelPosition() called");
    ScopeTimer timer("PositionEstimationCommand::UpdateVoxelPosition");
    auto kd_tree_root{ KDTreeAlgorithm<VoxelNode>::BuildKDTree(voxel_list, 0) };
    auto * kd_tree_root_ptr{ kd_tree_root.get() };

#ifdef USE_OPENMP
#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < voxel_list.size(); i++)
    {
        auto & voxel{ voxel_list[i] };
        auto voxel_position_origin{ m_selected_voxel_list[i].GetPosition() };
        std::array<float, 3> voxel_position_update{ 0.0f, 0.0f, 0.0f };
        float weight_sum{ 0.0f };
        auto knn_list{
            KDTreeAlgorithm<VoxelNode>::KNearestNeighbors(
                kd_tree_root_ptr, &voxel, static_cast<size_t>(m_options.knn_size))
        };

        for (auto neighbor_voxel : knn_list)
        {
            float K{
                std::exp(static_cast<float>(m_options.alpha) * std::log(neighbor_voxel->GetValue()))
            };
            voxel_position_update[0] += K * voxel_position_origin[0];
            voxel_position_update[1] += K * voxel_position_origin[1];
            voxel_position_update[2] += K * voxel_position_origin[2];
            weight_sum += K;
        }
        voxel_position_update[0] /= weight_sum;
        voxel_position_update[1] /= weight_sum;
        voxel_position_update[2] /= weight_sum;

        voxel.SetPosition(voxel_position_update);
    }
#else
    for (size_t i = 0; i < voxel_list.size(); i++)
    {
        auto & voxel{ voxel_list[i] };
        auto voxel_position_origin{ m_selected_voxel_list[i].GetPosition() };
        std::array<float, 3> voxel_position_update{ 0.0f, 0.0f, 0.0f };
        float weight_sum{ 0.0f };
        auto knn_list{
            KDTreeAlgorithm<VoxelNode>::KNearestNeighbors(
                kd_tree_root_ptr, &voxel, static_cast<size_t>(m_options.knn_size))
        };

        for (auto neighbor_voxel : knn_list)
        {
            float K{
                std::exp(static_cast<float>(m_options.alpha) * std::log(neighbor_voxel->GetValue()))
            };
            voxel_position_update[0] += K * voxel_position_origin[0];
            voxel_position_update[1] += K * voxel_position_origin[1];
            voxel_position_update[2] += K * voxel_position_origin[2];
            weight_sum += K;
        }
        voxel_position_update[0] /= weight_sum;
        voxel_position_update[1] /= weight_sum;
        voxel_position_update[2] /= weight_sum;

        voxel.SetPosition(voxel_position_update);
    }
#endif
}

void PositionEstimationCommand::DegenerateVoxelList(std::vector<VoxelNode> & voxel_list)
{
    Logger::Log(LogLevel::Debug, "PositionEstimationCommand::DegenerateVoxelList() called");
    ScopeTimer timer("PositionEstimationCommand::DegenerateVoxelList");
    if (voxel_list.empty())
    {
        Logger::Log(LogLevel::Warning, "Voxel list is empty. Nothing to degenerate.");
        return;
    }
    auto voxel_size_origin{ voxel_list.size() };
    auto tolerance{ 1.0e-2f };
    std::size_t removed_count{ 0 };
    for (std::size_t i = 0; i < voxel_list.size(); i++)
    {
        auto position_i{ voxel_list[i].GetPosition() };
        for (std::size_t j = i + 1; j < voxel_list.size(); j++)
        {
            auto position_j{ voxel_list[j].GetPosition() };
            if (ArrayStats<float>::ComputeNorm(position_i, position_j) < tolerance)
            {
                voxel_list.erase(voxel_list.begin() + j);
                removed_count++;
                j--; // Adjust index after removal
            }
        }
    }

    Logger::Log(LogLevel::Info,
        "Number of removed voxels = "
        + std::to_string(removed_count) +" / "+ std::to_string(voxel_size_origin) +
        ", remaining voxels = " + std::to_string(voxel_list.size()));
}

void PositionEstimationCommand::DisplayVoxelList(const std::vector<VoxelNode> & voxel_list) const
{
    Logger::Log(LogLevel::Debug, "PositionEstimationCommand::DisplayVoxelList() called");
    if (voxel_list.empty())
    {
        Logger::Log(LogLevel::Warning, "No voxels to display.");
        return;
    }
    Logger::Log(LogLevel::Info,
        "Displaying voxel position list: " + std::to_string(voxel_list.size()) + " voxels.");
    for (const auto & voxel : voxel_list)
    {
        const auto & position{ voxel.GetPosition() };
        Logger::Log(LogLevel::Info,
            "Voxel Position: ["
            + std::to_string(position[0]) + ", "
            + std::to_string(position[1]) + ", "
            + std::to_string(position[2]) + "]"
        );
    }
}

void PositionEstimationCommand::DisplayVoxelPositionChange(
    const std::vector<VoxelNode> & voxel_list) const
{
    Logger::Log(LogLevel::Debug, "PositionEstimationCommand::DisplayVoxelPositionChange() called");
    if (voxel_list.size() != m_selected_voxel_list.size())
    {
        Logger::Log(LogLevel::Error,
            "Voxel list size mismatch: "
            + std::to_string(voxel_list.size()) + " vs "
            + std::to_string(m_selected_voxel_list.size()));
        return;
    }

    for (size_t i = 0; i < voxel_list.size(); i++)
    {
        const auto & new_position{ voxel_list[i].GetPosition() };
        const auto & old_position{ m_selected_voxel_list[i].GetPosition() };
        Logger::Log(LogLevel::Info,
            "Voxel " + std::to_string(i) + " Position Change: ["
            + std::to_string(new_position[0] - old_position[0]) + ", "
            + std::to_string(new_position[1] - old_position[1]) + ", "
            + std::to_string(new_position[2] - old_position[2]) + "]"
        );
    }
}
