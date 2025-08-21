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

    if (m_options.iteration_count <= 0)
    {
        Logger::Log(LogLevel::Error, "Iteration count must be positive");
        return false;
    }
    if (m_options.knn_size == 0)
    {
        Logger::Log(LogLevel::Error, "KNN size must be positive");
        return false;
    }
    if (m_options.alpha <= 0.0f)
    {
        Logger::Log(LogLevel::Error, "Alpha must be positive");
        return false;
    }
    if (m_options.threshold_ratio <= 0.0f || m_options.threshold_ratio > 1.0f)
    {
        Logger::Log(LogLevel::Error, "Threshold ratio must be in (0, 1]");
        return false;
    }
    if (m_options.thread_size <= 0)
    {
        Logger::Log(LogLevel::Error, "Thread size must be positive");
        return false;
    }
    if (!FilePathHelper::EnsureFileExists(m_options.map_file_path, "Map file"))
    {
        return false;
    }

    return true;
}

void PositionEstimationCommand::RunMapValueConvergence(const MapObject * map_object)
{
    Logger::Log(LogLevel::Debug, "PositionEstimationCommand::RunMapValueConvergence() called");
    ScopeTimer timer("PositionEstimationCommand::RunMapValueConvergence");

    if (!BuildVoxelList(map_object))
    {
        Logger::Log(LogLevel::Warning, "No voxels selected. Skip map value convergence.");
        return;
    }

    if (m_selected_voxel_list.size() < m_options.knn_size)
    {
        if (m_selected_voxel_list.empty())
        {
            Logger::Log(LogLevel::Error, "No voxels available for position estimation.");
            return;
        }
        Logger::Log(LogLevel::Warning,
            "Selected voxel count (" + std::to_string(m_selected_voxel_list.size()) +
            ") is less than knn_size (" + std::to_string(m_options.knn_size) +
            "). Adjusting knn_size accordingly.");
        m_options.knn_size = m_selected_voxel_list.size();
    }

    std::vector<VoxelNode> query_point_list(m_selected_voxel_list);
    for (int t = 0; t < m_options.iteration_count; t++)
    {
        Logger::ProgressBar(static_cast<std::size_t>(t + 1),
            static_cast<std::size_t>(m_options.iteration_count), 50);
        UpdatePointList(query_point_list);
    }
    
    BuildUniquePointList(query_point_list);
}

bool PositionEstimationCommand::BuildVoxelList(const MapObject * map_object)
{
    Logger::Log(LogLevel::Debug, "PositionEstimationCommand::BuildVoxelList() called");
    ScopeTimer timer("PositionEstimationCommand::BuildVoxelList");
    m_selected_voxel_list.clear();
    auto array_size{ map_object->GetMapValueArraySize() };
    m_selected_voxel_list.reserve(array_size);
    auto threshold{ map_object->GetMapValueMax() * m_options.threshold_ratio };

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

    if (m_selected_voxel_list.empty())
    {
        Logger::Log(LogLevel::Warning, "No voxels were selected from the map.");
        return false;
    }

    m_kd_tree_root = KDTreeAlgorithm<VoxelNode>::BuildKDTree(m_selected_voxel_list);

    Logger::Log(LogLevel::Info,
        "Number of selected voxels to be estimated from map = "
        + std::to_string(m_selected_voxel_list.size()) + " / "
        + std::to_string(array_size) + " voxels."
    );
    return true;
}

void PositionEstimationCommand::UpdatePointList(std::vector<VoxelNode> & query_point_list)
{
    Logger::Log(LogLevel::Debug, "PositionEstimationCommand::UpdatePointList() called");
    auto update_point_position = [&](size_t index)
    {
        auto knn_list{
            KDTreeAlgorithm<VoxelNode>::KNearestNeighbors(
                m_kd_tree_root.get(), &query_point_list[index], m_options.knn_size)
        };
        size_t knn_count{ knn_list.size() };
        if (knn_count < m_options.knn_size)
        {
            Logger::Log(LogLevel::Warning,
                "KNN search returned " + std::to_string(knn_count) +
                " neighbors, less than requested " +
                std::to_string(m_options.knn_size) + ".");
        }
        if (knn_count == 0)
        {
            return;
        }

        const auto alpha{ m_options.alpha };
        float weight_sum{ 0.0f };
        std::array<float, 3> point_position_update{ 0.0f, 0.0f, 0.0f };
        for (size_t j = 0; j < knn_count; j++)
        {
            //auto w{ std::exp(alpha * std::log(knn_list[j]->GetValue())) };
            if (knn_list[j]->GetValue() <= 0.0f) continue;
            auto w{ std::pow(knn_list[j]->GetValue(), alpha) };
            auto & query_point_position{ knn_list[j]->GetPosition() };
            point_position_update[0] += w * query_point_position[0];
            point_position_update[1] += w * query_point_position[1];
            point_position_update[2] += w * query_point_position[2];
            weight_sum += w;
        }
        if (weight_sum == 0.0f)
        {
            Logger::Log(LogLevel::Warning,
                "Weight sum is non-positive for point index "
                + std::to_string(index) + ". Skipping update to avoid division by zero.");
            return;
        }
        point_position_update[0] /= weight_sum;
        point_position_update[1] /= weight_sum;
        point_position_update[2] /= weight_sum;
        query_point_list[index].SetPosition(point_position_update);
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

    struct Int3
    {
        int x, y, z;
        bool operator==(const Int3 & other) const noexcept
        {
            return x == other.x && y == other.y && z == other.z;
        }
    };
    struct Int3Hash
    {
        std::size_t operator()(const Int3 & p) const noexcept
        {
            std::size_t h1{ std::hash<int>{}(p.x) };
            std::size_t h2{ std::hash<int>{}(p.y) };
            std::size_t h3{ std::hash<int>{}(p.z) };
            return h1 ^ (h2 << 1) ^ (h3 << 2);
        }
    };

    auto point_size_origin{ point_list.size() };
    auto tolerance{ 1.0e-2f };
    auto inv_tolerance{ 1.0f / tolerance };
    std::unordered_set<Int3, Int3Hash> seen;
    seen.reserve(point_list.size());
    for (const auto & point : point_list)
    {
        const auto & position{ point.GetPosition() };
        Int3 key{
            static_cast<int>(std::floor(position[0] * inv_tolerance)),
            static_cast<int>(std::floor(position[1] * inv_tolerance)),
            static_cast<int>(std::floor(position[2] * inv_tolerance))
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

    auto map_file_name{
        "point_list_" + FilePathHelper::GetFileName(m_options.map_file_path, false) + ".cmm" };
    auto output_file{ m_options.folder_path / map_file_name };
    ChimeraXHelper::WriteCMMPoints(m_position_list, output_file, 0.05f);
    Logger::Log(LogLevel::Info, "Output file: " + output_file.string());
}
