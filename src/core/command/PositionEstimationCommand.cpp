#include "internal/command/PositionEstimationCommand.hpp"
#include <rhbm_gem/core/command/CommandApi.hpp>
#include "internal/command/CommandDataSupport.hpp"
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/data/io/DataObjectManager.hpp>
#include <rhbm_gem/utils/math/KDTreeAlgorithm.hpp>
#include <rhbm_gem/utils/domain/ScopeTimer.hpp>
#include <rhbm_gem/utils/math/ArrayStats.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/domain/ChimeraXHelper.hpp>
#include <rhbm_gem/utils/domain/FilePathHelper.hpp>

#include <random>
#include <memory>
#include <vector>
#include <array>
#include <algorithm>
#include <unordered_set>
#include <cmath>
#include <cstdint>

#ifdef USE_OPENMP
#include <omp.h>
#endif

namespace {
constexpr std::string_view kMapKey{ "map" };
constexpr std::string_view kMapOption{ "--map" };
constexpr std::string_view kIterationOption{ "--iter" };
constexpr std::string_view kKnnOption{ "--knn" };
constexpr std::string_view kAlphaOption{ "--alpha" };
constexpr std::string_view kThresholdOption{ "--threshold" };
constexpr std::string_view kDedupToleranceOption{ "--dedup-tolerance" };

struct QuantizedPointHash
{
    std::size_t operator()(const std::array<int64_t, 3> & a) const noexcept
    {
        std::size_t h1{ std::hash<int64_t>{}(a[0]) };
        std::size_t h2{ std::hash<int64_t>{}(a[1]) };
        std::size_t h3{ std::hash<int64_t>{}(a[2]) };
        return h1 ^ (h2 << 1) ^ (h3 << 2);
    }
};

} // namespace

namespace rhbm_gem {

PositionEstimationCommand::PositionEstimationCommand(CommonOptionProfile profile) :
    CommandWithOptions<PositionEstimationCommandOptions>{
        CommonOptionMaskForProfile(profile) },
    m_selected_voxel_list{}, m_query_point_list{}, m_position_list{},
    m_kd_tree_root{ nullptr }, m_map_object{ nullptr }
{
}

PositionEstimationCommand::~PositionEstimationCommand() = default;

void PositionEstimationCommand::ApplyRequest(const PositionEstimationRequest & request)
{
    ApplyCommonRequest(request.common);
    SetMapFilePath(request.map_file_path);
    SetIterationCount(request.iteration_count);
    SetKNNSize(request.knn_size);
    SetAlpha(request.alpha);
    SetThresholdRatio(request.threshold_ratio);
    SetDedupTolerance(request.dedup_tolerance);
}

bool PositionEstimationCommand::ExecuteImpl()
{
    if (BuildDataObject() == false) return false;
    if (BuildVoxelList() == false) return false;
    RunMapValueConvergence();
    RunUniquePointList(m_options.dedup_tolerance);
    OutputPointList();
    return true;
}

void PositionEstimationCommand::ResetRuntimeState()
{
    m_selected_voxel_list.clear();
    m_query_point_list.clear();
    m_position_list.clear();
    m_kd_tree_root.reset();
    m_map_object.reset();
}

void PositionEstimationCommand::SetMapFilePath(const std::filesystem::path & path)
{
    SetRequiredExistingPathOption(m_options.map_file_path, path, kMapOption, "Map file");
}

void PositionEstimationCommand::SetIterationCount(int value)
{
    SetNormalizedScalarOption(
        m_options.iteration_count,
        value,
        kIterationOption,
        [](int candidate) { return candidate > 0; },
        15,
        "Iteration count must be positive, reset to default 15");
}

void PositionEstimationCommand::SetKNNSize(int value)
{
    SetNormalizedScalarOption(
        m_options.knn_size,
        value,
        kKnnOption,
        [](int candidate) { return candidate > 0; },
        static_cast<size_t>(20),
        "KNN size must be positive, reset to default 20");
}

void PositionEstimationCommand::SetAlpha(double value)
{
    SetNormalizedScalarOption(
        m_options.alpha,
        value,
        kAlphaOption,
        [](double candidate) { return candidate > 0.0; },
        2.0f,
        "Alpha must be positive, reset to default 2.0");
}

void PositionEstimationCommand::SetThresholdRatio(double value)
{
    SetNormalizedScalarOption(
        m_options.threshold_ratio,
        value,
        kThresholdOption,
        [](double candidate) { return candidate > 0.0 && candidate <= 1.0; },
        0.01f,
        "Threshold ratio must be in (0, 1], reset to default 0.01");
}

void PositionEstimationCommand::SetDedupTolerance(double value)
{
    SetNormalizedScalarOption(
        m_options.dedup_tolerance,
        value,
        kDedupToleranceOption,
        [](double candidate) { return candidate > 0.0; },
        1.0e-2f,
        "Dedup tolerance must be positive, reset to default 0.01");
}

bool PositionEstimationCommand::BuildDataObject()
{
    ScopeTimer timer("PositionEstimationCommand::BuildDataObject");
    try
    {
        m_map_object = command_data_loader::ProcessMapFile(
            m_data_manager, m_options.map_file_path, kMapKey, "map file");
        NormalizeMapObject(*m_map_object);
    }
    catch (const std::exception & e)
    {
        Logger::Log(LogLevel::Error,
            "PositionEstimationCommand::BuildDataObject : " + std::string(e.what()));
        return false;
    }
    return true;
}

bool PositionEstimationCommand::BuildVoxelList()
{
    ScopeTimer timer("PositionEstimationCommand::BuildVoxelList");
    m_selected_voxel_list.clear();
    auto array_size{ m_map_object->GetMapValueArraySize() };
    auto threshold{ m_map_object->GetMapValueMax() * m_options.threshold_ratio };

    const float * map_values{ m_map_object->GetMapValueArray() };
    auto selected_count{
        static_cast<size_t>(std::count_if(map_values, map_values + array_size,
            [threshold](float v) { return v > threshold; }))
    };
    m_selected_voxel_list.reserve(selected_count);

    auto process_voxel = [&](size_t index, std::vector<VoxelNode> & list)
    {
        auto position{ m_map_object->GetGridPosition(index) };
        float value{ m_map_object->GetMapValue(index) };
        if (value <= threshold) return;
        list.emplace_back(position, value);
    };

#ifdef USE_OPENMP
    #pragma omp parallel num_threads(m_options.thread_size)
    {
        std::vector<VoxelNode> thread_local_list;
        thread_local_list.reserve(selected_count / static_cast<size_t>(m_options.thread_size) + 1);

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

    if (m_selected_voxel_list.empty())
    {
        Logger::Log(LogLevel::Warning, "No voxels were selected from the map.");
        return false;
    }

    Logger::Log(LogLevel::Info,
        " /- Building KD-Tree from "+ std::to_string(m_selected_voxel_list.size()) + " voxels...");
    m_kd_tree_root = KDTreeAlgorithm<VoxelNode>::BuildKDTree(
        m_selected_voxel_list, 0, m_options.thread_size, true);

    Logger::Log(LogLevel::Info,
        "Number of selected voxels to be estimated from map = "
        + std::to_string(m_selected_voxel_list.size()) + " / "
        + std::to_string(array_size) + " voxels."
    );
    return true;
}

void PositionEstimationCommand::RunMapValueConvergence()
{
    ScopeTimer timer("PositionEstimationCommand::RunMapValueConvergence");

    auto knn_size{ m_options.knn_size };
    if (m_selected_voxel_list.size() < knn_size)
    {
        if (m_selected_voxel_list.empty())
        {
            Logger::Log(LogLevel::Error, "No voxels available for position estimation.");
            return;
        }
        Logger::Log(LogLevel::Warning,
            "Selected voxel count (" + std::to_string(m_selected_voxel_list.size()) +
            ") is less than knn_size (" + std::to_string(knn_size) +
            "). Adjusting knn_size accordingly.");
        knn_size = m_selected_voxel_list.size();
    }

    m_query_point_list = m_selected_voxel_list;
    auto iteration_size{ static_cast<std::size_t>(m_options.iteration_count) };
    Logger::Log(LogLevel::Info, " /- Running map value convergence iteration...");
    for (size_t t = 1; t <= iteration_size; t++)
    {
        auto progress_message{
            " (" + std::to_string(t) + " / " + std::to_string(iteration_size) +")"
        };
        auto point_size{ m_query_point_list.size() };
        size_t point_count{ 0 };
#ifdef USE_OPENMP
        #pragma omp parallel for num_threads(m_options.thread_size)
#endif
        for (size_t i = 0; i < m_query_point_list.size(); i++)
        {
            UpdatePointPosition(i, knn_size);
#ifdef USE_OPENMP
            #pragma omp critical
#endif
            {
                point_count++;
                Logger::ProgressPercent(point_count, point_size, 50, progress_message);
            }
        }
    }
}

void PositionEstimationCommand::UpdatePointPosition(size_t index, size_t knn_size)
{
    static thread_local std::vector<VoxelNode *> knn_list;
    if (knn_list.capacity() < knn_size)
    {
        knn_list.reserve(knn_size);
    }
    KDTreeAlgorithm<VoxelNode>::KNearestNeighbors(
        m_kd_tree_root.get(), &m_query_point_list[index], knn_size, knn_list);
    size_t knn_count{ knn_list.size() };
    if (knn_count < knn_size)
    {
        Logger::Log(LogLevel::Warning,
            "KNN search returned " + std::to_string(knn_count) +
            " neighbors, less than requested " + std::to_string(knn_size) + ".");
    }
    if (knn_count == 0)
    {
        Logger::Log(LogLevel::Warning,
            "KNN search returned zero neighbors for point index ["
            + std::to_string(index) + "]. Skipping update.");
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
    m_query_point_list[index].SetPosition(point_position_update);
}

void PositionEstimationCommand::RunUniquePointList(float tolerance)
{
    ScopeTimer timer("PositionEstimationCommand::RunUniquePointList");
    if (m_query_point_list.empty())
    {
        Logger::Log(LogLevel::Warning, "Point list is empty. Nothing to degenerate.");
        return;
    }
    m_position_list.clear();
    m_position_list.reserve(m_query_point_list.size());

    auto point_size_origin{ m_query_point_list.size() };
    auto inv_tolerance{ 1.0f / tolerance };
    std::unordered_set<std::array<int64_t, 3>, QuantizedPointHash> unique_points;
    unique_points.reserve(m_query_point_list.size());
    for (const auto & point : m_query_point_list)
    {
        const auto & position{ point.GetPosition() };
        std::array<int64_t, 3> q{
            static_cast<int64_t>(std::floor(position[0] * inv_tolerance)),
            static_cast<int64_t>(std::floor(position[1] * inv_tolerance)),
            static_cast<int64_t>(std::floor(position[2] * inv_tolerance))
        };
        auto result{ unique_points.insert(q) };
        if (result.second)
        {
            m_position_list.emplace_back(std::array<float, 3>{
                static_cast<float>(q[0]) * tolerance,
                static_cast<float>(q[1]) * tolerance,
                static_cast<float>(q[2]) * tolerance
            });
        }
    }

    auto removed_count{ point_size_origin - m_position_list.size() };
    Logger::Log(LogLevel::Info,
        "Number of removed points = "
        + std::to_string(removed_count) +" / "+ std::to_string(point_size_origin) +
        ", remaining number of unique points = " + std::to_string(m_position_list.size())
    );
    
    m_selected_voxel_list.clear();
    m_selected_voxel_list.shrink_to_fit();
    m_query_point_list.clear();
    m_query_point_list.shrink_to_fit();
}

void PositionEstimationCommand::OutputPointList() const
{
    ScopeTimer timer("PositionEstimationCommand::OutputPointList");
    if (m_position_list.empty())
    {
        Logger::Log(LogLevel::Warning, "No points to output.");
        return;
    }
    Logger::Log(LogLevel::Info,
        "Outputting point position list: " + std::to_string(m_position_list.size()) + " points.");

    auto map_file_name{
        "point_list_" + FilePathHelper::GetFileName(m_options.map_file_path, false) };
    auto output_file{ BuildOutputPath(map_file_name, ".cmm") };
    ChimeraXHelper::WriteCMMPoints(m_position_list, output_file, 0.05f);
    Logger::Log(LogLevel::Info, "Output file: " + output_file.string());
}

} // namespace rhbm_gem
