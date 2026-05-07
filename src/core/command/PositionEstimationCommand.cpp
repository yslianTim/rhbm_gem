#include "detail/CommandExecutor.hpp"
#include <rhbm_gem/data/io/ModelMapFileIO.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/utils/math/KDTreeAlgorithm.hpp>
#include <rhbm_gem/utils/domain/ScopeTimer.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/domain/ChimeraXHelper.hpp>
#include <rhbm_gem/utils/domain/FilePathHelper.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <unordered_set>
#include <vector>

#ifdef USE_OPENMP
#include <omp.h>
#endif

namespace rhbm_gem {

class PositionEstimationCommand final : public CommandBase<PositionEstimationRequest>
{
public:
    PositionEstimationCommand();

private:
    void NormalizeAndValidateRequest() override;
    bool ExecuteImpl() override;
};

} // namespace rhbm_gem

namespace {

constexpr std::string_view kMapKey{ "map" };

struct VoxelNode
{
    std::array<float, 3> position;
    float value;

    VoxelNode(const std::array<float, 3> & position_in, float value_in) :
        position{ position_in },
        value{ value_in }
    {
    }

    const std::array<float, 3> & GetPosition() const { return position; }
    float GetValue() const { return value; }
    void SetPosition(const std::array<float, 3> & updated_position)
    {
        position = updated_position;
    }
};

struct PositionEstimationState
{
    std::vector<VoxelNode> selected_voxel_list;
    std::vector<VoxelNode> query_point_list;
    std::vector<std::array<float, 3>> position_list;
    std::unique_ptr<::KDNode<VoxelNode>> kd_tree_root;
};

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

std::filesystem::path BuildPositionEstimationOutputPath(
    const std::filesystem::path & output_dir,
    std::string_view stem,
    std::string_view extension)
{
    std::string normalized_extension{ extension };
    if (!normalized_extension.empty() && normalized_extension.front() != '.')
    {
        normalized_extension.insert(normalized_extension.begin(), '.');
    }
    return output_dir / (std::string(stem) + normalized_extension);
}

std::optional<std::unique_ptr<rhbm_gem::MapObject>> LoadPositionEstimationMap(
    const rhbm_gem::PositionEstimationRequest & request)
{
    ScopeTimer timer("PositionEstimationCommand::BuildDataObject");
    try
    {
        auto map_object{ rhbm_gem::ReadMap(request.map_file_path) };
        map_object->SetKeyTag(std::string(kMapKey));
        map_object->MapValueArrayNormalization();
        return std::optional<std::unique_ptr<rhbm_gem::MapObject>>{ std::move(map_object) };
    }
    catch (const std::exception & e)
    {
        Logger::Log(LogLevel::Error,
            "PositionEstimationCommand::BuildDataObject : Failed to load map file '"
                + request.map_file_path.string() + "' as MapObject: " + std::string(e.what()));
        return std::nullopt;
    }
}

bool BuildVoxelList(
    const rhbm_gem::PositionEstimationRequest & request,
    rhbm_gem::MapObject & map_object,
    PositionEstimationState & state,
    int thread_size)
{
    ScopeTimer timer("PositionEstimationCommand::BuildVoxelList");
    auto array_size{ map_object.GetMapValueArraySize() };
    auto threshold{ map_object.GetMapValueMax() * static_cast<float>(request.threshold_ratio) };

    const float * map_values{ map_object.GetMapValueArray() };
    auto selected_count{
        static_cast<size_t>(std::count_if(map_values, map_values + array_size,
            [threshold](float v) { return v > threshold; }))
    };
    state.selected_voxel_list.reserve(selected_count);

    auto process_voxel = [&](size_t index, std::vector<VoxelNode> & list)
    {
        auto position{ map_object.GetGridPosition(index) };
        float value{ map_object.GetMapValue(index) };
        if (value <= threshold) return;
        list.emplace_back(position, value);
    };

#ifdef USE_OPENMP
    #pragma omp parallel num_threads(thread_size)
    {
        std::vector<VoxelNode> thread_local_list;
        thread_local_list.reserve(selected_count / static_cast<size_t>(thread_size) + 1);

        #pragma omp for schedule(static)
        for (size_t i = 0; i < array_size; i++)
        {
            process_voxel(i, thread_local_list);
        }

        #pragma omp critical
        {
            state.selected_voxel_list.insert(
                state.selected_voxel_list.end(),
                thread_local_list.begin(), thread_local_list.end());
        }
    }
#else
    for (size_t i = 0; i < array_size; i++)
    {
        process_voxel(i, state.selected_voxel_list);
    }
#endif

    if (state.selected_voxel_list.empty())
    {
        Logger::Log(LogLevel::Warning, "No voxels were selected from the map.");
        return false;
    }

    Logger::Log(LogLevel::Info,
        " /- Building KD-Tree from "+ std::to_string(state.selected_voxel_list.size()) + " voxels...");
    state.kd_tree_root = KDTreeAlgorithm<VoxelNode>::BuildKDTree(
        state.selected_voxel_list, 0, thread_size, true);

    Logger::Log(LogLevel::Info,
        "Number of selected voxels to be estimated from map = "
        + std::to_string(state.selected_voxel_list.size()) + " / "
        + std::to_string(array_size) + " voxels."
    );
    return true;
}

void UpdatePointPosition(
    PositionEstimationState & state,
    size_t index,
    size_t knn_size,
    float alpha)
{
    static thread_local std::vector<VoxelNode *> knn_list;
    if (knn_list.capacity() < knn_size)
    {
        knn_list.reserve(knn_size);
    }
    KDTreeAlgorithm<VoxelNode>::KNearestNeighbors(
        state.kd_tree_root.get(), &state.query_point_list[index], knn_size, knn_list);
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
    state.query_point_list[index].SetPosition(point_position_update);
}

void RunMapValueConvergence(
    const rhbm_gem::PositionEstimationRequest & request,
    PositionEstimationState & state,
    int thread_size)
{
    ScopeTimer timer("PositionEstimationCommand::RunMapValueConvergence");

    auto knn_size{ request.knn_size };
    if (state.selected_voxel_list.size() < knn_size)
    {
        if (state.selected_voxel_list.empty())
        {
            Logger::Log(
                LogLevel::Error,
                "No voxels available for position estimation.");
            return;
        }
        Logger::Log(LogLevel::Warning,
            "Selected voxel count (" + std::to_string(state.selected_voxel_list.size()) +
            ") is less than knn_size (" + std::to_string(knn_size) +
            "). Adjusting knn_size accordingly.");
        knn_size = state.selected_voxel_list.size();
    }

    state.query_point_list = state.selected_voxel_list;
    auto iteration_size{ static_cast<std::size_t>(request.iteration_count) };
    const auto alpha{ static_cast<float>(request.alpha) };
    Logger::Log(LogLevel::Info, " /- Running map value convergence iteration...");
    for (size_t t = 1; t <= iteration_size; t++)
    {
        auto progress_message{
            " (" + std::to_string(t) + " / " + std::to_string(iteration_size) +")"
        };
        auto point_size{ state.query_point_list.size() };
        size_t point_count{ 0 };
#ifdef USE_OPENMP
        #pragma omp parallel for num_threads(thread_size)
#endif
        for (size_t i = 0; i < state.query_point_list.size(); i++)
        {
            UpdatePointPosition(state, i, knn_size, alpha);
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

void RunUniquePointList(PositionEstimationState & state, float tolerance)
{
    ScopeTimer timer("PositionEstimationCommand::RunUniquePointList");
    if (state.query_point_list.empty())
    {
        Logger::Log(LogLevel::Warning, "Point list is empty. Nothing to degenerate.");
        return;
    }
    state.position_list.clear();
    state.position_list.reserve(state.query_point_list.size());

    auto point_size_origin{ state.query_point_list.size() };
    auto inv_tolerance{ 1.0f / tolerance };
    std::unordered_set<std::array<int64_t, 3>, QuantizedPointHash> unique_points;
    unique_points.reserve(state.query_point_list.size());
    for (const auto & point : state.query_point_list)
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
            state.position_list.emplace_back(std::array<float, 3>{
                static_cast<float>(q[0]) * tolerance,
                static_cast<float>(q[1]) * tolerance,
                static_cast<float>(q[2]) * tolerance
            });
        }
    }

    auto removed_count{ point_size_origin - state.position_list.size() };
    Logger::Log(LogLevel::Info,
        "Number of removed points = "
        + std::to_string(removed_count) +" / "+ std::to_string(point_size_origin) +
        ", remaining number of unique points = " + std::to_string(state.position_list.size())
    );
}

void OutputPointList(
    const rhbm_gem::PositionEstimationRequest & request,
    const PositionEstimationState & state,
    const std::filesystem::path & output_dir)
{
    ScopeTimer timer("PositionEstimationCommand::OutputPointList");
    if (state.position_list.empty())
    {
        Logger::Log(LogLevel::Warning, "No points to output.");
        return;
    }
    Logger::Log(LogLevel::Info,
        "Outputting point position list: " + std::to_string(state.position_list.size()) + " points.");

    auto map_file_name{
        "point_list_" + rhbm_gem::path_helper::GetFileName(request.map_file_path, false) };
    auto output_file{ BuildPositionEstimationOutputPath(output_dir, map_file_name, ".cmm") };
    ChimeraXHelper::WriteCMMPoints(state.position_list, output_file, 0.05f);
    Logger::Log(LogLevel::Info, "Output file: " + output_file.string());
}

} // namespace

namespace rhbm_gem {

PositionEstimationCommand::PositionEstimationCommand() :
    CommandBase<PositionEstimationRequest>{}
{
}

void PositionEstimationCommand::NormalizeAndValidateRequest()
{
    auto & request{ MutableRequest() };
    ValidateRequiredPath(request.map_file_path, "--map", "Map file");
    CoercePositiveScalar(
        request.iteration_count,
        "--iter",
        15,
        LogLevel::Warning,
        "Iteration count");
    CoercePositiveScalar(
        request.knn_size,
        "--knn",
        static_cast<std::size_t>(20),
        LogLevel::Warning,
        "KNN size");
    CoerceFinitePositiveScalar(
        request.alpha,
        "--alpha",
        2.0,
        LogLevel::Warning,
        "Alpha");
    CoerceFiniteExclusiveInclusiveRangeScalar(
        request.threshold_ratio,
        "--threshold",
        0.0,
        1.0,
        0.01,
        LogLevel::Warning,
        "Threshold ratio");
    CoerceFinitePositiveScalar(
        request.dedup_tolerance,
        "--dedup-tolerance",
        1.0e-2,
        LogLevel::Warning,
        "Dedup tolerance");
}

bool PositionEstimationCommand::ExecuteImpl()
{
    const auto & request{ RequestOptions() };
    auto map_object{ LoadPositionEstimationMap(request) };
    if (!map_object.has_value()) return false;

    PositionEstimationState state;
    if (!BuildVoxelList(request, **map_object, state, ThreadSize())) return false;
    RunMapValueConvergence(request, state, ThreadSize());
    RunUniquePointList(state, static_cast<float>(request.dedup_tolerance));
    OutputPointList(request, state, OutputFolder());
    return true;
}

namespace command_internal {

CommandResult ExecutePositionEstimationCommand(const PositionEstimationRequest & request)
{
    return ExecuteCommandInstance<PositionEstimationCommand>(request);
}

} // namespace command_internal

} // namespace rhbm_gem
