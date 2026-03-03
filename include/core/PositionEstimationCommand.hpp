#pragma once

#include <array>
#include <filesystem>
#include <memory>
#include <string>
#include <vector>

#include <CLI/CLI.hpp>

#include "CommandBase.hpp"
#include "MapObject.hpp"

template <typename T> struct KDNode;

namespace rhbm_gem {

class PositionEstimationCommand : public CommandBase
{
public:
    struct Options : public CommandOptions
    {
        int iteration_count{ 15 };
        size_t knn_size{ 20 };
        float alpha{ 2.0 };
        float threshold_ratio{ 0.01f };
        float dedup_tolerance{ 1.0e-2f };
        std::filesystem::path map_file_path;
    };

private:
    Options m_options;
    std::vector<VoxelNode> m_selected_voxel_list;
    std::vector<VoxelNode> m_query_point_list;
    std::vector<std::array<float, 3>> m_position_list;
    std::unique_ptr<::KDNode<VoxelNode>> m_kd_tree_root;
    std::shared_ptr<MapObject> m_map_object;

public:
    PositionEstimationCommand();
    ~PositionEstimationCommand() override;
    void RegisterCLIOptionsExtend(::CLI::App * cmd) override;
    CommandId GetCommandId() const override { return CommandId::PositionEstimation; }
    void ResetRuntimeState() override;
    const CommandOptions & GetOptions() const override { return m_options; }
    CommandOptions & GetOptions() override { return m_options; }

    void SetMapFilePath(const std::filesystem::path & path);
    void SetIterationCount(int value);
    void SetKNNSize(int value);
    void SetAlpha(double value);
    void SetThresholdRatio(double value);
    void SetDedupTolerance(double value);

private:
    bool ExecuteImpl() override;
    bool BuildDataObject();
    bool BuildVoxelList();
    void RunMapValueConvergence();
    void UpdatePointPosition(size_t index, size_t knn_size);
    void RunUniquePointList(float tolerance);
    void OutputPointList() const;
};

} // namespace rhbm_gem
