#pragma once

#include <memory>
#include <string>
#include <vector>
#include <array>
#include <filesystem>
#include <CLI/CLI.hpp>

#include "CommandBase.hpp"
#include "MapObject.hpp"

template <typename T> struct KDNode;

class PositionEstimationCommand : public CommandBase
{
public:
    struct Options : public CommandOptions
    {
        int iteration_count{ 15 };
        size_t knn_size{ 20 };
        float alpha{ 2.0 };
        float threshold_ratio{ 0.01f };
        std::filesystem::path map_file_path;
    };

private:
    Options m_options;
    std::vector<VoxelNode> m_selected_voxel_list;
    std::vector<VoxelNode> m_query_point_list;
    std::vector<std::array<float, 3>> m_position_list;
    std::unique_ptr<KDNode<VoxelNode>> m_kd_tree_root;
    std::shared_ptr<MapObject> m_map_object;

public:
    PositionEstimationCommand(void);
    ~PositionEstimationCommand() = default;
    bool Execute(void) override;
    bool ValidateOptions(void) const override;
    void RegisterCLIOptionsExtend(CLI::App * cmd) override;
    const CommandOptions & GetOptions(void) const override { return m_options; }
    CommandOptions & GetOptions(void) override { return m_options; }

    void SetMapFilePath(const std::filesystem::path & path) { m_options.map_file_path = path; }
    void SetIterationCount(int value) { m_options.iteration_count = value; }
    void SetKNNSize(int value) { m_options.knn_size = static_cast<size_t>(value); }
    void SetAlpha(double value) { m_options.alpha = static_cast<float>(value); }
    void SetThresholdRatio(double value) { m_options.threshold_ratio = static_cast<float>(value); }

private:
    bool BuildDataObject(void);
    bool BuildVoxelList(void);
    void RunMapValueConvergence(void);
    void UpdatePointList(void);
    void RunUniquePointList(void);
    void OutputPointList(void) const;

};
