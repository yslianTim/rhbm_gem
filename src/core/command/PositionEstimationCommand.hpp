#pragma once

#include <array>
#include <filesystem>
#include <memory>
#include <string>
#include <vector>

#include "detail/CommandBase.hpp"
#include <rhbm_gem/data/object/MapObject.hpp>

template <typename T> struct KDNode;

namespace rhbm_gem {

class PositionEstimationCommand : public CommandWithRequest<PositionEstimationRequest>
{
private:
    std::vector<VoxelNode> m_selected_voxel_list;
    std::vector<VoxelNode> m_query_point_list;
    std::vector<std::array<float, 3>> m_position_list;
    std::unique_ptr<::KDNode<VoxelNode>> m_kd_tree_root;
    std::shared_ptr<MapObject> m_map_object;

public:
    PositionEstimationCommand();
    ~PositionEstimationCommand() override;

private:
    void NormalizeRequest() override;
    void ResetRuntimeState() override;
    bool ExecuteImpl() override;
    bool BuildDataObject();
    bool BuildVoxelList();
    void RunMapValueConvergence();
    void UpdatePointPosition(size_t index, size_t knn_size);
    void RunUniquePointList(float tolerance);
    void OutputPointList() const;
};

} // namespace rhbm_gem
