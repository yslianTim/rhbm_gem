#pragma once

#include <array>
#include <filesystem>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include "detail/CommandBase.hpp"

template <typename T> struct KDNode;

namespace rhbm_gem {

class ModelObject;
class MapObject;
class AtomObject;

class MapSimulationCommand : public CommandWithRequest<MapSimulationRequest>
{
private:
    std::vector<AtomObject *> m_selected_atom_list;
    std::unordered_map<int, double> m_atom_charge_map;
    std::shared_ptr<ModelObject> m_model_object;
    std::array<float, 3> m_atom_range_minimum, m_atom_range_maximum;

public:
    MapSimulationCommand();
    ~MapSimulationCommand() override = default;

private:
    void NormalizeRequest() override;
    void ValidateOptions() override;
    void ResetRuntimeState() override;
    bool ExecuteImpl() override;
    bool BuildDataObject();
    void RunMapSimulation();
    void BuildAtomList(ModelObject * model_object);
    std::unique_ptr<MapObject> CreateMapObject();
    void PopulateMapValueArray(MapObject * map_object, double blurring_width);
    std::array<int, 3> CalculateGridSize(const std::array<float, 3> & grid_spacing) const;
    std::array<float, 3> CalculateOrigin(const std::array<float, 3> & grid_spacing) const;
};

} // namespace rhbm_gem
