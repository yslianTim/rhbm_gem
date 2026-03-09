#pragma once

#include <array>
#include <filesystem>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include <rhbm_gem/core/command/CommandBase.hpp>
#include <rhbm_gem/core/command/OptionEnumClass.hpp>

namespace CLI
{
    class App;
}

template <typename T> struct KDNode;

namespace rhbm_gem {

class ModelObject;
class MapObject;
class AtomObject;

struct MapSimulationCommandOptions : public CommandOptions
{
    std::filesystem::path model_file_path;
    std::string map_file_name{"sim_map"};
    PotentialModel potential_model_choice{ PotentialModel::FIVE_GAUS_CHARGE };
    PartialCharge partial_charge_choice{ PartialCharge::PARTIAL };
    double cutoff_distance{ 5.0 };
    double grid_spacing{ 0.5 };
    std::vector<double> blurring_width_list{};
};

class MapSimulationCommand
    : public CommandWithProfileOptions<
          MapSimulationCommandOptions,
          CommandId::MapSimulation>
{
public:
    using Options = MapSimulationCommandOptions;

private:
    std::vector<AtomObject *> m_selected_atom_list;
    std::unordered_map<int, double> m_atom_charge_map;
    std::shared_ptr<ModelObject> m_model_object;
    std::array<float, 3> m_atom_range_minimum, m_atom_range_maximum;

public:
    MapSimulationCommand();
    ~MapSimulationCommand() override = default;
    void SetPotentialModelChoice(PotentialModel value);
    void SetPartialChargeChoice(PartialCharge value);
    void SetCutoffDistance(double value);
    void SetModelFilePath(const std::filesystem::path & value);
    void SetMapFileName(const std::string & value);
    void SetGridSpacing(double value);
    void SetBlurringWidthList(const std::string & value);

private:
    void RegisterCLIOptionsExtend(::CLI::App * cmd) override;
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
