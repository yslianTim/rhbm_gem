#pragma once

#include <memory>
#include <string>
#include <vector>
#include <array>
#include <unordered_map>
#include <filesystem>
#include <CLI/CLI.hpp>

#include "CommandBase.hpp"
#include "OptionEnumClass.hpp"

class ModelObject;
class MapObject;
class AtomObject;
template <typename T> struct KDNode;

class MapSimulationCommand : public CommandBase
{
public:
    struct Options : public CommandOptions
    {
        std::filesystem::path model_file_path;
        std::string map_file_name{"sim_map"};
        PotentialModel potential_model_choice{ PotentialModel::FIVE_GAUS_CHARGE };
        PartialCharge partial_charge_choice{ PartialCharge::PARTIAL };
        double cutoff_distance{ 5.0 };
        double grid_spacing{ 0.5 };
        std::vector<double> blurring_width_list{ 0.50 };
    };

private:
    Options m_options;
    std::vector<AtomObject *> m_selected_atom_list;
    std::unordered_map<int, double> m_atom_charge_map;
    std::unique_ptr<KDNode<AtomObject>> m_kd_tree_root;
    std::shared_ptr<ModelObject> m_model_object;
    std::array<float, 3> m_atom_range_minimum, m_atom_range_maximum;

public:
    MapSimulationCommand(void);
    ~MapSimulationCommand();
    bool Execute(void) override;
    bool ValidateOptions(void) const override;
    void RegisterCLIOptionsExtend(CLI::App * cmd) override;
    const CommandOptions & GetOptions(void) const override { return m_options; }
    CommandOptions & GetOptions(void) override { return m_options; }

    void SetThreadSize(int value) { m_options.thread_size = value; }
    void SetPotentialModelChoice(PotentialModel value) { m_options.potential_model_choice = value; }
    void SetPartialChargeChoice(PartialCharge value) { m_options.partial_charge_choice = value; }
    void SetCutoffDistance(double value) { m_options.cutoff_distance = value; }
    void SetModelFilePath(const std::filesystem::path & value) { m_options.model_file_path = value; }
    void SetFolderPath(const std::filesystem::path & path) { m_options.SetFolderPath(path); }
    void SetMapFileName(const std::string & value) { m_options.map_file_name = value; }
    void SetGridSpacing(double value) { m_options.grid_spacing = value; }
    void SetBlurringWidthList(const std::string & value);

private:
    bool BuildDataObject(void);
    void RunMapSimulation(void);
    void BuildAtomList(ModelObject * model_object);
    double CalculateAtomCharge(AtomObject * atom) const;
    void CalculateAtomRange(void);
    std::unique_ptr<MapObject> CreateSimulatedMapObject(double blurring_width);
    std::array<int, 3> CalculateGridSize(const std::array<float, 3> & grid_spacing) const;
    std::array<float, 3> CalculateOrigin(const std::array<float, 3> & grid_spacing) const;

};
