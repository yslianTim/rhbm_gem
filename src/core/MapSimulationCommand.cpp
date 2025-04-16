#include "MapSimulationCommand.hpp"
#include "AtomSelector.hpp"
#include "DataObjectManager.hpp"
#include "MapSimulationVisitor.hpp"

#include <sstream>

MapSimulationCommand::MapSimulationCommand(void) :
    m_atom_selector{ std::make_unique<AtomSelector>() }
{

}

MapSimulationCommand::~MapSimulationCommand()
{

}

void MapSimulationCommand::Execute(void)
{
    std::cout << "MapSimulationCommand::Execute() called." << std::endl;
    std::cout << "Total number of blurring width sets to be simulated: "
              << m_blurring_width_list.size() << std::endl;

    auto data_manager{ std::make_unique<DataObjectManager>() };
    data_manager->ProcessFile(m_model_file_path, "model");

    m_atom_selector->Print();

    auto analyzer{ std::make_unique<MapSimulationVisitor>(m_atom_selector.get()) };
    analyzer->SetModelObjectKeyTag("model");
    analyzer->SetFolderPath(m_folder_path);
    analyzer->SetThreadSize(static_cast<unsigned int>(m_thread_size));
    analyzer->SetPotentialModelChoice(m_potential_model_choice);
    analyzer->SetPartialChargeChoice(m_partial_charge_choice);
    analyzer->SetCutoffDistance(m_cutoff_distance);
    analyzer->SetGridSpacing(m_grid_spacing);
    analyzer->SetBlurringWidthList(m_blurring_width_list);

    data_manager->Accept(analyzer.get());
}

void MapSimulationCommand::SetBlurringWidthList(const std::string & value)
{
    std::stringstream ss(value);
    std::string segment;
    while (std::getline(ss, segment, ','))
    {
        if (segment == "") continue;
        m_blurring_width_list.emplace_back(std::stod(segment));
    }
}

void MapSimulationCommand::SetPickChainID(const std::string & value)
{
    m_atom_selector->PickChainID(value);
}

void MapSimulationCommand::SetPickIndicator(const std::string & value)
{
    m_atom_selector->PickIndicator(value);
}

void MapSimulationCommand::SetPickResidueType(const std::string & value)
{
    m_atom_selector->PickResidueType(value);
}

void MapSimulationCommand::SetPickElementType(const std::string & value)
{
    m_atom_selector->PickElementType(value);
}

void MapSimulationCommand::SetPickRemotenessType(const std::string & value)
{
    m_atom_selector->PickRemotenessType(value);
}

void MapSimulationCommand::SetPickBranchType(const std::string & value)
{
    m_atom_selector->PickBranchType(value);
}

void MapSimulationCommand::SetVetoChainID(const std::string & value)
{
    m_atom_selector->VetoChainID(value);
}

void MapSimulationCommand::SetVetoIndicator(const std::string & value)
{
    m_atom_selector->VetoIndicator(value);
}

void MapSimulationCommand::SetVetoResidueType(const std::string & value)
{
    m_atom_selector->VetoResidueType(value);
}

void MapSimulationCommand::SetVetoElementType(const std::string & value)
{
    m_atom_selector->VetoElementType(value);
}

void MapSimulationCommand::SetVetoRemotenessType(const std::string & value)
{
    m_atom_selector->VetoRemotenessType(value);
}

void MapSimulationCommand::SetVetoBranchType(const std::string & value)
{
    m_atom_selector->VetoBranchType(value);
}