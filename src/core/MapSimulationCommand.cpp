#include "MapSimulationCommand.hpp"
#include "AtomSelector.hpp"

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