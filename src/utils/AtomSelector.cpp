#include <rhbm_gem/utils/AtomSelector.hpp>
#include <rhbm_gem/utils/ChemicalDataHelper.hpp>
#include <rhbm_gem/utils/Logger.hpp>

#include <sstream>

AtomSelector::AtomSelector()
{

}

std::string AtomSelector::Describe() const
{
    std::ostringstream oss;
    oss<<"Atomic Picking List:\n";
    oss <<" - Chain set: ";
    for (auto & chain : pick_chain_set) oss << chain <<", ";
    oss <<"\n - Residue set: ";
    for (auto & residue : pick_residue_set) oss << ChemicalDataHelper::GetLabel(residue) <<", ";
    oss <<"\n - Element set: ";
    for (auto & element : pick_element_set) oss << ChemicalDataHelper::GetLabel(element) <<", ";
    oss <<"\n";

    oss <<"Atomic Vetoing List:\n";
    oss <<" - Chain set: ";
    for (auto & chain : veto_chain_set) oss << chain <<", ";
    oss <<"\n - Residue set: ";
    for (auto & residue : veto_residue_set) oss << ChemicalDataHelper::GetLabel(residue) <<", ";
    oss <<"\n - Element set: ";
    for (auto & element : veto_element_set) oss << ChemicalDataHelper::GetLabel(element) <<", ";
    oss <<"\n";

    return oss.str();
}

void AtomSelector::Print() const
{
    Logger::Log(LogLevel::Info, Describe());
}

bool AtomSelector::GetSelectionFlag(
    const std::string & chain_id,
    Residue residue,
    Element element) const
{
    auto selected{ true };
    if (!pick_chain_set.empty() && pick_chain_set.find(chain_id) == pick_chain_set.end()) selected = false;
    if (!pick_residue_set.empty() && pick_residue_set.find(residue) == pick_residue_set.end()) selected = false;
    if (!pick_element_set.empty() && pick_element_set.find(element) == pick_element_set.end()) selected = false;
    
    if (veto_chain_set.find(chain_id) != veto_chain_set.end()) selected = false;
    if (veto_residue_set.find(residue) != veto_residue_set.end()) selected = false;
    if (veto_element_set.find(element) != veto_element_set.end()) selected = false;
    
    return selected;
}

void AtomSelector::VetoChainID(const std::string & name)
{
    std::stringstream ss(name);
    std::string segment;
    while (std::getline(ss, segment, ',')) veto_chain_set.insert(segment);
}

void AtomSelector::VetoResidueType(const std::string & name)
{
    std::stringstream ss(name);
    std::string segment;
    while (std::getline(ss, segment, ','))
    {
        veto_residue_set.insert(ChemicalDataHelper::GetResidueFromString(segment));
    }
}

void AtomSelector::VetoElementType(const std::string & name)
{
    std::stringstream ss(name);
    std::string segment;
    while (std::getline(ss, segment, ','))
    {
        veto_element_set.insert(ChemicalDataHelper::GetElementFromString(segment));
    }
}

void AtomSelector::PickChainID(const std::string & name)
{
    std::stringstream ss(name);
    std::string segment;
    while (std::getline(ss, segment, ',')) pick_chain_set.insert(segment);
}

void AtomSelector::PickResidueType(const std::string & name)
{
    std::stringstream ss(name);
    std::string segment;
    while (std::getline(ss, segment, ','))
    {
        pick_residue_set.insert(ChemicalDataHelper::GetResidueFromString(segment));
    }
}

void AtomSelector::PickElementType(const std::string & name)
{
    std::stringstream ss(name);
    std::string segment;
    while (std::getline(ss, segment, ','))
    {
        pick_element_set.insert(ChemicalDataHelper::GetElementFromString(segment));
    }
}
