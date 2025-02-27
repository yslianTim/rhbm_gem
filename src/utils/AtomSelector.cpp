#include "AtomSelector.hpp"
#include "AtomicInfoHelper.hpp"

#include <sstream>

AtomSelector::AtomSelector(void)
{

}

AtomSelector::~AtomSelector()
{

}

bool AtomSelector::GetSelectionFlag(
    const std::string & chain_id, const std::string & indicator,
    int residue, int element, int remoteness, int branch)
{
    auto selected{ true };
    if (!pick_chain_set.empty() && pick_chain_set.find(chain_id) == pick_chain_set.end()) selected = false;
    if (!pick_indicator_set.empty() && pick_indicator_set.find(indicator) == pick_indicator_set.end()) selected = false;
    if (!pick_residue_set.empty() && pick_residue_set.find(residue) == pick_residue_set.end()) selected = false;
    if (!pick_element_set.empty() && pick_element_set.find(element) == pick_element_set.end()) selected = false;
    if (!pick_remoteness_set.empty() && pick_remoteness_set.find(remoteness) == pick_remoteness_set.end()) selected = false;
    if (!pick_branch_set.empty() && pick_branch_set.find(branch) == pick_branch_set.end()) selected = false;
    
    if (veto_chain_set.find(chain_id) != veto_chain_set.end()) selected = false;
    if (veto_indicator_set.find(indicator) != veto_indicator_set.end()) selected = false;
    if (veto_residue_set.find(residue) != veto_residue_set.end()) selected = false;
    if (veto_element_set.find(element) != veto_element_set.end()) selected = false;
    if (veto_remoteness_set.find(remoteness) != veto_remoteness_set.end()) selected = false;
    if (veto_branch_set.find(branch) != veto_branch_set.end()) selected = false;
    
    return selected;
}

void AtomSelector::VetoChainID(const std::string & name)
{
    std::stringstream ss(name);
    std::string segment;
    while (std::getline(ss, segment, ',')) veto_chain_set.insert(segment);
    is_up_to_date = false;
}

void AtomSelector::VetoIndicator(const std::string & name)
{
    std::stringstream ss(name);
    std::string segment;
    while (std::getline(ss, segment, ',')) veto_indicator_set.insert(segment);
    is_up_to_date = false;
}

void AtomSelector::VetoResidueType(const std::string & name)
{
    std::stringstream ss(name);
    std::string segment;
    while (std::getline(ss, segment, ','))
    {
        veto_residue_set.insert(AtomicInfoHelper::AtomInfoMapping<AtomicInfoHelper::ResidueTag>(segment));
    }
    is_up_to_date = false;
}

void AtomSelector::VetoElementType(const std::string & name)
{
    std::stringstream ss(name);
    std::string segment;
    while (std::getline(ss, segment, ','))
    {
        veto_element_set.insert(AtomicInfoHelper::AtomInfoMapping<AtomicInfoHelper::ElementTag>(segment));
    }
    is_up_to_date = false;
}

void AtomSelector::VetoRemotenessType(const std::string & name)
{
    std::stringstream ss(name);
    std::string segment;
    while (std::getline(ss, segment, ','))
    {
        veto_remoteness_set.insert(AtomicInfoHelper::AtomInfoMapping<AtomicInfoHelper::RemotenessTag>(segment));
    }
    is_up_to_date = false;
}

void AtomSelector::VetoBranchType(const std::string & name)
{
    std::stringstream ss(name);
    std::string segment;
    while (std::getline(ss, segment, ','))
    {
        veto_branch_set.insert(AtomicInfoHelper::AtomInfoMapping<AtomicInfoHelper::BranchTag>(segment));
    }
    is_up_to_date = false;
}

void AtomSelector::PickChainID(const std::string & name)
{
    std::stringstream ss(name);
    std::string segment;
    while (std::getline(ss, segment, ',')) pick_chain_set.insert(segment);
    is_up_to_date = false;
}

void AtomSelector::PickIndicator(const std::string & name)
{
    std::stringstream ss(name);
    std::string segment;
    while (std::getline(ss, segment, ',')) pick_indicator_set.insert(segment);
    is_up_to_date = false;
}

void AtomSelector::PickResidueType(const std::string & name)
{
    std::stringstream ss(name);
    std::string segment;
    while (std::getline(ss, segment, ','))
    {
        pick_residue_set.insert(AtomicInfoHelper::AtomInfoMapping<AtomicInfoHelper::ResidueTag>(segment));
    }
    is_up_to_date = false;
}

void AtomSelector::PickElementType(const std::string & name)
{
    std::stringstream ss(name);
    std::string segment;
    while (std::getline(ss, segment, ','))
    {
        pick_element_set.insert(AtomicInfoHelper::AtomInfoMapping<AtomicInfoHelper::ElementTag>(segment));
    }
    is_up_to_date = false;
}

void AtomSelector::PickRemotenessType(const std::string & name)
{
    std::stringstream ss(name);
    std::string segment;
    while (std::getline(ss, segment, ','))
    {
        pick_remoteness_set.insert(AtomicInfoHelper::AtomInfoMapping<AtomicInfoHelper::RemotenessTag>(segment));
    }
    is_up_to_date = false;
}

void AtomSelector::PickBranchType(const std::string & name)
{
    std::stringstream ss(name);
    std::string segment;
    while (std::getline(ss, segment, ','))
    {
        pick_branch_set.insert(AtomicInfoHelper::AtomInfoMapping<AtomicInfoHelper::BranchTag>(segment));
    }
    is_up_to_date = false;
}