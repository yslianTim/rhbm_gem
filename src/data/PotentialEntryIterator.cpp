#include "PotentialEntryIterator.hpp"
#include "ModelObject.hpp"
#include "GroupPotentialEntry.hpp"

PotentialEntryIterator::PotentialEntryIterator(ModelObject * model_object) :
    m_model_object{ model_object }
{
    auto group_entry1{ m_model_object->GetGroupPotentialEntry(AtomicInfoHelper::GetElementClassKey()) };
    m_element_class_group_entry = dynamic_cast<GroupPotentialEntry<ElementKeyType> *>(group_entry1);
    auto group_entry2{ m_model_object->GetGroupPotentialEntry(AtomicInfoHelper::GetResidueClassKey()) };
    m_residue_class_group_entry = dynamic_cast<GroupPotentialEntry<ResidueKeyType> *>(group_entry2);
}

PotentialEntryIterator::~PotentialEntryIterator()
{

}

bool PotentialEntryIterator::IsAvailableGroupKey(ElementKeyType & group_key) const
{
    return CheckGroupKey(group_key, false);
}

bool PotentialEntryIterator::IsAvailableGroupKey(ResidueKeyType & group_key) const
{
    return CheckGroupKey(group_key, false);
}

std::tuple<double, double> PotentialEntryIterator::GetGausEstimatePrior(ElementKeyType & group_key) const
{
    if (CheckGroupKey(group_key) == false)
    {
        return std::make_tuple(0.0, 0.0);
    }
    return m_element_class_group_entry->GetGausEstimatePrior(&group_key);
}

std::tuple<double, double> PotentialEntryIterator::GetGausEstimatePrior(ResidueKeyType & group_key) const
{
    if (CheckGroupKey(group_key) == false)
    {
        return std::make_tuple(0.0, 0.0);
    }
    return m_residue_class_group_entry->GetGausEstimatePrior(&group_key);
}

std::tuple<double, double> PotentialEntryIterator::GetGausVariancePrior(ElementKeyType & group_key) const
{
    if (CheckGroupKey(group_key) == false)
    {
        return std::make_tuple(0.0, 0.0);
    }
    return m_element_class_group_entry->GetGausVariancePrior(&group_key);
}

std::tuple<double, double> PotentialEntryIterator::GetGausVariancePrior(ResidueKeyType & group_key) const
{
    if (CheckGroupKey(group_key) == false)
    {
        return std::make_tuple(0.0, 0.0);
    }
    return m_residue_class_group_entry->GetGausVariancePrior(&group_key);
}

bool PotentialEntryIterator::CheckGroupKey(ElementKeyType & group_key, bool verbose) const
{
    const auto & group_key_set{ m_element_class_group_entry->GetGroupKeySet() };
    if (group_key_set.find(group_key) == group_key_set.end())
    {
        if (verbose == true)
        {
            std::cout <<"Elelemt class group key :"
                      << std::get<0>(group_key) <<", "
                      << std::get<1>(group_key) <<", "
                      << std::get<2>(group_key) << std::boolalpha <<" not found." << std::endl;
        }
        return false;
    }
    return true;
}

bool PotentialEntryIterator::CheckGroupKey(ResidueKeyType & group_key, bool verbose) const
{
    const auto & group_key_set{ m_residue_class_group_entry->GetGroupKeySet() };
    if (group_key_set.find(group_key) == group_key_set.end())
    {
        if (verbose == true)
        {
            std::cout <<"Residue class group key : tuple<"
                      << std::get<0>(group_key) <<", "
                      << std::get<1>(group_key) <<", "
                      << std::get<2>(group_key) <<", "
                      << std::get<3>(group_key) <<", "
                      << std::get<4>(group_key) << std::boolalpha <<"> not found." << std::endl;
        }
        return false;
    }
    return true;
}