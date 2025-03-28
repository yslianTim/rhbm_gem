#pragma once

#include <tuple>
#include <vector>
#include <unordered_map>
#include "AtomicInfoHelper.hpp"
#include "GroupPotentialEntry.hpp"

class AtomObject;
class ModelObject;
class AtomicPotentialEntry;

class PotentialEntryIterator
{
    using ElementKeyType = GroupKeyMapping<ElementGroupClassifierTag>::type;
    using ResidueKeyType = GroupKeyMapping<ResidueGroupClassifierTag>::type;

    AtomObject * m_atom_object;
    ModelObject * m_model_object;
    AtomicPotentialEntry * m_atomic_entry;
    GroupPotentialEntry<ElementKeyType> * m_element_class_group_entry;
    GroupPotentialEntry<ResidueKeyType> * m_residue_class_group_entry;

public:
    PotentialEntryIterator(ModelObject * model_object);
    PotentialEntryIterator(AtomObject * atom_object);
    ~PotentialEntryIterator();
    bool IsAvailableGroupKey(ElementKeyType & group_key) const;
    bool IsAvailableGroupKey(ResidueKeyType & group_key) const;
    std::tuple<double, double> GetGausEstimatePrior(ElementKeyType & group_key) const;
    std::tuple<double, double> GetGausEstimatePrior(ResidueKeyType & group_key) const;
    std::tuple<double, double> GetGausVariancePrior(ElementKeyType & group_key) const;
    std::tuple<double, double> GetGausVariancePrior(ResidueKeyType & group_key) const;

    const std::vector<std::tuple<float, float>> & GetDistanceAndMapValueList(void) const;

private:
    bool CheckGroupKey(ElementKeyType & group_key, bool verbose=true) const;
    bool CheckGroupKey(ResidueKeyType & group_key, bool verbose=true) const;

};