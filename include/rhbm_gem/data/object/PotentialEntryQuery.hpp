#pragma once

#include <cstddef>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

#include <rhbm_gem/utils/domain/GlobalEnumClass.hpp>

namespace rhbm_gem {

class AtomObject;
class BondObject;
class LocalPotentialEntry;
class LocalPotentialView;
class ModelPotentialView;
class ModelObject;

class PotentialEntryQuery
{
    AtomObject * m_atom_object;
    BondObject * m_bond_object;
    ModelObject * m_model_object;
    LocalPotentialEntry * m_atom_local_entry;
    LocalPotentialEntry * m_bond_local_entry;

public:
    explicit PotentialEntryQuery(ModelObject * model_object);
    explicit PotentialEntryQuery(AtomObject * atom_object);
    explicit PotentialEntryQuery(BondObject * bond_object);
    ~PotentialEntryQuery();

    double GetAtomGausEstimateMinimum(int par_id, Element element) const;
    double GetBondGausEstimateMinimum(int par_id) const;
    bool IsOutlierAtom(const std::string & class_key) const;
    bool IsOutlierBond(const std::string & class_key) const;
    bool IsAvailableAtomGroupKey(GroupKey group_key, const std::string & class_key, bool varbose=false) const;
    bool IsAvailableBondGroupKey(GroupKey group_key, const std::string & class_key, bool varbose=false) const;
    double GetAtomGausEstimateMean(GroupKey group_key, const std::string & class_key, int par_id) const;
    double GetAtomGausEstimatePrior(GroupKey group_key, const std::string & class_key, int par_id) const;
    double GetBondGausEstimatePrior(GroupKey group_key, const std::string & class_key, int par_id) const;
    double GetAtomGausVariancePrior(GroupKey group_key, const std::string & class_key, int par_id) const;
    double GetBondGausVariancePrior(GroupKey group_key, const std::string & class_key, int par_id) const;
    const std::vector<AtomObject *> & GetAtomObjectList(GroupKey group_key, const std::string & class_key) const;
    const std::vector<BondObject *> & GetBondObjectList(GroupKey group_key, const std::string & class_key) const;
    std::vector<AtomObject *> GetOutlierAtomObjectList(GroupKey group_key, const std::string & class_key) const;
    std::unordered_map<int, AtomObject *> GetAtomObjectMap(GroupKey group_key, const std::string & class_key) const;

    std::vector<std::tuple<float, float>> GetLinearModelDistanceAndMapValueList() const;
    const std::vector<std::tuple<float, float>> & GetDistanceAndMapValueList() const;
    std::vector<std::tuple<float, float>> GetBinnedDistanceAndMapValueList(int bin_size=15, double x_min=0.0, double x_max=1.5) const;
    std::tuple<float, float> GetDistanceRange(double margin_rate=0.0) const;
    std::tuple<float, float> GetMapValueRange(double margin_rate=0.0) const;
    std::tuple<float, float> GetBinnedMapValueRange(int bin_size=15, double x_min=0.0, double x_max=1.5, double margin_rate=0.0) const;

    double GetAmplitudeEstimateMDPDE() const;
    double GetAmplitudeEstimatePosterior(const std::string & key) const;
    double GetAmplitudeVariancePosterior(const std::string & key) const;
    double GetWidthEstimateMDPDE() const;
    double GetWidthEstimatePosterior(const std::string & key) const;
    double GetWidthVariancePosterior(const std::string & key) const;
    double GetAlphaR() const;
    double GetAtomAlphaR(GroupKey group_key, const std::string & class_key) const;
    double GetBondAlphaR(GroupKey group_key, const std::string & class_key) const;
    double GetAtomAlphaG(GroupKey group_key, const std::string & class_key) const;
    double GetBondAlphaG(GroupKey group_key, const std::string & class_key) const;
    Residue GetResidueFromAtomGroupKey(GroupKey group_key, const std::string & class_key) const;
    Residue GetResidueFromBondGroupKey(GroupKey group_key, const std::string & class_key) const;

    AtomObject * GetAtomObject() const;
    BondObject * GetBondObject() const;
    ModelObject * GetModelObject() const;
    LocalPotentialEntry * GetAtomLocalEntry() const;
    LocalPotentialEntry * GetBondLocalEntry() const;

private:
    size_t GetAtomResidueCount(const std::string & class_key, Residue residue, Structure structure=static_cast<Structure>(0)) const;
    size_t GetBondResidueCount(const std::string & class_key, Residue residue) const;
    ModelPotentialView GetModelView() const;
    LocalPotentialView GetLocalView() const;
};

} // namespace rhbm_gem
