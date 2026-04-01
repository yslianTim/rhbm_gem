#include <rhbm_gem/data/object/PotentialEntryQuery.hpp>

#include <stdexcept>

#include <rhbm_gem/data/object/AtomClassifier.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/BondClassifier.hpp>
#include <rhbm_gem/data/object/BondObject.hpp>
#include <rhbm_gem/data/object/LocalPotentialEntry.hpp>
#include <rhbm_gem/data/object/LocalPotentialView.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/data/object/ModelPotentialView.hpp>
#include <rhbm_gem/utils/domain/ChemicalDataHelper.hpp>

namespace rhbm_gem {

PotentialEntryQuery::PotentialEntryQuery(ModelObject * model_object) :
    m_atom_object{ nullptr }, m_bond_object{ nullptr }, m_model_object{ model_object },
    m_atom_local_entry{ nullptr }, m_bond_local_entry{ nullptr }
{
}

PotentialEntryQuery::PotentialEntryQuery(AtomObject * atom_object) :
    m_atom_object{ atom_object }, m_bond_object{ nullptr }, m_model_object{ nullptr },
    m_atom_local_entry{ atom_object != nullptr ? atom_object->GetLocalPotentialEntry() : nullptr },
    m_bond_local_entry{ nullptr }
{
}

PotentialEntryQuery::PotentialEntryQuery(BondObject * bond_object) :
    m_atom_object{ nullptr }, m_bond_object{ bond_object }, m_model_object{ nullptr },
    m_atom_local_entry{ nullptr },
    m_bond_local_entry{ bond_object != nullptr ? bond_object->GetLocalPotentialEntry() : nullptr }
{
}

PotentialEntryQuery::~PotentialEntryQuery() = default;

double PotentialEntryQuery::GetAtomGausEstimateMinimum(int par_id, Element element) const
{
    return GetModelView().GetAtomGausEstimateMinimum(par_id, element);
}

double PotentialEntryQuery::GetBondGausEstimateMinimum(int par_id) const
{
    return GetModelView().GetBondGausEstimateMinimum(par_id);
}

bool PotentialEntryQuery::IsOutlierAtom(const std::string & class_key) const
{
    if (m_atom_local_entry == nullptr)
    {
        throw std::runtime_error("Atom local entry is not available.");
    }
    return m_atom_local_entry->GetOutlierTag(class_key);
}

bool PotentialEntryQuery::IsOutlierBond(const std::string & class_key) const
{
    if (m_bond_local_entry == nullptr)
    {
        throw std::runtime_error("Bond local entry is not available.");
    }
    return m_bond_local_entry->GetOutlierTag(class_key);
}

bool PotentialEntryQuery::IsAvailableAtomGroupKey(
    GroupKey group_key, const std::string & class_key, bool verbose) const
{
    return GetModelView().HasAtomGroup(group_key, class_key, verbose);
}

bool PotentialEntryQuery::IsAvailableBondGroupKey(
    GroupKey group_key, const std::string & class_key, bool verbose) const
{
    return GetModelView().HasBondGroup(group_key, class_key, verbose);
}

size_t PotentialEntryQuery::GetAtomResidueCount(
    const std::string & class_key, Residue residue, Structure structure) const
{
    GroupKey group_key{ 0 };
    AtomClassifier classifier;
    if (class_key == ChemicalDataHelper::GetComponentAtomClassKey())
    {
        group_key = classifier.GetMainChainComponentAtomClassGroupKey(0, residue);
    }
    else if (class_key == ChemicalDataHelper::GetStructureAtomClassKey())
    {
        group_key = classifier.GetMainChainStructureAtomClassGroupKey(0, structure, residue);
    }
    if (!IsAvailableAtomGroupKey(group_key, class_key))
    {
        return 0;
    }
    return GetAtomObjectList(group_key, class_key).size();
}

size_t PotentialEntryQuery::GetBondResidueCount(
    const std::string & class_key, Residue residue) const
{
    GroupKey group_key{ 0 };
    BondClassifier classifier;
    if (class_key == ChemicalDataHelper::GetComponentBondClassKey())
    {
        group_key = classifier.GetMainChainComponentBondClassGroupKey(0, residue);
    }
    if (!IsAvailableBondGroupKey(group_key, class_key))
    {
        return 0;
    }
    return GetBondObjectList(group_key, class_key).size();
}

double PotentialEntryQuery::GetAtomGausEstimateMean(
    GroupKey group_key, const std::string & class_key, int par_id) const
{
    return GetModelView().GetAtomGausEstimateMean(group_key, class_key, par_id);
}

double PotentialEntryQuery::GetAtomGausEstimatePrior(
    GroupKey group_key, const std::string & class_key, int par_id) const
{
    return GetModelView().GetAtomGausEstimatePrior(group_key, class_key, par_id);
}

double PotentialEntryQuery::GetBondGausEstimatePrior(
    GroupKey group_key, const std::string & class_key, int par_id) const
{
    return GetModelView().GetBondGausEstimatePrior(group_key, class_key, par_id);
}

double PotentialEntryQuery::GetAtomGausVariancePrior(
    GroupKey group_key, const std::string & class_key, int par_id) const
{
    return GetModelView().GetAtomGausVariancePrior(group_key, class_key, par_id);
}

double PotentialEntryQuery::GetBondGausVariancePrior(
    GroupKey group_key, const std::string & class_key, int par_id) const
{
    return GetModelView().GetBondGausVariancePrior(group_key, class_key, par_id);
}

const std::vector<AtomObject *> & PotentialEntryQuery::GetAtomObjectList(
    GroupKey group_key, const std::string & class_key) const
{
    return GetModelView().GetAtomObjectList(group_key, class_key);
}

const std::vector<BondObject *> & PotentialEntryQuery::GetBondObjectList(
    GroupKey group_key, const std::string & class_key) const
{
    return GetModelView().GetBondObjectList(group_key, class_key);
}

std::vector<AtomObject *> PotentialEntryQuery::GetOutlierAtomObjectList(
    GroupKey group_key, const std::string & class_key) const
{
    return GetModelView().GetOutlierAtomObjectList(group_key, class_key);
}

std::unordered_map<int, AtomObject *> PotentialEntryQuery::GetAtomObjectMap(
    GroupKey group_key, const std::string & class_key) const
{
    return GetModelView().GetAtomObjectMap(group_key, class_key);
}

std::vector<std::tuple<float, float>>
PotentialEntryQuery::GetLinearModelDistanceAndMapValueList() const
{
    return GetLocalView().GetLinearModelDistanceAndMapValueList();
}

const std::vector<std::tuple<float, float>> &
PotentialEntryQuery::GetDistanceAndMapValueList() const
{
    return GetLocalView().GetDistanceAndMapValueList();
}

std::vector<std::tuple<float, float>> PotentialEntryQuery::GetBinnedDistanceAndMapValueList(
    int bin_size, double x_min, double x_max) const
{
    return GetLocalView().GetBinnedDistanceAndMapValueList(bin_size, x_min, x_max);
}

std::tuple<float, float> PotentialEntryQuery::GetDistanceRange(double margin_rate) const
{
    return GetLocalView().GetDistanceRange(margin_rate);
}

std::tuple<float, float> PotentialEntryQuery::GetMapValueRange(double margin_rate) const
{
    return GetLocalView().GetMapValueRange(margin_rate);
}

std::tuple<float, float> PotentialEntryQuery::GetBinnedMapValueRange(
    int bin_size, double x_min, double x_max, double margin_rate) const
{
    return GetLocalView().GetBinnedMapValueRange(bin_size, x_min, x_max, margin_rate);
}

double PotentialEntryQuery::GetAmplitudeEstimateMDPDE() const
{
    return GetLocalView().GetEstimateMDPDE().amplitude;
}

double PotentialEntryQuery::GetAmplitudeEstimatePosterior(const std::string & class_key) const
{
    return GetLocalView().GetPosterior(class_key).estimate.amplitude;
}

double PotentialEntryQuery::GetAmplitudeVariancePosterior(const std::string & class_key) const
{
    return GetLocalView().GetPosterior(class_key).variance.amplitude;
}

double PotentialEntryQuery::GetWidthEstimateMDPDE() const
{
    return GetLocalView().GetEstimateMDPDE().width;
}

double PotentialEntryQuery::GetWidthEstimatePosterior(const std::string & class_key) const
{
    return GetLocalView().GetPosterior(class_key).estimate.width;
}

double PotentialEntryQuery::GetWidthVariancePosterior(const std::string & class_key) const
{
    return GetLocalView().GetPosterior(class_key).variance.width;
}

double PotentialEntryQuery::GetAlphaR() const
{
    return GetLocalView().GetAlphaR();
}

double PotentialEntryQuery::GetAtomAlphaR(
    GroupKey group_key, const std::string & class_key) const
{
    return GetModelView().GetAtomAlphaR(group_key, class_key);
}

double PotentialEntryQuery::GetBondAlphaR(
    GroupKey group_key, const std::string & class_key) const
{
    return GetModelView().GetBondAlphaR(group_key, class_key);
}

double PotentialEntryQuery::GetAtomAlphaG(
    GroupKey group_key, const std::string & class_key) const
{
    return GetModelView().GetAtomAlphaG(group_key, class_key);
}

double PotentialEntryQuery::GetBondAlphaG(
    GroupKey group_key, const std::string & class_key) const
{
    return GetModelView().GetBondAlphaG(group_key, class_key);
}

Residue PotentialEntryQuery::GetResidueFromAtomGroupKey(
    GroupKey group_key, const std::string & class_key) const
{
    return GetModelView().GetResidueFromAtomGroupKey(group_key, class_key);
}

Residue PotentialEntryQuery::GetResidueFromBondGroupKey(
    GroupKey group_key, const std::string & class_key) const
{
    return GetModelView().GetResidueFromBondGroupKey(group_key, class_key);
}

AtomObject * PotentialEntryQuery::GetAtomObject() const
{
    return m_atom_object;
}

BondObject * PotentialEntryQuery::GetBondObject() const
{
    return m_bond_object;
}

ModelObject * PotentialEntryQuery::GetModelObject() const
{
    return m_model_object;
}

LocalPotentialEntry * PotentialEntryQuery::GetAtomLocalEntry() const
{
    return m_atom_local_entry;
}

LocalPotentialEntry * PotentialEntryQuery::GetBondLocalEntry() const
{
    return m_bond_local_entry;
}

ModelPotentialView PotentialEntryQuery::GetModelView() const
{
    if (m_model_object == nullptr)
    {
        throw std::runtime_error("Model object is not available.");
    }
    return ModelPotentialView(*m_model_object);
}

LocalPotentialView PotentialEntryQuery::GetLocalView() const
{
    if (m_atom_local_entry != nullptr)
    {
        return LocalPotentialView(*m_atom_local_entry);
    }
    if (m_bond_local_entry != nullptr)
    {
        return LocalPotentialView(*m_bond_local_entry);
    }
    throw std::runtime_error("Local entry is not available.");
}

} // namespace rhbm_gem
