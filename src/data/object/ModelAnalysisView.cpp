#include <rhbm_gem/data/object/ModelAnalysisView.hpp>

#include "data/detail/ModelAnalysisData.hpp"
#include "data/detail/LocalPotentialEntry.hpp"

#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/BondObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/domain/ChemicalDataHelper.hpp>
#include <rhbm_gem/utils/domain/KeyPacker.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/math/ArrayStats.hpp>

#include <sstream>
#include <stdexcept>

namespace rhbm_gem {

namespace {

const LocalPotentialEntry * FindResolvedLocalEntry(const LocalPotentialView & view)
{
    if (view.GetAtomObjectPtr() != nullptr)
    {
        return ModelAnalysisData::FindLocalEntry(*view.GetAtomObjectPtr());
    }
    if (view.GetBondObjectPtr() != nullptr)
    {
        return ModelAnalysisData::FindLocalEntry(*view.GetBondObjectPtr());
    }
    return nullptr;
}

const LocalPotentialEntry & RequireResolvedLocalEntry(
    const LocalPotentialView & view,
    const char * context)
{
    return ModelAnalysisData::RequireLocalEntry(FindResolvedLocalEntry(view), context);
}

template <typename EntryT>
std::vector<GroupKey> CollectGroupKeys(const EntryT * entry)
{
    return entry == nullptr ? std::vector<GroupKey>{} : entry->CollectGroupKeys();
}

template <typename EntryT>
const EntryT * RequireGroupEntry(
    const EntryT * entry,
    const char * context)
{
    if (entry == nullptr)
    {
        throw std::runtime_error(std::string(context) + " is not available.");
    }
    return entry;
}

} // namespace

LocalPotentialView::LocalPotentialView(const AtomObject * atom_object) :
    m_atom_object{ atom_object }
{
}

LocalPotentialView::LocalPotentialView(const BondObject * bond_object) :
    m_bond_object{ bond_object }
{
}

LocalPotentialView LocalPotentialView::For(const AtomObject & atom_object)
{
    return LocalPotentialView(&atom_object);
}

LocalPotentialView LocalPotentialView::For(const BondObject & bond_object)
{
    return LocalPotentialView(&bond_object);
}

bool LocalPotentialView::IsAvailable() const
{
    return FindResolvedLocalEntry(*this) != nullptr;
}

const GaussianEstimate & LocalPotentialView::GetEstimateOLS() const
{
    return RequireResolvedLocalEntry(*this, "Local estimate OLS").GetEstimateOLS();
}

const GaussianEstimate & LocalPotentialView::GetEstimateMDPDE() const
{
    return RequireResolvedLocalEntry(*this, "Local estimate MDPDE").GetEstimateMDPDE();
}

const std::vector<std::tuple<float, float>> & LocalPotentialView::GetDistanceAndMapValueList() const
{
    return RequireResolvedLocalEntry(*this, "Local distance-map series").GetDistanceAndMapValueList();
}

int LocalPotentialView::GetDistanceAndMapValueListSize() const
{
    return RequireResolvedLocalEntry(*this, "Local distance-map size").GetDistanceAndMapValueListSize();
}

double LocalPotentialView::GetAlphaR() const
{
    return RequireResolvedLocalEntry(*this, "Local alpha-r").GetAlphaR();
}

std::optional<LocalPotentialAnnotationView> LocalPotentialView::FindAnnotation(
    const std::string & key) const
{
    const auto * annotation{ RequireResolvedLocalEntry(*this, "Local annotation").FindAnnotation(key) };
    if (annotation == nullptr)
    {
        return std::nullopt;
    }
    return LocalPotentialAnnotationView{
        annotation->posterior,
        annotation->is_outlier,
        annotation->statistical_distance
    };
}

double LocalPotentialView::GetMapValueNearCenter() const
{
    return RequireResolvedLocalEntry(*this, "Local center map value").GetMapValueNearCenter();
}

double LocalPotentialView::GetMomentZeroEstimate() const
{
    return RequireResolvedLocalEntry(*this, "Local moment-zero estimate").GetMomentZeroEstimate();
}

double LocalPotentialView::GetMomentTwoEstimate() const
{
    return RequireResolvedLocalEntry(*this, "Local moment-two estimate").GetMomentTwoEstimate();
}

double LocalPotentialView::CalculateQScore(int par_choice) const
{
    return RequireResolvedLocalEntry(*this, "Local q-score").CalculateQScore(par_choice);
}

ModelAnalysisView::ModelAnalysisView(const ModelObject & model_object) :
    m_model_object(model_object)
{
}

ModelAnalysisView ModelAnalysisView::Of(const ModelObject & model_object)
{
    return ModelAnalysisView(model_object);
}

bool ModelAnalysisView::HasLocalAnalysis(const AtomObject & atom_object)
{
    return ModelAnalysisData::FindLocalEntry(atom_object) != nullptr;
}

bool ModelAnalysisView::HasLocalAnalysis(const BondObject & bond_object)
{
    return ModelAnalysisData::FindLocalEntry(bond_object) != nullptr;
}

LocalPotentialView ModelAnalysisView::FindLocalPotential(const AtomObject & atom_object)
{
    return LocalPotentialView::For(atom_object);
}

LocalPotentialView ModelAnalysisView::FindLocalPotential(const BondObject & bond_object)
{
    return LocalPotentialView::For(bond_object);
}

LocalPotentialView ModelAnalysisView::RequireLocalPotential(const AtomObject & atom_object)
{
    auto view{ LocalPotentialView::For(atom_object) };
    (void)RequireResolvedLocalEntry(view, "Atom local analysis");
    return view;
}

LocalPotentialView ModelAnalysisView::RequireLocalPotential(const BondObject & bond_object)
{
    auto view{ LocalPotentialView::For(bond_object) };
    (void)RequireResolvedLocalEntry(view, "Bond local analysis");
    return view;
}

bool ModelAnalysisView::HasGroupedAnalysisData() const
{
    const auto & analysis_data{ ModelAnalysisData::Of(m_model_object) };
    return !analysis_data.AtomGroupEntries().empty() || !analysis_data.BondGroupEntries().empty();
}

double ModelAnalysisView::GetAtomGausEstimateMinimum(int par_id, Element element) const
{
    std::vector<double> gaus_estimate_list;
    gaus_estimate_list.reserve(m_model_object.GetSelectedAtomCount());
    for (const auto * atom : m_model_object.GetSelectedAtoms())
    {
        if (atom->GetElement() != element) continue;
        gaus_estimate_list.emplace_back(RequireLocalPotential(*atom).GetEstimateMDPDE().GetParameter(par_id));
    }
    return ArrayStats<double>::ComputeMin(gaus_estimate_list.data(), gaus_estimate_list.size());
}

double ModelAnalysisView::GetBondGausEstimateMinimum(int par_id) const
{
    std::vector<double> gaus_estimate_list;
    gaus_estimate_list.reserve(m_model_object.GetSelectedBondCount());
    for (const auto * bond : m_model_object.GetSelectedBonds())
    {
        gaus_estimate_list.emplace_back(RequireLocalPotential(*bond).GetEstimateMDPDE().GetParameter(par_id));
    }
    return ArrayStats<double>::ComputeMin(gaus_estimate_list.data(), gaus_estimate_list.size());
}

bool ModelAnalysisView::HasAtomGroup(
    GroupKey group_key,
    const std::string & class_key,
    bool verbose) const
{
    const auto * entry{ ModelAnalysisData::Of(m_model_object).FindAtomGroupEntry(class_key) };
    if (entry == nullptr)
    {
        if (verbose)
        {
            Logger::Log(LogLevel::Error, "Atom class key not found: " + class_key);
        }
        return false;
    }
    if (!entry->HasGroup(group_key))
    {
        if (verbose)
        {
            std::ostringstream oss;
            if (class_key == ChemicalDataHelper::GetSimpleAtomClassKey())
            {
                oss << "Atom class group key : "
                    << KeyPackerSimpleAtomClass::GetKeyString(group_key)
                    << " not found.";
            }
            else if (class_key == ChemicalDataHelper::GetComponentAtomClassKey())
            {
                oss << "Atom class group key : "
                    << KeyPackerComponentAtomClass::GetKeyString(group_key)
                    << " not found.";
            }
            else if (class_key == ChemicalDataHelper::GetStructureAtomClassKey())
            {
                oss << "Atom class group key : "
                    << KeyPackerStructureAtomClass::GetKeyString(group_key)
                    << " not found.";
            }
            Logger::Log(LogLevel::Error, oss.str());
        }
        return false;
    }
    return true;
}

bool ModelAnalysisView::HasBondGroup(
    GroupKey group_key,
    const std::string & class_key,
    bool verbose) const
{
    const auto * entry{ ModelAnalysisData::Of(m_model_object).FindBondGroupEntry(class_key) };
    if (entry == nullptr)
    {
        if (verbose)
        {
            Logger::Log(LogLevel::Error, "Bond class key not found: " + class_key);
        }
        return false;
    }
    if (!entry->HasGroup(group_key))
    {
        if (verbose)
        {
            std::ostringstream oss;
            if (class_key == ChemicalDataHelper::GetSimpleBondClassKey())
            {
                oss << "Bond class group key : "
                    << KeyPackerSimpleBondClass::GetKeyString(group_key)
                    << " not found.";
            }
            else if (class_key == ChemicalDataHelper::GetComponentBondClassKey())
            {
                oss << "Bond class group key : "
                    << KeyPackerComponentBondClass::GetKeyString(group_key)
                    << " not found.";
            }
            Logger::Log(LogLevel::Error, oss.str());
        }
        return false;
    }
    return true;
}

const GaussianEstimate & ModelAnalysisView::GetAtomGroupMean(
    GroupKey group_key,
    const std::string & class_key) const
{
    const auto * entry{
        RequireGroupEntry(
            ModelAnalysisData::Of(m_model_object).FindAtomGroupEntry(class_key),
            "Atom group entry")
    };
    return entry->GetMean(group_key);
}

const GaussianEstimate & ModelAnalysisView::GetBondGroupMean(
    GroupKey group_key,
    const std::string & class_key) const
{
    const auto * entry{
        RequireGroupEntry(
            ModelAnalysisData::Of(m_model_object).FindBondGroupEntry(class_key),
            "Bond group entry")
    };
    return entry->GetMean(group_key);
}

const GaussianEstimate & ModelAnalysisView::GetAtomGroupMDPDE(
    GroupKey group_key,
    const std::string & class_key) const
{
    const auto * entry{
        RequireGroupEntry(
            ModelAnalysisData::Of(m_model_object).FindAtomGroupEntry(class_key),
            "Atom group entry")
    };
    return entry->GetMDPDE(group_key);
}

const GaussianEstimate & ModelAnalysisView::GetBondGroupMDPDE(
    GroupKey group_key,
    const std::string & class_key) const
{
    const auto * entry{
        RequireGroupEntry(
            ModelAnalysisData::Of(m_model_object).FindBondGroupEntry(class_key),
            "Bond group entry")
    };
    return entry->GetMDPDE(group_key);
}

const GaussianEstimate & ModelAnalysisView::GetAtomGroupPrior(
    GroupKey group_key,
    const std::string & class_key) const
{
    const auto * entry{
        RequireGroupEntry(
            ModelAnalysisData::Of(m_model_object).FindAtomGroupEntry(class_key),
            "Atom group entry")
    };
    return entry->GetPrior(group_key);
}

const GaussianEstimate & ModelAnalysisView::GetBondGroupPrior(
    GroupKey group_key,
    const std::string & class_key) const
{
    const auto * entry{
        RequireGroupEntry(
            ModelAnalysisData::Of(m_model_object).FindBondGroupEntry(class_key),
            "Bond group entry")
    };
    return entry->GetPrior(group_key);
}

GaussianPosterior ModelAnalysisView::GetAtomGroupPriorPosterior(
    GroupKey group_key,
    const std::string & class_key) const
{
    const auto * entry{
        RequireGroupEntry(
            ModelAnalysisData::Of(m_model_object).FindAtomGroupEntry(class_key),
            "Atom group entry")
    };
    return entry->BuildPriorPosterior(group_key);
}

GaussianPosterior ModelAnalysisView::GetBondGroupPriorPosterior(
    GroupKey group_key,
    const std::string & class_key) const
{
    const auto * entry{
        RequireGroupEntry(
            ModelAnalysisData::Of(m_model_object).FindBondGroupEntry(class_key),
            "Bond group entry")
    };
    return entry->BuildPriorPosterior(group_key);
}

double ModelAnalysisView::GetAtomGausEstimatePrior(
    GroupKey group_key,
    const std::string & class_key,
    int par_id) const
{
    return GetAtomGroupPrior(group_key, class_key).GetParameter(par_id);
}

double ModelAnalysisView::GetBondGausEstimatePrior(
    GroupKey group_key,
    const std::string & class_key,
    int par_id) const
{
    return GetBondGroupPrior(group_key, class_key).GetParameter(par_id);
}

double ModelAnalysisView::GetAtomGausVariancePrior(
    GroupKey group_key,
    const std::string & class_key,
    int par_id) const
{
    return GetAtomGroupPriorPosterior(group_key, class_key).GetVariance(par_id);
}

double ModelAnalysisView::GetBondGausVariancePrior(
    GroupKey group_key,
    const std::string & class_key,
    int par_id) const
{
    return GetBondGroupPriorPosterior(group_key, class_key).GetVariance(par_id);
}

const std::vector<AtomObject *> & ModelAnalysisView::GetAtomObjectList(
    GroupKey group_key,
    const std::string & class_key) const
{
    const auto * entry{
        RequireGroupEntry(
            ModelAnalysisData::Of(m_model_object).FindAtomGroupEntry(class_key),
            "Atom group entry")
    };
    return entry->GetMembers(group_key);
}

const std::vector<BondObject *> & ModelAnalysisView::GetBondObjectList(
    GroupKey group_key,
    const std::string & class_key) const
{
    const auto * entry{
        RequireGroupEntry(
            ModelAnalysisData::Of(m_model_object).FindBondGroupEntry(class_key),
            "Bond group entry")
    };
    return entry->GetMembers(group_key);
}

std::vector<AtomObject *> ModelAnalysisView::GetOutlierAtomObjectList(
    GroupKey group_key,
    const std::string & class_key) const
{
    std::vector<AtomObject *> outlier_atom_list;
    const auto & atom_list{ GetAtomObjectList(group_key, class_key) };
    outlier_atom_list.reserve(atom_list.size());
    for (auto * atom : atom_list)
    {
        const auto annotation{ RequireLocalPotential(*atom).FindAnnotation(class_key) };
        if (annotation.has_value() && annotation->is_outlier)
        {
            outlier_atom_list.emplace_back(atom);
        }
    }
    return outlier_atom_list;
}

std::unordered_map<int, AtomObject *> ModelAnalysisView::GetAtomObjectMap(
    GroupKey group_key,
    const std::string & class_key) const
{
    std::unordered_map<int, AtomObject *> atom_object_map;
    for (auto * atom_object : GetAtomObjectList(group_key, class_key))
    {
        atom_object_map[atom_object->GetSerialID()] = atom_object;
    }
    return atom_object_map;
}

double ModelAnalysisView::GetAtomAlphaR(
    GroupKey group_key,
    const std::string & class_key) const
{
    const auto & atom_list{ GetAtomObjectList(group_key, class_key) };
    if (atom_list.empty())
    {
        throw std::runtime_error("Atom group has no members.");
    }
    return RequireLocalPotential(*atom_list.front()).GetAlphaR();
}

double ModelAnalysisView::GetAtomAlphaG(
    GroupKey group_key,
    const std::string & class_key) const
{
    const auto * entry{
        RequireGroupEntry(
            ModelAnalysisData::Of(m_model_object).FindAtomGroupEntry(class_key),
            "Atom group entry")
    };
    return entry->GetAlphaG(group_key);
}

double ModelAnalysisView::GetBondAlphaG(
    GroupKey group_key,
    const std::string & class_key) const
{
    const auto * entry{
        RequireGroupEntry(
            ModelAnalysisData::Of(m_model_object).FindBondGroupEntry(class_key),
            "Bond group entry")
    };
    return entry->GetAlphaG(group_key);
}

Residue ModelAnalysisView::DecodeResidueFromAtomGroupKey(
    GroupKey group_key,
    const std::string & class_key) const
{
    if (class_key == ChemicalDataHelper::GetSimpleAtomClassKey())
    {
        Logger::Log(LogLevel::Error, "Simple class group key is not recording Residue info.");
        return Residue::UNK;
    }
    if (class_key == ChemicalDataHelper::GetComponentAtomClassKey())
    {
        auto unpack_key{ KeyPackerComponentAtomClass::Unpack(group_key) };
        return static_cast<Residue>(std::get<0>(unpack_key));
    }
    if (class_key == ChemicalDataHelper::GetStructureAtomClassKey())
    {
        auto unpack_key{ KeyPackerStructureAtomClass::Unpack(group_key) };
        return static_cast<Residue>(std::get<1>(unpack_key));
    }
    return Residue::UNK;
}

Residue ModelAnalysisView::DecodeResidueFromBondGroupKey(
    GroupKey group_key,
    const std::string & class_key) const
{
    if (class_key == ChemicalDataHelper::GetSimpleBondClassKey())
    {
        Logger::Log(LogLevel::Error, "Simple class group key is not recording Residue info.");
        return Residue::UNK;
    }
    if (class_key == ChemicalDataHelper::GetComponentBondClassKey())
    {
        auto unpack_key{ KeyPackerComponentBondClass::Unpack(group_key) };
        return static_cast<Residue>(std::get<0>(unpack_key));
    }
    return Residue::UNK;
}

Residue ModelAnalysisView::GetResidueFromAtomGroupKey(
    GroupKey group_key,
    const std::string & class_key) const
{
    return DecodeResidueFromAtomGroupKey(group_key, class_key);
}

Residue ModelAnalysisView::GetResidueFromBondGroupKey(
    GroupKey group_key,
    const std::string & class_key) const
{
    return DecodeResidueFromBondGroupKey(group_key, class_key);
}

std::vector<GroupKey> ModelAnalysisView::CollectAtomGroupKeys(const std::string & class_key) const
{
    return CollectGroupKeys(ModelAnalysisData::Of(m_model_object).FindAtomGroupEntry(class_key));
}

std::vector<GroupKey> ModelAnalysisView::CollectBondGroupKeys(const std::string & class_key) const
{
    return CollectGroupKeys(ModelAnalysisData::Of(m_model_object).FindBondGroupEntry(class_key));
}

} // namespace rhbm_gem
