#include <rhbm_gem/data/object/ModelAnalysisView.hpp>

#include "data/detail/ModelAnalysisData.hpp"
#include "data/detail/LocalPotentialEntry.hpp"

#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/BondObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/domain/ChemicalDataHelper.hpp>
#include <rhbm_gem/utils/domain/KeyPacker.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/math/ArrayHelper.hpp>

#include <sstream>
#include <stdexcept>

namespace rhbm_gem {

namespace {

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

template <typename ClassKeyFn, typename GroupCountFn>
std::string DescribeGrouping(
    const std::string & title,
    size_t class_count,
    ClassKeyFn get_class_key,
    GroupCountFn get_group_count)
{
    std::ostringstream oss;
    oss << title;
    for (size_t i = 0; i < class_count; i++)
    {
        const auto & class_key{ get_class_key(i) };
        oss << '\n'
            << " - Class type: " << class_key
            << " include " << get_group_count(class_key) << " groups.";
    }
    return oss.str();
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

LocalPotentialView LocalPotentialView::RequireFor(const AtomObject & atom_object)
{
    auto view{ LocalPotentialView::For(atom_object) };
    (void)view.RequireEntry("Atom local analysis");
    return view;
}

LocalPotentialView LocalPotentialView::RequireFor(const BondObject & bond_object)
{
    auto view{ LocalPotentialView::For(bond_object) };
    (void)view.RequireEntry("Bond local analysis");
    return view;
}

bool LocalPotentialView::IsAvailable() const
{
    return FindEntry() != nullptr;
}

const LocalPotentialEntry * LocalPotentialView::FindEntry() const
{
    if (m_atom_object != nullptr)
    {
        return ModelAnalysisData::FindLocalEntry(*m_atom_object);
    }
    if (m_bond_object != nullptr)
    {
        return ModelAnalysisData::FindLocalEntry(*m_bond_object);
    }
    return nullptr;
}

const LocalPotentialEntry & LocalPotentialView::RequireEntry(const char * context) const
{
    return ModelAnalysisData::RequireLocalEntry(FindEntry(), context);
}

const LocalGaussianResult & LocalPotentialView::GetGaussianResult() const
{
    return RequireEntry("Local Gaussian result").GetGaussianResult();
}

const GaussianModel3D & LocalPotentialView::GetEstimateOLS() const
{
    return RequireEntry("Local estimate OLS").GetEstimateOLS();
}

const GaussianModel3D & LocalPotentialView::GetEstimateMDPDE() const
{
    return RequireEntry("Local estimate MDPDE").GetEstimateMDPDE();
}

const LocalPotentialSampleList & LocalPotentialView::GetSamplingEntries() const
{
    return RequireEntry("Local sampling entries").GetSamplingEntries();
}

std::tuple<float, float> LocalPotentialView::GetDistanceRange(double margin_rate) const
{
    return RequireEntry("Local distance range").GetDistanceRange(margin_rate);
}

std::tuple<float, float> LocalPotentialView::GetResponseRange(double margin_rate) const
{
    return RequireEntry("Local response range").GetResponseRange(margin_rate);
}

SeriesPointList LocalPotentialView::GetBinnedDistanceResponseSeries(
    int bin_size,
    double x_min,
    double x_max) const
{
    return RequireEntry("Local binned distance-response series")
        .GetBinnedDistanceResponseSeries(bin_size, x_min, x_max);
}

double LocalPotentialView::GetAlphaR() const
{
    return RequireEntry("Local alpha-r").GetAlphaR();
}

std::optional<LocalPotentialAnnotation> LocalPotentialView::FindAnnotation(
    const std::string & key) const
{
    const auto * annotation{ RequireEntry("Local annotation").FindAnnotation(key) };
    if (annotation == nullptr)
    {
        return std::nullopt;
    }
    return LocalPotentialAnnotation{
        annotation->gaussian,
        annotation->is_outlier,
        annotation->statistical_distance
    };
}

double LocalPotentialView::GetMapValueNearCenter() const
{
    return RequireEntry("Local center map value").GetMapValueNearCenter();
}

double LocalPotentialView::CalculateQScore(int par_choice) const
{
    return RequireEntry("Local q-score").CalculateQScore(par_choice);
}

ModelAnalysisView::ModelAnalysisView(const ModelObject & model_object) :
    m_model_object(model_object)
{
}

bool ModelAnalysisView::HasGroupedAnalysisData() const
{
    const auto & analysis_data{ ModelAnalysisData::Of(m_model_object) };
    return !analysis_data.AtomGroupEntries().empty() || !analysis_data.BondGroupEntries().empty();
}

std::string ModelAnalysisView::DescribeAtomGrouping() const
{
    return DescribeGrouping(
        "Atom Grouping Summary:",
        ChemicalDataHelper::GetGroupAtomClassCount(),
        ChemicalDataHelper::GetGroupAtomClassKey,
        [this](const std::string & class_key)
        {
            return CollectAtomGroupKeys(class_key).size();
        });
}

double ModelAnalysisView::GetAtomGausEstimateMinimum(int par_id, Element element) const
{
    std::vector<double> gaus_estimate_list;
    gaus_estimate_list.reserve(m_model_object.GetSelectedAtomCount());
    for (const auto * atom : m_model_object.GetSelectedAtoms())
    {
        if (atom->GetElement() != element) continue;
        gaus_estimate_list.emplace_back(
            LocalPotentialView::RequireFor(*atom).GetEstimateMDPDE().GetDisplayParameter(par_id));
    }
    return array_helper::ComputeMin(gaus_estimate_list.data(), gaus_estimate_list.size());
}

double ModelAnalysisView::GetBondGausEstimateMinimum(int par_id) const
{
    std::vector<double> gaus_estimate_list;
    gaus_estimate_list.reserve(m_model_object.GetSelectedBondCount());
    for (const auto * bond : m_model_object.GetSelectedBonds())
    {
        gaus_estimate_list.emplace_back(
            LocalPotentialView::RequireFor(*bond).GetEstimateMDPDE().GetDisplayParameter(par_id));
    }
    return array_helper::ComputeMin(gaus_estimate_list.data(), gaus_estimate_list.size());
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

const GaussianModel3D & ModelAnalysisView::GetAtomGroupMean(
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

const GaussianModel3D & ModelAnalysisView::GetBondGroupMean(
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

const GaussianModel3D & ModelAnalysisView::GetAtomGroupMDPDE(
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

const GaussianModel3D & ModelAnalysisView::GetBondGroupMDPDE(
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

const GaussianModel3D & ModelAnalysisView::GetAtomGroupPrior(
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

const GaussianModel3D & ModelAnalysisView::GetBondGroupPrior(
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

GaussianModel3DWithUncertainty ModelAnalysisView::GetAtomGroupPriorWithUncertainty(
    GroupKey group_key,
    const std::string & class_key) const
{
    const auto * entry{
        RequireGroupEntry(
            ModelAnalysisData::Of(m_model_object).FindAtomGroupEntry(class_key),
            "Atom group entry")
    };
    return entry->GetPriorWithUncertainty(group_key);
}

GaussianModel3DWithUncertainty ModelAnalysisView::GetBondGroupPriorWithUncertainty(
    GroupKey group_key,
    const std::string & class_key) const
{
    const auto * entry{
        RequireGroupEntry(
            ModelAnalysisData::Of(m_model_object).FindBondGroupEntry(class_key),
            "Bond group entry")
    };
    return entry->GetPriorWithUncertainty(group_key);
}

double ModelAnalysisView::GetAtomGausEstimatePrior(
    GroupKey group_key,
    const std::string & class_key,
    int par_id) const
{
    return GetAtomGroupPrior(group_key, class_key).GetDisplayParameter(par_id);
}

double ModelAnalysisView::GetBondGausEstimatePrior(
    GroupKey group_key,
    const std::string & class_key,
    int par_id) const
{
    return GetBondGroupPrior(group_key, class_key).GetDisplayParameter(par_id);
}

double ModelAnalysisView::GetAtomGausPriorStandardDeviation(
    GroupKey group_key,
    const std::string & class_key,
    int par_id) const
{
    return GetAtomGroupPriorWithUncertainty(group_key, class_key).GetDisplayStandardDeviation(par_id);
}

double ModelAnalysisView::GetBondGausPriorStandardDeviation(
    GroupKey group_key,
    const std::string & class_key,
    int par_id) const
{
    return GetBondGroupPriorWithUncertainty(group_key, class_key).GetDisplayStandardDeviation(par_id);
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
        const auto annotation{ LocalPotentialView::RequireFor(*atom).FindAnnotation(class_key) };
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
    return LocalPotentialView::RequireFor(*atom_list.front()).GetAlphaR();
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
