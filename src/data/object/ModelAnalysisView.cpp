#include <rhbm_gem/data/object/ModelAnalysisView.hpp>

#include "data/detail/ModelAnalysisData.hpp"
#include "data/detail/LocalPotentialEntry.hpp"

#include <rhbm_gem/data/object/AtomObject.hpp>
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

const LocalPotentialEntry & RequireLocalEntry(
    const LocalPotentialEntry * entry,
    const char * context)
{
    if (entry == nullptr)
    {
        throw std::runtime_error(std::string(context) + " is not available.");
    }
    return *entry;
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

AtomLocalPotentialView::AtomLocalPotentialView(const AtomObject * atom_object) :
    m_atom_object{ atom_object }
{
}

AtomLocalPotentialView AtomLocalPotentialView::For(const AtomObject & atom_object)
{
    return AtomLocalPotentialView(&atom_object);
}

AtomLocalPotentialView AtomLocalPotentialView::RequireFor(const AtomObject & atom_object)
{
    auto view{ AtomLocalPotentialView::For(atom_object) };
    (void)view.RequireEntry("Atom local analysis");
    return view;
}

bool AtomLocalPotentialView::IsAvailable() const
{
    return FindEntry() != nullptr;
}

const LocalPotentialEntry * AtomLocalPotentialView::FindEntry() const
{
    if (m_atom_object != nullptr && m_atom_object->m_owner_model != nullptr)
    {
        return ModelAnalysisData::Of(*m_atom_object->m_owner_model).FindAtomLocalEntry(*m_atom_object);
    }
    return nullptr;
}

const LocalPotentialEntry & AtomLocalPotentialView::RequireEntry(const char * context) const
{
    return RequireLocalEntry(FindEntry(), context);
}

const LocalGaussianResult & AtomLocalPotentialView::GetGaussianResult() const
{
    return RequireEntry("Local Gaussian result").GetGaussianResult();
}

const GaussianModel3D & AtomLocalPotentialView::GetEstimateOLS() const
{
    return RequireEntry("Local estimate OLS").GetEstimateOLS();
}

const GaussianModel3D & AtomLocalPotentialView::GetEstimateMDPDE() const
{
    return RequireEntry("Local estimate MDPDE").GetEstimateMDPDE();
}

const LocalPotentialSampleList & AtomLocalPotentialView::GetSamplingEntries() const
{
    return RequireEntry("Local sampling entries").GetSamplingEntries();
}

std::tuple<float, float> AtomLocalPotentialView::GetDistanceRange(double margin_rate) const
{
    return RequireEntry("Local distance range").GetDistanceRange(margin_rate);
}

std::tuple<float, float> AtomLocalPotentialView::GetResponseRange(double margin_rate) const
{
    return RequireEntry("Local response range").GetResponseRange(margin_rate);
}

SeriesPointList AtomLocalPotentialView::GetBinnedDistanceResponseSeries(
    int bin_size,
    double x_min,
    double x_max) const
{
    return RequireEntry("Local binned distance-response series")
        .GetBinnedDistanceResponseSeries(bin_size, x_min, x_max);
}

double AtomLocalPotentialView::GetAlphaR() const
{
    return RequireEntry("Local alpha-r").GetAlphaR();
}

std::optional<LocalPotentialAnnotation> AtomLocalPotentialView::FindAnnotation(
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

double AtomLocalPotentialView::GetMapValueNearCenter() const
{
    return RequireEntry("Local center map value").GetMapValueNearCenter();
}

double AtomLocalPotentialView::CalculateQScore(int par_choice) const
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
    return !analysis_data.AtomGroupEntries().empty();
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
            AtomLocalPotentialView::RequireFor(*atom).GetEstimateMDPDE().GetDisplayParameter(par_id));
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

double ModelAnalysisView::GetAtomGausEstimatePrior(
    GroupKey group_key,
    const std::string & class_key,
    int par_id) const
{
    return GetAtomGroupPrior(group_key, class_key).GetDisplayParameter(par_id);
}

double ModelAnalysisView::GetAtomGausPriorStandardDeviation(
    GroupKey group_key,
    const std::string & class_key,
    int par_id) const
{
    return GetAtomGroupPriorWithUncertainty(group_key, class_key).GetDisplayStandardDeviation(par_id);
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

std::vector<AtomObject *> ModelAnalysisView::GetOutlierAtomObjectList(
    GroupKey group_key,
    const std::string & class_key) const
{
    std::vector<AtomObject *> outlier_atom_list;
    const auto & atom_list{ GetAtomObjectList(group_key, class_key) };
    outlier_atom_list.reserve(atom_list.size());
    for (auto * atom : atom_list)
    {
        const auto annotation{ AtomLocalPotentialView::RequireFor(*atom).FindAnnotation(class_key) };
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
    return AtomLocalPotentialView::RequireFor(*atom_list.front()).GetAlphaR();
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

Residue ModelAnalysisView::GetResidueFromAtomGroupKey(
    GroupKey group_key,
    const std::string & class_key) const
{
    return DecodeResidueFromAtomGroupKey(group_key, class_key);
}

std::vector<GroupKey> ModelAnalysisView::CollectAtomGroupKeys(const std::string & class_key) const
{
    return CollectGroupKeys(ModelAnalysisData::Of(m_model_object).FindAtomGroupEntry(class_key));
}

} // namespace rhbm_gem
