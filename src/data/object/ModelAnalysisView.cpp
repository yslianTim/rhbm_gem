#include <rhbm_gem/data/object/ModelAnalysisView.hpp>

#include "data/detail/ModelAnalysisData.hpp"
#include "data/detail/LocalPotentialEntry.hpp"

#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/domain/ChemicalDataHelper.hpp>
#include <rhbm_gem/utils/domain/KeyPacker.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/math/ArrayHelper.hpp>

#include <cmath>
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
    return RequireEntry("Local Gaussian result").GaussianResult();
}

const GaussianModel3D & AtomLocalPotentialView::GetEstimateOLS() const
{
    return RequireEntry("Local estimate OLS").GaussianResult().ols.GetModel();
}

const GaussianModel3D & AtomLocalPotentialView::GetEstimateMDPDE() const
{
    return RequireEntry("Local estimate MDPDE").GaussianResult().mdpde.GetModel();
}

const LocalPotentialSampleList & AtomLocalPotentialView::GetSamplingEntries() const
{
    return RequireEntry("Local sampling entries").SamplingEntries();
}

std::tuple<float, float> AtomLocalPotentialView::GetDistanceRange(double margin_rate) const
{
    const auto & sampling_entries{ RequireEntry("Local distance range").SamplingEntries() };
    std::vector<float> distance_array;
    distance_array.reserve(sampling_entries.size());
    for (const auto & sample : sampling_entries)
    {
        distance_array.emplace_back(sample.distance);
    }
    return array_helper::ComputeScalingRangeTuple(
        distance_array, static_cast<float>(margin_rate));
}

std::tuple<float, float> AtomLocalPotentialView::GetResponseRange(double margin_rate) const
{
    const auto & sampling_entries{ RequireEntry("Local response range").SamplingEntries() };
    std::vector<float> map_value_array;
    map_value_array.reserve(sampling_entries.size());
    for (const auto & sample : sampling_entries)
    {
        map_value_array.emplace_back(sample.response);
    }
    return array_helper::ComputeScalingRangeTuple(
        map_value_array, static_cast<float>(margin_rate));
}

SeriesPointList AtomLocalPotentialView::GetBinnedDistanceResponseSeries(
    int bin_size,
    double x_min,
    double x_max) const
{
    const auto & sampling_entries{
        RequireEntry("Local binned distance-response series").SamplingEntries() };
    const auto bin_spacing{ (x_max - x_min) / static_cast<double>(bin_size) };
    std::vector<std::vector<float>> bin_map(static_cast<size_t>(bin_size));
    for (const auto & sample : sampling_entries)
    {
        const auto distance{ sample.distance };
        if (distance < x_min || distance >= x_max)
        {
            continue;
        }
        const auto shifted_distance{ static_cast<double>(distance) - x_min };
        const auto bin_index{ static_cast<int>(std::floor(shifted_distance / bin_spacing)) };
        if (bin_index < 0 || bin_index >= bin_size)
        {
            continue;
        }
        bin_map.at(static_cast<size_t>(bin_index)).emplace_back(sample.response);
    }

    SeriesPointList binned_distance_response_series;
    binned_distance_response_series.reserve(static_cast<size_t>(bin_size));
    for (int i = 0; i < bin_size; i++)
    {
        const auto x_value{ static_cast<float>(x_min + (i + 0.5) * bin_spacing) };
        const auto & bin_values{ bin_map.at(static_cast<size_t>(i)) };
        const auto y_value{ bin_values.empty() ? 0.0f : array_helper::ComputeMedian(bin_values) };
        binned_distance_response_series.emplace_back(SeriesPoint{ { x_value }, y_value });
    }
    return binned_distance_response_series;
}

double AtomLocalPotentialView::GetAlphaR() const
{
    return RequireEntry("Local alpha-r").GaussianResult().alpha_r;
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
    const auto & data_list{ RequireEntry("Local center map value").SamplingEntries() };
    if (data_list.empty()) return 0.0;
    int count{ 0 };
    double sum{ 0.0 };
    for (const auto & sample : data_list)
    {
        if (sample.distance > 0.05f) continue;
        sum += sample.response;
        count++;
    }
    if (count == 0) return 0.0;
    return sum/static_cast<double>(count);
}

double AtomLocalPotentialView::CalculateQScore(int par_choice) const
{
    const auto & entry{ RequireEntry("Local q-score") };
    const auto & sampling_entries{ entry.SamplingEntries() };
    if (sampling_entries.empty())
    {
        return 0.0;
    }

    auto q_score{ 0.0 };
    auto amplitude{ 0.0 };
    auto width{ 0.0 };
    auto intersect{ 0.0 };
    if (par_choice == 0)
    {
        amplitude = 0.05; // TODO
        width = 0.6; // TODO
        intersect = -0.005; // TODO
    }
    else if (par_choice == 1)
    {
        const auto & estimate{ entry.GaussianResult().mdpde.GetModel() };
        amplitude = estimate.Intensity();
        width = estimate.GetWidth();
        intersect = 0.0;
    }

    if (std::isfinite(amplitude) == false || std::isfinite(width) == false ||
        std::isfinite(intersect) == false || width <= 0.0)
    {
        return 0.0;
    }

    std::vector<float> map_value_list;
    std::vector<float> reference_value_list;
    map_value_list.reserve(sampling_entries.size());
    reference_value_list.reserve(sampling_entries.size());
    for (const auto & sample : sampling_entries)
    {
        map_value_list.emplace_back(sample.response);
        reference_value_list.emplace_back(
            amplitude * std::exp(-0.5 * std::pow(sample.distance / width, 2)) + intersect);
    }
    if (map_value_list.empty())
    {
        return 0.0;
    }

    float map_value_mean{ array_helper::ComputeMean(map_value_list.data(), map_value_list.size()) };
    float reference_value_mean{ array_helper::ComputeMean(reference_value_list.data(), reference_value_list.size()) };

    auto numerator{ 0.0 };
    auto denominator{ 0.0 };
    auto map_value_norm_square{ 0.0 };
    auto reference_value_norm_square{ 0.0 };
    for (size_t i = 0; i < map_value_list.size(); i++)
    {
        auto map_value_diff{ map_value_list.at(i) - map_value_mean };
        auto reference_value_diff{ reference_value_list.at(i) - reference_value_mean };
        numerator += map_value_diff * reference_value_diff;
        map_value_norm_square += map_value_diff * map_value_diff;
        reference_value_norm_square += reference_value_diff * reference_value_diff;
    }
    denominator = std::sqrt(map_value_norm_square) * std::sqrt(reference_value_norm_square);
    if (denominator == 0.0 || std::isfinite(denominator) == false || std::isfinite(numerator) == false)
    {
        return 0.0;
    }
    q_score = numerator/denominator;

    if (std::isfinite(q_score) == false)
    {
        return 0.0;
    }

    return q_score;
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
