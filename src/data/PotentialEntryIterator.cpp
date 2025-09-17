#include "PotentialEntryIterator.hpp"
#include "AtomObject.hpp"
#include "ModelObject.hpp"
#include "LocalPotentialEntry.hpp"
#include "GroupPotentialEntry.hpp"
#include "AtomicInfoHelper.hpp"
#include "ROOTHelper.hpp"
#include "ArrayStats.hpp"
#include "KeyPacker.hpp"
#include "Constants.hpp"
#include "GlobalEnumClass.hpp"
#include "KDTreeAlgorithm.hpp"
#include "AtomClassifier.hpp"
#include "Logger.hpp"

#include <sstream>
#include <cmath>

#ifdef HAVE_ROOT
#include <TGraphErrors.h>
#include <TGraph2DErrors.h>
#include <TF1.h>
#include <TH1.h>
#endif

PotentialEntryIterator::PotentialEntryIterator(ModelObject * model_object) :
    m_atom_object{ nullptr }, m_model_object{ model_object }, m_local_entry{ nullptr }
{
}

PotentialEntryIterator::PotentialEntryIterator(AtomObject * atom_object) :
    m_atom_object{ atom_object }, m_model_object{ nullptr }
{
    m_local_entry = atom_object->GetLocalPotentialEntry();
}

PotentialEntryIterator::~PotentialEntryIterator()
{

}

double PotentialEntryIterator::GetGausEstimateMinimum(int par_id, Element element) const
{
    if (m_model_object == nullptr)
    {
        throw std::runtime_error("Model object is not available.");
    }
    std::vector<double> gaus_estimate_list;
    gaus_estimate_list.reserve(m_model_object->GetNumberOfSelectedAtom());
    for (auto & atom : m_model_object->GetSelectedAtomList())
    {
        if (atom->GetElement() != element) continue;
        auto entry{ atom->GetLocalPotentialEntry() };
        gaus_estimate_list.emplace_back(entry->GetGausEstimateMDPDE(par_id));
    }
    return ArrayStats<double>::ComputeMin(gaus_estimate_list.data(), gaus_estimate_list.size());
}

bool PotentialEntryIterator::IsOutlierAtom(const std::string & class_key) const
{
    if (IsAtomicEntryAvailable() == false)
    {
        throw std::runtime_error("Atomic entry is not available.");
    }
    return m_local_entry->GetOutlierTag(class_key);
}

bool PotentialEntryIterator::IsAvailableGroupKey(
    GroupKey group_key, const std::string & class_key, bool varbose) const
{
    return CheckGroupKey(group_key, class_key, varbose);
}

size_t PotentialEntryIterator::GetResidueCount(
    const std::string & class_key, Residue residue, Structure structure) const
{
    GroupKey group_key{ 0 };
    AtomClassifier classifier;
    if (class_key == AtomicInfoHelper::GetComponentAtomClassKey())
    {
        group_key = classifier.GetMainChainComponentAtomClassGroupKey(0, residue);
    }
    else if (class_key == AtomicInfoHelper::GetStructureAtomClassKey())
    {
        group_key = classifier.GetMainChainStructureAtomClassGroupKey(0, structure, residue);
    }
    if (IsAvailableGroupKey(group_key, class_key) == false)
    {
        return 0;
    }
    return GetAtomObjectList(group_key, class_key).size();
}

double PotentialEntryIterator::GetGausEstimatePrior(
    GroupKey group_key, const std::string & class_key, int par_id) const
{
    if (CheckGroupKey(group_key, class_key) == false)
    {
        throw std::runtime_error("Group key is not available.");
    }
    return m_model_object->GetAtomGroupPotentialEntry(class_key)->GetGausEstimatePrior(group_key, par_id);
}

double PotentialEntryIterator::GetGausVariancePrior(
    GroupKey group_key, const std::string & class_key, int par_id) const
{
    if (CheckGroupKey(group_key, class_key) == false)
    {
        throw std::runtime_error("Group key is not available.");
    }
    return m_model_object->GetAtomGroupPotentialEntry(class_key)->GetGausVariancePrior(group_key, par_id);
}

const std::vector<AtomObject *> & PotentialEntryIterator::GetAtomObjectList(
    GroupKey group_key, const std::string & class_key) const
{
    if (CheckGroupKey(group_key, class_key) == false)
    {
        throw std::runtime_error("Group key is not available.");
    }
    return m_model_object->GetAtomGroupPotentialEntry(class_key)->GetAtomObjectPtrList(group_key);
}

std::unordered_map<int, AtomObject *> PotentialEntryIterator::GetAtomObjectMap(
    GroupKey group_key, const std::string & class_key) const
{
    std::unordered_map<int, AtomObject *> atom_object_map;
    for (auto & atom_object : GetAtomObjectList(group_key, class_key))
    {
        atom_object_map[atom_object->GetSerialID()] = atom_object;
    }
    return atom_object_map;
}

const std::vector<std::tuple<float, float>> & PotentialEntryIterator::GetDistanceAndMapValueList(void) const
{
    if (IsAtomicEntryAvailable() == false)
    {
        throw std::runtime_error("Atomic entry is not available.");
    }
    return m_local_entry->GetDistanceAndMapValueList();
}

std::vector<std::tuple<float, float>> PotentialEntryIterator::GetBinnedDistanceAndMapValueList(
    int bin_size, double x_min, double x_max) const
{
    if (IsAtomicEntryAvailable() == false)
    {
        throw std::runtime_error("Atomic entry is not available.");
    }
    
    auto bin_spacing{ (x_max - x_min) / static_cast<double>(bin_size) };
    std::unordered_map<int, std::vector<float>> bin_map;
    for (auto & [distance, map_value] : m_local_entry->GetDistanceAndMapValueList())
    {
        auto bin_index{ static_cast<int>(std::floor(distance / bin_spacing)) };
        bin_map[bin_index].emplace_back(map_value);
    }

    std::vector<std::tuple<float, float>> binned_distance_and_map_value_list;
    binned_distance_and_map_value_list.reserve(static_cast<size_t>(bin_size));
    for (int i = 0; i < bin_size; i++)
    {
        auto x_value{ static_cast<float>(x_min + (i + 0.5) * bin_spacing) };
        auto y_value{ ArrayStats<float>::ComputeMedian(bin_map.at(i)) };
        binned_distance_and_map_value_list.emplace_back(std::make_tuple(x_value, y_value));
    }
    return binned_distance_and_map_value_list;
}

std::tuple<float, float> PotentialEntryIterator::GetDistanceRange(double margin_rate) const
{
    if (IsAtomicEntryAvailable() == false)
    {
        throw std::runtime_error("Atomic entry is not available.");
    }
    std::vector<float> distance_array;
    distance_array.reserve(static_cast<size_t>(m_local_entry->GetDistanceAndMapValueListSize()));
    for (auto & [distance, map_value] : m_local_entry->GetDistanceAndMapValueList())
    {
        distance_array.emplace_back(distance);
    }
    return ArrayStats<float>::ComputeScalingRangeTuple(distance_array, static_cast<float>(margin_rate));
}

std::tuple<float, float> PotentialEntryIterator::GetMapValueRange(double margin_rate) const
{
    if (IsAtomicEntryAvailable() == false)
    {
        throw std::runtime_error("Atomic entry is not available.");
    }
    std::vector<float> map_value_array;
    map_value_array.reserve(static_cast<size_t>(m_local_entry->GetDistanceAndMapValueListSize()));
    for (auto & [distance, map_value] : m_local_entry->GetDistanceAndMapValueList())
    {
        map_value_array.emplace_back(map_value);
    }
    return ArrayStats<float>::ComputeScalingRangeTuple(map_value_array, static_cast<float>(margin_rate));
}

double PotentialEntryIterator::GetAmplitudeEstimateMDPDE(void) const
{
    if (IsAtomicEntryAvailable() == false)
    {
        return 0.0;
    }
    return m_local_entry->GetAmplitudeEstimateMDPDE();
}

double PotentialEntryIterator::GetWidthEstimateMDPDE(void) const
{
    if (IsAtomicEntryAvailable() == false)
    {
        return 0.0;
    }
    return m_local_entry->GetWidthEstimateMDPDE();
}

bool PotentialEntryIterator::IsAtomObjectAvailable(void) const
{
    if (m_atom_object == nullptr)
    {
        Logger::Log(LogLevel::Error, "Atom object is not available.");
        return false;
    }
    return true;
}

bool PotentialEntryIterator::IsAtomicEntryAvailable(void) const
{
    if (m_local_entry == nullptr)
    {
        Logger::Log(LogLevel::Error, "Atomic entry is not available.");
        return false;
    }
    return true;
}

bool PotentialEntryIterator::IsModelObjectAvailable(void) const
{
    if (m_model_object == nullptr)
    {
        Logger::Log(LogLevel::Error, "Model object is not available.");
        return false;
    }
    return true;
}

bool PotentialEntryIterator::CheckGroupKey(GroupKey group_key, const std::string & class_key, bool verbose) const
{
    const auto & group_key_set{ m_model_object->GetAtomGroupPotentialEntry(class_key)->GetGroupKeySet() };
    if (group_key_set.find(group_key) == group_key_set.end())
    {
        std::ostringstream oss;
        if (verbose == true)
        {
            if (class_key == AtomicInfoHelper::GetSimpleAtomClassKey())
            {
                oss <<"Atom class group key : "
                    << KeyPackerSimpleAtomClass::GetKeyString(group_key)
                    <<" not found." << std::endl;
            }
            else if (class_key == AtomicInfoHelper::GetComponentAtomClassKey())
            {
                oss <<"Atom class group key : "
                    << KeyPackerComponentAtomClass::GetKeyString(group_key)
                    <<" not found." << std::endl;
            }
            else if (class_key == AtomicInfoHelper::GetStructureAtomClassKey())
            {
                oss <<"Atom class group key : "
                    << KeyPackerStructureAtomClass::GetKeyString(group_key)
                    <<" not found." << std::endl;
            }
            Logger::Log(LogLevel::Error, oss.str());
        }
        return false;
    }
    return true;
}

Residue PotentialEntryIterator::GetResidueFromGroupKey(
    GroupKey group_key, const std::string & class_key) const
{
    if (class_key == AtomicInfoHelper::GetSimpleAtomClassKey())
    {
        Logger::Log(LogLevel::Error, "Element class group key is not recording Residue info.");
        return Residue::UNK;
    }
    else if (class_key == AtomicInfoHelper::GetComponentAtomClassKey())
    {
        auto unpack_key{ KeyPackerComponentAtomClass::Unpack(group_key) };
        return static_cast<Residue>(std::get<0>(unpack_key));
    }
    else if (class_key == AtomicInfoHelper::GetStructureAtomClassKey())
    {
        auto unpack_key{ KeyPackerStructureAtomClass::Unpack(group_key) };
        return static_cast<Residue>(std::get<1>(unpack_key));
    }
    return Residue::UNK;
}

#ifdef HAVE_ROOT
std::unique_ptr<TH1D> PotentialEntryIterator::CreateResidueCountHistogram(
    const std::string & class_key, Structure structure)
{
    if (IsModelObjectAvailable() == false)
    {
        return nullptr;
    }
    auto hist{ ROOTHelper::CreateHist1D("hist", "Residue Count", 20, -0.5, 19.5) };
    for (auto & residue : AtomicInfoHelper::GetStandardResidueList())
    {
        auto count{ GetResidueCount(class_key, residue, structure) };
        hist->SetBinContent(static_cast<int>(residue), static_cast<double>(count));
    }
    return hist;
}

std::vector<std::unique_ptr<TH1D>> PotentialEntryIterator::CreateMainChainRankHistogram(
    int par_id, int & chain_size, Residue residue,
    size_t extra_id, std::vector<Residue> veto_residues_list)
{
    if (IsModelObjectAvailable() == false)
    {
        return {};
    }
    
    std::unordered_map<std::string, std::unordered_map<int, std::array<double, 4>>> values_map;
    for (auto & atom : m_model_object->GetSelectedAtomList())
    {
        size_t id;
        if (AtomClassifier::IsMainChainMember(atom->GetSpot(), id) == false) continue;
        if (residue != Residue::UNK && atom->GetResidue() != residue) continue;
        if (atom->GetSpecialAtomFlag() == true) continue;
        bool is_veto_residue{ false };
        for (auto & veto_residue : veto_residues_list)
        {
            if (atom->GetResidue() == veto_residue) is_veto_residue = true;
        }
        if (is_veto_residue == true) continue;
        auto residue_id{ atom->GetResidueID() };
        auto chain_id{ atom->GetChainID() };
        auto entry{ atom->GetLocalPotentialEntry() };
        auto gaus_value{ entry->GetGausEstimateMDPDE(par_id) };
        values_map[chain_id][residue_id].at(id) = gaus_value;
    }
    chain_size = static_cast<int>(values_map.size());

    std::vector<std::unique_ptr<TH1D>> hist_list;
    for (size_t i = 0; i < 4; i++)
    {
        auto name{ Form("h%d_%d_%d_%d", static_cast<int>(extra_id), static_cast<int>(i), static_cast<int>(residue), par_id) };
        auto hist{ ROOTHelper::CreateHist1D(name, "", 4,  0.5, 4.5) };
        for (auto & [chain_id, values_map_tmp] : values_map)
        {
            for (auto & [residue_id, values] : values_map_tmp)
            {
                hist->Fill(ArrayStats<double>::ComputeRank(values, i));
            }
        }
        hist_list.emplace_back(std::move(hist));
    }
    return hist_list;
}

std::unique_ptr<TGraphErrors> PotentialEntryIterator::CreateAmplitudeRatioToWidthScatterGraph(
    size_t target_id, size_t reference_id, Residue residue)
{
    if (IsModelObjectAvailable() == false)
    {
        return nullptr;
    }
    auto graph{ ROOTHelper::CreateGraphErrors() };
    auto element_size{ AtomClassifier::GetMainChainMemberCount() };
    if (target_id >= element_size || reference_id >= element_size)
    {
        Logger::Log(LogLevel::Error, "Target or reference ID exceeds the number of main chain elements.");
        return nullptr;
    }
    std::unordered_map<int, std::tuple<double, double>> gaus_estimate_map[4];
    size_t current_id{ 0 };
    for (auto & atom : m_model_object->GetSelectedAtomList())
    {
        if (atom->GetSpecialAtomFlag() == true) continue;
        if (AtomClassifier::IsMainChainMember(atom->GetSpot(), current_id) == false) continue;
        if (residue != Residue::UNK && atom->GetResidue() != residue) continue;
        auto entry{ atom->GetLocalPotentialEntry() };
        auto residue_id{ atom->GetResidueID() };
        auto amplitude_estimate{ entry->GetAmplitudeEstimateMDPDE() };
        auto width_estimate{ entry->GetWidthEstimateMDPDE() };
        (gaus_estimate_map[current_id])[residue_id] = std::make_tuple(amplitude_estimate, width_estimate);
    }

    auto count{ 0 };
    for (auto & [residue_id, gaus_estimate] : gaus_estimate_map[target_id])
    {
        if (gaus_estimate_map[reference_id].find(residue_id) == gaus_estimate_map[reference_id].end()) continue;
        auto target_amplitude{ std::get<0>(gaus_estimate) };
        auto target_width{ std::get<1>(gaus_estimate) };
        auto reference_amplitude{ std::get<0>(gaus_estimate_map[reference_id].at(residue_id)) };
        if (reference_amplitude <= 0.0) continue;
        auto ratio{ target_amplitude / reference_amplitude };
        graph->SetPoint(count, target_width, ratio);
        graph->SetPointError(count, 0.0, 0.0);
        count++;
    }
    return graph;
}

std::unique_ptr<TGraphErrors> PotentialEntryIterator::CreateNormalizedGausEstimateScatterGraph(
    Element element, double reference_amplitude, bool reverse)
{
    if (IsModelObjectAvailable() == false)
    {
        return nullptr;
    }
    auto graph{ ROOTHelper::CreateGraphErrors() };
    std::unordered_map<int, double> amplitude_diff_to_carbonyl_oxygen_map;
    for (auto & atom : m_model_object->GetSelectedAtomList())
    {
        if (atom->GetSpot() != Spot::O) continue;
        if (atom->GetSpecialAtomFlag() == false)
        {
            auto entry{ atom->GetLocalPotentialEntry() };
            auto residue_id{ atom->GetResidueID() };
            auto amplitude_estimate{ entry->GetAmplitudeEstimateMDPDE() };
            amplitude_diff_to_carbonyl_oxygen_map[residue_id] = amplitude_estimate - reference_amplitude;
        }
    }
    auto count{ 0 };
    for (auto & atom : m_model_object->GetSelectedAtomList())
    {
        if (atom->GetElement() != element) continue;
        auto residue_id{ atom->GetResidueID() };
        auto entry{ atom->GetLocalPotentialEntry() };
        auto normalized_amplitude{ entry->GetAmplitudeEstimateMDPDE() };
        if (amplitude_diff_to_carbonyl_oxygen_map.find(residue_id) != amplitude_diff_to_carbonyl_oxygen_map.end())
        {
            normalized_amplitude -= amplitude_diff_to_carbonyl_oxygen_map.at(residue_id);
        }
        if (reverse == false)
        {
            graph->SetPoint(count, normalized_amplitude, entry->GetWidthEstimateMDPDE());
        }
        else
        {
            graph->SetPoint(count, entry->GetWidthEstimateMDPDE(), normalized_amplitude);
        }
        count++;
    }
    return graph;
}

std::unique_ptr<TGraphErrors> PotentialEntryIterator::CreateBfactorToWidthScatterGraph(
    GroupKey group_key, const std::string & class_key)
{
    if (IsModelObjectAvailable() == false)
    {
        return nullptr;
    }
    if (IsAvailableGroupKey(group_key, class_key) == false)
    {
        Logger::Log(LogLevel::Error, "Group key is not available.");
        return nullptr;
    }
    auto graph{ ROOTHelper::CreateGraphErrors() };
    const auto & atom_list{ GetAtomObjectList(group_key, class_key) };

    auto count{ 0 };
    for (auto & atom : atom_list)
    {
        auto entry{ atom->GetLocalPotentialEntry() };
        graph->SetPoint(count, atom->GetTemperature()/(8.0*Constants::pi_square), entry->GetWidthEstimateMDPDE());
        count++;
    }

    return graph;
}

std::unordered_map<std::string, std::unique_ptr<TGraphErrors>>
PotentialEntryIterator::CreateGausEstimateToResidueIDGraphMap(
    size_t main_chain_element_id, const int par_id, Residue residue)
{
    if (IsModelObjectAvailable() == false)
    {
        return {};
    }
    
    std::unordered_map<std::string, std::unique_ptr<TGraphErrors>> graph_map;
    std::unordered_map<std::string, int> count_map;
    for (auto & atom : m_model_object->GetSelectedAtomList())
    {
        if (atom->GetElement() != AtomClassifier::GetMainChainElement(main_chain_element_id)) continue;
        if (atom->GetSpot() != AtomClassifier::GetMainChainSpot(main_chain_element_id)) continue;
        if (residue != Residue::UNK && atom->GetResidue() != residue) continue;
        auto entry{ atom->GetLocalPotentialEntry() };
        auto residue_id{ atom->GetResidueID() };
        auto chain_id{ atom->GetChainID() };
        if (residue_id < 0) continue;
        if (graph_map.find(chain_id) == graph_map.end())
        {
            graph_map[chain_id] = ROOTHelper::CreateGraphErrors();
            count_map[chain_id] = 0;
        }
        auto x_value{ static_cast<double>(residue_id) };
        graph_map[chain_id]->SetPoint(count_map[chain_id], x_value, entry->GetGausEstimateMDPDE(par_id));
        count_map[chain_id]++;
    }
    return graph_map;
}

std::unique_ptr<TGraphErrors> PotentialEntryIterator::CreateGausEstimateToResidueGraph(
    std::vector<GroupKey> & group_key_list, const std::string & class_key, const int par_id)
{
    if (IsModelObjectAvailable() == false)
    {
        return nullptr;
    }
    auto graph{ ROOTHelper::CreateGraphErrors() };
    
    auto count{ 0 };
    for (auto & group_key : group_key_list)
    {
        if (IsAvailableGroupKey(group_key, class_key) == false) continue;
        auto x_value{ static_cast<int>(GetResidueFromGroupKey(group_key, class_key)) - 1 };
        auto y_value{ m_model_object->GetAtomGroupPotentialEntry(class_key)->GetGausEstimatePrior(group_key, par_id) };
        auto y_error{ m_model_object->GetAtomGroupPotentialEntry(class_key)->GetGausVariancePrior(group_key, par_id) };
        graph->SetPoint(count, x_value, y_value);
        graph->SetPointError(count, 0.0, y_error);
        count++;
    }
    return graph;
}

std::unique_ptr<TGraphErrors> PotentialEntryIterator::CreateGausEstimateScatterGraph(
    std::vector<GroupKey> & group_key_list, const std::string & class_key,
    int par1_id, int par2_id)
{
    if (IsModelObjectAvailable() == false)
    {
        return nullptr;
    }
    auto graph{ ROOTHelper::CreateGraphErrors() };
    
    auto count{ 0 };
    for (auto & group_key : group_key_list)
    {
        if (IsAvailableGroupKey(group_key, class_key) == false) continue;
        graph->SetPoint(count,
            GetGausEstimatePrior(group_key, class_key, par1_id),
            GetGausEstimatePrior(group_key, class_key, par2_id));
        graph->SetPointError(count,
            GetGausVariancePrior(group_key, class_key, par1_id),
            GetGausVariancePrior(group_key, class_key, par2_id));
        count++;
    }
    return graph;
}

std::unique_ptr<TGraphErrors> PotentialEntryIterator::CreateGausEstimateScatterGraph(
    Element element, bool reverse)
{
    if (IsModelObjectAvailable() == false)
    {
        return nullptr;
    }
    
    auto graph{ ROOTHelper::CreateGraphErrors() };
    auto count{ 0 };
    for (auto & atom : m_model_object->GetSelectedAtomList())
    {
        if (atom->GetElement() != element) continue;
        auto entry{ atom->GetLocalPotentialEntry() };
        if (reverse == false)
        {
            graph->SetPoint(count, entry->GetAmplitudeEstimateMDPDE(), entry->GetWidthEstimateMDPDE());
        }
        else
        {
            graph->SetPoint(count, entry->GetWidthEstimateMDPDE(), entry->GetAmplitudeEstimateMDPDE());
        }
        count++;
    }
    return graph;
}

std::unique_ptr<TGraphErrors> PotentialEntryIterator::CreateGausEstimateScatterGraph(
    GroupKey group_key1, GroupKey group_key2, const std::string & class_key, const int par_id)
{
    if (IsModelObjectAvailable() == false)
    {
        return nullptr;
    }
    if (IsAvailableGroupKey(group_key1, class_key) == false || IsAvailableGroupKey(group_key2, class_key) == false)
    {
        Logger::Log(LogLevel::Error, "Group key is not available.");
        return nullptr;
    }
    auto graph{ ROOTHelper::CreateGraphErrors() };

    const auto & atom_list1{ GetAtomObjectList(group_key1, class_key) };
    const auto & atom_list2{ GetAtomObjectList(group_key2, class_key) };

    auto count{ 0 };
    for (auto atom1 : atom_list1)
    {
        auto entry1{ atom1->GetLocalPotentialEntry() };
        for (auto atom2 : atom_list2)
        {
            auto entry2{ atom2->GetLocalPotentialEntry() };
            if (atom1->GetResidueID() == atom2->GetResidueID() &&
                atom1->GetChainID() == atom2->GetChainID())
            {
                graph->SetPoint(count,
                    entry1->GetGausEstimateMDPDE(par_id),
                    entry2->GetGausEstimateMDPDE(par_id));
                count++;
                break;
            } 
        }
    }
    return graph;
}

std::unique_ptr<TGraphErrors> PotentialEntryIterator::CreateDistanceToMapValueGraph(void)
{
    if (IsAtomicEntryAvailable() == false)
    {
        return nullptr;
    }
    auto graph{ ROOTHelper::CreateGraphErrors() };

    auto count{ 0 };
    for (auto & [distance, map_value] : m_local_entry->GetDistanceAndMapValueList())
    {
        graph->SetPoint(count, distance, map_value);
        count++;
    }
    return graph;
}

std::unique_ptr<TGraphErrors> PotentialEntryIterator::CreateBinnedDistanceToMapValueGraph(
    int bin_size, double x_min, double x_max)
{
    auto data_array{ GetBinnedDistanceAndMapValueList(bin_size, x_min, x_max) };
    auto graph{ ROOTHelper::CreateGraphErrors(bin_size) };
    auto count{ 0 };
    for (auto & [distance, map_value] : data_array)
    {
        graph->SetPoint(count, distance, map_value);
        count++;
    }
    return graph;
}

std::unique_ptr<TGraphErrors> PotentialEntryIterator::CreateInRangeAtomsToGausEstimateGraph(
    GroupKey group_key, const std::string & class_key, double range, int par_id)
{
    if (IsModelObjectAvailable() == false)
    {
        return nullptr;
    }
    auto graph{ ROOTHelper::CreateGraphErrors() };
    auto kd_tree_root{ m_model_object->GetKDTreeRoot() };
    if (kd_tree_root == nullptr)
    {
        Logger::Log(LogLevel::Error, "KDTree is not available for the model object.");
        return nullptr;
    }

    auto count{ 0 };
    for (auto & atom : GetAtomObjectList(group_key, class_key))
    {
        auto in_range_atom_list{ KDTreeAlgorithm<AtomObject>::RangeSearch(kd_tree_root, atom, range) };
        auto atom_entry{ atom->GetLocalPotentialEntry() };
        graph->SetPoint(count,
            static_cast<double>(in_range_atom_list.size()),
            atom_entry->GetGausEstimateMDPDE(par_id));
        count++;
    }
    return graph;
}

std::unique_ptr<TGraphErrors> PotentialEntryIterator::CreateCOMDistanceToGausEstimateGraph(
    GroupKey group_key, const std::string & class_key, int par_id)
{
    if (IsModelObjectAvailable() == false)
    {
        return nullptr;
    }
    auto graph{ ROOTHelper::CreateGraphErrors() };
    auto center_of_mass_pos{ m_model_object->GetCenterOfMassPosition() };

    auto count{ 0 };
    for (auto & atom : GetAtomObjectList(group_key, class_key))
    {
        const auto & atom_pos{ atom->GetPositionRef() };
        auto distance{ ArrayStats<float>::ComputeNorm(atom_pos, center_of_mass_pos) };
        auto atom_entry{ atom->GetLocalPotentialEntry() };
        graph->SetPoint(count,
            static_cast<double>(distance), atom_entry->GetGausEstimateMDPDE(par_id));
        count++;
    }
    return graph;
}

std::unique_ptr<TGraphErrors> PotentialEntryIterator::CreateXYPositionTomographyGraph(
    double normalized_z_pos, double z_ratio_window, bool com_center)
{
    if (IsModelObjectAvailable() == false)
    {
        return nullptr;
    }
    auto com_pos{ com_center ? m_model_object->GetCenterOfMassPosition() : std::array<float, 3>{0.0, 0.0, 0.0} };
    auto graph{ ROOTHelper::CreateGraphErrors() };
    auto z_pos{ m_model_object->GetModelPosition(2, normalized_z_pos) };
    auto window_width{ 0.5 * m_model_object->GetModelLength(2) * z_ratio_window };
    auto z_window_min{ z_pos - window_width };
    auto z_window_max{ z_pos + window_width };
    
    auto count{ 0 };
    for (auto & atom : m_model_object->GetAtomList())
    {
        auto position{ atom->GetPosition() };
        if (position.at(2) < z_window_min || position.at(2) >= z_window_max)
        {
            continue;
        }
        graph->SetPoint(count, position.at(0) - com_pos.at(0), position.at(1) - com_pos.at(1));
        count++;
    }
    return graph;
}

std::unordered_map<size_t, std::unique_ptr<TGraphErrors>>
PotentialEntryIterator::CreateXYPositionTomographyGraphMap(
    double normalized_z_pos, double z_ratio_window, bool com_center)
{
    if (IsModelObjectAvailable() == false)
    {
        return {};
    }

    auto com_pos{ com_center ? m_model_object->GetCenterOfMassPosition() : std::array<float, 3>{0.0, 0.0, 0.0} };
    auto z_pos{ m_model_object->GetModelPosition(2, normalized_z_pos) };
    auto window_width{ 0.5 * m_model_object->GetModelLength(2) * z_ratio_window };
    auto z_window_min{ z_pos - window_width };
    auto z_window_max{ z_pos + window_width };
    
    std::unordered_map<size_t, std::unique_ptr<TGraphErrors>> graph_map;
    std::unordered_map<size_t, int> count_map;
    size_t current_id{ 0 };
    for (auto & atom : m_model_object->GetAtomList())
    {
        if (atom->GetSpecialAtomFlag() == true) continue;
        if (AtomClassifier::IsMainChainMember(atom->GetSpot(), current_id) == false) continue;

        auto position{ atom->GetPosition() };
        if (position.at(2) < z_window_min || position.at(2) >= z_window_max) continue;
        if (graph_map.find(current_id) == graph_map.end())
        {
            graph_map[current_id] = ROOTHelper::CreateGraphErrors();
            count_map[current_id] = 0;
        }
        graph_map[current_id]->SetPoint(
            count_map[current_id],
            position.at(0) - com_pos.at(0),
            position.at(1) - com_pos.at(1));
        count_map[current_id]++;
    }
    return graph_map;
}

std::unique_ptr<TGraphErrors> PotentialEntryIterator::CreateXYPositionTomographyGraph(
    std::vector<GroupKey> & group_key_list, const std::string & class_key, double normalized_z_pos, double z_ratio_window)
{
    if (IsModelObjectAvailable() == false)
    {
        return nullptr;
    }
    auto graph{ ROOTHelper::CreateGraphErrors() };
    auto z_range_tuple{ m_model_object->GetModelPositionRange(2) };
    auto z_pos_min{ std::get<0>(z_range_tuple) };
    auto z_pos_max{ std::get<1>(z_range_tuple) };
    auto z_range{ z_pos_max - z_pos_min };
    auto window_width{ 0.5 * z_range * z_ratio_window };
    auto z_window_min{ z_pos_min + normalized_z_pos * z_range - window_width };
    auto z_window_max{ z_pos_min + normalized_z_pos * z_range + window_width };
    
    auto count{ 0 };
    for (auto & group_key : group_key_list)
    {
        if (IsAvailableGroupKey(group_key, class_key) == false) continue;
        for (auto & atom : GetAtomObjectList(group_key, class_key))
        {
            auto position{ atom->GetPosition() };
            if (position.at(2) < z_window_min || position.at(2) >= z_window_max)
            {
                continue;
            }
            graph->SetPoint(count, position.at(0), position.at(1));
            count++;
        }
    }
    return graph;
}

std::unique_ptr<TGraph2DErrors> PotentialEntryIterator::CreateXYPositionTomographyToGausEstimateGraph2D(
    std::vector<GroupKey> & group_key_list, const std::string & class_key,
    double normalized_z_pos, double z_ratio_window, int par_id)
{
    if (IsModelObjectAvailable() == false)
    {
        return nullptr;
    }
    auto graph{ ROOTHelper::CreateGraph2DErrors() };
    auto z_range_tuple{ m_model_object->GetModelPositionRange(2) };
    auto z_pos_min{ std::get<0>(z_range_tuple) };
    auto z_pos_max{ std::get<1>(z_range_tuple) };
    auto z_range{ z_pos_max - z_pos_min };
    auto window_width{ 0.5 * z_range * z_ratio_window };
    auto z_window_min{ z_pos_min + normalized_z_pos * z_range - window_width };
    auto z_window_max{ z_pos_min + normalized_z_pos * z_range + window_width };
    
    auto count{ 0 };
    for (auto & group_key : group_key_list)
    {
        if (IsAvailableGroupKey(group_key, class_key) == false) continue;
        for (auto & atom : GetAtomObjectList(group_key, class_key))
        {
            auto position{ atom->GetPosition() };
            if (position.at(2) < z_window_min || position.at(2) >= z_window_max)
            {
                continue;
            }
            auto entry_iter{ atom->GetLocalPotentialEntry() };
            graph->SetPoint(count,
                position.at(0), position.at(1),
                entry_iter->GetGausEstimatePosterior(class_key, par_id));
            count++;
        }
    }
    return graph;
}

std::unordered_map<size_t, std::unique_ptr<TGraph2DErrors>>
PotentialEntryIterator::CreateXYPositionTomographyToGausEstimateGraph2DMap(
    double normalized_z_pos, double z_ratio_window, int par_id, bool com_center)
{
    if (IsModelObjectAvailable() == false)
    {
        return {};
    }

    auto com_pos{ com_center ? m_model_object->GetCenterOfMassPosition() : std::array<float, 3>{0.0, 0.0, 0.0} };
    auto z_range_tuple{ m_model_object->GetModelPositionRange(2) };
    auto z_pos_min{ std::get<0>(z_range_tuple) };
    auto z_pos_max{ std::get<1>(z_range_tuple) };
    auto z_range{ z_pos_max - z_pos_min };
    auto window_width{ 0.5 * z_range * z_ratio_window };
    auto z_window_min{ z_pos_min + normalized_z_pos * z_range - window_width };
    auto z_window_max{ z_pos_min + normalized_z_pos * z_range + window_width };
    
    std::unordered_map<size_t, std::unique_ptr<TGraph2DErrors>> graph_map;
    std::unordered_map<size_t, int> count_map;
    size_t current_id{ 0 };
    for (auto & atom : m_model_object->GetSelectedAtomList())
    {
        if (atom->GetSpecialAtomFlag() == true) continue;
        if (AtomClassifier::IsMainChainMember(atom->GetSpot(), current_id) == false) continue;
        auto position{ atom->GetPosition() };
        if (position.at(2) < z_window_min || position.at(2) >= z_window_max) continue;
        if (graph_map.find(current_id) == graph_map.end())
        {
            graph_map[current_id] = ROOTHelper::CreateGraph2DErrors();
            count_map[current_id] = 0;
        }
        auto entry{ atom->GetLocalPotentialEntry() };
        graph_map[current_id]->SetPoint(
            count_map[current_id],
            position.at(0) - com_pos.at(0),
            position.at(1) - com_pos.at(1),
            entry->GetGausEstimateMDPDE(par_id));
        count_map[current_id]++;
    }
    return graph_map;
}

std::unique_ptr<TF1> PotentialEntryIterator::CreateGroupGausFunctionPrior(
    GroupKey group_key, const std::string & class_key) const
{
    auto amplitude
    {
        m_model_object->GetAtomGroupPotentialEntry(class_key)->GetGausEstimatePrior(group_key, 0)
    };
    auto width
    {
        m_model_object->GetAtomGroupPotentialEntry(class_key)->GetGausEstimatePrior(group_key, 1)
    };
    return ROOTHelper::CreateGausFunction1D("gaus", amplitude, width);
}
#endif
