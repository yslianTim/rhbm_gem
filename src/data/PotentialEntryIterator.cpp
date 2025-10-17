#include "PotentialEntryIterator.hpp"
#include "AtomObject.hpp"
#include "BondObject.hpp"
#include "ModelObject.hpp"
#include "LocalPotentialEntry.hpp"
#include "GroupPotentialEntry.hpp"
#include "ChemicalDataHelper.hpp"
#include "ROOTHelper.hpp"
#include "ArrayStats.hpp"
#include "KeyPacker.hpp"
#include "Constants.hpp"
#include "GlobalEnumClass.hpp"
#include "KDTreeAlgorithm.hpp"
#include "AtomClassifier.hpp"
#include "BondClassifier.hpp"
#include "Logger.hpp"

#include <sstream>
#include <cmath>
#include <map>

#ifdef HAVE_ROOT
#include <TGraphErrors.h>
#include <TGraph2DErrors.h>
#include <TF1.h>
#include <TH1.h>
#endif

PotentialEntryIterator::PotentialEntryIterator(ModelObject * model_object) :
    m_atom_object{ nullptr }, m_bond_object{ nullptr }, m_model_object{ model_object },
    m_atom_local_entry{ nullptr }, m_bond_local_entry{ nullptr }
{
}

PotentialEntryIterator::PotentialEntryIterator(AtomObject * atom_object) :
    m_atom_object{ atom_object }, m_bond_object{ nullptr }, m_model_object{ nullptr },
    m_bond_local_entry{ nullptr }
{
    m_atom_local_entry = atom_object->GetLocalPotentialEntry();
}

PotentialEntryIterator::PotentialEntryIterator(BondObject * bond_object) :
    m_atom_object{ nullptr }, m_bond_object{ bond_object }, m_model_object{ nullptr },
    m_atom_local_entry{ nullptr }
{
    m_bond_local_entry = bond_object->GetLocalPotentialEntry();
}

PotentialEntryIterator::~PotentialEntryIterator()
{

}

double PotentialEntryIterator::GetAtomGausEstimateMinimum(int par_id, Element element) const
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

double PotentialEntryIterator::GetBondGausEstimateMinimum(int par_id) const
{
    if (m_model_object == nullptr)
    {
        throw std::runtime_error("Model object is not available.");
    }
    std::vector<double> gaus_estimate_list;
    gaus_estimate_list.reserve(m_model_object->GetNumberOfSelectedBond());
    for (auto & bond : m_model_object->GetSelectedBondList())
    {
        auto entry{ bond->GetLocalPotentialEntry() };
        gaus_estimate_list.emplace_back(entry->GetGausEstimateMDPDE(par_id));
    }
    return ArrayStats<double>::ComputeMin(gaus_estimate_list.data(), gaus_estimate_list.size());
}

bool PotentialEntryIterator::IsOutlierAtom(const std::string & class_key) const
{
    if (IsAtomLocalEntryAvailable() == false)
    {
        throw std::runtime_error("Atom local entry is not available.");
    }
    return m_atom_local_entry->GetOutlierTag(class_key);
}

bool PotentialEntryIterator::IsOutlierBond(const std::string & class_key) const
{
    if (IsBondLocalEntryAvailable() == false)
    {
        throw std::runtime_error("Bond local entry is not available.");
    }
    return m_bond_local_entry->GetOutlierTag(class_key);
}

bool PotentialEntryIterator::IsAvailableAtomGroupKey(
    GroupKey group_key, const std::string & class_key, bool varbose) const
{
    return CheckAtomGroupKey(group_key, class_key, varbose);
}

bool PotentialEntryIterator::IsAvailableBondGroupKey(
    GroupKey group_key, const std::string & class_key, bool varbose) const
{
    return CheckBondGroupKey(group_key, class_key, varbose);
}

size_t PotentialEntryIterator::GetAtomResidueCount(
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
    if (IsAvailableAtomGroupKey(group_key, class_key) == false)
    {
        return 0;
    }
    return GetAtomObjectList(group_key, class_key).size();
}

size_t PotentialEntryIterator::GetBondResidueCount(
    const std::string & class_key, Residue residue) const
{
    GroupKey group_key{ 0 };
    BondClassifier classifier;
    if (class_key == ChemicalDataHelper::GetComponentBondClassKey())
    {
        group_key = classifier.GetMainChainComponentBondClassGroupKey(0, residue);
    }
    if (IsAvailableBondGroupKey(group_key, class_key) == false)
    {
        return 0;
    }
    return GetBondObjectList(group_key, class_key).size();
}

double PotentialEntryIterator::GetAtomGausEstimatePrior(
    GroupKey group_key, const std::string & class_key, int par_id) const
{
    if (CheckAtomGroupKey(group_key, class_key) == false)
    {
        throw std::runtime_error("Group key is not available.");
    }
    return m_model_object->GetAtomGroupPotentialEntry(class_key)->GetGausEstimatePrior(group_key, par_id);
}

double PotentialEntryIterator::GetBondGausEstimatePrior(
    GroupKey group_key, const std::string & class_key, int par_id) const
{
    if (CheckBondGroupKey(group_key, class_key) == false)
    {
        throw std::runtime_error("Group key is not available.");
    }
    return m_model_object->GetBondGroupPotentialEntry(class_key)->GetGausEstimatePrior(group_key, par_id);
}

double PotentialEntryIterator::GetAtomGausVariancePrior(
    GroupKey group_key, const std::string & class_key, int par_id) const
{
    if (CheckAtomGroupKey(group_key, class_key) == false)
    {
        throw std::runtime_error("Group key is not available.");
    }
    return m_model_object->GetAtomGroupPotentialEntry(class_key)->GetGausVariancePrior(group_key, par_id);
}

double PotentialEntryIterator::GetBondGausVariancePrior(
    GroupKey group_key, const std::string & class_key, int par_id) const
{
    if (CheckBondGroupKey(group_key, class_key) == false)
    {
        throw std::runtime_error("Group key is not available.");
    }
    return m_model_object->GetBondGroupPotentialEntry(class_key)->GetGausVariancePrior(group_key, par_id);
}

const std::vector<AtomObject *> & PotentialEntryIterator::GetAtomObjectList(
    GroupKey group_key, const std::string & class_key) const
{
    if (CheckAtomGroupKey(group_key, class_key) == false)
    {
        throw std::runtime_error("Atom group key is not available.");
    }
    return m_model_object->GetAtomGroupPotentialEntry(class_key)->GetAtomObjectPtrList(group_key);
}

const std::vector<BondObject *> & PotentialEntryIterator::GetBondObjectList(
    GroupKey group_key, const std::string & class_key) const
{
    if (CheckBondGroupKey(group_key, class_key) == false)
    {
        throw std::runtime_error("Bond group key is not available.");
    }
    return m_model_object->GetBondGroupPotentialEntry(class_key)->GetBondObjectPtrList(group_key);
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
    if (IsAtomLocalEntryAvailable() == true)
    {
        return m_atom_local_entry->GetDistanceAndMapValueList();
    }
    else if (IsBondLocalEntryAvailable() == true)
    {
        return m_bond_local_entry->GetDistanceAndMapValueList();
    }
    else
    {
        throw std::runtime_error("Local entry is not available.");
    }
}

std::vector<std::tuple<float, float>> PotentialEntryIterator::GetBinnedDistanceAndMapValueList(
    int bin_size, double x_min, double x_max) const
{
    auto bin_spacing{ (x_max - x_min) / static_cast<double>(bin_size) };
    std::map<int, std::vector<float>> bin_map;

    if (IsAtomLocalEntryAvailable() == true)
    {
        for (auto & [distance, map_value] : m_atom_local_entry->GetDistanceAndMapValueList())
        {
            auto bin_index{ static_cast<int>(std::floor(distance / bin_spacing)) };
            bin_map[bin_index].emplace_back(map_value);
        }
    }
    else if (IsBondLocalEntryAvailable() == true)
    {
        for (auto & [distance, map_value] : m_bond_local_entry->GetDistanceAndMapValueList())
        {
            auto bin_index{ static_cast<int>(std::floor(distance / bin_spacing)) };
            bin_map[bin_index].emplace_back(map_value);
        }
    }
    else
    {
        throw std::runtime_error("Local entry is not available.");
    }

    std::vector<std::tuple<float, float>> binned_distance_and_map_value_list;
    binned_distance_and_map_value_list.reserve(static_cast<size_t>(bin_size));
    for (int i = 0; i < bin_size; i++)
    {
        auto x_value{ static_cast<float>(x_min + (i + 0.5) * bin_spacing) };
        auto y_value{ (bin_map.find(i) == bin_map.end()) ?
            0.0f : ArrayStats<float>::ComputeMedian(bin_map.at(i))
        };
        binned_distance_and_map_value_list.emplace_back(std::make_tuple(x_value, y_value));
    }
    return binned_distance_and_map_value_list;
}

std::tuple<float, float> PotentialEntryIterator::GetDistanceRange(double margin_rate) const
{
    std::vector<float> distance_array;
    if (IsAtomLocalEntryAvailable() == true)
    {
        distance_array.reserve(static_cast<size_t>(m_atom_local_entry->GetDistanceAndMapValueListSize()));
        for (auto & [distance, map_value] : m_atom_local_entry->GetDistanceAndMapValueList())
        {
            distance_array.emplace_back(distance);
        }
    }
    else if (IsBondLocalEntryAvailable() == true)
    {
        distance_array.reserve(static_cast<size_t>(m_bond_local_entry->GetDistanceAndMapValueListSize()));
        for (auto & [distance, map_value] : m_bond_local_entry->GetDistanceAndMapValueList())
        {
            distance_array.emplace_back(distance);
        }
    }
    else
    {
        throw std::runtime_error("Local entry is not available.");
    }
    return ArrayStats<float>::ComputeScalingRangeTuple(distance_array, static_cast<float>(margin_rate));
}

std::tuple<float, float> PotentialEntryIterator::GetMapValueRange(double margin_rate) const
{
    std::vector<float> map_value_array;
    if (IsAtomLocalEntryAvailable() == true)
    {
        map_value_array.reserve(static_cast<size_t>(m_atom_local_entry->GetDistanceAndMapValueListSize()));
        for (auto & [distance, map_value] : m_atom_local_entry->GetDistanceAndMapValueList())
        {
            map_value_array.emplace_back(map_value);
        }
    }
    else if (IsBondLocalEntryAvailable() == true)
    {
        map_value_array.reserve(static_cast<size_t>(m_bond_local_entry->GetDistanceAndMapValueListSize()));
        for (auto & [distance, map_value] : m_bond_local_entry->GetDistanceAndMapValueList())
        {
            map_value_array.emplace_back(map_value);
        }
    }
    else
    {
        throw std::runtime_error("Local entry is not available.");
    }
    return ArrayStats<float>::ComputeScalingRangeTuple(map_value_array, static_cast<float>(margin_rate));
}

double PotentialEntryIterator::GetAmplitudeEstimateMDPDE(void) const
{
    if (IsAtomLocalEntryAvailable() == true)
    {
        return m_atom_local_entry->GetAmplitudeEstimateMDPDE();
    }
    else if (IsBondLocalEntryAvailable() == true)
    {
        return m_bond_local_entry->GetAmplitudeEstimateMDPDE();
    }
    return 0.0;
}

double PotentialEntryIterator::GetWidthEstimateMDPDE(void) const
{
    if (IsAtomLocalEntryAvailable() == true)
    {
        return m_atom_local_entry->GetWidthEstimateMDPDE();
    }
    else if (IsBondLocalEntryAvailable() == true)
    {
        return m_bond_local_entry->GetWidthEstimateMDPDE();
    }
    return 0.0;
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

bool PotentialEntryIterator::IsBondObjectAvailable(void) const
{
    if (m_bond_object == nullptr)
    {
        Logger::Log(LogLevel::Error, "Bond object is not available.");
        return false;
    }
    return true;
}

bool PotentialEntryIterator::IsAtomLocalEntryAvailable(void) const
{
    if (m_atom_local_entry == nullptr)
    {
        return false;
    }
    return true;
}

bool PotentialEntryIterator::IsBondLocalEntryAvailable(void) const
{
    if (m_bond_local_entry == nullptr)
    {
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

bool PotentialEntryIterator::CheckAtomGroupKey(
    GroupKey group_key, const std::string & class_key, bool verbose) const
{
    const auto & group_key_set{ m_model_object->GetAtomGroupPotentialEntry(class_key)->GetGroupKeySet() };
    if (group_key_set.find(group_key) == group_key_set.end())
    {
        std::ostringstream oss;
        if (verbose == true)
        {
            if (class_key == ChemicalDataHelper::GetSimpleAtomClassKey())
            {
                oss <<"Atom class group key : "
                    << KeyPackerSimpleAtomClass::GetKeyString(group_key)
                    <<" not found." << std::endl;
            }
            else if (class_key == ChemicalDataHelper::GetComponentAtomClassKey())
            {
                oss <<"Atom class group key : "
                    << KeyPackerComponentAtomClass::GetKeyString(group_key)
                    <<" not found." << std::endl;
            }
            else if (class_key == ChemicalDataHelper::GetStructureAtomClassKey())
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

bool PotentialEntryIterator::CheckBondGroupKey(
    GroupKey group_key, const std::string & class_key, bool verbose) const
{
    const auto & group_key_set{ m_model_object->GetBondGroupPotentialEntry(class_key)->GetGroupKeySet() };
    if (group_key_set.find(group_key) == group_key_set.end())
    {
        std::ostringstream oss;
        if (verbose == true)
        {
            if (class_key == ChemicalDataHelper::GetSimpleBondClassKey())
            {
                oss <<"Bond class group key : "
                    << KeyPackerSimpleBondClass::GetKeyString(group_key)
                    <<" not found." << std::endl;
            }
            else if (class_key == ChemicalDataHelper::GetComponentBondClassKey())
            {
                oss <<"Bond class group key : "
                    << KeyPackerComponentBondClass::GetKeyString(group_key)
                    <<" not found." << std::endl;
            }
            Logger::Log(LogLevel::Error, oss.str());
        }
        return false;
    }
    return true;
}

Residue PotentialEntryIterator::GetResidueFromAtomGroupKey(
    GroupKey group_key, const std::string & class_key) const
{
    if (class_key == ChemicalDataHelper::GetSimpleAtomClassKey())
    {
        Logger::Log(LogLevel::Error, "Simple class group key is not recording Residue info.");
        return Residue::UNK;
    }
    else if (class_key == ChemicalDataHelper::GetComponentAtomClassKey())
    {
        auto unpack_key{ KeyPackerComponentAtomClass::Unpack(group_key) };
        return static_cast<Residue>(std::get<0>(unpack_key));
    }
    else if (class_key == ChemicalDataHelper::GetStructureAtomClassKey())
    {
        auto unpack_key{ KeyPackerStructureAtomClass::Unpack(group_key) };
        return static_cast<Residue>(std::get<1>(unpack_key));
    }
    return Residue::UNK;
}

Residue PotentialEntryIterator::GetResidueFromBondGroupKey(
    GroupKey group_key, const std::string & class_key) const
{
    if (class_key == ChemicalDataHelper::GetSimpleBondClassKey())
    {
        Logger::Log(LogLevel::Error, "Simple class group key is not recording Residue info.");
        return Residue::UNK;
    }
    else if (class_key == ChemicalDataHelper::GetComponentBondClassKey())
    {
        auto unpack_key{ KeyPackerComponentBondClass::Unpack(group_key) };
        return static_cast<Residue>(std::get<0>(unpack_key));
    }
    return Residue::UNK;
}

#ifdef HAVE_ROOT
std::unique_ptr<TH1D> PotentialEntryIterator::CreateAtomResidueCountHistogram(
    const std::string & class_key, Structure structure)
{
    if (IsModelObjectAvailable() == false)
    {
        return nullptr;
    }
    auto hist{ ROOTHelper::CreateHist1D("hist", "Residue Count", 20, -0.5, 19.5) };
    for (auto & residue : ChemicalDataHelper::GetStandardAminoAcidList())
    {
        auto count{ GetAtomResidueCount(class_key, residue, structure) };
        hist->SetBinContent(static_cast<int>(residue), static_cast<double>(count));
    }
    return hist;
}

std::unique_ptr<TH1D> PotentialEntryIterator::CreateBondResidueCountHistogram(
    const std::string & class_key)
{
    if (IsModelObjectAvailable() == false)
    {
        return nullptr;
    }
    auto hist{ ROOTHelper::CreateHist1D("hist", "Residue Count", 20, -0.5, 19.5) };
    for (auto & residue : ChemicalDataHelper::GetStandardAminoAcidList())
    {
        auto count{ GetBondResidueCount(class_key, residue) };
        hist->SetBinContent(static_cast<int>(residue), static_cast<double>(count));
    }
    return hist;
}

std::vector<std::unique_ptr<TH1D>> PotentialEntryIterator::CreateMainChainAtomGausRankHistogram(
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
        auto sequence_id{ atom->GetSequenceID() };
        auto chain_id{ atom->GetChainID() };
        auto entry{ atom->GetLocalPotentialEntry() };
        auto gaus_value{ entry->GetGausEstimateMDPDE(par_id) };
        values_map[chain_id][sequence_id].at(id) = gaus_value;
    }
    chain_size = static_cast<int>(values_map.size());

    std::vector<std::unique_ptr<TH1D>> hist_list;
    for (size_t i = 0; i < 4; i++)
    {
        auto name{ Form("h%d_%d_%d_%d", static_cast<int>(extra_id), static_cast<int>(i), static_cast<int>(residue), par_id) };
        auto hist{ ROOTHelper::CreateHist1D(name, "", 4,  0.5, 4.5) };
        for (auto & [chain_id, values_map_tmp] : values_map)
        {
            for (auto & [sequence_id, values] : values_map_tmp)
            {
                hist->Fill(ArrayStats<double>::ComputeRank(values, i));
            }
        }
        hist_list.emplace_back(std::move(hist));
    }
    return hist_list;
}

std::unique_ptr<TGraphErrors> PotentialEntryIterator::CreateNormalizedAtomGausEstimateScatterGraph(
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
            auto sequence_id{ atom->GetSequenceID() };
            auto amplitude_estimate{ entry->GetAmplitudeEstimateMDPDE() };
            amplitude_diff_to_carbonyl_oxygen_map[sequence_id] = amplitude_estimate - reference_amplitude;
        }
    }
    auto count{ 0 };
    for (auto & atom : m_model_object->GetSelectedAtomList())
    {
        if (atom->GetElement() != element) continue;
        auto sequence_id{ atom->GetSequenceID() };
        auto entry{ atom->GetLocalPotentialEntry() };
        auto normalized_amplitude{ entry->GetAmplitudeEstimateMDPDE() };
        if (amplitude_diff_to_carbonyl_oxygen_map.find(sequence_id) != amplitude_diff_to_carbonyl_oxygen_map.end())
        {
            normalized_amplitude -= amplitude_diff_to_carbonyl_oxygen_map.at(sequence_id);
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

std::unique_ptr<TGraphErrors> PotentialEntryIterator::CreateNormalizedBondGausEstimateScatterGraph(
    Element element, double reference_amplitude, bool reverse)
{
    if (IsModelObjectAvailable() == false)
    {
        return nullptr;
    }
    auto graph{ ROOTHelper::CreateGraphErrors() };
    std::unordered_map<int, double> amplitude_diff_to_carbonyl_oxygen_map;
    for (auto & bond : m_model_object->GetSelectedBondList())
    {
        if (bond->GetSpecialBondFlag() == false)
        {
            auto entry{ bond->GetLocalPotentialEntry() };
            auto sequence_id{ bond->GetAtomObject1()->GetSequenceID() };
            auto amplitude_estimate{ entry->GetAmplitudeEstimateMDPDE() };
            amplitude_diff_to_carbonyl_oxygen_map[sequence_id] = amplitude_estimate - reference_amplitude;
        }
    }
    auto count{ 0 };
    for (auto & bond : m_model_object->GetSelectedBondList())
    {
        if (bond->GetAtomObject1()->GetElement() != element) continue;
        auto sequence_id{ bond->GetAtomObject1()->GetSequenceID() };
        auto entry{ bond->GetLocalPotentialEntry() };
        auto normalized_amplitude{ entry->GetAmplitudeEstimateMDPDE() };
        if (amplitude_diff_to_carbonyl_oxygen_map.find(sequence_id) != amplitude_diff_to_carbonyl_oxygen_map.end())
        {
            normalized_amplitude -= amplitude_diff_to_carbonyl_oxygen_map.at(sequence_id);
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

std::unordered_map<std::string, std::unique_ptr<TGraphErrors>>
PotentialEntryIterator::CreateAtomGausEstimateToSequenceIDGraphMap(
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
        auto sequence_id{ atom->GetSequenceID() };
        auto chain_id{ atom->GetChainID() };
        if (sequence_id < 0) continue;
        if (graph_map.find(chain_id) == graph_map.end())
        {
            graph_map[chain_id] = ROOTHelper::CreateGraphErrors();
            count_map[chain_id] = 0;
        }
        auto x_value{ static_cast<double>(sequence_id) };
        graph_map[chain_id]->SetPoint(count_map[chain_id], x_value, entry->GetGausEstimateMDPDE(par_id));
        count_map[chain_id]++;
    }
    return graph_map;
}

std::unordered_map<std::string, std::unique_ptr<TGraphErrors>>
PotentialEntryIterator::CreateBondGausEstimateToSequenceIDGraphMap(
    size_t main_chain_element_id, const int par_id, Residue residue)
{
    if (IsModelObjectAvailable() == false)
    {
        return {};
    }
    
    std::unordered_map<std::string, std::unique_ptr<TGraphErrors>> graph_map;
    std::unordered_map<std::string, int> count_map;
    auto class_key{ ChemicalDataHelper::GetSimpleBondClassKey() };
    auto group_key{ BondClassifier::GetMainChainSimpleBondClassGroupKey(main_chain_element_id) };
    for (auto & bond : GetBondObjectList(group_key, class_key))
    {
        if (residue != Residue::UNK && bond->GetAtomObject1()->GetResidue() != residue) continue;
        auto entry{ bond->GetLocalPotentialEntry() };
        auto sequence_id{ bond->GetAtomObject1()->GetSequenceID() };
        auto chain_id{ bond->GetAtomObject1()->GetChainID() };
        if (sequence_id < 0) continue;
        if (graph_map.find(chain_id) == graph_map.end())
        {
            graph_map[chain_id] = ROOTHelper::CreateGraphErrors();
            count_map[chain_id] = 0;
        }
        auto x_value{ static_cast<double>(sequence_id) };
        auto y_value{ entry->GetGausEstimateMDPDE(par_id) };
        graph_map[chain_id]->SetPoint(count_map[chain_id], x_value, y_value);
        count_map[chain_id]++;
    }

    for (auto & [chain, graph] : graph_map)
    {
        graph->Sort(); // sort points along x axis
    }
    return graph_map;
}

std::unique_ptr<TGraphErrors> PotentialEntryIterator::CreateAtomGausEstimateToResidueGraph(
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
        if (IsAvailableAtomGroupKey(group_key, class_key) == false) continue;
        auto x_value{ static_cast<int>(GetResidueFromAtomGroupKey(group_key, class_key)) - 1 };
        auto y_value{ m_model_object->GetAtomGroupPotentialEntry(class_key)->GetGausEstimatePrior(group_key, par_id) };
        auto y_error{ m_model_object->GetAtomGroupPotentialEntry(class_key)->GetGausVariancePrior(group_key, par_id) };
        graph->SetPoint(count, x_value, y_value);
        graph->SetPointError(count, 0.0, y_error);
        count++;
    }
    return graph;
}

std::unique_ptr<TGraphErrors> PotentialEntryIterator::CreateBondGausEstimateToResidueGraph(
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
        if (IsAvailableBondGroupKey(group_key, class_key) == false) continue;
        auto x_value{ static_cast<int>(GetResidueFromBondGroupKey(group_key, class_key)) - 1 };
        auto y_value{ m_model_object->GetBondGroupPotentialEntry(class_key)->GetGausEstimatePrior(group_key, par_id) };
        auto y_error{ m_model_object->GetBondGroupPotentialEntry(class_key)->GetGausVariancePrior(group_key, par_id) };
        graph->SetPoint(count, x_value, y_value);
        graph->SetPointError(count, 0.0, y_error);
        count++;
    }
    return graph;
}

std::unique_ptr<TGraphErrors> PotentialEntryIterator::CreateAtomGausEstimateScatterGraph(
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
        if (IsAvailableAtomGroupKey(group_key, class_key) == false) continue;
        graph->SetPoint(count,
            GetAtomGausEstimatePrior(group_key, class_key, par1_id),
            GetAtomGausEstimatePrior(group_key, class_key, par2_id));
        graph->SetPointError(count,
            GetAtomGausVariancePrior(group_key, class_key, par1_id),
            GetAtomGausVariancePrior(group_key, class_key, par2_id));
        count++;
    }
    return graph;
}

std::unique_ptr<TGraphErrors> PotentialEntryIterator::CreateBondGausEstimateScatterGraph(
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
        if (IsAvailableBondGroupKey(group_key, class_key) == false) continue;
        graph->SetPoint(count,
            GetBondGausEstimatePrior(group_key, class_key, par1_id),
            GetBondGausEstimatePrior(group_key, class_key, par2_id));
        graph->SetPointError(count,
            GetBondGausVariancePrior(group_key, class_key, par1_id),
            GetBondGausVariancePrior(group_key, class_key, par2_id));
        count++;
    }
    return graph;
}

std::unique_ptr<TGraphErrors> PotentialEntryIterator::CreateAtomGausEstimateScatterGraph(
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

std::unique_ptr<TGraphErrors> PotentialEntryIterator::CreateBondGausEstimateScatterGraph(
    Element element, bool reverse)
{
    if (IsModelObjectAvailable() == false)
    {
        return nullptr;
    }

    auto graph{ ROOTHelper::CreateGraphErrors() };
    auto count{ 0 };
    for (auto & bond : m_model_object->GetSelectedBondList())
    {
        if (bond->GetAtomObject1()->GetElement() != element) continue;
        auto entry{ bond->GetLocalPotentialEntry() };
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
    if (IsAvailableAtomGroupKey(group_key1, class_key) == false || IsAvailableAtomGroupKey(group_key2, class_key) == false)
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
            if (atom1->GetSequenceID() == atom2->GetSequenceID() &&
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
    auto graph{ ROOTHelper::CreateGraphErrors() };
    auto count{ 0 };
    if (IsAtomLocalEntryAvailable() == true)
    {
        for (auto & [distance, map_value] : m_atom_local_entry->GetDistanceAndMapValueList())
        {
            graph->SetPoint(count, distance, map_value);
            count++;
        }
    }
    else if (IsBondLocalEntryAvailable() == true)
    {
        for (auto & [distance, map_value] : m_bond_local_entry->GetDistanceAndMapValueList())
        {
            graph->SetPoint(count, distance, map_value);
            count++;
        }
    }
    else return nullptr;
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

std::unique_ptr<TGraphErrors> PotentialEntryIterator::CreateAtomXYPositionTomographyGraph(
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

std::unique_ptr<TF1> PotentialEntryIterator::CreateAtomGroupGausFunctionPrior(
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
    return ROOTHelper::CreateGaus3DFunctionIn1D("gaus", amplitude, width);
}

std::unique_ptr<TF1> PotentialEntryIterator::CreateBondGroupGausFunctionPrior(
    GroupKey group_key, const std::string & class_key) const
{
    auto amplitude
    {
        m_model_object->GetBondGroupPotentialEntry(class_key)->GetGausEstimatePrior(group_key, 0)
    };
    auto width
    {
        m_model_object->GetBondGroupPotentialEntry(class_key)->GetGausEstimatePrior(group_key, 1)
    };
    return ROOTHelper::CreateGaus2DFunctionIn1D("gaus", amplitude, width);
}
#endif
