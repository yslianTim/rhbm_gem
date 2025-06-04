#include "ChargeEntryIterator.hpp"
#include "AtomObject.hpp"
#include "ModelObject.hpp"
#include "AtomicChargeEntry.hpp"
#include "GroupChargeEntry.hpp"
#include "AtomicInfoHelper.hpp"
#include "ROOTHelper.hpp"
#include "ArrayStats.hpp"
#include "KeyPacker.hpp"
#include "Constants.hpp"
#include "GlobalEnumClass.hpp"
#include "AtomClassifier.hpp"

#include <cmath>
#include <iostream>

#ifdef HAVE_ROOT
#include <TGraphErrors.h>
#include <TGraph2DErrors.h>
#include <TF1.h>
#include <TH1.h>
#endif

ChargeEntryIterator::ChargeEntryIterator(ModelObject * model_object) :
    m_atom_object{ nullptr }, m_model_object{ model_object }, m_atomic_entry{ nullptr }
{
}

ChargeEntryIterator::ChargeEntryIterator(AtomObject * atom_object) :
    m_atom_object{ atom_object }, m_model_object{ nullptr }
{
    m_atomic_entry = atom_object->GetAtomicChargeEntry();
}

ChargeEntryIterator::~ChargeEntryIterator()
{

}

double ChargeEntryIterator::GetModelEstimateMinimum(int par_id, Element element) const
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
        auto entry{ atom->GetAtomicChargeEntry() };
        gaus_estimate_list.emplace_back(entry->GetModelEstimateMDPDE(par_id));
    }
    return ArrayStats<double>::ComputeMin(gaus_estimate_list.data(), gaus_estimate_list.size());
}

bool ChargeEntryIterator::IsOutlierAtom(const std::string & class_key) const
{
    if (IsAtomicEntryAvailable() == false)
    {
        throw std::runtime_error("Atomic entry is not available.");
    }
    return m_atomic_entry->GetOutlierTag(class_key);
}

bool ChargeEntryIterator::IsAvailableGroupKey(uint64_t group_key, const std::string & class_key) const
{
    return CheckGroupKey(group_key, class_key, false);
}

size_t ChargeEntryIterator::GetResidueCount(
    const std::string & class_key, Residue residue, Structure structure) const
{
    uint64_t group_key{ 0 };
    if (class_key == AtomicInfoHelper::GetResidueClassKey())
    {
        group_key = KeyPackerResidueClass::Pack(residue, Element::CARBON, Remoteness::ALPHA, Branch::NONE, false);
    }
    else if (class_key == AtomicInfoHelper::GetStructureClassKey())
    {
        group_key = KeyPackerStructureClass::Pack(structure, residue, Element::CARBON, Remoteness::ALPHA, Branch::NONE, false);
    }
    if (IsAvailableGroupKey(group_key, class_key) == false)
    {
        return 0;
    }
    return GetAtomObjectList(group_key, class_key).size();
}

double ChargeEntryIterator::GetModelEstimatePrior(
    uint64_t group_key, const std::string & class_key, int par_id) const
{
    if (CheckGroupKey(group_key, class_key) == false)
    {
        throw std::runtime_error("Group key is not available.");
    }
    return m_model_object->GetGroupChargeEntry(class_key)->GetModelEstimatePrior(group_key, par_id);
}

double ChargeEntryIterator::GetModelVariancePrior(
    uint64_t group_key, const std::string & class_key, int par_id) const
{
    if (CheckGroupKey(group_key, class_key) == false)
    {
        throw std::runtime_error("Group key is not available.");
    }
    return m_model_object->GetGroupChargeEntry(class_key)->GetModelVariancePrior(group_key, par_id);
}

const std::vector<AtomObject *> & ChargeEntryIterator::GetAtomObjectList(
    uint64_t group_key, const std::string & class_key) const
{
    if (CheckGroupKey(group_key, class_key) == false)
    {
        throw std::runtime_error("Group key is not available.");
    }
    return m_model_object->GetGroupChargeEntry(class_key)->GetAtomObjectPtrList(group_key);
}

std::unordered_map<int, AtomObject *> ChargeEntryIterator::GetAtomObjectMap(
    uint64_t group_key, const std::string & class_key) const
{
    std::unordered_map<int, AtomObject *> atom_object_map;
    for (auto & atom_object : GetAtomObjectList(group_key, class_key))
    {
        atom_object_map[atom_object->GetSerialID()] = atom_object;
    }
    return atom_object_map;
}

const std::vector<std::tuple<float, float>> & ChargeEntryIterator::GetDistanceAndMapValueList(void) const
{
    if (IsAtomicEntryAvailable() == false)
    {
        throw std::runtime_error("Atomic entry is not available.");
    }
    return m_atomic_entry->GetDistanceAndMapValueList();
}

std::tuple<float, float> ChargeEntryIterator::GetDistanceRange(double margin_rate) const
{
    if (IsAtomicEntryAvailable() == false)
    {
        throw std::runtime_error("Atomic entry is not available.");
    }
    std::vector<float> distance_array;
    distance_array.reserve(static_cast<size_t>(m_atomic_entry->GetDistanceAndMapValueListSize()));
    for (auto & [distance, map_value] : m_atomic_entry->GetDistanceAndMapValueList())
    {
        distance_array.emplace_back(distance);
    }
    return ArrayStats<float>::ComputeScalingRangeTuple(distance_array, static_cast<float>(margin_rate));
}

std::tuple<float, float> ChargeEntryIterator::GetMapValueRange(double margin_rate) const
{
    if (IsAtomicEntryAvailable() == false)
    {
        throw std::runtime_error("Atomic entry is not available.");
    }
    std::vector<float> map_value_array;
    map_value_array.reserve(static_cast<size_t>(m_atomic_entry->GetDistanceAndMapValueListSize()));
    for (auto & [distance, map_value] : m_atomic_entry->GetDistanceAndMapValueList())
    {
        map_value_array.emplace_back(map_value);
    }
    return ArrayStats<float>::ComputeScalingRangeTuple(map_value_array, static_cast<float>(margin_rate));
}

double ChargeEntryIterator::GetModelEstimateMDPDE(int par_id) const
{
    if (IsAtomicEntryAvailable() == false)
    {
        throw std::runtime_error("Atomic entry is not available.");
    }
    return m_atomic_entry->GetModelEstimateMDPDE(par_id);
}

double ChargeEntryIterator::GetInterceptEstimateMDPDE(void) const
{
    if (IsAtomicEntryAvailable() == false)
    {
        return 0.0;
    }
    return m_atomic_entry->GetModelEstimateMDPDE(0);
}

double ChargeEntryIterator::GetScaleEstimateMDPDE(void) const
{
    if (IsAtomicEntryAvailable() == false)
    {
        return 0.0;
    }
    return m_atomic_entry->GetModelEstimateMDPDE(1);
}

double ChargeEntryIterator::GetChargeEstimateMDPDE(void) const
{
    if (IsAtomicEntryAvailable() == false)
    {
        return 0.0;
    }
    return m_atomic_entry->GetModelEstimateMDPDE(2);
}

bool ChargeEntryIterator::IsAtomObjectAvailable(void) const
{
    if (m_atom_object == nullptr)
    {
        std::cerr << "Atom object is not available." << std::endl;
        return false;
    }
    return true;
}

bool ChargeEntryIterator::IsAtomicEntryAvailable(void) const
{
    if (m_atomic_entry == nullptr)
    {
        std::cerr << "Atomic entry is not available." << std::endl;
        return false;
    }
    return true;
}

bool ChargeEntryIterator::IsModelObjectAvailable(void) const
{
    if (m_model_object == nullptr)
    {
        std::cerr << "Model object is not available." << std::endl;
        return false;
    }
    return true;
}

bool ChargeEntryIterator::CheckGroupKey(uint64_t group_key, const std::string & class_key, bool verbose) const
{
    const auto & group_key_set{ m_model_object->GetGroupChargeEntry(class_key)->GetGroupKeySet() };
    if (group_key_set.find(group_key) == group_key_set.end())
    {
        if (verbose == true)
        {
            if (class_key == AtomicInfoHelper::GetElementClassKey())
            {
                auto unpack_key{ KeyPackerElementClass::Unpack(group_key) };
                std::cout <<"Elelemt class group key :" << std::boolalpha
                          << static_cast<int>(std::get<0>(unpack_key)) <<", "
                          << static_cast<int>(std::get<1>(unpack_key)) <<", "
                          << std::get<2>(unpack_key) <<" not found." << std::endl;
            }
            else if (class_key == AtomicInfoHelper::GetResidueClassKey())
            {
                auto unpack_key{ KeyPackerResidueClass::Unpack(group_key) };
                std::cout <<"Residue class group key : tuple<" << std::boolalpha
                          << static_cast<int>(std::get<0>(unpack_key)) <<", "
                          << static_cast<int>(std::get<1>(unpack_key)) <<", "
                          << static_cast<int>(std::get<2>(unpack_key)) <<", "
                          << static_cast<int>(std::get<3>(unpack_key)) <<", "
                          << std::get<4>(unpack_key) <<"> not found." << std::endl;
            }
            else if (class_key == AtomicInfoHelper::GetStructureClassKey())
            {
                auto unpack_key{ KeyPackerStructureClass::Unpack(group_key) };
                std::cout <<"Structure class group key : tuple<" << std::boolalpha
                          << static_cast<int>(std::get<0>(unpack_key)) <<", "
                          << static_cast<int>(std::get<1>(unpack_key)) <<", "
                          << static_cast<int>(std::get<2>(unpack_key)) <<", "
                          << static_cast<int>(std::get<3>(unpack_key)) <<", "
                          << static_cast<int>(std::get<4>(unpack_key)) <<", "
                          << std::get<5>(unpack_key) <<"> not found." << std::endl;
            }
        }
        return false;
    }
    return true;
}

Residue ChargeEntryIterator::GetResidueFromGroupKey(
    uint64_t group_key, const std::string & class_key) const
{
    if (class_key == AtomicInfoHelper::GetElementClassKey())
    {
        std::cout <<"Element class group key is not recording Residue info."<< std::endl;
        return Residue::UNK;
    }
    else if (class_key == AtomicInfoHelper::GetResidueClassKey())
    {
        auto unpack_key{ KeyPackerResidueClass::Unpack(group_key) };
        return static_cast<Residue>(std::get<0>(unpack_key));
    }
    else if (class_key == AtomicInfoHelper::GetStructureClassKey())
    {
        auto unpack_key{ KeyPackerStructureClass::Unpack(group_key) };
        return static_cast<Residue>(std::get<1>(unpack_key));
    }
    return Residue::UNK;
}

#ifdef HAVE_ROOT
std::unique_ptr<TH1D> ChargeEntryIterator::CreateResidueCountHistogram(
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
        hist->SetBinContent(static_cast<int>(residue) + 1, static_cast<double>(count));
    }
    return hist;
}

std::unique_ptr<TGraphErrors> ChargeEntryIterator::CreateModelEstimateToResidueGraph(
    std::vector<uint64_t> & group_key_list, const std::string & class_key, const int par_id)
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
        auto x_value{ static_cast<int>(GetResidueFromGroupKey(group_key, class_key)) };
        auto y_value{ m_model_object->GetGroupChargeEntry(class_key)->GetModelEstimatePrior(group_key, par_id) };
        auto y_error{ m_model_object->GetGroupChargeEntry(class_key)->GetModelVariancePrior(group_key, par_id) };
        //if (par_id == 2)
        //{
            //auto beta_1{ m_model_object->GetGroupChargeEntry(class_key)->GetModelEstimatePrior(group_key, 1) };
            //auto beta_2{ m_model_object->GetGroupChargeEntry(class_key)->GetModelEstimatePrior(group_key, 2) };
            //y_value = beta_2 / (beta_1 + beta_2);
            //y_value = beta_2 / beta_1;
        //}
        graph->SetPoint(count, x_value, y_value);
        graph->SetPointError(count, 0.0, y_error);
        count++;
    }
    return graph;
}

std::unique_ptr<TGraphErrors> ChargeEntryIterator::CreateModelEstimateScatterGraph(
    std::vector<uint64_t> & group_key_list, const std::string & class_key,
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
            GetModelEstimatePrior(group_key, class_key, par1_id),
            GetModelEstimatePrior(group_key, class_key, par2_id));
        graph->SetPointError(count,
            GetModelVariancePrior(group_key, class_key, par1_id),
            GetModelVariancePrior(group_key, class_key, par2_id));
        count++;
    }
    return graph;
}

std::unordered_map<std::string, std::unique_ptr<TGraphErrors>>
ChargeEntryIterator::CreateModelEstimateToResidueIDGraphMap(
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
        if (atom->GetRemoteness() != AtomClassifier::GetMainChainRemoteness(main_chain_element_id)) continue;
        if (residue != Residue::UNK && atom->GetResidue() != residue) continue;
        auto entry{ atom->GetAtomicChargeEntry() };
        auto residue_id{ atom->GetResidueID() };
        auto chain_id{ atom->GetChainID() };
        if (residue_id < 0) continue;
        if (graph_map.find(chain_id) == graph_map.end())
        {
            graph_map[chain_id] = ROOTHelper::CreateGraphErrors();
            count_map[chain_id] = 0;
        }
        auto x_value{ static_cast<double>(residue_id) };
        graph_map[chain_id]->SetPoint(count_map[chain_id], x_value, entry->GetModelEstimateMDPDE(par_id));
        count_map[chain_id]++;
    }
    return graph_map;
}

std::unique_ptr<TGraphErrors> ChargeEntryIterator::CreateModelBasisToMapValueGraph(int basis_id)
{
    if (IsAtomObjectAvailable() == false)
    {
        return nullptr;
    }
    auto data_list{ m_atomic_entry->GetDistanceAndMapValueList() };
    auto map_0_list{ m_atomic_entry->GetDistanceAndNeutralMapValueList() };
    auto map_pos_list{ m_atomic_entry->GetDistanceAndPositiveMapValueList() };
    auto map_neg_list{ m_atomic_entry->GetDistanceAndNegativeMapValueList() };
    auto graph{ ROOTHelper::CreateGraphErrors() };
    auto x{ 0.0 };
    auto y{ 0.0 };
    for (size_t p = 0; p < data_list.size(); p++)
    {
        auto distance_data{ std::get<0>(data_list.at(p)) };
        auto distance_map_0{ std::get<0>(map_0_list.at(p)) };
        auto distance_map_pos{ std::get<0>(map_pos_list.at(p)) };
        auto distance_map_neg{ std::get<0>(map_neg_list.at(p)) };
        if (distance_data != distance_map_0 ||
            distance_data != distance_map_pos ||
            distance_data != distance_map_neg)
        {
            std::cerr << "Distance mismatch in ChargeEntryIterator::CreateModelBasisToMapValueGraph." << std::endl;
            continue;
        }
        if      (basis_id == 0) x = std::get<1>(map_0_list.at(p));
        //else if (basis_id == 1) x = std::get<1>(map_pos_list.at(p)) - std::get<1>(map_0_list.at(p));
        //else if (basis_id == 2) x = std::get<1>(map_0_list.at(p)) - std::get<1>(map_neg_list.at(p));
        else if (basis_id == 1) x = std::get<1>(map_pos_list.at(p));
        else if (basis_id == 2) x = std::get<1>(map_neg_list.at(p));
        y = std::get<1>(data_list.at(p));
        graph->SetPoint(static_cast<int>(p), x, y);
    }
    return graph;
}

std::unique_ptr<TGraphErrors> ChargeEntryIterator::CreateRegressionDataGraph(int basis_id)
{
    if (IsAtomObjectAvailable() == false)
    {
        return nullptr;
    }
    auto data_list{ m_atomic_entry->GetRegressionDataList() };
    auto graph{ ROOTHelper::CreateGraphErrors() };
    auto x{ 0.0 };
    auto count{ 0 };
    for (auto & [x1, x2, y] : data_list)
    {
        if      (basis_id == 0) x = x1;
        else if (basis_id == 1) x = x2;
        graph->SetPoint(count, x, y);
        count++;
    }
    return graph;
}
#endif