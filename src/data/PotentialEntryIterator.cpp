#include "PotentialEntryIterator.hpp"
#include "AtomObject.hpp"
#include "ModelObject.hpp"
#include "AtomicPotentialEntry.hpp"
#include "GroupPotentialEntry.hpp"
#include "ROOTHelper.hpp"
#include "ArrayStats.hpp"

#include <cmath>

#ifdef HAVE_ROOT
#include <TGraphErrors.h>
#include <TF1.h>
#endif

PotentialEntryIterator::PotentialEntryIterator(ModelObject * model_object) :
    m_atom_object{ nullptr }, m_model_object{ model_object }, m_atomic_entry{ nullptr }
{
    auto group_entry1{ m_model_object->GetGroupPotentialEntry(AtomicInfoHelper::GetElementClassKey()) };
    m_element_class_group_entry = dynamic_cast<GroupPotentialEntry<ElementKeyType> *>(group_entry1);
    auto group_entry2{ m_model_object->GetGroupPotentialEntry(AtomicInfoHelper::GetResidueClassKey()) };
    m_residue_class_group_entry = dynamic_cast<GroupPotentialEntry<ResidueKeyType> *>(group_entry2);
}

PotentialEntryIterator::PotentialEntryIterator(AtomObject * atom_object) :
    m_atom_object{ atom_object }, m_model_object{ nullptr },
    m_element_class_group_entry{ nullptr }, m_residue_class_group_entry{ nullptr }
{
    m_atomic_entry = atom_object->GetAtomicPotentialEntry();
}

PotentialEntryIterator::~PotentialEntryIterator()
{

}

bool PotentialEntryIterator::IsAvailableGroupKey(ElementKeyType & group_key) const
{
    return CheckGroupKey(group_key, false);
}

bool PotentialEntryIterator::IsAvailableGroupKey(ResidueKeyType & group_key) const
{
    return CheckGroupKey(group_key, false);
}

std::tuple<double, double> PotentialEntryIterator::GetGausEstimatePrior(ElementKeyType & group_key) const
{
    if (CheckGroupKey(group_key) == false)
    {
        return std::make_tuple(0.0, 0.0);
    }
    return m_element_class_group_entry->GetGausEstimatePrior(&group_key);
}

std::tuple<double, double> PotentialEntryIterator::GetGausEstimatePrior(ResidueKeyType & group_key) const
{
    if (CheckGroupKey(group_key) == false)
    {
        return std::make_tuple(0.0, 0.0);
    }
    return m_residue_class_group_entry->GetGausEstimatePrior(&group_key);
}

std::tuple<double, double> PotentialEntryIterator::GetGausVariancePrior(ElementKeyType & group_key) const
{
    if (CheckGroupKey(group_key) == false)
    {
        return std::make_tuple(0.0, 0.0);
    }
    return m_element_class_group_entry->GetGausVariancePrior(&group_key);
}

std::tuple<double, double> PotentialEntryIterator::GetGausVariancePrior(ResidueKeyType & group_key) const
{
    if (CheckGroupKey(group_key) == false)
    {
        return std::make_tuple(0.0, 0.0);
    }
    return m_residue_class_group_entry->GetGausVariancePrior(&group_key);
}

const std::vector<AtomObject *> & PotentialEntryIterator::GetAtomObjectList(
    ElementKeyType & group_key) const
{
    if (CheckGroupKey(group_key) == false)
    {
        throw std::runtime_error("Group key is not available.");
    }
    return m_element_class_group_entry->GetAtomObjectPtrList(&group_key);
}

const std::vector<AtomObject *> & PotentialEntryIterator::GetAtomObjectList(
    ResidueKeyType & group_key) const
{
    if (CheckGroupKey(group_key) == false)
    {
        throw std::runtime_error("Group key is not available.");
    }
    return m_residue_class_group_entry->GetAtomObjectPtrList(&group_key);
}

const std::vector<std::tuple<float, float>> & PotentialEntryIterator::GetDistanceAndMapValueList(void) const
{
    if (IsAtomicEntryAvailable() == false)
    {
        throw std::runtime_error("Atomic entry is not available.");
    }
    return m_atomic_entry->GetDistanceAndMapValueList();
}

std::tuple<float, float> PotentialEntryIterator::GetDistanceRange(double margin_rate) const
{
    if (IsAtomicEntryAvailable() == false)
    {
        throw std::runtime_error("Atomic entry is not available.");
    }
    std::vector<float> distance_array;
    distance_array.reserve(m_atomic_entry->GetDistanceAndMapValueListSize());
    for (auto & [distance, map_value] : m_atomic_entry->GetDistanceAndMapValueList())
    {
        distance_array.emplace_back(distance);
    }
    return ArrayStats<float>::ComputeScalingRangeTuple(distance_array, margin_rate);
}

std::tuple<float, float> PotentialEntryIterator::GetMapValueRange(double margin_rate) const
{
    if (IsAtomicEntryAvailable() == false)
    {
        throw std::runtime_error("Atomic entry is not available.");
    }
    std::vector<float> map_value_array;
    map_value_array.reserve(m_atomic_entry->GetDistanceAndMapValueListSize());
    for (auto & [distance, map_value] : m_atomic_entry->GetDistanceAndMapValueList())
    {
        map_value_array.emplace_back(map_value);
    }
    return ArrayStats<float>::ComputeScalingRangeTuple(map_value_array, margin_rate);
}

double PotentialEntryIterator::GetAmplitudeEstimateMDPDE(void) const
{
    if (IsAtomicEntryAvailable() == false)
    {
        return 0.0;
    }
    return m_atomic_entry->GetAmplitudeEstimateMDPDE();
}

double PotentialEntryIterator::GetWidthEstimateMDPDE(void) const
{
    if (IsAtomicEntryAvailable() == false)
    {
        return 0.0;
    }
    return m_atomic_entry->GetWidthEstimateMDPDE();
}

bool PotentialEntryIterator::IsAtomObjectAvailable(void) const
{
    if (m_atom_object == nullptr)
    {
        std::cerr << "Atom object is not available." << std::endl;
        return false;
    }
    return true;
}

bool PotentialEntryIterator::IsAtomicEntryAvailable(void) const
{
    if (m_atomic_entry == nullptr)
    {
        std::cerr << "Atomic entry is not available." << std::endl;
        return false;
    }
    return true;
}

bool PotentialEntryIterator::IsModelObjectAvailable(void) const
{
    if (m_model_object == nullptr)
    {
        std::cerr << "Model object is not available." << std::endl;
        return false;
    }
    return true;
}

bool PotentialEntryIterator::CheckGroupKey(ElementKeyType & group_key, bool verbose) const
{
    const auto & group_key_set{ m_element_class_group_entry->GetGroupKeySet() };
    if (group_key_set.find(group_key) == group_key_set.end())
    {
        if (verbose == true)
        {
            std::cout <<"Elelemt class group key :"
                      << std::get<0>(group_key) <<", "
                      << std::get<1>(group_key) <<", "
                      << std::get<2>(group_key) << std::boolalpha <<" not found." << std::endl;
        }
        return false;
    }
    return true;
}

bool PotentialEntryIterator::CheckGroupKey(ResidueKeyType & group_key, bool verbose) const
{
    const auto & group_key_set{ m_residue_class_group_entry->GetGroupKeySet() };
    if (group_key_set.find(group_key) == group_key_set.end())
    {
        if (verbose == true)
        {
            std::cout <<"Residue class group key : tuple<"
                      << std::get<0>(group_key) <<", "
                      << std::get<1>(group_key) <<", "
                      << std::get<2>(group_key) <<", "
                      << std::get<3>(group_key) <<", "
                      << std::get<4>(group_key) << std::boolalpha <<"> not found." << std::endl;
        }
        return false;
    }
    return true;
}

#ifdef HAVE_ROOT
std::unique_ptr<TGraphErrors> PotentialEntryIterator::CreateDistanceToMapValueGraph(void)
{
    if (IsAtomicEntryAvailable() == false)
    {
        return nullptr;
    }
    auto graph{ ROOTHelper::CreateGraphErrors() };

    auto count{ 0 };
    for (auto & [distance, map_value] : m_atomic_entry->GetDistanceAndMapValueList())
    {
        graph->SetPoint(count, distance, map_value);
        count++;
    }
    return graph;
}

std::unique_ptr<TGraphErrors> PotentialEntryIterator::CreateBinnedDistanceToMapValueGraph(
    int bin_size, double x_min, double x_max)
{
    if (IsAtomicEntryAvailable() == false)
    {
        return nullptr;
    }
    auto bin_spacing{ (x_max - x_min) / static_cast<double>(bin_size) };
    auto graph{ ROOTHelper::CreateGraphErrors(bin_size) };
    std::unordered_map<int, std::vector<float>> bin_map;

    for (auto & [distance, map_value] : m_atomic_entry->GetDistanceAndMapValueList())
    {
        auto bin_index{ static_cast<int>(std::floor(distance / bin_spacing)) };
        bin_map[bin_index].emplace_back(map_value);
    }

    for (int i = 0; i < bin_size; i++)
    {
        auto x_value{ x_min + (i + 0.5) * bin_spacing };
        if (bin_map.find(i) == bin_map.end())
        {
            graph->SetPoint(i, x_value, 0.0);
            continue;
        }
        graph->SetPoint(i, x_value, ArrayStats<float>::ComputeMedian(bin_map.at(i)));
    }
    return graph;
}

std::unique_ptr<TF1> PotentialEntryIterator::CreateGroupGausFunctionPrior(
    ElementKeyType & group_key) const
{
    auto gaus_estimate{ GetGausEstimatePrior(group_key) };
    auto amplitude{ std::get<0>(gaus_estimate) };
    auto width{  std::get<1>(gaus_estimate) };
    return ROOTHelper::CreateGausFunction1D("gaus", amplitude, width);
}

std::unique_ptr<TF1> PotentialEntryIterator::CreateGroupGausFunctionPrior(
    ResidueKeyType & group_key) const
{
    auto gaus_estimate{ GetGausEstimatePrior(group_key) };
    auto amplitude{ std::get<0>(gaus_estimate) };
    auto width{  std::get<1>(gaus_estimate) };
    return ROOTHelper::CreateGausFunction1D("gaus", amplitude, width);
}
#endif