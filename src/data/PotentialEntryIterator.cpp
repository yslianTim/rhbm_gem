#include "PotentialEntryIterator.hpp"
#include "AtomObject.hpp"
#include "ModelObject.hpp"
#include "AtomicPotentialEntry.hpp"
#include "GroupPotentialEntry.hpp"
#include "AtomicInfoHelper.hpp"
#include "ROOTHelper.hpp"
#include "ArrayStats.hpp"
#include "KeyPacker.hpp"
#include "Constants.hpp"

#include <cmath>
#include <iostream>

#ifdef HAVE_ROOT
#include <TGraphErrors.h>
#include <TGraph2DErrors.h>
#include <TF1.h>
#endif

PotentialEntryIterator::PotentialEntryIterator(ModelObject * model_object) :
    m_atom_object{ nullptr }, m_model_object{ model_object }, m_atomic_entry{ nullptr }
{
}

PotentialEntryIterator::PotentialEntryIterator(AtomObject * atom_object) :
    m_atom_object{ atom_object }, m_model_object{ nullptr }
{
    m_atomic_entry = atom_object->GetAtomicPotentialEntry();
}

PotentialEntryIterator::~PotentialEntryIterator()
{

}

bool PotentialEntryIterator::IsAvailableGroupKey(uint64_t group_key, const std::string & class_key) const
{
    return CheckGroupKey(group_key, class_key, false);
}

std::tuple<double, double> PotentialEntryIterator::GetGausEstimatePrior(
    uint64_t group_key, const std::string & class_key) const
{
    if (CheckGroupKey(group_key, class_key) == false)
    {
        return std::make_tuple(0.0, 0.0);
    }
    return m_model_object->GetGroupPotentialEntry(class_key)->GetGausEstimatePrior(group_key);
}

std::tuple<double, double> PotentialEntryIterator::GetGausVariancePrior(
    uint64_t group_key, const std::string & class_key) const
{
    if (CheckGroupKey(group_key, class_key) == false)
    {
        return std::make_tuple(0.0, 0.0);
    }
    return m_model_object->GetGroupPotentialEntry(class_key)->GetGausVariancePrior(group_key);
}

double PotentialEntryIterator::GetGausEstimatePrior(uint64_t group_key, const std::string & class_key, int par_id) const
{
    if (CheckParameterIndex(par_id) == false) return 0.0;
    auto gaus_estimate{ GetGausEstimatePrior(group_key, class_key) };
    return (par_id == 0) ? std::get<0>(gaus_estimate) : std::get<1>(gaus_estimate);
}

double PotentialEntryIterator::GetGausVariancePrior(uint64_t group_key, const std::string & class_key, int par_id) const
{
    if (CheckParameterIndex(par_id) == false) return 0.0;
    auto gaus_variance{ GetGausVariancePrior(group_key, class_key) };
    return (par_id == 0) ? std::get<0>(gaus_variance) : std::get<1>(gaus_variance);
}

const std::vector<AtomObject *> & PotentialEntryIterator::GetAtomObjectList(
    uint64_t group_key, const std::string & class_key) const
{
    if (CheckGroupKey(group_key, class_key) == false)
    {
        throw std::runtime_error("Group key is not available.");
    }
    return m_model_object->GetGroupPotentialEntry(class_key)->GetAtomObjectPtrList(group_key);
}

std::unordered_map<int, AtomObject *> PotentialEntryIterator::GetAtomObjectMap(
    uint64_t group_key, const std::string & class_key) const
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
    return m_atomic_entry->GetDistanceAndMapValueList();
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
    for (auto & [distance, map_value] : m_atomic_entry->GetDistanceAndMapValueList())
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
    distance_array.reserve(static_cast<size_t>(m_atomic_entry->GetDistanceAndMapValueListSize()));
    for (auto & [distance, map_value] : m_atomic_entry->GetDistanceAndMapValueList())
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
    map_value_array.reserve(static_cast<size_t>(m_atomic_entry->GetDistanceAndMapValueListSize()));
    for (auto & [distance, map_value] : m_atomic_entry->GetDistanceAndMapValueList())
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

bool PotentialEntryIterator::CheckParameterIndex(int par_id) const
{
    if (par_id < 0 || par_id > 1)
    {
        std::cerr << "Invalid parameter index: " << par_id << std::endl;
        return false;
    }
    return true;
}

bool PotentialEntryIterator::CheckGroupKey(uint64_t group_key, const std::string & class_key, bool verbose) const
{
    const auto & group_key_set{ m_model_object->GetGroupPotentialEntry(class_key)->GetGroupKeySet() };
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
        }
        return false;
    }
    return true;
}

#ifdef HAVE_ROOT
std::unique_ptr<TGraphErrors> PotentialEntryIterator::CreateBfactorToWidthScatterGraph(
    uint64_t group_key, const std::string & class_key)
{
    if (IsModelObjectAvailable() == false)
    {
        return nullptr;
    }
    if (IsAvailableGroupKey(group_key, class_key) == false)
    {
        std::cerr << "Group key is not available." << std::endl;
        return nullptr;
    }
    auto graph{ ROOTHelper::CreateGraphErrors() };
    const auto & atom_list{ GetAtomObjectList(group_key, class_key) };

    auto count{ 0 };
    for (auto & atom : atom_list)
    {
        auto entry{ atom->GetAtomicPotentialEntry() };
        graph->SetPoint(count, atom->GetTemperature()/(8.0*Constants::pi_square), entry->GetWidthEstimateMDPDE());
        count++;
    }

    return graph;
}

std::unique_ptr<TGraphErrors> PotentialEntryIterator::CreateGausEstimateToResidueGraph(
    std::vector<uint64_t> & group_key_list, const std::string & class_key, const int par_id)
{
    if (IsModelObjectAvailable() == false || CheckParameterIndex(par_id) == false)
    {
        return nullptr;
    }
    auto graph{ ROOTHelper::CreateGraphErrors() };
    
    auto count{ 0 };
    for (auto & group_key : group_key_list)
    {
        if (IsAvailableGroupKey(group_key, class_key) == false) continue;
        auto unpack_key{ KeyPackerResidueClass::Unpack(group_key) };
        auto residue_id{ std::get<0>(unpack_key) };
        auto gaus_estimate{ GetGausEstimatePrior(group_key, class_key) };
        auto gaus_variance{ GetGausVariancePrior(group_key, class_key) };
        auto y_value{ (par_id == 0) ? std::get<0>(gaus_estimate) : std::get<1>(gaus_estimate) };
        auto y_error{ (par_id == 0) ? std::get<0>(gaus_variance) : std::get<1>(gaus_variance) };
        auto x_value{ static_cast<int>(residue_id) };
        graph->SetPoint(count, x_value, y_value);
        graph->SetPointError(count, 0.0, y_error);
        count++;
    }
    return graph;
}

std::unique_ptr<TGraphErrors> PotentialEntryIterator::CreateGausEstimateScatterGraph(
    std::vector<uint64_t> & group_key_list, const std::string & class_key, bool reverse)
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
        auto gaus_estimate{ GetGausEstimatePrior(group_key, class_key) };
        auto gaus_variance{ GetGausVariancePrior(group_key, class_key) };
        if (reverse == false)
        {
            graph->SetPoint(count, std::get<0>(gaus_estimate), std::get<1>(gaus_estimate));
            graph->SetPointError(count, std::get<0>(gaus_variance), std::get<1>(gaus_variance));
        }
        else
        {
            graph->SetPoint(count, std::get<1>(gaus_estimate), std::get<0>(gaus_estimate));
            graph->SetPointError(count, std::get<1>(gaus_variance), std::get<0>(gaus_variance));
        }
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
    for (auto & atom : m_model_object->GetComponentsList())
    {
        if (atom->GetElement() == element && atom->GetSelectedFlag() == true)
        {
            auto entry{ atom->GetAtomicPotentialEntry() };
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
    }
    return graph;
}

std::unique_ptr<TGraphErrors> PotentialEntryIterator::CreateGausEstimateScatterGraph(
    uint64_t group_key1, uint64_t group_key2, const std::string & class_key, const int par_id)
{
    if (IsModelObjectAvailable() == false)
    {
        return nullptr;
    }
    if (IsAvailableGroupKey(group_key1, class_key) == false || IsAvailableGroupKey(group_key2, class_key) == false)
    {
        std::cerr << "Group key is not available." << std::endl;
        return nullptr;
    }
    auto graph{ ROOTHelper::CreateGraphErrors() };

    const auto & atom_list1{ GetAtomObjectList(group_key1, class_key) };
    const auto & atom_list2{ GetAtomObjectList(group_key2, class_key) };

    auto count{ 0 };
    for (auto atom1 : atom_list1)
    {
        auto entry1{ atom1->GetAtomicPotentialEntry() };
        for (auto atom2 : atom_list2)
        {
            auto entry2{ atom2->GetAtomicPotentialEntry() };
            if (atom1->GetResidueID() == atom2->GetResidueID() &&
                atom1->GetChainID() == atom2->GetChainID())
            {
                auto gaus_estimate1{ (par_id == 0) ? entry1->GetAmplitudeEstimateMDPDE() : entry1->GetWidthEstimateMDPDE() };
                auto gaus_estimate2{ (par_id == 0) ? entry2->GetAmplitudeEstimateMDPDE() : entry2->GetWidthEstimateMDPDE() };
                graph->SetPoint(count, gaus_estimate1, gaus_estimate2);
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
    uint64_t group_key, const std::string & class_key, double range, int par_id)
{
    if (IsModelObjectAvailable() == false || CheckParameterIndex(par_id) == false)
    {
        return nullptr;
    }
    auto graph{ ROOTHelper::CreateGraphErrors() };
    auto kd_tree_root{ m_model_object->GetKDTreeRoot() };

    auto count{ 0 };
    for (auto & atom : GetAtomObjectList(group_key, class_key))
    {
        auto in_range_atom_list{ KDTreeAlgorithm<AtomObject>::RangeSearch(kd_tree_root, atom, range) };
        auto atom_entry{ atom->GetAtomicPotentialEntry() };
        auto gaus_estimate{ (par_id == 0) ? atom_entry->GetAmplitudeEstimateMDPDE() : atom_entry->GetWidthEstimateMDPDE() };
        graph->SetPoint(count, static_cast<double>(in_range_atom_list.size()), gaus_estimate);
        count++;
    }
    return graph;
}

std::unique_ptr<TGraphErrors> PotentialEntryIterator::CreateXYPositionTomographyGraph(
    double normalized_z_pos, double z_ratio_window)
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
    for (auto & atom : m_model_object->GetComponentsList())
    {
        auto position{ atom->GetPosition() };
        if (position.at(2) < z_window_min || position.at(2) >= z_window_max)
        {
            continue;
        }
        graph->SetPoint(count, position.at(0), position.at(1));
        count++;
    }
    return graph;
}

std::unique_ptr<TGraphErrors> PotentialEntryIterator::CreateXYPositionTomographyGraph(
    std::vector<uint64_t> & group_key_list, const std::string & class_key, double normalized_z_pos, double z_ratio_window)
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
    std::vector<uint64_t> & group_key_list, const std::string & class_key, double normalized_z_pos, double z_ratio_window, int par_id)
{
    if (IsModelObjectAvailable() == false || CheckParameterIndex(par_id) == false)
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
            auto entry_iter{ atom->GetAtomicPotentialEntry() };
            auto gaus_estimate{ 0.0 };
            switch(par_id)
            {
                case 0:
                    gaus_estimate = entry_iter->GetAmplitudeEstimatePosterior(AtomicInfoHelper::GetResidueClassKey());
                    break;
                case 1:
                    gaus_estimate = entry_iter->GetWidthEstimatePosterior(AtomicInfoHelper::GetResidueClassKey());
                    break;
                default:
                    std::cerr << "Invalid parameter id." << std::endl;
                    return nullptr;
            }
            graph->SetPoint(count, position.at(0), position.at(1), gaus_estimate);
            count++;
        }
    }
    return graph;
}

std::unique_ptr<TF1> PotentialEntryIterator::CreateGroupGausFunctionPrior(
    uint64_t group_key, const std::string & class_key) const
{
    auto gaus_estimate{ GetGausEstimatePrior(group_key, class_key) };
    auto amplitude{ std::get<0>(gaus_estimate) };
    auto width{  std::get<1>(gaus_estimate) };
    return ROOTHelper::CreateGausFunction1D("gaus", amplitude, width);
}
#endif