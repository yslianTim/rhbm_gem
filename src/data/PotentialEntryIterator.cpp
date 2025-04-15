#include "PotentialEntryIterator.hpp"
#include "AtomObject.hpp"
#include "ModelObject.hpp"
#include "AtomicPotentialEntry.hpp"
#include "GroupPotentialEntry.hpp"
#include "AtomicInfoHelper.hpp"
#include "ROOTHelper.hpp"
#include "ArrayStats.hpp"

#include <cmath>

#ifdef HAVE_ROOT
#include <TGraphErrors.h>
#include <TGraph2DErrors.h>
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

double PotentialEntryIterator::GetGausEstimatePrior(ElementKeyType & group_key, int par_id) const
{
    if (CheckParameterIndex(par_id) == false) return 0.0;
    auto gaus_estimate{ GetGausEstimatePrior(group_key) };
    return (par_id == 0) ? std::get<0>(gaus_estimate) : std::get<1>(gaus_estimate);
}

double PotentialEntryIterator::GetGausEstimatePrior(ResidueKeyType & group_key, int par_id) const
{
    if (CheckParameterIndex(par_id) == false) return 0.0;
    auto gaus_estimate{ GetGausEstimatePrior(group_key) };
    return (par_id == 0) ? std::get<0>(gaus_estimate) : std::get<1>(gaus_estimate);
}

double PotentialEntryIterator::GetGausVariancePrior(ElementKeyType & group_key, int par_id) const
{
    if (CheckParameterIndex(par_id) == false) return 0.0;
    auto gaus_variance{ GetGausVariancePrior(group_key) };
    return (par_id == 0) ? std::get<0>(gaus_variance) : std::get<1>(gaus_variance);
}

double PotentialEntryIterator::GetGausVariancePrior(ResidueKeyType & group_key, int par_id) const
{
    if (CheckParameterIndex(par_id) == false) return 0.0;
    auto gaus_variance{ GetGausVariancePrior(group_key) };
    return (par_id == 0) ? std::get<0>(gaus_variance) : std::get<1>(gaus_variance);
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

std::unordered_map<int, AtomObject *> PotentialEntryIterator::GetAtomObjectMap(
    ElementKeyType & group_key) const
{
    std::unordered_map<int, AtomObject *> atom_object_map;
    for (auto & atom_object : GetAtomObjectList(group_key))
    {
        atom_object_map[atom_object->GetSerialID()] = atom_object;
    }
    return atom_object_map;
}

std::unordered_map<int, AtomObject *> PotentialEntryIterator::GetAtomObjectMap(
    ResidueKeyType & group_key) const
{
    std::unordered_map<int, AtomObject *> atom_object_map;
    for (auto & atom_object : GetAtomObjectList(group_key))
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
    binned_distance_and_map_value_list.reserve(bin_size);
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

bool PotentialEntryIterator::CheckParameterIndex(int par_id) const
{
    if (par_id < 0 || par_id > 1)
    {
        std::cerr << "Invalid parameter index: " << par_id << std::endl;
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
std::unique_ptr<TGraphErrors> PotentialEntryIterator::CreateGausEstimateToResidueGraph(
    std::vector<ResidueKeyType> & group_key_list, const int par_id)
{
    if (IsModelObjectAvailable() == false || CheckParameterIndex(par_id) == false)
    {
        return nullptr;
    }
    auto graph{ ROOTHelper::CreateGraphErrors() };
    
    auto count{ 0 };
    for (auto & group_key : group_key_list)
    {
        if (IsAvailableGroupKey(group_key) == false) continue;
        auto gaus_estimate{ GetGausEstimatePrior(group_key) };
        auto gaus_variance{ GetGausVariancePrior(group_key) };
        auto y_value{ (par_id == 0) ? std::get<0>(gaus_estimate) : std::get<1>(gaus_estimate) };
        auto y_error{ (par_id == 0) ? std::get<0>(gaus_variance) : std::get<1>(gaus_variance) };
        auto x_value{ std::get<0>(group_key) };
        graph->SetPoint(count, x_value, y_value);
        graph->SetPointError(count, 0.0, y_error);
        count++;
    }
    return graph;
}

std::unique_ptr<TGraphErrors> PotentialEntryIterator::CreateGausEstimateScatterGraph(
    std::vector<ResidueKeyType> & group_key_list, bool reverse)
{
    if (IsModelObjectAvailable() == false)
    {
        return nullptr;
    }
    auto graph{ ROOTHelper::CreateGraphErrors() };
    
    auto count{ 0 };
    for (auto & group_key : group_key_list)
    {
        if (IsAvailableGroupKey(group_key) == false) continue;
        auto gaus_estimate{ GetGausEstimatePrior(group_key) };
        auto gaus_variance{ GetGausVariancePrior(group_key) };
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
    ElementKeyType & group_key, double range, int par_id)
{
    if (IsModelObjectAvailable() == false || CheckParameterIndex(par_id) == false)
    {
        return nullptr;
    }
    auto graph{ ROOTHelper::CreateGraphErrors() };
    auto kd_tree_root{ m_model_object->GetKDTreeRoot() };

    auto count{ 0 };
    for (auto & atom : GetAtomObjectList(group_key))
    {
        auto in_range_atom_list{ KDTreeAlgorithm<AtomObject>::RangeSearch(kd_tree_root, atom, range) };
        auto atom_entry{ atom->GetAtomicPotentialEntry() };
        auto gaus_estimate{ (par_id == 0) ? atom_entry->GetAmplitudeEstimateMDPDE() : atom_entry->GetWidthEstimateMDPDE() };
        graph->SetPoint(count, in_range_atom_list.size(), gaus_estimate);
        count++;
    }
    return graph;
}

std::unique_ptr<TGraphErrors> PotentialEntryIterator::CreateInRangeAtomsToGausEstimateGraph(
    ResidueKeyType & group_key, double range, int par_id)
{
    if (IsModelObjectAvailable() == false || CheckParameterIndex(par_id) == false)
    {
        return nullptr;
    }
    auto graph{ ROOTHelper::CreateGraphErrors() };
    auto kd_tree_root{ m_model_object->GetKDTreeRoot() };

    auto count{ 0 };
    for (auto & atom : GetAtomObjectList(group_key))
    {
        auto in_range_atom_list{ KDTreeAlgorithm<AtomObject>::RangeSearch(kd_tree_root, atom, range) };
        auto atom_entry{ atom->GetAtomicPotentialEntry() };
        auto gaus_estimate{ (par_id == 0) ? atom_entry->GetAmplitudeEstimateMDPDE() : atom_entry->GetWidthEstimateMDPDE() };
        graph->SetPoint(count, in_range_atom_list.size(), gaus_estimate);
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
    std::vector<ResidueKeyType> & group_key_list, double normalized_z_pos, double z_ratio_window)
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
        if (IsAvailableGroupKey(group_key) == false) continue;
        for (auto & atom : GetAtomObjectList(group_key))
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
    std::vector<ResidueKeyType> & group_key_list, double normalized_z_pos, double z_ratio_window, int par_id)
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
        if (IsAvailableGroupKey(group_key) == false) continue;
        for (auto & atom : GetAtomObjectList(group_key))
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