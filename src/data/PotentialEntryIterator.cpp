#include <rhbm_gem/data/PotentialEntryIterator.hpp>
#include <rhbm_gem/data/AtomObject.hpp>
#include <rhbm_gem/data/BondObject.hpp>
#include <rhbm_gem/data/ModelObject.hpp>
#include <rhbm_gem/data/LocalPotentialEntry.hpp>
#include <rhbm_gem/data/GroupPotentialEntry.hpp>
#include <rhbm_gem/utils/ChemicalDataHelper.hpp>
#include <rhbm_gem/utils/ROOTHelper.hpp>
#include <rhbm_gem/utils/ArrayStats.hpp>
#include <rhbm_gem/utils/KeyPacker.hpp>
#include <rhbm_gem/utils/Constants.hpp>
#include <rhbm_gem/utils/GlobalEnumClass.hpp>
#include <rhbm_gem/utils/KDTreeAlgorithm.hpp>
#include <rhbm_gem/data/AtomClassifier.hpp>
#include <rhbm_gem/data/BondClassifier.hpp>
#include <rhbm_gem/utils/GausLinearTransformHelper.hpp>
#include <rhbm_gem/utils/Logger.hpp>

#include <Eigen/Dense>
#include <sstream>
#include <cmath>
#include <map>

#ifdef HAVE_ROOT
#include <TGraphErrors.h>
#include <TGraph2DErrors.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#endif

namespace rhbm_gem {

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

double PotentialEntryIterator::GetAtomGausEstimateMean(
    GroupKey group_key, const std::string & class_key, int par_id) const
{
    if (CheckAtomGroupKey(group_key, class_key) == false)
    {
        throw std::runtime_error("Group key is not available.");
    }
    return m_model_object->GetAtomGroupPotentialEntry(class_key)->GetGausEstimateMean(group_key, par_id);
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

std::vector<AtomObject *> PotentialEntryIterator::GetOutlierAtomObjectList(
    GroupKey group_key, const std::string & class_key) const
{
    auto atom_list{ GetAtomObjectList(group_key, class_key) };
    std::vector<AtomObject *> outlier_atom_list;
    outlier_atom_list.reserve(atom_list.size());
    for (auto & atom : atom_list)
    {
        if (atom->GetLocalPotentialEntry()->GetOutlierTag(class_key) == true)
        {
            outlier_atom_list.emplace_back(atom);
        }
    }
    return outlier_atom_list;
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

std::vector<std::tuple<float, float>>
PotentialEntryIterator::GetLinearModelDistanceAndMapValueList() const
{
    Eigen::VectorXd model_par_init{ Eigen::VectorXd::Zero(3) };
    const auto & data_array{ GetDistanceAndMapValueList() };
    std::vector<std::tuple<float, float>> linear_model_distance_and_map_value_list;
    linear_model_distance_and_map_value_list.reserve(data_array.size());
    for (auto & [distance, map_value] : data_array)
    {
        if (map_value <= 0.0f) continue;
        auto data_vector{
            GausLinearTransformHelper::BuildLinearModelDataVector(distance, map_value, model_par_init)
        };
        linear_model_distance_and_map_value_list.emplace_back(
            std::make_tuple(data_vector(1), data_vector(2)))
        ;
    }
    return linear_model_distance_and_map_value_list;
}

const std::vector<std::tuple<float, float>> &
PotentialEntryIterator::GetDistanceAndMapValueList() const
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
    for (auto & [distance, map_value] : GetDistanceAndMapValueList())
    {
        auto bin_index{ static_cast<int>(std::floor(distance / bin_spacing)) };
        bin_map[bin_index].emplace_back(map_value);
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
    const auto & data_array{ GetDistanceAndMapValueList() };
    std::vector<float> distance_array;
    distance_array.reserve(static_cast<size_t>(data_array.size()));
    for (auto & [distance, map_value] : data_array)
    {
        distance_array.emplace_back(distance);
    }
    return ArrayStats<float>::ComputeScalingRangeTuple(distance_array, static_cast<float>(margin_rate));
}

std::tuple<float, float> PotentialEntryIterator::GetMapValueRange(double margin_rate) const
{
    const auto & data_array{ GetDistanceAndMapValueList() };
    std::vector<float> map_value_array;
    map_value_array.reserve(data_array.size());
    for (auto & [distance, map_value] : data_array)
    {
        map_value_array.emplace_back(map_value);
    }
    return ArrayStats<float>::ComputeScalingRangeTuple(map_value_array, static_cast<float>(margin_rate));
}

std::tuple<float, float> PotentialEntryIterator::GetBinnedMapValueRange(
    int bin_size, double x_min, double x_max, double margin_rate) const
{
    auto data_array{ GetBinnedDistanceAndMapValueList(bin_size, x_min, x_max) };
    std::vector<float> map_value_array;
    map_value_array.reserve(data_array.size());
    for (auto & [distance, map_value] : data_array)
    {
        map_value_array.emplace_back(map_value);
    }
    return ArrayStats<float>::ComputeScalingRangeTuple(map_value_array, static_cast<float>(margin_rate));
}

double PotentialEntryIterator::GetAmplitudeEstimateMDPDE() const
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

double PotentialEntryIterator::GetAmplitudeEstimatePosterior(const std::string & class_key) const
{
    if (IsAtomLocalEntryAvailable() == true)
    {
        return m_atom_local_entry->GetAmplitudeEstimatePosterior(class_key);
    }
    else if (IsBondLocalEntryAvailable() == true)
    {
        return m_bond_local_entry->GetAmplitudeEstimatePosterior(class_key);
    }
    return 0.0;
}

double PotentialEntryIterator::GetAmplitudeVariancePosterior(const std::string & class_key) const
{
    if (IsAtomLocalEntryAvailable() == true)
    {
        return m_atom_local_entry->GetAmplitudeVariancePosterior(class_key);
    }
    else if (IsBondLocalEntryAvailable() == true)
    {
        return m_bond_local_entry->GetAmplitudeVariancePosterior(class_key);
    }
    return 0.0;
}

double PotentialEntryIterator::GetWidthEstimateMDPDE() const
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

double PotentialEntryIterator::GetWidthEstimatePosterior(const std::string & class_key) const
{
    if (IsAtomLocalEntryAvailable() == true)
    {
        return m_atom_local_entry->GetWidthEstimatePosterior(class_key);
    }
    else if (IsBondLocalEntryAvailable() == true)
    {
        return m_bond_local_entry->GetWidthEstimatePosterior(class_key);
    }
    return 0.0;
}

double PotentialEntryIterator::GetWidthVariancePosterior(const std::string & class_key) const
{
    if (IsAtomLocalEntryAvailable() == true)
    {
        return m_atom_local_entry->GetWidthVariancePosterior(class_key);
    }
    else if (IsBondLocalEntryAvailable() == true)
    {
        return m_bond_local_entry->GetWidthVariancePosterior(class_key);
    }
    return 0.0;
}

double PotentialEntryIterator::GetAlphaR() const
{
    if (IsAtomLocalEntryAvailable() == true)
    {
        return m_atom_local_entry->GetAlphaR();
    }
    else if (IsBondLocalEntryAvailable() == true)
    {
        return m_bond_local_entry->GetAlphaR();
    }
    throw std::runtime_error(
        "Local entry is not available. Use GetAtomAlphaR/GetBondAlphaR for model-level access.");
}

double PotentialEntryIterator::GetAtomAlphaR(
    GroupKey group_key, const std::string & class_key) const
{
    const auto & atom_list{ GetAtomObjectList(group_key, class_key) };
    if (atom_list.empty())
    {
        throw std::runtime_error("Atom group has no members.");
    }
    auto local_entry{ atom_list.front()->GetLocalPotentialEntry() };
    if (local_entry == nullptr)
    {
        throw std::runtime_error("Local entry of the first atom member is not available.");
    }
    return local_entry->GetAlphaR();
}

double PotentialEntryIterator::GetBondAlphaR(
    GroupKey group_key, const std::string & class_key) const
{
    const auto & bond_list{ GetBondObjectList(group_key, class_key) };
    if (bond_list.empty())
    {
        throw std::runtime_error("Bond group has no members.");
    }
    auto local_entry{ bond_list.front()->GetLocalPotentialEntry() };
    if (local_entry == nullptr)
    {
        throw std::runtime_error("Local entry of the first bond member is not available.");
    }
    return local_entry->GetAlphaR();
}

double PotentialEntryIterator::GetAtomAlphaG(
    GroupKey group_key, const std::string & class_key) const
{
    if (CheckAtomGroupKey(group_key, class_key) == false)
    {
        throw std::runtime_error("Group key is not available.");
    }
    return m_model_object->GetAtomGroupPotentialEntry(class_key)->GetAlphaG(group_key);
}

double PotentialEntryIterator::GetBondAlphaG(
    GroupKey group_key, const std::string & class_key) const
{
    if (CheckBondGroupKey(group_key, class_key) == false)
    {
        throw std::runtime_error("Group key is not available.");
    }
    return m_model_object->GetBondGroupPotentialEntry(class_key)->GetAlphaG(group_key);
}

bool PotentialEntryIterator::IsAtomObjectAvailable() const
{
    if (m_atom_object == nullptr)
    {
        Logger::Log(LogLevel::Error, "Atom object is not available.");
        return false;
    }
    return true;
}

bool PotentialEntryIterator::IsBondObjectAvailable() const
{
    if (m_bond_object == nullptr)
    {
        Logger::Log(LogLevel::Error, "Bond object is not available.");
        return false;
    }
    return true;
}

bool PotentialEntryIterator::IsAtomLocalEntryAvailable() const
{
    if (m_atom_local_entry == nullptr)
    {
        return false;
    }
    return true;
}

bool PotentialEntryIterator::IsBondLocalEntryAvailable() const
{
    if (m_bond_local_entry == nullptr)
    {
        return false;
    }
    return true;
}

bool PotentialEntryIterator::IsModelObjectAvailable() const
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
std::unique_ptr<TH1D> PotentialEntryIterator::CreateComponentCountHistogram(
    std::vector<GroupKey> & group_key_list, const std::string & class_key) const
{
    if (IsModelObjectAvailable() == false)
    {
        return nullptr;
    }
    auto component_size{ static_cast<int>(group_key_list.size()) };
    auto hist{
        ROOTHelper::CreateHist1D("hist_" + class_key, "Component Count",
            component_size, -0.5, component_size-0.5)
    };
    for (size_t i = 0; i < group_key_list.size(); i++)
    {
        auto count{ GetAtomObjectList(group_key_list.at(i), class_key).size() };
        hist->SetBinContent(static_cast<int>(i+1), static_cast<double>(count));
    }
    return hist;
}

std::unique_ptr<TH2D> PotentialEntryIterator::CreateAtomResidueCountHistogram2D(
    const std::string & class_key)
{
    if (IsModelObjectAvailable() == false)
    {
        return nullptr;
    }
    auto hist{ ROOTHelper::CreateHist2D("hist", "Residue Count", 20, -0.5, 19.5, 3, -0.5, 2.5) };
    std::vector<Structure> structure_list{ Structure::FREE, Structure::HELX_P, Structure::SHEET };
    for (size_t i = 0; i < structure_list.size(); i++)
    {
        auto structure{ structure_list.at(i) };
        for (auto & residue : ChemicalDataHelper::GetStandardAminoAcidList())
        {
            auto count{ GetAtomResidueCount(class_key, residue, structure) };
            hist->SetBinContent(
                static_cast<int>(residue), static_cast<int>(i+1), static_cast<double>(count)
            );
        }
    }
    return hist;
}

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

std::unique_ptr<TH1D> PotentialEntryIterator::CreateAtomGausEstimateHistogram(
    GroupKey group_key, const std::string & class_key, int par_id) const
{
    if (IsModelObjectAvailable() == false)
    {
        return nullptr;
    }
    if (IsAvailableAtomGroupKey(group_key, class_key) == false)
    {
        Logger::Log(LogLevel::Error, "Group key is not available.");
        return nullptr;
    }

    const auto & atom_list{ GetAtomObjectList(group_key, class_key) };
    std::vector<double> gaus_estimate_list;
    gaus_estimate_list.reserve(atom_list.size());
    for (auto atom : atom_list)
    {
        auto local_entry{ atom->GetLocalPotentialEntry() };
        gaus_estimate_list.emplace_back(local_entry->GetGausEstimateMDPDE(par_id));
    }

    double x_min{ 0.0 };
    double x_max{ 1.0 };

    auto estimate_range{
        ArrayStats<double>::ComputeScalingPercentileRangeTuple(gaus_estimate_list, 0.1, 0.05, 0.95)
    };

    if (gaus_estimate_list.size() > 1)
    {
        x_min = std::get<0>(estimate_range);
        x_max = std::get<1>(estimate_range);
    }
    else
    {
        x_max = 2.0 * std::ceil(gaus_estimate_list.at(0));
    }
    if (x_max == 0.0) x_max = 1.0;

    auto estimate_average{ 0.5 * (x_max + x_min) };
    if (par_id == 0 && (x_max - x_min) < 5.0)
    {
        x_max = estimate_average + 2.5;
        x_min = estimate_average - 2.5;
    }
    else if (par_id == 1 && (x_max - x_min) < 0.1)
    {
        x_max = estimate_average + 0.05;
        x_min = estimate_average - 0.05;
    }

    auto hist_name{ std::to_string(group_key) + "_" + class_key + "_par" + std::to_string(par_id) };
    auto hist{ ROOTHelper::CreateHist1D(hist_name.data(), "", 25, x_min, x_max) };
    for (auto & value : gaus_estimate_list) hist->Fill(value);
    return hist;
}

std::unique_ptr<TH1D> PotentialEntryIterator::CreateLinearModelDataHistogram(int dimension_id) const
{
    auto data_array{ GetLinearModelDistanceAndMapValueList() };
    std::vector<float> data_list;
    data_list.reserve(data_array.size());
    for (auto & [distance, map_value] : data_array)
    {
        switch (dimension_id)
        {
            case 0:
                data_list.emplace_back(distance);
                break;
            case 1:
                data_list.emplace_back(map_value);
                break;
            default:
                throw std::runtime_error("Dimension id is invalid.");
        }
    }

    auto data_range{ ArrayStats<float>::ComputeScalingRangeTuple(data_list, 0.1f) };
    double x_min{ std::get<0>(data_range) };
    double x_max{ std::get<1>(data_range) };

    auto hist_name{ "data_dim" + std::to_string(dimension_id) };
    auto hist{ ROOTHelper::CreateHist1D(hist_name.data(), "", 25, x_min, x_max) };
    for (auto & value : data_list) hist->Fill(value);
    return hist;
}

std::unique_ptr<TH2D> PotentialEntryIterator::CreateDistanceToMapValueHistogram(
    int x_bin_size, int y_bin_size) const
{
    auto distance_range{ GetDistanceRange(0.0) };
    auto map_value_range{ GetMapValueRange(0.1) };
    auto hist{
        ROOTHelper::CreateHist2D(
            "hist_distance_mapvalue", "Distance vs Map Value",
            x_bin_size, std::get<0>(distance_range), std::get<1>(distance_range),
            y_bin_size, std::get<0>(map_value_range), std::get<1>(map_value_range))
    };
    for (auto & [distance, map_value] : GetDistanceAndMapValueList())
    {
        hist->Fill(distance, map_value);
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
    
    //auto class_key{ ChemicalDataHelper::GetComponentAtomClassKey() };
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
        //auto gaus_value{ entry->GetGausEstimatePosterior(class_key, par_id) };
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
PotentialEntryIterator::CreateAtomMapValueToSequenceIDGraphMap(
    size_t main_chain_element_id, Residue residue)
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
        graph_map[chain_id]->SetPoint(count_map[chain_id], x_value, entry->GetMapValueNearCenter());
        count_map[chain_id]++;
    }
    return graph_map;
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

std::unordered_map<std::string, std::unique_ptr<TGraphErrors>>
PotentialEntryIterator::CreateAtomGausEstimatePosteriorToSequenceIDGraphMap(
    size_t main_chain_element_id, const std::string & class_key, const int par_id, Residue residue)
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
        graph_map[chain_id]->SetPoint(
            count_map[chain_id], x_value, entry->GetGausEstimatePosterior(class_key, par_id)
        );
        graph_map[chain_id]->SetPointError(
            count_map[chain_id], 0.0, entry->GetGausVariancePosterior(class_key, par_id)
        );
        count_map[chain_id]++;
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

std::unique_ptr<TGraphErrors> PotentialEntryIterator::CreateAtomGausEstimateToSpotGraph(
    std::vector<GroupKey> & group_key_list, const std::string & class_key, const int par_id)
{
    if (IsModelObjectAvailable() == false)
    {
        return nullptr;
    }
    auto graph{ ROOTHelper::CreateGraphErrors() };
    
    auto count{ 0 };
    for (size_t i = 0; i < group_key_list.size(); i++)
    {
        auto & group_key{ group_key_list[i] };
        if (IsAvailableAtomGroupKey(group_key, class_key) == false) continue;
        auto x_value{ static_cast<double>(i) };
        auto y_value{ m_model_object->GetAtomGroupPotentialEntry(class_key)->GetGausEstimatePrior(group_key, par_id) };
        auto y_error{ m_model_object->GetAtomGroupPotentialEntry(class_key)->GetGausVariancePrior(group_key, par_id) };
        graph->SetPoint(count, x_value, y_value);
        graph->SetPointError(count, 0.0, y_error);
        count++;
    }
    return graph;
}

std::unique_ptr<TGraphErrors> PotentialEntryIterator::CreateAtomGausEstimateToAtomIdGraph(
    const std::map<std::string, GroupKey> & group_key_map,
    const std::vector<std::string> & atom_id_list,
    const std::string & class_key, const int par_id)
{
    if (IsModelObjectAvailable() == false)
    {
        return nullptr;
    }
    auto graph{ ROOTHelper::CreateGraphErrors() };
    auto count{ 0 };
    for (size_t i = 0; i < atom_id_list.size(); i++)
    {
        auto & atom_id{ atom_id_list[i] };
        if (group_key_map.find(atom_id) == group_key_map.end()) continue;
        auto & group_key{ group_key_map.at(atom_id) };
        if (IsAvailableAtomGroupKey(group_key, class_key) == false) continue;
        auto x_value{ static_cast<double>(i) };
        auto y_value{ m_model_object->GetAtomGroupPotentialEntry(class_key)->GetGausEstimatePrior(group_key, par_id) };
        auto y_error{ m_model_object->GetAtomGroupPotentialEntry(class_key)->GetGausVariancePrior(group_key, par_id) };
        graph->SetPoint(count, x_value, y_value);
        graph->SetPointError(count, 0.0, y_error);
        count++;
    }
    return graph;
}

std::unique_ptr<TGraphErrors> PotentialEntryIterator::CreateAtomGausEstimateScatterGraph(
    GroupKey group_key, const std::string & class_key, int par1_id, int par2_id, bool select_outliers)
{
    if (IsModelObjectAvailable() == false)
    {
        return nullptr;
    }
    auto atom_list{ GetAtomObjectList(group_key, class_key) };
    auto graph{ ROOTHelper::CreateGraphErrors() };
    auto count{ 0 };
    for (auto atom : atom_list)
    {
        auto entry{ atom->GetLocalPotentialEntry() };
        auto is_outlier{ entry->GetOutlierTag(class_key) };
        if (select_outliers == true && is_outlier == false) continue;
        graph->SetPoint(count, entry->GetGausEstimateMDPDE(par1_id), entry->GetGausEstimateMDPDE(par2_id));
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

std::unique_ptr<TGraphErrors> PotentialEntryIterator::CreateDistanceToMapValueGraph()
{
    auto graph{ ROOTHelper::CreateGraphErrors() };
    auto count{ 0 };
    for (auto & [distance, map_value] : GetDistanceAndMapValueList())
    {
        graph->SetPoint(count, distance, map_value);
        count++;
    }
    return graph;
}

std::unique_ptr<TGraphErrors> PotentialEntryIterator::CreateLinearModelDistanceToMapValueGraph()
{
    auto graph{ ROOTHelper::CreateGraphErrors() };
    auto count{ 0 };
    for (auto & [x, y] : GetLinearModelDistanceAndMapValueList())
    {
        graph->SetPoint(count, x, y);
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

std::unique_ptr<TF1> PotentialEntryIterator::CreateAtomLocalLinearModelFunctionOLS() const
{
    if (IsAtomLocalEntryAvailable() == false)
    {
        return nullptr;
    }
    auto beta_0{ m_atom_local_entry->GetBetaEstimateOLS(0) };
    auto beta_1{ m_atom_local_entry->GetBetaEstimateOLS(1) };
    return ROOTHelper::CreateLinearModelFunction("linear", beta_0, beta_1);
}

std::unique_ptr<TF1> PotentialEntryIterator::CreateAtomLocalLinearModelFunctionMDPDE() const
{
    if (IsAtomLocalEntryAvailable() == false)
    {
        return nullptr;
    }
    auto beta_0{ m_atom_local_entry->GetBetaEstimateMDPDE(0) };
    auto beta_1{ m_atom_local_entry->GetBetaEstimateMDPDE(1) };
    return ROOTHelper::CreateLinearModelFunction("linear", beta_0, beta_1);
}

std::unique_ptr<TF1> PotentialEntryIterator::CreateAtomLocalGausFunctionOLS() const
{
    if (IsAtomLocalEntryAvailable() == false)
    {
        return nullptr;
    }
    auto amplitude{ m_atom_local_entry->GetGausEstimateOLS(0) };
    auto width{ m_atom_local_entry->GetGausEstimateOLS(1) };
    return ROOTHelper::CreateGaus3DFunctionIn1D("gaus", amplitude, width);
}

std::unique_ptr<TF1> PotentialEntryIterator::CreateAtomLocalGausFunctionMDPDE() const
{
    if (IsAtomLocalEntryAvailable() == false)
    {
        return nullptr;
    }
    auto amplitude{ m_atom_local_entry->GetGausEstimateMDPDE(0) };
    auto width{ m_atom_local_entry->GetGausEstimateMDPDE(1) };
    return ROOTHelper::CreateGaus3DFunctionIn1D("gaus", amplitude, width);
}

std::unique_ptr<TF1> PotentialEntryIterator::CreateAtomGroupLinearModelFunctionMean(
    GroupKey group_key, const std::string & class_key, double x_min, double x_max) const
{
    auto mu_0
    {
        m_model_object->GetAtomGroupPotentialEntry(class_key)->GetMuEstimateMean(group_key, 0)
    };
    auto mu_1
    {
        m_model_object->GetAtomGroupPotentialEntry(class_key)->GetMuEstimateMean(group_key, 1)
    };
    return ROOTHelper::CreateLinearModelFunction("linear_mean", mu_0, mu_1, x_min, x_max);
}

std::unique_ptr<TF1> PotentialEntryIterator::CreateAtomGroupLinearModelFunctionPrior(
    GroupKey group_key, const std::string & class_key, double x_min, double x_max) const
{
    auto mu_0
    {
        m_model_object->GetAtomGroupPotentialEntry(class_key)->GetMuEstimatePrior(group_key, 0)
    };
    auto mu_1
    {
        m_model_object->GetAtomGroupPotentialEntry(class_key)->GetMuEstimatePrior(group_key, 1)
    };
    return ROOTHelper::CreateLinearModelFunction("linear_prior", mu_0, mu_1, x_min, x_max);
}

std::unique_ptr<TF1> PotentialEntryIterator::CreateAtomGroupGausFunctionMean(
    GroupKey group_key, const std::string & class_key) const
{
    auto amplitude
    {
        m_model_object->GetAtomGroupPotentialEntry(class_key)->GetGausEstimateMean(group_key, 0)
    };
    auto width
    {
        m_model_object->GetAtomGroupPotentialEntry(class_key)->GetGausEstimateMean(group_key, 1)
    };
    return ROOTHelper::CreateGaus3DFunctionIn1D("group_gaus_mean", amplitude, width);
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
    return ROOTHelper::CreateGaus3DFunctionIn1D("group_gaus_prior", amplitude, width);
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
    return ROOTHelper::CreateGaus2DFunctionIn1D("group_gaus_prior", amplitude, width);
}
#endif
} // namespace rhbm_gem
