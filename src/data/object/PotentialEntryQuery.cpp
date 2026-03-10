#include <rhbm_gem/data/object/PotentialEntryQuery.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/BondObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/data/object/LocalPotentialEntry.hpp>
#include <rhbm_gem/data/object/GroupPotentialEntry.hpp>
#include <rhbm_gem/utils/domain/ChemicalDataHelper.hpp>
#include <rhbm_gem/utils/math/ArrayStats.hpp>
#include <rhbm_gem/utils/domain/KeyPacker.hpp>
#include <rhbm_gem/utils/domain/Constants.hpp>
#include <rhbm_gem/utils/domain/GlobalEnumClass.hpp>
#include <rhbm_gem/data/object/AtomClassifier.hpp>
#include <rhbm_gem/data/object/BondClassifier.hpp>
#include <rhbm_gem/utils/math/GausLinearTransformHelper.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>

#include <Eigen/Dense>
#include <stdexcept>
#include <sstream>
#include <cmath>
#include <map>


namespace rhbm_gem {

PotentialEntryQuery::PotentialEntryQuery(ModelObject * model_object) :
    m_atom_object{ nullptr }, m_bond_object{ nullptr }, m_model_object{ model_object },
    m_atom_local_entry{ nullptr }, m_bond_local_entry{ nullptr }
{
}

PotentialEntryQuery::PotentialEntryQuery(AtomObject * atom_object) :
    m_atom_object{ atom_object }, m_bond_object{ nullptr }, m_model_object{ nullptr },
    m_bond_local_entry{ nullptr }
{
    m_atom_local_entry = atom_object->GetLocalPotentialEntry();
}

PotentialEntryQuery::PotentialEntryQuery(BondObject * bond_object) :
    m_atom_object{ nullptr }, m_bond_object{ bond_object }, m_model_object{ nullptr },
    m_atom_local_entry{ nullptr }
{
    m_bond_local_entry = bond_object->GetLocalPotentialEntry();
}

PotentialEntryQuery::~PotentialEntryQuery()
{

}

double PotentialEntryQuery::GetAtomGausEstimateMinimum(int par_id, Element element) const
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

double PotentialEntryQuery::GetBondGausEstimateMinimum(int par_id) const
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

bool PotentialEntryQuery::IsOutlierAtom(const std::string & class_key) const
{
    if (IsAtomLocalEntryAvailable() == false)
    {
        throw std::runtime_error("Atom local entry is not available.");
    }
    return m_atom_local_entry->GetOutlierTag(class_key);
}

bool PotentialEntryQuery::IsOutlierBond(const std::string & class_key) const
{
    if (IsBondLocalEntryAvailable() == false)
    {
        throw std::runtime_error("Bond local entry is not available.");
    }
    return m_bond_local_entry->GetOutlierTag(class_key);
}

bool PotentialEntryQuery::IsAvailableAtomGroupKey(
    GroupKey group_key, const std::string & class_key, bool varbose) const
{
    return CheckAtomGroupKey(group_key, class_key, varbose);
}

bool PotentialEntryQuery::IsAvailableBondGroupKey(
    GroupKey group_key, const std::string & class_key, bool varbose) const
{
    return CheckBondGroupKey(group_key, class_key, varbose);
}

size_t PotentialEntryQuery::GetAtomResidueCount(
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

size_t PotentialEntryQuery::GetBondResidueCount(
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

double PotentialEntryQuery::GetAtomGausEstimateMean(
    GroupKey group_key, const std::string & class_key, int par_id) const
{
    if (CheckAtomGroupKey(group_key, class_key) == false)
    {
        throw std::runtime_error("Group key is not available.");
    }
    return m_model_object->GetAtomGroupPotentialEntry(class_key)->GetGausEstimateMean(group_key, par_id);
}

double PotentialEntryQuery::GetAtomGausEstimatePrior(
    GroupKey group_key, const std::string & class_key, int par_id) const
{
    if (CheckAtomGroupKey(group_key, class_key) == false)
    {
        throw std::runtime_error("Group key is not available.");
    }
    return m_model_object->GetAtomGroupPotentialEntry(class_key)->GetGausEstimatePrior(group_key, par_id);
}

double PotentialEntryQuery::GetBondGausEstimatePrior(
    GroupKey group_key, const std::string & class_key, int par_id) const
{
    if (CheckBondGroupKey(group_key, class_key) == false)
    {
        throw std::runtime_error("Group key is not available.");
    }
    return m_model_object->GetBondGroupPotentialEntry(class_key)->GetGausEstimatePrior(group_key, par_id);
}

double PotentialEntryQuery::GetAtomGausVariancePrior(
    GroupKey group_key, const std::string & class_key, int par_id) const
{
    if (CheckAtomGroupKey(group_key, class_key) == false)
    {
        throw std::runtime_error("Group key is not available.");
    }
    return m_model_object->GetAtomGroupPotentialEntry(class_key)->GetGausVariancePrior(group_key, par_id);
}

double PotentialEntryQuery::GetBondGausVariancePrior(
    GroupKey group_key, const std::string & class_key, int par_id) const
{
    if (CheckBondGroupKey(group_key, class_key) == false)
    {
        throw std::runtime_error("Group key is not available.");
    }
    return m_model_object->GetBondGroupPotentialEntry(class_key)->GetGausVariancePrior(group_key, par_id);
}

const std::vector<AtomObject *> & PotentialEntryQuery::GetAtomObjectList(
    GroupKey group_key, const std::string & class_key) const
{
    if (CheckAtomGroupKey(group_key, class_key) == false)
    {
        throw std::runtime_error("Atom group key is not available.");
    }
    return m_model_object->GetAtomGroupPotentialEntry(class_key)->GetAtomObjectPtrList(group_key);
}

const std::vector<BondObject *> & PotentialEntryQuery::GetBondObjectList(
    GroupKey group_key, const std::string & class_key) const
{
    if (CheckBondGroupKey(group_key, class_key) == false)
    {
        throw std::runtime_error("Bond group key is not available.");
    }
    return m_model_object->GetBondGroupPotentialEntry(class_key)->GetBondObjectPtrList(group_key);
}

std::vector<AtomObject *> PotentialEntryQuery::GetOutlierAtomObjectList(
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

std::unordered_map<int, AtomObject *> PotentialEntryQuery::GetAtomObjectMap(
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
PotentialEntryQuery::GetLinearModelDistanceAndMapValueList() const
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
PotentialEntryQuery::GetDistanceAndMapValueList() const
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

std::vector<std::tuple<float, float>> PotentialEntryQuery::GetBinnedDistanceAndMapValueList(
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

std::tuple<float, float> PotentialEntryQuery::GetDistanceRange(double margin_rate) const
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

std::tuple<float, float> PotentialEntryQuery::GetMapValueRange(double margin_rate) const
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

std::tuple<float, float> PotentialEntryQuery::GetBinnedMapValueRange(
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

double PotentialEntryQuery::GetAmplitudeEstimateMDPDE() const
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

double PotentialEntryQuery::GetAmplitudeEstimatePosterior(const std::string & class_key) const
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

double PotentialEntryQuery::GetAmplitudeVariancePosterior(const std::string & class_key) const
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

double PotentialEntryQuery::GetWidthEstimateMDPDE() const
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

double PotentialEntryQuery::GetWidthEstimatePosterior(const std::string & class_key) const
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

double PotentialEntryQuery::GetWidthVariancePosterior(const std::string & class_key) const
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

double PotentialEntryQuery::GetAlphaR() const
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

double PotentialEntryQuery::GetAtomAlphaR(
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

double PotentialEntryQuery::GetBondAlphaR(
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

double PotentialEntryQuery::GetAtomAlphaG(
    GroupKey group_key, const std::string & class_key) const
{
    if (CheckAtomGroupKey(group_key, class_key) == false)
    {
        throw std::runtime_error("Group key is not available.");
    }
    return m_model_object->GetAtomGroupPotentialEntry(class_key)->GetAlphaG(group_key);
}

double PotentialEntryQuery::GetBondAlphaG(
    GroupKey group_key, const std::string & class_key) const
{
    if (CheckBondGroupKey(group_key, class_key) == false)
    {
        throw std::runtime_error("Group key is not available.");
    }
    return m_model_object->GetBondGroupPotentialEntry(class_key)->GetAlphaG(group_key);
}

bool PotentialEntryQuery::IsAtomObjectAvailable() const
{
    if (m_atom_object == nullptr)
    {
        Logger::Log(LogLevel::Error, "Atom object is not available.");
        return false;
    }
    return true;
}

bool PotentialEntryQuery::IsBondObjectAvailable() const
{
    if (m_bond_object == nullptr)
    {
        Logger::Log(LogLevel::Error, "Bond object is not available.");
        return false;
    }
    return true;
}

bool PotentialEntryQuery::IsAtomLocalEntryAvailable() const
{
    if (m_atom_local_entry == nullptr)
    {
        return false;
    }
    return true;
}

bool PotentialEntryQuery::IsBondLocalEntryAvailable() const
{
    if (m_bond_local_entry == nullptr)
    {
        return false;
    }
    return true;
}

bool PotentialEntryQuery::IsModelObjectAvailable() const
{
    if (m_model_object == nullptr)
    {
        Logger::Log(LogLevel::Error, "Model object is not available.");
        return false;
    }
    return true;
}

bool PotentialEntryQuery::CheckAtomGroupKey(
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

bool PotentialEntryQuery::CheckBondGroupKey(
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

Residue PotentialEntryQuery::GetResidueFromAtomGroupKey(
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

Residue PotentialEntryQuery::GetResidueFromBondGroupKey(
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

AtomObject * PotentialEntryQuery::GetAtomObject() const
{
    return m_atom_object;
}

BondObject * PotentialEntryQuery::GetBondObject() const
{
    return m_bond_object;
}

ModelObject * PotentialEntryQuery::GetModelObject() const
{
    return m_model_object;
}

LocalPotentialEntry * PotentialEntryQuery::GetAtomLocalEntry() const
{
    return m_atom_local_entry;
}

LocalPotentialEntry * PotentialEntryQuery::GetBondLocalEntry() const
{
    return m_bond_local_entry;
}

} // namespace rhbm_gem
