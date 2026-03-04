#include "LegacyModelObjectReader.hpp"

#include "AtomClassifier.hpp"
#include "AtomObject.hpp"
#include "BondClassifier.hpp"
#include "BondObject.hpp"
#include "ChemicalComponentEntry.hpp"
#include "ChemicalDataHelper.hpp"
#include "GroupPotentialEntry.hpp"
#include "LocalPotentialEntry.hpp"
#include "ModelObject.hpp"
#include "SQLiteWrapper.hpp"

#include <cctype>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>

namespace
{
    std::string FormatSQL(std::string_view templates, const std::string & table)
    {
        auto pos{ templates.find("{}") };
        if (pos == std::string_view::npos)
        {
            return std::string(templates);
        }
        std::string result;
        result.reserve(templates.size() - 2 + table.size());
        result.append(templates.substr(0, pos));
        result.append(table);
        result.append(templates.substr(pos + 2));
        return result;
    }

    constexpr std::string_view SELECT_MODEL_LIST_SQL = R"sql(
        SELECT
            key_tag, atom_size, pdb_id, emd_id, map_resolution, resolution_method
        FROM {} WHERE key_tag = ? LIMIT 1;
    )sql";

    constexpr std::string_view SELECT_ATOM_LIST_SQL = R"sql(
        SELECT
            serial_id, sequence_id, component_id, atom_id, chain_id, indicator,
            occupancy, temperature, element, structure, is_special_atom,
            position_x, position_y, position_z, component_key, atom_key
        FROM {};
    )sql";

    constexpr std::string_view SELECT_BOND_LIST_SQL = R"sql(
        SELECT
            atom_serial_id_1, atom_serial_id_2,
            bond_key, bond_type, bond_order, is_special_bond
        FROM {};
    )sql";

    constexpr std::string_view SELECT_COMPONENT_ENTRY_LIST_SQL = R"sql(
        SELECT
            component_key, id, name, type, formula,
            molecular_weight, is_standard_monomer
        FROM {};
    )sql";

    constexpr std::string_view SELECT_COMPONENT_ATOM_ENTRY_LIST_SQL = R"sql(
        SELECT
            component_key, atom_key, atom_id,
            element_type, aromatic_atom_flag, stereo_config
        FROM {};
    )sql";

    constexpr std::string_view SELECT_COMPONENT_BOND_ENTRY_LIST_SQL = R"sql(
        SELECT
            component_key, bond_key, bond_id,
            bond_type, bond_order, aromatic_atom_flag, stereo_config
        FROM {};
    )sql";

    constexpr std::string_view SELECT_ATOM_LOCAL_ENTRY_SQL = R"sql(
        SELECT
            serial_id, sampling_size, distance_and_map_value_list,
            amplitude_estimate_ols, width_estimate_ols,
            amplitude_estimate_mdpde, width_estimate_mdpde,
            alpha_r
        FROM {};
    )sql";

    constexpr std::string_view SELECT_ATOM_LOCAL_SUBLIST_SQL = R"sql(
        SELECT
            serial_id, amplitude_estimate_posterior, width_estimate_posterior,
            amplitude_variance_posterior, width_variance_posterior,
            outlier_tag, statistical_distance
        FROM {};
    )sql";

    constexpr std::string_view SELECT_BOND_LOCAL_ENTRY_SQL = R"sql(
        SELECT
            atom_serial_id_1, atom_serial_id_2,
            sampling_size, distance_and_map_value_list,
            amplitude_estimate_ols, width_estimate_ols,
            amplitude_estimate_mdpde, width_estimate_mdpde,
            alpha_r
        FROM {};
    )sql";

    constexpr std::string_view SELECT_BOND_LOCAL_SUBLIST_SQL = R"sql(
        SELECT
            atom_serial_id_1, atom_serial_id_2,
            amplitude_estimate_posterior, width_estimate_posterior,
            amplitude_variance_posterior, width_variance_posterior,
            outlier_tag, statistical_distance
        FROM {};
    )sql";

    constexpr std::string_view SELECT_GROUP_ENTRY_SQL = R"sql(
        SELECT
            key_id, member_size,
            amplitude_estimate_mean, width_estimate_mean,
            amplitude_estimate_mdpde, width_estimate_mdpde,
            amplitude_estimate_prior, width_estimate_prior,
            amplitude_variance_prior, width_variance_prior,
            alpha_g
        FROM {};
    )sql";

    constexpr std::string_view SELECT_ROW_COUNT_SQL = "SELECT COUNT(*) FROM {};";
    constexpr std::string_view TABLE_EXISTS_SQL =
        "SELECT name FROM sqlite_master WHERE type = 'table' AND name = ? LIMIT 1;";
}

namespace rhbm_gem {

LegacyModelObjectReader::LegacyModelObjectReader(SQLiteWrapper * database) :
    m_database{ database }, m_table_cache{}
{
    if (m_database == nullptr)
    {
        throw std::runtime_error("LegacyModelObjectReader requires a valid database handle.");
    }
}

std::unique_ptr<ModelObject> LegacyModelObjectReader::Load(const std::string & key_tag)
{
    auto sanitized_key_tag{ SanitizeTableName(key_tag) };
    auto model_object{ std::make_unique<ModelObject>() };

    LoadChemicalComponentEntryList(
        model_object.get(), "chemical_component_entry_in_" + sanitized_key_tag);
    LoadComponentAtomEntryList(
        model_object.get(), "component_atom_entry_in_" + sanitized_key_tag);
    LoadComponentBondEntryList(
        model_object.get(), "component_bond_entry_in_" + sanitized_key_tag);
    model_object->SetAtomList(LoadAtomObjectList(sanitized_key_tag));
    model_object->SetBondList(LoadBondObjectList(sanitized_key_tag, model_object.get()));

    m_database->Prepare(FormatSQL(SELECT_MODEL_LIST_SQL, "model_list"));
    SQLiteWrapper::StatementGuard guard(*m_database);
    m_database->Bind<std::string>(1, key_tag);
    const auto return_code{ m_database->StepNext() };
    int atom_size{ 0 };
    if (return_code == SQLiteWrapper::StepRow())
    {
        model_object->SetKeyTag(m_database->GetColumn<std::string>(0));
        model_object->SetPdbID(m_database->GetColumn<std::string>(2));
        model_object->SetEmdID(m_database->GetColumn<std::string>(3));
        model_object->SetResolution(m_database->GetColumn<double>(4));
        model_object->SetResolutionMethod(m_database->GetColumn<std::string>(5));
        atom_size = m_database->GetColumn<int>(1);
    }
    else if (return_code == SQLiteWrapper::StepDone())
    {
        throw std::runtime_error("Cannot find the row with key_tag = " + key_tag);
    }
    else
    {
        throw std::runtime_error("Step failed: " + m_database->ErrorMessage());
    }

    if (atom_size != static_cast<int>(model_object->GetNumberOfAtom()))
    {
        throw std::runtime_error(
            "The number of atoms in the model object does not match the database record.");
    }

    for (size_t i = 0; i < ChemicalDataHelper::GetGroupAtomClassCount(); i++)
    {
        const auto group_class_key{ ChemicalDataHelper::GetGroupAtomClassKey(i) };
        auto group_potential_entry{ std::make_unique<GroupPotentialEntry>() };
        model_object->AddAtomGroupPotentialEntry(group_class_key, group_potential_entry);
        LoadAtomGroupPotentialEntryList(
            model_object.get(),
            group_class_key,
            group_class_key + "_atom_group_potential_entry_in_" + sanitized_key_tag);
    }

    for (size_t i = 0; i < ChemicalDataHelper::GetGroupBondClassCount(); i++)
    {
        const auto group_class_key{ ChemicalDataHelper::GetGroupBondClassKey(i) };
        auto group_potential_entry{ std::make_unique<GroupPotentialEntry>() };
        model_object->AddBondGroupPotentialEntry(group_class_key, group_potential_entry);
        LoadBondGroupPotentialEntryList(
            model_object.get(),
            group_class_key,
            group_class_key + "_bond_group_potential_entry_in_" + sanitized_key_tag);
    }
    return model_object;
}

std::vector<std::unique_ptr<AtomObject>> LegacyModelObjectReader::LoadAtomObjectList(
    const std::string & key_tag)
{
    const auto sanitized_key_tag{ SanitizeTableName(key_tag) };
    const auto atom_list_table_name{ "atom_list_in_" + sanitized_key_tag };
    if (!TableExists(atom_list_table_name))
    {
        throw std::runtime_error("Required table not found: " + atom_list_table_name);
    }

    const auto count_row_vec{
        m_database->Query<int>(FormatSQL(SELECT_ROW_COUNT_SQL, atom_list_table_name)) };
    const auto atom_count{
        count_row_vec.empty() ? 0 : std::get<0>(count_row_vec.front()) };

    std::unordered_map<int, std::unique_ptr<LocalPotentialEntry>> local_potential_entry_map;
    const auto local_potential_entry_table_name{
        "atom_local_potential_entry_in_" + sanitized_key_tag };
    if (TableExists(local_potential_entry_table_name))
    {
        local_potential_entry_map =
            LoadAtomLocalPotentialEntryMap(local_potential_entry_table_name);
    }

    std::vector<std::unique_ptr<AtomObject>> atom_object_list;
    atom_object_list.reserve(static_cast<size_t>(atom_count));
    auto iter{
        m_database->IterateQuery<
            int, int, std::string, std::string, std::string, std::string,
            double, double, int, int, int, double, double, double, ComponentKey, AtomKey>(
            FormatSQL(SELECT_ATOM_LIST_SQL, atom_list_table_name)) };
    std::tuple<int, int, std::string, std::string, std::string, std::string,
               double, double, int, int, int, double, double, double, ComponentKey, AtomKey> row;
    while (iter.Next(row))
    {
        auto atom_object{ std::make_unique<AtomObject>() };
        atom_object->SetSerialID(std::get<0>(row));
        atom_object->SetSequenceID(std::get<1>(row));
        atom_object->SetComponentID(std::get<2>(row));
        atom_object->SetAtomID(std::get<3>(row));
        atom_object->SetChainID(std::get<4>(row));
        atom_object->SetIndicator(std::get<5>(row));
        atom_object->SetOccupancy(static_cast<float>(std::get<6>(row)));
        atom_object->SetTemperature(static_cast<float>(std::get<7>(row)));
        atom_object->SetElement(static_cast<Element>(std::get<8>(row)));
        atom_object->SetStructure(static_cast<Structure>(std::get<9>(row)));
        atom_object->SetSpecialAtomFlag(static_cast<bool>(std::get<10>(row)));
        atom_object->SetPosition(
            static_cast<float>(std::get<11>(row)),
            static_cast<float>(std::get<12>(row)),
            static_cast<float>(std::get<13>(row)));
        atom_object->SetComponentKey(std::get<14>(row));
        atom_object->SetAtomKey(std::get<15>(row));

        const auto serial_id{ atom_object->GetSerialID() };
        if (!local_potential_entry_map.empty())
        {
            auto entry_iter{ local_potential_entry_map.find(serial_id) };
            atom_object->SetSelectedFlag(entry_iter != local_potential_entry_map.end());
            if (entry_iter != local_potential_entry_map.end())
            {
                atom_object->AddLocalPotentialEntry(std::move(entry_iter->second));
            }
        }

        atom_object_list.emplace_back(std::move(atom_object));
    }
    return atom_object_list;
}

std::vector<std::unique_ptr<BondObject>> LegacyModelObjectReader::LoadBondObjectList(
    const std::string & key_tag, const ModelObject * model_object)
{
    const auto sanitized_key_tag{ SanitizeTableName(key_tag) };
    const auto bond_list_table_name{ "bond_list_in_" + sanitized_key_tag };
    if (!TableExists(bond_list_table_name))
    {
        throw std::runtime_error("Required table not found: " + bond_list_table_name);
    }

    const auto count_row_vec{
        m_database->Query<int>(FormatSQL(SELECT_ROW_COUNT_SQL, bond_list_table_name)) };
    const auto bond_count{
        count_row_vec.empty() ? 0 : std::get<0>(count_row_vec.front()) };

    std::map<std::pair<int, int>, std::unique_ptr<LocalPotentialEntry>> local_potential_entry_map;
    const auto local_potential_entry_table_name{
        "bond_local_potential_entry_in_" + sanitized_key_tag };
    if (TableExists(local_potential_entry_table_name))
    {
        local_potential_entry_map =
            LoadBondLocalPotentialEntryMap(local_potential_entry_table_name);
    }

    std::vector<std::unique_ptr<BondObject>> bond_object_list;
    bond_object_list.reserve(static_cast<size_t>(bond_count));
    auto iter{
        m_database->IterateQuery<int, int, BondKey, int, int, int>(
            FormatSQL(SELECT_BOND_LIST_SQL, bond_list_table_name)) };
    std::tuple<int, int, BondKey, int, int, int> row;
    while (iter.Next(row))
    {
        auto atom_object_1{ model_object->GetAtomPtr(std::get<0>(row)) };
        auto atom_object_2{ model_object->GetAtomPtr(std::get<1>(row)) };
        auto bond_object{ std::make_unique<BondObject>(atom_object_1, atom_object_2) };
        bond_object->SetBondKey(std::get<2>(row));
        bond_object->SetBondType(static_cast<BondType>(std::get<3>(row)));
        bond_object->SetBondOrder(static_cast<BondOrder>(std::get<4>(row)));
        bond_object->SetSpecialBondFlag(static_cast<bool>(std::get<5>(row)));

        const auto atom_serial_id_pair{
            std::make_pair(atom_object_1->GetSerialID(), atom_object_2->GetSerialID()) };
        if (!local_potential_entry_map.empty())
        {
            auto entry_iter{ local_potential_entry_map.find(atom_serial_id_pair) };
            bond_object->SetSelectedFlag(entry_iter != local_potential_entry_map.end());
            if (entry_iter != local_potential_entry_map.end())
            {
                bond_object->AddLocalPotentialEntry(std::move(entry_iter->second));
            }
        }

        bond_object_list.emplace_back(std::move(bond_object));
    }
    return bond_object_list;
}

void LegacyModelObjectReader::LoadChemicalComponentEntryList(
    ModelObject * model_obj, const std::string & table_name)
{
    if (!TableExists(table_name))
    {
        return;
    }

    auto iter{
        m_database->IterateQuery<int, std::string, std::string, std::string, std::string, double, int>(
            FormatSQL(SELECT_COMPONENT_ENTRY_LIST_SQL, table_name)) };
    std::tuple<int, std::string, std::string, std::string, std::string, double, int> row;
    while (iter.Next(row))
    {
        auto component_entry{ std::make_unique<ChemicalComponentEntry>() };
        const auto component_key{ static_cast<ComponentKey>(std::get<0>(row)) };
        const auto component_id{ std::get<1>(row) };
        component_entry->SetComponentId(component_id);
        component_entry->SetComponentName(std::get<2>(row));
        component_entry->SetComponentType(std::get<3>(row));
        component_entry->SetComponentFormula(std::get<4>(row));
        component_entry->SetComponentMolecularWeight(static_cast<float>(std::get<5>(row)));
        component_entry->SetStandardMonomerFlag(static_cast<bool>(std::get<6>(row)));
        model_obj->AddChemicalComponentEntry(component_key, std::move(component_entry));
        model_obj->GetComponentKeySystemPtr()->RegisterComponent(component_id, component_key);
    }
}

void LegacyModelObjectReader::LoadComponentAtomEntryList(
    ModelObject * model_obj, const std::string & table_name)
{
    if (!TableExists(table_name))
    {
        return;
    }

    auto iter{
        m_database->IterateQuery<ComponentKey, AtomKey, std::string, int, int, int>(
            FormatSQL(SELECT_COMPONENT_ATOM_ENTRY_LIST_SQL, table_name)) };
    std::tuple<ComponentKey, AtomKey, std::string, int, int, int> row;
    while (iter.Next(row))
    {
        const auto component_key{ std::get<0>(row) };
        if (model_obj->GetChemicalComponentEntryMap().find(component_key)
            == model_obj->GetChemicalComponentEntryMap().end())
        {
            continue;
        }

        auto & component_entry{ model_obj->GetChemicalComponentEntryMap().at(component_key) };
        const auto atom_key{ std::get<1>(row) };
        const auto atom_id{ std::get<2>(row) };
        ComponentAtomEntry atom_entry;
        atom_entry.atom_id = atom_id;
        atom_entry.element_type = static_cast<Element>(std::get<3>(row));
        atom_entry.aromatic_atom_flag = static_cast<bool>(std::get<4>(row));
        atom_entry.stereo_config = static_cast<StereoChemistry>(std::get<5>(row));
        component_entry->AddComponentAtomEntry(atom_key, atom_entry);
        model_obj->GetAtomKeySystemPtr()->RegisterAtom(atom_id, atom_key);
    }
}

void LegacyModelObjectReader::LoadComponentBondEntryList(
    ModelObject * model_obj, const std::string & table_name)
{
    if (!TableExists(table_name))
    {
        return;
    }

    auto iter{
        m_database->IterateQuery<ComponentKey, BondKey, std::string, int, int, int, int>(
            FormatSQL(SELECT_COMPONENT_BOND_ENTRY_LIST_SQL, table_name)) };
    std::tuple<ComponentKey, BondKey, std::string, int, int, int, int> row;
    while (iter.Next(row))
    {
        const auto component_key{ std::get<0>(row) };
        if (model_obj->GetChemicalComponentEntryMap().find(component_key)
            == model_obj->GetChemicalComponentEntryMap().end())
        {
            continue;
        }

        auto & component_entry{ model_obj->GetChemicalComponentEntryMap().at(component_key) };
        const auto bond_key{ std::get<1>(row) };
        const auto bond_id{ std::get<2>(row) };
        ComponentBondEntry bond_entry;
        bond_entry.bond_id = bond_id;
        bond_entry.bond_type = static_cast<BondType>(std::get<3>(row));
        bond_entry.bond_order = static_cast<BondOrder>(std::get<4>(row));
        bond_entry.aromatic_atom_flag = static_cast<bool>(std::get<5>(row));
        bond_entry.stereo_config = static_cast<StereoChemistry>(std::get<6>(row));
        component_entry->AddComponentBondEntry(bond_key, bond_entry);
        model_obj->GetBondKeySystemPtr()->RegisterBond(bond_id, bond_key);
    }
}

std::unordered_map<int, std::unique_ptr<LocalPotentialEntry>>
LegacyModelObjectReader::LoadAtomLocalPotentialEntryMap(const std::string & table_name)
{
    if (!TableExists(table_name))
    {
        return {};
    }

    const auto count_row_vec{
        m_database->Query<int>(FormatSQL(SELECT_ROW_COUNT_SQL, table_name)) };
    const auto entry_count{
        count_row_vec.empty() ? 0 : std::get<0>(count_row_vec.front()) };

    std::unordered_map<int, std::unique_ptr<LocalPotentialEntry>> local_potential_entry_map;
    local_potential_entry_map.reserve(static_cast<size_t>(entry_count));
    auto iter{
        m_database->IterateQuery<
            int, int, std::vector<std::tuple<float, float>>, double, double, double, double, double>(
            FormatSQL(SELECT_ATOM_LOCAL_ENTRY_SQL, table_name)) };
    std::tuple<int, int, std::vector<std::tuple<float, float>>,
               double, double, double, double, double> row;
    while (iter.Next(row))
    {
        auto local_potential_entry{ std::make_unique<LocalPotentialEntry>() };
        const auto serial_id{ std::get<0>(row) };
        local_potential_entry->AddDistanceAndMapValueList(std::move(std::get<2>(row)));
        local_potential_entry->AddGausEstimateOLS(std::get<3>(row), std::get<4>(row));
        local_potential_entry->AddGausEstimateMDPDE(std::get<5>(row), std::get<6>(row));
        local_potential_entry->SetAlphaR(std::get<7>(row));
        local_potential_entry_map[serial_id] = std::move(local_potential_entry);
    }

    for (size_t i = 0; i < ChemicalDataHelper::GetGroupAtomClassCount(); i++)
    {
        const auto class_key{ ChemicalDataHelper::GetGroupAtomClassKey(i) };
        LoadAtomLocalPotentialEntrySubList(class_key + "_" + table_name, class_key, local_potential_entry_map);
    }
    return local_potential_entry_map;
}

void LegacyModelObjectReader::LoadAtomLocalPotentialEntrySubList(
    const std::string & table_name,
    const std::string & class_key,
    std::unordered_map<int, std::unique_ptr<LocalPotentialEntry>> & entry_map)
{
    if (!TableExists(table_name))
    {
        return;
    }

    auto iter{
        m_database->IterateQuery<int, double, double, double, double, int, double>(
            FormatSQL(SELECT_ATOM_LOCAL_SUBLIST_SQL, table_name)) };
    std::tuple<int, double, double, double, double, int, double> row;
    while (iter.Next(row))
    {
        const auto serial_id{ std::get<0>(row) };
        auto entry_iter{ entry_map.find(serial_id) };
        if (entry_iter == entry_map.end())
        {
            continue;
        }
        auto & entry{ entry_iter->second };
        entry->AddGausEstimatePosterior(class_key, std::get<1>(row), std::get<2>(row));
        entry->AddGausVariancePosterior(class_key, std::get<3>(row), std::get<4>(row));
        entry->AddOutlierTag(class_key, static_cast<bool>(std::get<5>(row)));
        entry->AddStatisticalDistance(class_key, std::get<6>(row));
    }
}

std::map<std::pair<int, int>, std::unique_ptr<LocalPotentialEntry>>
LegacyModelObjectReader::LoadBondLocalPotentialEntryMap(const std::string & table_name)
{
    if (!TableExists(table_name))
    {
        return {};
    }

    std::map<std::pair<int, int>, std::unique_ptr<LocalPotentialEntry>> local_potential_entry_map;
    auto iter{
        m_database->IterateQuery<
            int, int, int, std::vector<std::tuple<float, float>>, double, double, double, double, double>(
            FormatSQL(SELECT_BOND_LOCAL_ENTRY_SQL, table_name)) };
    std::tuple<int, int, int,
               std::vector<std::tuple<float, float>>, double, double, double, double, double> row;
    while (iter.Next(row))
    {
        auto local_potential_entry{ std::make_unique<LocalPotentialEntry>() };
        const auto atom_serial_id_pair{ std::make_pair(std::get<0>(row), std::get<1>(row)) };
        local_potential_entry->AddDistanceAndMapValueList(std::move(std::get<3>(row)));
        local_potential_entry->AddGausEstimateOLS(std::get<4>(row), std::get<5>(row));
        local_potential_entry->AddGausEstimateMDPDE(std::get<6>(row), std::get<7>(row));
        local_potential_entry->SetAlphaR(std::get<8>(row));
        local_potential_entry_map[atom_serial_id_pair] = std::move(local_potential_entry);
    }

    for (size_t i = 0; i < ChemicalDataHelper::GetGroupBondClassCount(); i++)
    {
        const auto class_key{ ChemicalDataHelper::GetGroupBondClassKey(i) };
        LoadBondLocalPotentialEntrySubList(class_key + "_" + table_name, class_key, local_potential_entry_map);
    }
    return local_potential_entry_map;
}

void LegacyModelObjectReader::LoadBondLocalPotentialEntrySubList(
    const std::string & table_name,
    const std::string & class_key,
    std::map<std::pair<int, int>, std::unique_ptr<LocalPotentialEntry>> & entry_map)
{
    if (!TableExists(table_name))
    {
        return;
    }

    auto iter{
        m_database->IterateQuery<int, int, double, double, double, double, int, double>(
            FormatSQL(SELECT_BOND_LOCAL_SUBLIST_SQL, table_name)) };
    std::tuple<int, int, double, double, double, double, int, double> row;
    while (iter.Next(row))
    {
        const auto atom_serial_id_pair{ std::make_pair(std::get<0>(row), std::get<1>(row)) };
        auto entry_iter{ entry_map.find(atom_serial_id_pair) };
        if (entry_iter == entry_map.end())
        {
            continue;
        }
        auto & entry{ entry_iter->second };
        entry->AddGausEstimatePosterior(class_key, std::get<2>(row), std::get<3>(row));
        entry->AddGausVariancePosterior(class_key, std::get<4>(row), std::get<5>(row));
        entry->AddOutlierTag(class_key, static_cast<bool>(std::get<6>(row)));
        entry->AddStatisticalDistance(class_key, std::get<7>(row));
    }
}

void LegacyModelObjectReader::LoadAtomGroupPotentialEntryList(
    ModelObject * model_obj, const std::string & class_key, const std::string & table_name)
{
    if (!TableExists(table_name))
    {
        return;
    }

    auto group_entry{ model_obj->GetAtomGroupPotentialEntry(class_key) };
    auto iter{
        m_database->IterateQuery<GroupKey, int, double, double, double, double, double, double, double, double, double>(
            FormatSQL(SELECT_GROUP_ENTRY_SQL, table_name)) };
    std::tuple<GroupKey, int, double, double, double, double, double, double, double, double, double> row;
    while (iter.Next(row))
    {
        const auto group_key{ std::get<0>(row) };
        group_entry->InsertGroupKey(group_key);
        group_entry->ReserveAtomObjectPtrList(group_key, std::get<1>(row));
        group_entry->AddGausEstimateMean(group_key, std::get<2>(row), std::get<3>(row));
        group_entry->AddGausEstimateMDPDE(group_key, std::get<4>(row), std::get<5>(row));
        group_entry->AddGausEstimatePrior(group_key, std::get<6>(row), std::get<7>(row));
        group_entry->AddGausVariancePrior(group_key, std::get<8>(row), std::get<9>(row));
        group_entry->AddAlphaG(group_key, std::get<10>(row));
    }

    for (auto & atom : model_obj->GetSelectedAtomList())
    {
        group_entry->AddAtomObjectPtr(AtomClassifier::GetGroupKeyInClass(atom, class_key), atom);
    }
}

void LegacyModelObjectReader::LoadBondGroupPotentialEntryList(
    ModelObject * model_obj, const std::string & class_key, const std::string & table_name)
{
    if (!TableExists(table_name))
    {
        return;
    }

    auto group_entry{ model_obj->GetBondGroupPotentialEntry(class_key) };
    auto iter{
        m_database->IterateQuery<GroupKey, int, double, double, double, double, double, double, double, double, double>(
            FormatSQL(SELECT_GROUP_ENTRY_SQL, table_name)) };
    std::tuple<GroupKey, int, double, double, double, double, double, double, double, double, double> row;
    while (iter.Next(row))
    {
        const auto group_key{ std::get<0>(row) };
        group_entry->InsertGroupKey(group_key);
        group_entry->ReserveBondObjectPtrList(group_key, std::get<1>(row));
        group_entry->AddGausEstimateMean(group_key, std::get<2>(row), std::get<3>(row));
        group_entry->AddGausEstimateMDPDE(group_key, std::get<4>(row), std::get<5>(row));
        group_entry->AddGausEstimatePrior(group_key, std::get<6>(row), std::get<7>(row));
        group_entry->AddGausVariancePrior(group_key, std::get<8>(row), std::get<9>(row));
        group_entry->AddAlphaG(group_key, std::get<10>(row));
    }

    for (auto & bond : model_obj->GetSelectedBondList())
    {
        group_entry->AddBondObjectPtr(BondClassifier::GetGroupKeyInClass(bond, class_key), bond);
    }
}

bool LegacyModelObjectReader::TableExists(const std::string & table_name) const
{
    if (m_table_cache.find(table_name) != m_table_cache.end())
    {
        return true;
    }

    m_database->Prepare(TABLE_EXISTS_SQL.data());
    SQLiteWrapper::StatementGuard guard(*m_database);
    m_database->Bind<std::string>(1, table_name);
    const auto rc{ m_database->StepNext() };
    if (rc != SQLiteWrapper::StepRow() && rc != SQLiteWrapper::StepDone())
    {
        throw std::runtime_error(m_database->ErrorMessage());
    }
    const bool exists{ rc == SQLiteWrapper::StepRow() };
    if (exists)
    {
        m_table_cache.insert(table_name);
    }
    return exists;
}

std::string LegacyModelObjectReader::SanitizeTableName(const std::string & key_tag) const
{
    std::string sanitized_key_tag;
    sanitized_key_tag.reserve(key_tag.size());
    for (char ch : key_tag)
    {
        if (std::isalnum(static_cast<unsigned char>(ch)) || ch == '_')
        {
            sanitized_key_tag.push_back(ch);
        }
        else
        {
            sanitized_key_tag.push_back('_');
        }
    }
    return sanitized_key_tag;
}

} // namespace rhbm_gem
