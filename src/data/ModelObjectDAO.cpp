#include "ModelObjectDAO.hpp"
#include "SQLiteWrapper.hpp"
#include "DataObjectBase.hpp"
#include "ModelObject.hpp"
#include "AtomicInfoHelper.hpp"
#include "AtomObject.hpp"
#include "AtomicPotentialEntry.hpp"
#include "GroupPotentialEntry.hpp"
#include "KeyPacker.hpp"

#include <cctype>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>

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

    constexpr std::string_view INSERT_MODEL_LIST_SQL = R"sql(
        INSERT INTO {} (
            key_tag, atom_size, pdb_id, emd_id, map_resolution, resolution_method
        ) VALUES (?, ?, ?, ?, ?, ?)
        ON CONFLICT(key_tag)
        DO UPDATE SET
            key_tag = excluded.key_tag,
            atom_size = excluded.atom_size,
            pdb_id  = excluded.pdb_id,
            emd_id  = excluded.emd_id,
            map_resolution = excluded.map_resolution,
            resolution_method = excluded.resolution_method
        )sql";

    constexpr std::string_view SELECT_MODEL_LIST_SQL = R"sql(
        SELECT
            key_tag, atom_size, pdb_id, emd_id, map_resolution, resolution_method
        FROM {} WHERE key_tag = ? LIMIT 1; )sql";

    constexpr std::string_view CREATE_MODEL_TABLE_SQL = R"sql(
        CREATE TABLE IF NOT EXISTS {} (
            key_tag           TEXT PRIMARY KEY,
            atom_size         INTEGER,
            pdb_id            TEXT,
            emd_id            TEXT,
            map_resolution    DOUBLE,
            resolution_method TEXT
        ) )sql";

    constexpr std::string_view CREATE_ATOM_TABLE_SQL = R"sql(
        CREATE TABLE IF NOT EXISTS {} (
            serial_id INTEGER PRIMARY KEY,
            residue_id INTEGER,
            chain_id TEXT,
            indicator TEXT,
            occupancy DOUBLE,
            temperature DOUBLE,
            residue_type INTEGER,
            element_type INTEGER,
            remoteness_type INTEGER,
            branch_type INTEGER,
            structure INTEGER,
            is_special_atom INTEGER,
            position_x DOUBLE,
            position_y DOUBLE,
            position_z DOUBLE
        ) )sql";

    constexpr std::string_view CREATE_ATOMIC_ENTRY_TABLE_SQL = R"sql(
        CREATE TABLE IF NOT EXISTS {} (
            serial_id INTEGER PRIMARY KEY,
            sampling_size INTEGER,
            distance_and_map_value_list BLOB,
            amplitude_estimate_ols DOUBLE,
            width_estimate_ols DOUBLE,
            amplitude_estimate_mdpde DOUBLE,
            width_estimate_mdpde DOUBLE
        ) )sql";

    constexpr std::string_view CREATE_ATOMIC_SUBLIST_TABLE_SQL = R"sql(
        CREATE TABLE IF NOT EXISTS {} (
            serial_id INTEGER PRIMARY KEY,
            amplitude_estimate_posterior DOUBLE,
            width_estimate_posterior DOUBLE,
            amplitude_variance_posterior DOUBLE,
            width_variance_posterior DOUBLE,
            outlier_tag INTEGER,
            statistical_distance DOUBLE
        ) )sql";

    constexpr std::string_view CREATE_GROUP_ENTRY_TABLE_SQL = R"sql(
        CREATE TABLE IF NOT EXISTS {} (
            key_id INTEGER PRIMARY KEY,
            atom_size INTEGER,
            amplitude_estimate_mean DOUBLE,
            width_estimate_mean DOUBLE,
            amplitude_estimate_mdpde DOUBLE,
            width_estimate_mdpde DOUBLE,
            amplitude_estimate_prior DOUBLE,
            width_estimate_prior DOUBLE,
            amplitude_variance_prior DOUBLE,
            width_variance_prior DOUBLE
        ) )sql";

    constexpr std::string_view INSERT_ATOM_LIST_SQL = R"sql(
        INSERT OR REPLACE INTO {} (
            serial_id, residue_id, chain_id, indicator,
            occupancy, temperature, residue_type,
            element_type, remoteness_type, branch_type,
            structure, is_special_atom,
            position_x, position_y, position_z
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?); )sql";

    constexpr std::string_view INSERT_ATOMIC_ENTRY_SQL = R"sql(
        INSERT OR REPLACE INTO {} (
            serial_id, sampling_size, distance_and_map_value_list,
            amplitude_estimate_ols, width_estimate_ols,
            amplitude_estimate_mdpde, width_estimate_mdpde
        ) VALUES (?, ?, ?, ?, ?, ?, ?); )sql";

    constexpr std::string_view INSERT_ATOMIC_SUBLIST_SQL = R"sql(
        INSERT OR REPLACE INTO {} (
            serial_id, amplitude_estimate_posterior, width_estimate_posterior,
            amplitude_variance_posterior, width_variance_posterior,
            outlier_tag, statistical_distance
        ) VALUES (?, ?, ?, ?, ?, ?, ?); )sql";

    constexpr std::string_view INSERT_GROUP_ENTRY_SQL = R"sql(
        INSERT OR REPLACE INTO {} (
            key_id, atom_size,
            amplitude_estimate_mean, width_estimate_mean,
            amplitude_estimate_mdpde, width_estimate_mdpde,
            amplitude_estimate_prior, width_estimate_prior,
            amplitude_variance_prior, width_variance_prior
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?); )sql";

    constexpr std::string_view SELECT_ATOM_LIST_SQL = R"sql(
        SELECT
            serial_id, residue_id, chain_id, indicator,
            occupancy, temperature, residue_type,
            element_type, remoteness_type, branch_type,
            structure, is_special_atom,
            position_x, position_y, position_z
        FROM {}; )sql";

    constexpr std::string_view SELECT_ATOMIC_ENTRY_SQL = R"sql(
        SELECT
            serial_id, sampling_size, distance_and_map_value_list,
            amplitude_estimate_ols, width_estimate_ols,
            amplitude_estimate_mdpde, width_estimate_mdpde
        FROM {}; )sql";

    constexpr std::string_view SELECT_ATOMIC_SUBLIST_SQL = R"sql(
        SELECT
            serial_id, amplitude_estimate_posterior, width_estimate_posterior,
            amplitude_variance_posterior, width_variance_posterior,
            outlier_tag, statistical_distance
        FROM {}; )sql";

    constexpr std::string_view SELECT_GROUP_ENTRY_SQL = R"sql(
        SELECT
            key_id, atom_size,
            amplitude_estimate_mean, width_estimate_mean,
            amplitude_estimate_mdpde, width_estimate_mdpde,
            amplitude_estimate_prior, width_estimate_prior,
            amplitude_variance_prior, width_variance_prior
        FROM {}; )sql";

    constexpr std::string_view TABLE_EXISTS_SQL =
        "SELECT name FROM sqlite_master WHERE type = 'table' AND name = ? LIMIT 1;";
}

ModelObjectDAO::ModelObjectDAO(SQLiteWrapper * database) :
    m_database{ database }, m_table_cache{}
{

}

ModelObjectDAO::~ModelObjectDAO()
{

}

void ModelObjectDAO::Save(const DataObjectBase * obj)
{
    auto model_obj{ dynamic_cast<const ModelObject *>(obj) };
    if (!model_obj)
    {
        throw std::runtime_error("ModelObjectDAO::Save() failed: object is not a ModelObject instance.");
    }

    auto key_tag{ model_obj->GetKeyTag() };
    auto sanitized_key_tag{ SanitizeTableName(key_tag) };

    SQLiteWrapper::TransactionGuard transaction(*m_database);
    
    // Save model object list
    auto model_list_table_name{ "model_list" };
    m_database->Execute(FormatSQL(CREATE_MODEL_TABLE_SQL, model_list_table_name));
    m_database->Prepare(FormatSQL(INSERT_MODEL_LIST_SQL, model_list_table_name));
    SQLiteWrapper::StatementGuard guard(*m_database);
    m_database->Bind<std::string>(1, model_obj->GetKeyTag());
    m_database->Bind<int>(2, static_cast<int>(model_obj->GetNumberOfAtom()));
    m_database->Bind<std::string>(3, model_obj->GetPdbID());
    m_database->Bind<std::string>(4, model_obj->GetEmdID());
    m_database->Bind<double>(5, model_obj->GetResolution());
    m_database->Bind<std::string>(6, model_obj->GetResolutionMethod());
    m_database->StepOnce();

    // Save atom object list
    auto atom_list_table_name{ "atom_list_in_" + sanitized_key_tag };
    SaveAtomObjectList(model_obj, atom_list_table_name);

    // Save atomic potential entries
    if (model_obj->GetGroupPotentialEntryMap().empty() == false)
    {
        auto atomic_potential_entry_table_name{ "atomic_potential_entry_in_" + sanitized_key_tag };
        SaveAtomicPotentialEntryList(model_obj, atomic_potential_entry_table_name);

        auto group_potential_entry_table_name{ "group_potential_entry_in_" + sanitized_key_tag };
        for (auto & [class_key, group_entry] : model_obj->GetGroupPotentialEntryMap())
        {
            auto group_table_name{ class_key + "_" + group_potential_entry_table_name };
            auto atom_table_name{ class_key + "_" + atomic_potential_entry_table_name };
            SaveAtomicPotentialEntrySubList(model_obj, atom_table_name, class_key);
            SaveGroupPotentialEntryList(group_entry.get(), group_table_name);
        }
    }
}

std::unique_ptr<DataObjectBase> ModelObjectDAO::Load(const std::string & key_tag)
{
    auto sanitized_key_tag{ SanitizeTableName(key_tag) };
    auto atom_object_list{ LoadAtomObjectList(sanitized_key_tag) };
    auto model_object{ std::make_unique<ModelObject>(std::move(atom_object_list)) };

    std::string model_table_name{ "model_list" };
    m_database->Prepare(FormatSQL(SELECT_MODEL_LIST_SQL, model_table_name));
    SQLiteWrapper::StatementGuard guard(*m_database);
    m_database->Bind<std::string>(1, key_tag);
    auto return_code{ m_database->StepNext() };
    auto atom_size{ 0 };
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
        throw std::runtime_error("The number of atoms in the model object does not match the database record.");
    }

    for (size_t i = 0; i < AtomicInfoHelper::GetGroupClassCount(); i++)
    {
        auto group_class_key{ AtomicInfoHelper::GetGroupClassKey(i) };
        // Load potential entry
        auto potential_table_name{ group_class_key + "_group_potential_entry_in_" + sanitized_key_tag };
        auto group_potential_entry{ std::make_unique<GroupPotentialEntry>() };
        model_object->AddGroupPotentialEntry(group_class_key, group_potential_entry);
        LoadGroupPotentialEntryList(model_object.get(), group_class_key, potential_table_name);
    }
    return model_object;
}

void ModelObjectDAO::SaveAtomObjectList(
    const ModelObject * model_obj, const std::string & table_name)
{
    m_database->Execute(FormatSQL(CREATE_ATOM_TABLE_SQL, table_name));
    m_database->ClearTable(table_name);
    m_database->Prepare(FormatSQL(INSERT_ATOM_LIST_SQL, table_name));
    SQLiteWrapper::StatementGuard guard(*m_database);

    for (auto & atom_object : model_obj->GetComponentsList())
    {
        m_database->Bind<int>(1, atom_object->GetSerialID());
        m_database->Bind<int>(2, atom_object->GetResidueID());
        m_database->Bind<std::string>(3, atom_object->GetChainID());
        m_database->Bind<std::string>(4, atom_object->GetIndicator());
        m_database->Bind<double>(5, static_cast<double>(atom_object->GetOccupancy()));
        m_database->Bind<double>(6, static_cast<double>(atom_object->GetTemperature()));
        m_database->Bind<int>(7, static_cast<int>(atom_object->GetResidue()));
        m_database->Bind<int>(8, static_cast<int>(atom_object->GetElement()));
        m_database->Bind<int>(9, static_cast<int>(atom_object->GetRemoteness()));
        m_database->Bind<int>(10, static_cast<int>(atom_object->GetBranch()));
        m_database->Bind<int>(11, static_cast<int>(atom_object->GetStructure()));
        m_database->Bind<int>(12, static_cast<int>(atom_object->GetSpecialAtomFlag()));
        m_database->Bind<double>(13, static_cast<double>(atom_object->GetPosition().at(0)));
        m_database->Bind<double>(14, static_cast<double>(atom_object->GetPosition().at(1)));
        m_database->Bind<double>(15, static_cast<double>(atom_object->GetPosition().at(2)));

        m_database->StepOnce();
        m_database->Reset();
    }
}

void ModelObjectDAO::SaveAtomicPotentialEntryList(
    const ModelObject * model_obj, const std::string & table_name)
{
    m_database->Execute(FormatSQL(CREATE_ATOMIC_ENTRY_TABLE_SQL, table_name));
    m_database->ClearTable(table_name);
    m_database->Prepare(FormatSQL(INSERT_ATOMIC_ENTRY_SQL, table_name));
    SQLiteWrapper::StatementGuard guard(*m_database);

    for (auto & atom_object : model_obj->GetComponentsList())
    {
        auto entry{ atom_object->GetAtomicPotentialEntry() };
        if (entry == nullptr) continue;

        m_database->Bind<int>(1, atom_object->GetSerialID());
        m_database->Bind<int>(2, entry->GetDistanceAndMapValueListSize());
        m_database->Bind<std::vector<std::tuple<float, float>>>(3, entry->GetDistanceAndMapValueList());
        m_database->Bind<double>(4, entry->GetAmplitudeEstimateOLS());
        m_database->Bind<double>(5, entry->GetWidthEstimateOLS());
        m_database->Bind<double>(6, entry->GetAmplitudeEstimateMDPDE());
        m_database->Bind<double>(7, entry->GetWidthEstimateMDPDE());

        m_database->StepOnce();
        m_database->Reset();
    }
}

void ModelObjectDAO::SaveAtomicPotentialEntrySubList(
    const ModelObject * model_obj, const std::string & table_name, const std::string & class_key)
{
    m_database->Execute(FormatSQL(CREATE_ATOMIC_SUBLIST_TABLE_SQL, table_name));
    m_database->ClearTable(table_name);
    m_database->Prepare(FormatSQL(INSERT_ATOMIC_SUBLIST_SQL, table_name));
    SQLiteWrapper::StatementGuard guard(*m_database);
    for (auto & atom_object : model_obj->GetComponentsList())
    {
        auto entry{ atom_object->GetAtomicPotentialEntry() };
        if (entry == nullptr) continue;
        m_database->Bind<int>(1, atom_object->GetSerialID());
        m_database->Bind<double>(2, entry->GetAmplitudeEstimatePosterior(class_key));
        m_database->Bind<double>(3, entry->GetWidthEstimatePosterior(class_key));
        m_database->Bind<double>(4, entry->GetAmplitudeVariancePosterior(class_key));
        m_database->Bind<double>(5, entry->GetWidthVariancePosterior(class_key));
        m_database->Bind<int>(6, static_cast<int>(entry->GetOutlierTag(class_key)));
        m_database->Bind<double>(7, entry->GetStatisticalDistance(class_key));
        m_database->StepOnce();
        m_database->Reset();
    }
}

void ModelObjectDAO::SaveGroupPotentialEntryList(
    const GroupPotentialEntry * group_entry, const std::string & table_name)
{
    m_database->Execute(FormatSQL(CREATE_GROUP_ENTRY_TABLE_SQL, table_name));
    m_database->ClearTable(table_name);
    m_database->Prepare(FormatSQL(INSERT_GROUP_ENTRY_SQL, table_name));
    SQLiteWrapper::StatementGuard guard(*m_database);
    for (auto & group_key : group_entry->GetGroupKeySet())
    {
        m_database->Bind<int64_t>(1, static_cast<int64_t>(group_key));
        m_database->Bind<int>(2, group_entry->GetAtomObjectPtrListSize(group_key));
        m_database->Bind<double>(3, std::get<0>(group_entry->GetGausEstimateMean(group_key)));
        m_database->Bind<double>(4, std::get<1>(group_entry->GetGausEstimateMean(group_key)));
        m_database->Bind<double>(5, std::get<0>(group_entry->GetGausEstimateMDPDE(group_key)));
        m_database->Bind<double>(6, std::get<1>(group_entry->GetGausEstimateMDPDE(group_key)));
        m_database->Bind<double>(7, std::get<0>(group_entry->GetGausEstimatePrior(group_key)));
        m_database->Bind<double>(8, std::get<1>(group_entry->GetGausEstimatePrior(group_key)));
        m_database->Bind<double>(9, std::get<0>(group_entry->GetGausVariancePrior(group_key)));
        m_database->Bind<double>(10, std::get<1>(group_entry->GetGausVariancePrior(group_key)));
        m_database->StepOnce();
        m_database->Reset();
    }
}

std::vector<std::unique_ptr<AtomObject>> ModelObjectDAO::LoadAtomObjectList(
    const std::string & key_tag)
{
    auto sanitized_key_tag{ SanitizeTableName(key_tag) };
    auto atom_list_table_name{ "atom_list_in_" + sanitized_key_tag };
    if (TableExists(atom_list_table_name) == false)
    {
        throw std::runtime_error("Required table not found: " + atom_list_table_name);
    }

    auto atomic_potential_entry_table_name{ "atomic_potential_entry_in_" + sanitized_key_tag };
    if (TableExists(atomic_potential_entry_table_name) == false)
    {
        throw std::runtime_error("Required table not found: " + atomic_potential_entry_table_name);
    }

    auto atomic_potential_entry_map{ LoadAtomicPotentialEntryMap(atomic_potential_entry_table_name) };

    std::vector<std::unique_ptr<AtomObject>> atom_object_list;
    auto iter{
        m_database->IterateQuery<int, int, std::string, std::string,
                double, double, int, int, int, int, int, int, double, double, double>(
            FormatSQL(SELECT_ATOM_LIST_SQL, atom_list_table_name)) };
    std::tuple<int, int, std::string, std::string,
               double, double, int, int, int, int, int, int, double, double, double> row;
    while (iter.Next(row))
    {
        auto atom_object{ std::make_unique<AtomObject>() };
        atom_object->SetSerialID(std::get<0>(row));
        atom_object->SetResidueID(std::get<1>(row));
        atom_object->SetChainID(std::get<2>(row));
        atom_object->SetIndicator(std::get<3>(row));
        atom_object->SetOccupancy(static_cast<float>(std::get<4>(row)));
        atom_object->SetTemperature(static_cast<float>(std::get<5>(row)));
        atom_object->SetResidue(static_cast<Residue>(std::get<6>(row)));
        atom_object->SetElement(static_cast<Element>(std::get<7>(row)));
        atom_object->SetRemoteness(static_cast<Remoteness>(std::get<8>(row)));
        atom_object->SetBranch(static_cast<Branch>(std::get<9>(row)));
        atom_object->SetStructure(static_cast<Structure>(std::get<10>(row)));
        atom_object->SetSpecialAtomFlag(static_cast<bool>(std::get<11>(row)));
        atom_object->SetPosition(
            static_cast<float>(std::get<12>(row)),
            static_cast<float>(std::get<13>(row)),
            static_cast<float>(std::get<14>(row)) );

        auto serial_id{ atom_object->GetSerialID() };

        // Load potential entry
        if (atomic_potential_entry_map.empty() == false)
        {
            if (atomic_potential_entry_map.find(serial_id) != atomic_potential_entry_map.end())
            {
                atom_object->SetSelectedFlag(true);
                atom_object->AddAtomicPotentialEntry(std::move(atomic_potential_entry_map.at(serial_id)));
            }
            else
            {
                atom_object->SetSelectedFlag(false);
            }
        }

        atom_object_list.emplace_back(std::move(atom_object));
    }
    return atom_object_list;
}

std::unordered_map<int, std::unique_ptr<AtomicPotentialEntry>>
ModelObjectDAO::LoadAtomicPotentialEntryMap(const std::string & table_name)
{
    if (TableExists(table_name) == false) return {};
    auto serial_id{ 0 };
    std::unordered_map<int, std::unique_ptr<AtomicPotentialEntry>> atomic_potential_entry_map;
    auto iter{
        m_database->IterateQuery<int, int, std::vector<std::tuple<float, float>>, double, double, double, double>(
            FormatSQL(SELECT_ATOMIC_ENTRY_SQL, table_name)) };
    std::tuple<int, int, std::vector<std::tuple<float, float>>, double, double, double, double> row;
    while (iter.Next(row))
    {
        auto atomic_potential_entry{ std::make_unique<AtomicPotentialEntry>() };
        serial_id = std::get<0>(row);
        atomic_potential_entry->AddDistanceAndMapValueList(std::get<2>(row));
        atomic_potential_entry->AddGausEstimateOLS(std::get<3>(row), std::get<4>(row));
        atomic_potential_entry->AddGausEstimateMDPDE(std::get<5>(row), std::get<6>(row));
        atomic_potential_entry_map[serial_id] = std::move(atomic_potential_entry);
    }

    for (size_t i = 0; i < AtomicInfoHelper::GetGroupClassCount(); i++)
    {
        auto class_key{ AtomicInfoHelper::GetGroupClassKey(i) };
        auto table_name_with_class_key{ class_key + "_" + table_name };
        LoadAtomicPotentialEntrySubList(table_name_with_class_key, class_key, atomic_potential_entry_map);
    }
    return atomic_potential_entry_map;
}

void ModelObjectDAO::LoadAtomicPotentialEntrySubList(
    const std::string & table_name, const std::string & class_key,
    std::unordered_map<int, std::unique_ptr<AtomicPotentialEntry>> & entry_map)
{
    auto serial_id{ 0 };
    auto iter{ m_database->IterateQuery<int, double, double, double, double, int, double>(
            FormatSQL(SELECT_ATOMIC_SUBLIST_SQL, table_name)) };
    std::tuple<int, double, double, double, double, int, double> row;
    while (iter.Next(row))
    {
        serial_id = std::get<0>(row);
        if (entry_map.find(serial_id) == entry_map.end()) continue;
        auto & entry{ entry_map.at(serial_id) };
        entry->AddGausEstimatePosterior(class_key, std::get<1>(row), std::get<2>(row));
        entry->AddGausVariancePosterior(class_key, std::get<3>(row), std::get<4>(row));
        entry->AddOutlierTag(class_key, static_cast<bool>(std::get<5>(row)));
        entry->AddStatisticalDistance(class_key, std::get<6>(row));
    }
}

void ModelObjectDAO::LoadGroupPotentialEntryList(
    ModelObject * model_obj, const std::string & class_key, const std::string & table_name)
{
    if (TableExists(table_name) == false) return;
    auto group_entry{ model_obj->GetGroupPotentialEntry(class_key) };
    auto iter{
        m_database->IterateQuery<int64_t, int,
                                 double, double, double, double,
                                 double, double, double, double>(
            FormatSQL(SELECT_GROUP_ENTRY_SQL, table_name)) };
    std::tuple<int64_t, int,
               double, double, double, double,
               double, double, double, double> row;
    while (iter.Next(row))
    {
        auto group_key{ static_cast<uint64_t>(std::get<0>(row)) };
        group_entry->InsertGroupKey(group_key);
        group_entry->ReserveAtomObjectPtrList(group_key, std::get<1>(row));
        group_entry->AddGausEstimateMean(group_key, std::get<2>(row), std::get<3>(row));
        group_entry->AddGausEstimateMDPDE(group_key, std::get<4>(row), std::get<5>(row));
        group_entry->AddGausEstimatePrior(group_key, std::get<6>(row), std::get<7>(row));
        group_entry->AddGausVariancePrior(group_key, std::get<8>(row), std::get<9>(row));
    }
    LoadGroupPotentialEntrySubList(model_obj, class_key);
}

void ModelObjectDAO::LoadGroupPotentialEntrySubList(
    ModelObject * model_obj, const std::string & class_key)
{
    auto group_entry{ model_obj->GetGroupPotentialEntry(class_key) };
    if (class_key == AtomicInfoHelper::GetElementClassKey())
    {
        for (auto & atom : model_obj->GetSelectedAtomList())
        {
            auto group_key{ KeyPackerElementClass::Pack(
                atom->GetElement(), atom->GetRemoteness(), atom->GetSpecialAtomFlag()) };
            group_entry->AddAtomObjectPtr(group_key, atom);
        }
    }
    else if (class_key == AtomicInfoHelper::GetResidueClassKey())
    {
        for (auto & atom : model_obj->GetSelectedAtomList())
        {
            auto group_key{
                KeyPackerResidueClass::Pack(
                    atom->GetResidue(), atom->GetElement(),
                    atom->GetRemoteness(), atom->GetBranch(), atom->GetSpecialAtomFlag())
            };
            group_entry->AddAtomObjectPtr(group_key, atom);
        }
    }
    else if (class_key == AtomicInfoHelper::GetStructureClassKey())
    {
        for (auto & atom : model_obj->GetSelectedAtomList())
        {
            auto group_key{
                KeyPackerStructureClass::Pack(
                    atom->GetStructure(), atom->GetResidue(),
                    atom->GetElement(), atom->GetRemoteness(),
                    atom->GetBranch(), atom->GetSpecialAtomFlag())
            };
            group_entry->AddAtomObjectPtr(group_key, atom);
        }
    }
    else
    {
        throw std::runtime_error("ModelObjectDAO::LoadGroupPotentialEntrySubList() failed: "
                                 "Unsupported group class key: " + class_key);
    }
}

bool ModelObjectDAO::TableExists(const std::string & table_name) const
{
    if (m_table_cache.find(table_name) != m_table_cache.end())
    {
        return true;
    }

    m_database->Prepare(TABLE_EXISTS_SQL.data());
    SQLiteWrapper::StatementGuard guard(*m_database);
    m_database->Bind<std::string>(1, table_name);
    auto rc{ m_database->StepNext() };
    if (rc != SQLiteWrapper::StepRow() && rc != SQLiteWrapper::StepDone())
    {
        throw std::runtime_error(m_database->ErrorMessage());
    }
    bool exists{ rc == SQLiteWrapper::StepRow() };
    if (exists)
    {
        m_table_cache.insert(table_name);
    }
    return exists;
}

std::string ModelObjectDAO::SanitizeTableName(const std::string & key_tag) const
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
