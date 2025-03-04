#include "ModelObjectDAO.hpp"
#include "SQLiteWrapper.hpp"
#include "DataObjectBase.hpp"
#include "ModelObject.hpp"
#include "AtomObject.hpp"
#include "AtomicPotentialEntry.hpp"

#include <sstream>
#include <stdexcept>

ModelObjectDAO::ModelObjectDAO(SQLiteWrapper * database) :
    m_database{ database }
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

    auto model_list_table_name{ "model_list" };
    auto atom_list_table_name{ "atom_list_in_" + key_tag };
    auto potential_entry_table_name{ "atom_potential_entry_in_" + key_tag };
    CreateModelObjectListTable(model_list_table_name);
    CreateAtomObjectListTable(atom_list_table_name);
    CreateAtomicPotentialEntryListTable(potential_entry_table_name);

    std::stringstream sql;
    sql <<"INSERT INTO "<< model_list_table_name << R"(
        (
            key_tag, atom_size, pdb_id, emd_id
        ) VALUES (?, ?, ?, ?)
        ON CONFLICT(key_tag)
        DO UPDATE SET
            key_tag = excluded.key_tag,
            atom_size = excluded.atom_size,
            pdb_id  = excluded.pdb_id,
            emd_id  = excluded.emd_id
    )";

    m_database->Prepare(sql.str());
    m_database->Bind<std::string>(1, model_obj->GetKeyTag());
    m_database->Bind<int>(2, model_obj->GetNumberOfAtom());
    m_database->Bind<std::string>(3, model_obj->GetPdbID());
    m_database->Bind<std::string>(4, model_obj->GetEmdID());

    m_database->StepOnce();
    m_database->Finalize();

    SaveAtomObjectList(model_obj, atom_list_table_name);
    SaveAtomicPotentialEntryList(model_obj, potential_entry_table_name);
}

std::unique_ptr<DataObjectBase> ModelObjectDAO::Load(const std::string & key_tag)
{
    auto atom_object_list{ LoadAtomObjectList("atom_list_in_" + key_tag) };
    auto model_object{ std::make_unique<ModelObject>(std::move(atom_object_list)) };

    std::string model_table_name{ "model_list" };
    std::stringstream sql;
    sql <<"SELECT "<< R"(
            key_tag, atom_size, pdb_id, emd_id
    )"<<" FROM "<< model_table_name << " WHERE key_tag = ? LIMIT 1;";

    m_database->Prepare(sql.str());
    m_database->Bind<std::string>(1, key_tag);
    auto return_code{ m_database->StepNext() };
    auto atom_size{ 0 };
    if (return_code == SQLiteWrapper::StepRow())
    {
        model_object->SetKeyTag(m_database->GetColumn<std::string>(0));
        model_object->SetPdbID(m_database->GetColumn<std::string>(2));
        model_object->SetEmdID(m_database->GetColumn<std::string>(3));
        atom_size = m_database->GetColumn<int>(1);
        m_database->Finalize();
    }
    else if (return_code == SQLiteWrapper::StepDone())
    {
        m_database->Finalize();
        throw std::runtime_error("Cannot find the row with key_tag = " + key_tag);
    }
    else
    {
        m_database->Finalize();
        throw std::runtime_error("Step failed: " + m_database->ErrorMessage());
    }

    if (atom_size != model_object->GetNumberOfAtom())
    {
        throw std::runtime_error("The number of atoms in the model object does not match the database record.");
    }

    return model_object;
}

void ModelObjectDAO::CreateModelObjectListTable(const std::string & table_name)
{
    std::stringstream sql;
    sql <<"CREATE TABLE IF NOT EXISTS "<< table_name << R"(
        (
            key_tag           TEXT PRIMARY KEY,
            atom_size         INTEGER,
            pdb_id            TEXT,
            emd_id            TEXT
        )
    )";
    m_database->Execute(sql.str());
}

void ModelObjectDAO::CreateAtomObjectListTable(const std::string & table_name)
{
    std::stringstream sql;
    sql <<"CREATE TABLE IF NOT EXISTS "<< table_name << R"(
        (
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
            status INTEGER,
            is_special_atom INTEGER,
            position_x DOUBLE,
            position_y DOUBLE,
            position_z DOUBLE
        )
    )";
    m_database->Execute(sql.str());
}

void ModelObjectDAO::CreateAtomicPotentialEntryListTable(const std::string & table_name)
{
    std::stringstream sql;
    sql <<"CREATE TABLE IF NOT EXISTS "<< table_name << R"(
        (
            serial_id INTEGER PRIMARY KEY,
            sampling_size INTEGER,
            distance_and_map_value_list BLOB,
            amplitude_estimate_ols DOUBLE,
            width_estimate_ols DOUBLE,
            amplitude_estimate_mdpde DOUBLE,
            width_estimate_mdpde DOUBLE
        )
    )";
    m_database->Execute(sql.str());
}

void ModelObjectDAO::CreateAtomicPotentialEntrySubListTable(
    const std::string & table_name, const std::string & group_key)
{
    std::stringstream sql;
    sql <<"CREATE TABLE IF NOT EXISTS "<< group_key <<"_"<< table_name << R"(
        (
            serial_id INTEGER PRIMARY KEY,
            amplitude_estimate_posterior DOUBLE,
            width_estimate_posterior DOUBLE,
            amplitude_variance_posterior DOUBLE,
            width_variance_posterior DOUBLE,
            outlier_tag INTEGER,
            statistical_distance DOUBLE
        )
    )";
    m_database->Execute(sql.str());
}

void ModelObjectDAO::SaveAtomObjectList(
    const ModelObject * model_obj, const std::string & table_name)
{
    m_database->BeginTransaction();

    std::stringstream sql;
    sql <<"INSERT OR REPLACE INTO "<< table_name << R"(
        (
            serial_id, residue_id, chain_id, indicator,
            occupancy, temperature, residue_type,
            element_type, remoteness_type, branch_type,
            status, is_special_atom,
            position_x, position_y, position_z
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?);
    )";

    m_database->Prepare(sql.str());

    for (auto & atom_object : model_obj->GetComponentsList())
    {
        m_database->Bind<int>(1, atom_object->GetSerialID());
        m_database->Bind<int>(2, atom_object->GetResidueID());
        m_database->Bind<std::string>(3, atom_object->GetChainID());
        m_database->Bind<std::string>(4, atom_object->GetIndicator());
        m_database->Bind<double>(5, static_cast<double>(atom_object->GetOccupancy()));
        m_database->Bind<double>(6, static_cast<double>(atom_object->GetTemperature()));
        m_database->Bind<int>(7, atom_object->GetResidue());
        m_database->Bind<int>(8, atom_object->GetElement());
        m_database->Bind<int>(9, atom_object->GetRemoteness());
        m_database->Bind<int>(10, atom_object->GetBranch());
        m_database->Bind<int>(11, atom_object->GetStatus());
        m_database->Bind<int>(12, static_cast<int>(atom_object->GetSpecialAtomFlag()));
        m_database->Bind<double>(13, static_cast<double>(atom_object->GetPosition().at(0)));
        m_database->Bind<double>(14, static_cast<double>(atom_object->GetPosition().at(1)));
        m_database->Bind<double>(15, static_cast<double>(atom_object->GetPosition().at(2)));

        m_database->StepOnce();
        m_database->Reset();
    }

    m_database->Finalize();
    m_database->CommitTransaction();
}

void ModelObjectDAO::SaveAtomicPotentialEntryList(
    const ModelObject * model_obj, const std::string & table_name)
{
    m_database->BeginTransaction();

    std::stringstream sql;
    sql <<"INSERT OR REPLACE INTO "<< table_name << R"(
        (
            serial_id, sampling_size, distance_and_map_value_list,
            amplitude_estimate_ols, width_estimate_ols,
            amplitude_estimate_mdpde, width_estimate_mdpde
        ) VALUES (?, ?, ?, ?, ?, ?, ?);
    )";

    m_database->Prepare(sql.str());

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

    m_database->Finalize();
    m_database->CommitTransaction();

    std::vector<std::string> group_list;
    for (auto & [key, entry] : model_obj->GetGroupPotentialEntryMap())
    {
        group_list.emplace_back(key);
        CreateAtomicPotentialEntrySubListTable(table_name, key);
        SaveAtomicPotentialEntrySubList(model_obj, table_name, key);
    }

}

void ModelObjectDAO::SaveAtomicPotentialEntrySubList(
    const ModelObject * model_obj, const std::string & table_name, const std::string & group_key)
{
    m_database->BeginTransaction();

    std::stringstream sql;
    sql <<"INSERT OR REPLACE INTO "<< group_key <<"_"<< table_name << R"(
        (
            serial_id, amplitude_estimate_posterior, width_estimate_posterior,
            amplitude_variance_posterior, width_variance_posterior,
            outlier_tag, statistical_distance
        ) VALUES (?, ?, ?, ?, ?, ?, ?);
    )";

    m_database->Prepare(sql.str());

    for (auto & atom_object : model_obj->GetComponentsList())
    {
        auto entry{ atom_object->GetAtomicPotentialEntry() };
        if (entry == nullptr) continue;

        m_database->Bind<int>(1, atom_object->GetSerialID());
        m_database->Bind<double>(2, entry->GetAmplitudeEstimatePosterior(group_key));
        m_database->Bind<double>(3, entry->GetWidthEstimatePosterior(group_key));
        m_database->Bind<double>(4, entry->GetAmplitudeVariancePosterior(group_key));
        m_database->Bind<double>(5, entry->GetWidthVariancePosterior(group_key));
        m_database->Bind<int>(6, static_cast<int>(entry->GetOutlierTag(group_key)));
        m_database->Bind<double>(7, entry->GetStatisticalDistance(group_key));

        m_database->StepOnce();
        m_database->Reset();
    }

    m_database->Finalize();
    m_database->CommitTransaction();
}

std::vector<std::unique_ptr<AtomObject>> ModelObjectDAO::LoadAtomObjectList(
    const std::string & table_name)
{
    std::stringstream sql;
    sql <<"SELECT "<< R"(
            serial_id, residue_id, chain_id, indicator,
            occupancy, temperature, residue_type,
            element_type, remoteness_type, branch_type,
            status, is_special_atom,
            position_x, position_y, position_z
        )"<<" FROM "<< table_name <<";";
    
    auto row_list
    {
        m_database->Query<int, int, std::string, std::string,
                  double, double, int, int, int, int, int, int, double, double, double>(sql.str())
    };

    std::vector<std::unique_ptr<AtomObject>> atom_object_list;
    atom_object_list.reserve(row_list.size());
    for (auto & row : row_list)
    {
        auto atom_object{ std::make_unique<AtomObject>() };
        atom_object->SetSerialID(std::get<0>(row));
        atom_object->SetResidueID(std::get<1>(row));
        atom_object->SetChainID(std::get<2>(row));
        atom_object->SetIndicator(std::get<3>(row));
        atom_object->SetOccupancy(std::get<4>(row));
        atom_object->SetTemperature(std::get<5>(row));
        atom_object->SetResidue(std::get<6>(row));
        atom_object->SetElement(std::get<7>(row));
        atom_object->SetRemoteness(std::get<8>(row));
        atom_object->SetBranch(std::get<9>(row));
        atom_object->SetStatus(std::get<10>(row));
        atom_object->SetSpecialAtomFlag(static_cast<bool>(std::get<11>(row)));
        atom_object->SetPosition(std::get<12>(row), std::get<13>(row), std::get<14>(row));
        atom_object_list.emplace_back(std::move(atom_object));
    }
    return atom_object_list;
}