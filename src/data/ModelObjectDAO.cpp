#include "ModelObjectDAO.hpp"
#include "DatabaseManager.hpp"
#include "SQLiteWrapper.hpp"
#include "DataObjectBase.hpp"
#include "ModelObject.hpp"
#include "AtomObject.hpp"

#include <sstream>
#include <stdexcept>

ModelObjectDAO::ModelObjectDAO(DatabaseManager & db_manager) :
    m_db_manager{ db_manager }
{

}

ModelObjectDAO::~ModelObjectDAO()
{

}

void ModelObjectDAO::Save(const DataObjectBase * obj, const std::string & key_tag)
{
    auto model_obj{ dynamic_cast<const ModelObject *>(obj) };
    if (!model_obj)
    {
        throw std::runtime_error("ModelObjectDAO::Save() failed: object is not a ModelObject instance.");
    }
    auto db{ m_db_manager.GetDatabase() };
    CreateModelObjectListTable(db);

    std::string sql{ R"(
        INSERT INTO model_list (
            key_tag, atom_size, pdb_id, emd_id
        ) VALUES (?, ?, ?, ?)
        ON CONFLICT(key_tag)
        DO UPDATE SET
            key_tag = excluded.key_tag,
            atom_size = excluded.atom_size,
            pdb_id  = excluded.pdb_id,
            emd_id  = excluded.emd_id
    )" };

    db->Prepare(sql);
    db->Bind<std::string>(1, model_obj->GetKeyTag());
    db->Bind<int>(2, model_obj->GetNumberOfAtom());
    db->Bind<std::string>(3, model_obj->GetPdbID());
    db->Bind<std::string>(4, model_obj->GetEmdID());

    db->StepOnce();
    db->Finalize();

    SaveAtomObjectList(model_obj);
}

std::unique_ptr<DataObjectBase> ModelObjectDAO::Load(const std::string & key_tag)
{
    auto atom_object_list{ LoadAtomObjectList(key_tag) };
    auto model_object{ std::make_unique<ModelObject>(std::move(atom_object_list)) };
    auto db{ m_db_manager.GetDatabase() };

    std::string model_table_name{ "model_list" };
    std::stringstream sql;
    sql <<"SELECT "<< R"(
            key_tag, atom_size, pdb_id, emd_id
    )"<<" FROM "<< model_table_name << " WHERE key_tag = ? LIMIT 1;";

    db->Prepare(sql.str());
    db->Bind<std::string>(1, key_tag);
    auto return_code{ db->StepNext() };
    auto atom_size{ 0 };
    if (return_code == SQLiteWrapper::StepRow())
    {
        model_object->SetKeyTag(db->GetColumn<std::string>(0));
        model_object->SetPdbID(db->GetColumn<std::string>(2));
        model_object->SetEmdID(db->GetColumn<std::string>(3));
        atom_size = db->GetColumn<int>(1);
        db->Finalize();
    }
    else if (return_code == SQLiteWrapper::StepDone())
    {
        db->Finalize();
        throw std::runtime_error("Cannot find the row with key_tag = " + key_tag);
    }
    else
    {
        db->Finalize();
        throw std::runtime_error("Step failed: " + db->ErrorMessage());
    }

    if (atom_size != model_object->GetNumberOfAtom())
    {
        throw std::runtime_error("The number of atoms in the model object does not match the database record.");
    }

    return model_object;
}

void ModelObjectDAO::CreateModelObjectListTable(SQLiteWrapper * db)
{
    std::string sql_info{ R"(
        CREATE TABLE IF NOT EXISTS model_list
        (
            key_tag           TEXT PRIMARY KEY,
            atom_size         INTEGER,
            pdb_id            TEXT,
            emd_id            TEXT
        );
    )" };
    db->Execute(sql_info);
}

void ModelObjectDAO::CreateAtomObjectListTable(SQLiteWrapper * db, const std::string & table_name)
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
    db->Execute(sql.str());
}

void ModelObjectDAO::SaveAtomObjectList(const ModelObject * model_obj)
{
    auto key_tag{ model_obj->GetKeyTag() };
    auto table_name{ "atom_list_in_" + key_tag };
    const auto & atom_object_list{ model_obj->GetComponentsList() };
    auto db{ m_db_manager.GetDatabase() };
    CreateAtomObjectListTable(db, table_name);
    db->BeginTransaction();

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

    db->Prepare(sql.str());

    for (auto & atom_object : atom_object_list)
    {
        db->Bind<int>(1, atom_object->GetSerialID());
        db->Bind<int>(2, atom_object->GetResidueID());
        db->Bind<std::string>(3, atom_object->GetChainID());
        db->Bind<std::string>(4, atom_object->GetIndicator());
        db->Bind<double>(5, static_cast<double>(atom_object->GetOccupancy()));
        db->Bind<double>(6, static_cast<double>(atom_object->GetTemperature()));
        db->Bind<int>(7, atom_object->GetResidue());
        db->Bind<int>(8, atom_object->GetElement());
        db->Bind<int>(9, atom_object->GetRemoteness());
        db->Bind<int>(10, atom_object->GetBranch());
        db->Bind<int>(11, atom_object->GetStatus());
        db->Bind<int>(12, static_cast<int>(atom_object->GetSpecialAtomFlag()));
        db->Bind<double>(13, static_cast<double>(atom_object->GetPosition().at(0)));
        db->Bind<double>(14, static_cast<double>(atom_object->GetPosition().at(1)));
        db->Bind<double>(15, static_cast<double>(atom_object->GetPosition().at(2)));

        db->StepOnce();
        db->Reset();
    }

    db->Finalize();
    db->CommitTransaction();
}

std::vector<std::unique_ptr<AtomObject>> ModelObjectDAO::LoadAtomObjectList(const std::string & key_tag)
{
    auto table_name{ "atom_list_in_" + key_tag };
    auto db{ m_db_manager.GetDatabase() };
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
        db->Query<int, int, std::string, std::string,
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