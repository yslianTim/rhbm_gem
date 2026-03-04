#include "ModelStructurePersistence.hpp"

#include "ModelSchemaSql.hpp"
#include "SQLiteStatementBatch.hpp"

#include "AtomObject.hpp"
#include "BondObject.hpp"
#include "ChemicalComponentEntry.hpp"
#include "LocalPotentialEntry.hpp"
#include "ModelObject.hpp"
#include "SQLiteWrapper.hpp"

#include <stdexcept>
#include <string>
#include <unordered_map>

namespace rhbm_gem::model_io {

namespace {

void SaveModelObjectRow(
    SQLiteWrapper & database,
    const ModelObject & model_obj,
    const std::string & key_tag)
{
    SQLiteStatementBatch batch{ database, std::string(kUpsertModelObjectSql) };
    batch.Execute([&](SQLiteWrapper & statement_db)
    {
        statement_db.Bind<std::string>(1, key_tag);
        statement_db.Bind<int>(2, static_cast<int>(model_obj.GetNumberOfAtom()));
        statement_db.Bind<std::string>(3, model_obj.GetPdbID());
        statement_db.Bind<std::string>(4, model_obj.GetEmdID());
        statement_db.Bind<double>(5, model_obj.GetResolution());
        statement_db.Bind<std::string>(6, model_obj.GetResolutionMethod());
    });
}

void SaveChainMap(
    SQLiteWrapper & database,
    const ModelObject & model_obj,
    const std::string & key_tag)
{
    SQLiteStatementBatch batch{ database, std::string(kInsertModelChainMapSql) };
    for (const auto & [entity_id, chain_id_list] : model_obj.GetChainIDListMap())
    {
        for (size_t ordinal = 0; ordinal < chain_id_list.size(); ordinal++)
        {
            batch.Execute([&](SQLiteWrapper & statement_db)
            {
                statement_db.Bind<std::string>(1, key_tag);
                statement_db.Bind<std::string>(2, entity_id);
                statement_db.Bind<int>(3, static_cast<int>(ordinal));
                statement_db.Bind<std::string>(4, chain_id_list.at(ordinal));
            });
        }
    }
}

void SaveChemicalComponentEntryList(
    SQLiteWrapper & database,
    const ModelObject & model_obj,
    const std::string & key_tag)
{
    SQLiteStatementBatch batch{ database, std::string(kInsertModelComponentSql) };
    for (const auto & [component_key, component_entry] : model_obj.GetChemicalComponentEntryMap())
    {
        batch.Execute([&](SQLiteWrapper & statement_db)
        {
            statement_db.Bind<std::string>(1, key_tag);
            statement_db.Bind<ComponentKey>(2, component_key);
            statement_db.Bind<std::string>(3, component_entry->GetComponentId());
            statement_db.Bind<std::string>(4, component_entry->GetComponentName());
            statement_db.Bind<std::string>(5, component_entry->GetComponentType());
            statement_db.Bind<std::string>(6, component_entry->GetComponentFormula());
            statement_db.Bind<double>(7, static_cast<double>(component_entry->GetComponentMolecularWeight()));
            statement_db.Bind<int>(8, static_cast<int>(component_entry->IsStandardMonomer()));
        });
    }
}

void SaveComponentAtomEntryList(
    SQLiteWrapper & database,
    const ModelObject & model_obj,
    const std::string & key_tag)
{
    SQLiteStatementBatch batch{ database, std::string(kInsertModelComponentAtomSql) };
    for (const auto & [component_key, component_entry] : model_obj.GetChemicalComponentEntryMap())
    {
        for (const auto & [atom_key, atom_entry] : component_entry->GetComponentAtomEntryMap())
        {
            batch.Execute([&](SQLiteWrapper & statement_db)
            {
                statement_db.Bind<std::string>(1, key_tag);
                statement_db.Bind<ComponentKey>(2, component_key);
                statement_db.Bind<AtomKey>(3, atom_key);
                statement_db.Bind<std::string>(4, atom_entry.atom_id);
                statement_db.Bind<int>(5, static_cast<int>(atom_entry.element_type));
                statement_db.Bind<int>(6, static_cast<int>(atom_entry.aromatic_atom_flag));
                statement_db.Bind<int>(7, static_cast<int>(atom_entry.stereo_config));
            });
        }
    }
}

void SaveComponentBondEntryList(
    SQLiteWrapper & database,
    const ModelObject & model_obj,
    const std::string & key_tag)
{
    SQLiteStatementBatch batch{ database, std::string(kInsertModelComponentBondSql) };
    for (const auto & [component_key, component_entry] : model_obj.GetChemicalComponentEntryMap())
    {
        for (const auto & [bond_key, bond_entry] : component_entry->GetComponentBondEntryMap())
        {
            batch.Execute([&](SQLiteWrapper & statement_db)
            {
                statement_db.Bind<std::string>(1, key_tag);
                statement_db.Bind<ComponentKey>(2, component_key);
                statement_db.Bind<BondKey>(3, bond_key);
                statement_db.Bind<std::string>(4, bond_entry.bond_id);
                statement_db.Bind<int>(5, static_cast<int>(bond_entry.bond_type));
                statement_db.Bind<int>(6, static_cast<int>(bond_entry.bond_order));
                statement_db.Bind<int>(7, static_cast<int>(bond_entry.aromatic_atom_flag));
                statement_db.Bind<int>(8, static_cast<int>(bond_entry.stereo_config));
            });
        }
    }
}

void SaveAtomObjectList(
    SQLiteWrapper & database,
    const ModelObject & model_obj,
    const std::string & key_tag)
{
    SQLiteStatementBatch batch{ database, std::string(kInsertModelAtomSql) };
    for (const auto & atom_object : model_obj.GetAtomList())
    {
        batch.Execute([&](SQLiteWrapper & statement_db)
        {
            statement_db.Bind<std::string>(1, key_tag);
            statement_db.Bind<int>(2, atom_object->GetSerialID());
            statement_db.Bind<int>(3, atom_object->GetSequenceID());
            statement_db.Bind<std::string>(4, atom_object->GetComponentID());
            statement_db.Bind<std::string>(5, atom_object->GetAtomID());
            statement_db.Bind<std::string>(6, atom_object->GetChainID());
            statement_db.Bind<std::string>(7, atom_object->GetIndicator());
            statement_db.Bind<double>(8, static_cast<double>(atom_object->GetOccupancy()));
            statement_db.Bind<double>(9, static_cast<double>(atom_object->GetTemperature()));
            statement_db.Bind<int>(10, static_cast<int>(atom_object->GetElement()));
            statement_db.Bind<int>(11, static_cast<int>(atom_object->GetStructure()));
            statement_db.Bind<int>(12, static_cast<int>(atom_object->GetSpecialAtomFlag()));
            statement_db.Bind<double>(13, static_cast<double>(atom_object->GetPosition().at(0)));
            statement_db.Bind<double>(14, static_cast<double>(atom_object->GetPosition().at(1)));
            statement_db.Bind<double>(15, static_cast<double>(atom_object->GetPosition().at(2)));
            statement_db.Bind<int>(16, static_cast<int>(atom_object->GetComponentKey()));
            statement_db.Bind<int>(17, static_cast<int>(atom_object->GetAtomKey()));
        });
    }
}

void SaveBondObjectList(
    SQLiteWrapper & database,
    const ModelObject & model_obj,
    const std::string & key_tag)
{
    SQLiteStatementBatch batch{ database, std::string(kInsertModelBondSql) };
    for (const auto & bond_object : model_obj.GetBondList())
    {
        batch.Execute([&](SQLiteWrapper & statement_db)
        {
            statement_db.Bind<std::string>(1, key_tag);
            statement_db.Bind<int>(2, bond_object->GetAtomSerialID1());
            statement_db.Bind<int>(3, bond_object->GetAtomSerialID2());
            statement_db.Bind<int>(4, bond_object->GetBondKey());
            statement_db.Bind<int>(5, static_cast<int>(bond_object->GetBondType()));
            statement_db.Bind<int>(6, static_cast<int>(bond_object->GetBondOrder()));
            statement_db.Bind<int>(7, static_cast<int>(bond_object->GetSpecialBondFlag()));
        });
    }
}

void LoadModelObjectRow(SQLiteWrapper & database, ModelObject & model_obj, const std::string & key_tag)
{
    database.Prepare(std::string(kSelectModelObjectSql));
    SQLiteWrapper::StatementGuard guard(database);
    database.Bind<std::string>(1, key_tag);
    const auto rc{ database.StepNext() };
    if (rc == SQLiteWrapper::StepDone())
    {
        throw std::runtime_error("Cannot find the row with key_tag = " + key_tag);
    }
    if (rc != SQLiteWrapper::StepRow())
    {
        throw std::runtime_error("Step failed: " + database.ErrorMessage());
    }

    const auto atom_size{ database.GetColumn<int>(1) };
    model_obj.SetKeyTag(database.GetColumn<std::string>(0));
    model_obj.SetPdbID(database.GetColumn<std::string>(2));
    model_obj.SetEmdID(database.GetColumn<std::string>(3));
    model_obj.SetResolution(database.GetColumn<double>(4));
    model_obj.SetResolutionMethod(database.GetColumn<std::string>(5));
    if (atom_size != static_cast<int>(model_obj.GetNumberOfAtom()))
    {
        throw std::runtime_error(
            "The number of atoms in the model object does not match the database record.");
    }
}

void LoadChainMap(SQLiteWrapper & database, ModelObject & model_obj, const std::string & key_tag)
{
    std::unordered_map<std::string, std::vector<std::string>> chain_map;
    database.Prepare(std::string(kSelectModelChainMapSql));
    SQLiteWrapper::StatementGuard guard(database);
    database.Bind<std::string>(1, key_tag);
    while (true)
    {
        const auto rc{ database.StepNext() };
        if (rc == SQLiteWrapper::StepDone())
        {
            break;
        }
        if (rc != SQLiteWrapper::StepRow())
        {
            throw std::runtime_error("Step failed: " + database.ErrorMessage());
        }
        auto & chain_id_list{ chain_map[database.GetColumn<std::string>(0)] };
        const auto ordinal{ database.GetColumn<int>(1) };
        if (static_cast<int>(chain_id_list.size()) <= ordinal)
        {
            chain_id_list.resize(static_cast<size_t>(ordinal + 1));
        }
        chain_id_list.at(static_cast<size_t>(ordinal)) = database.GetColumn<std::string>(2);
    }
    model_obj.SetChainIDListMap(chain_map);
}

void LoadChemicalComponentEntryList(
    SQLiteWrapper & database,
    ModelObject & model_obj,
    const std::string & key_tag)
{
    database.Prepare(std::string(kSelectModelComponentSql));
    SQLiteWrapper::StatementGuard guard(database);
    database.Bind<std::string>(1, key_tag);
    while (true)
    {
        const auto rc{ database.StepNext() };
        if (rc == SQLiteWrapper::StepDone())
        {
            break;
        }
        if (rc != SQLiteWrapper::StepRow())
        {
            throw std::runtime_error("Step failed: " + database.ErrorMessage());
        }

        auto component_entry{ std::make_unique<ChemicalComponentEntry>() };
        const auto component_key{ database.GetColumn<ComponentKey>(0) };
        const auto component_id{ database.GetColumn<std::string>(1) };
        component_entry->SetComponentId(component_id);
        component_entry->SetComponentName(database.GetColumn<std::string>(2));
        component_entry->SetComponentType(database.GetColumn<std::string>(3));
        component_entry->SetComponentFormula(database.GetColumn<std::string>(4));
        component_entry->SetComponentMolecularWeight(
            static_cast<float>(database.GetColumn<double>(5)));
        component_entry->SetStandardMonomerFlag(static_cast<bool>(database.GetColumn<int>(6)));
        model_obj.AddChemicalComponentEntry(component_key, std::move(component_entry));
        model_obj.GetComponentKeySystemPtr()->RegisterComponent(component_id, component_key);
    }
}

void LoadComponentAtomEntryList(
    SQLiteWrapper & database,
    ModelObject & model_obj,
    const std::string & key_tag)
{
    database.Prepare(std::string(kSelectModelComponentAtomSql));
    SQLiteWrapper::StatementGuard guard(database);
    database.Bind<std::string>(1, key_tag);
    while (true)
    {
        const auto rc{ database.StepNext() };
        if (rc == SQLiteWrapper::StepDone())
        {
            break;
        }
        if (rc != SQLiteWrapper::StepRow())
        {
            throw std::runtime_error("Step failed: " + database.ErrorMessage());
        }

        const auto component_key{ database.GetColumn<ComponentKey>(0) };
        if (model_obj.GetChemicalComponentEntryMap().find(component_key)
            == model_obj.GetChemicalComponentEntryMap().end())
        {
            continue;
        }
        const ComponentAtomEntry atom_entry{
            database.GetColumn<std::string>(2),
            static_cast<Element>(database.GetColumn<int>(3)),
            static_cast<bool>(database.GetColumn<int>(4)),
            static_cast<StereoChemistry>(database.GetColumn<int>(5))
        };
        const auto atom_key{ database.GetColumn<AtomKey>(1) };
        model_obj.GetChemicalComponentEntryMap().at(component_key)->AddComponentAtomEntry(
            atom_key, atom_entry);
        model_obj.GetAtomKeySystemPtr()->RegisterAtom(atom_entry.atom_id, atom_key);
    }
}

void LoadComponentBondEntryList(
    SQLiteWrapper & database,
    ModelObject & model_obj,
    const std::string & key_tag)
{
    database.Prepare(std::string(kSelectModelComponentBondSql));
    SQLiteWrapper::StatementGuard guard(database);
    database.Bind<std::string>(1, key_tag);
    while (true)
    {
        const auto rc{ database.StepNext() };
        if (rc == SQLiteWrapper::StepDone())
        {
            break;
        }
        if (rc != SQLiteWrapper::StepRow())
        {
            throw std::runtime_error("Step failed: " + database.ErrorMessage());
        }

        const auto component_key{ database.GetColumn<ComponentKey>(0) };
        if (model_obj.GetChemicalComponentEntryMap().find(component_key)
            == model_obj.GetChemicalComponentEntryMap().end())
        {
            continue;
        }
        ComponentBondEntry bond_entry;
        const auto bond_key{ database.GetColumn<BondKey>(1) };
        bond_entry.bond_id = database.GetColumn<std::string>(2);
        bond_entry.bond_type = static_cast<BondType>(database.GetColumn<int>(3));
        bond_entry.bond_order = static_cast<BondOrder>(database.GetColumn<int>(4));
        bond_entry.aromatic_atom_flag = static_cast<bool>(database.GetColumn<int>(5));
        bond_entry.stereo_config = static_cast<StereoChemistry>(database.GetColumn<int>(6));
        model_obj.GetChemicalComponentEntryMap().at(component_key)->AddComponentBondEntry(
            bond_key, bond_entry);
        model_obj.GetBondKeySystemPtr()->RegisterBond(bond_entry.bond_id, bond_key);
    }
}

std::vector<std::unique_ptr<AtomObject>> LoadAtomObjectList(
    SQLiteWrapper & database,
    const std::string & key_tag)
{
    std::vector<std::unique_ptr<AtomObject>> atom_object_list;

    database.Prepare(std::string(kSelectModelAtomSql));
    SQLiteWrapper::StatementGuard guard(database);
    database.Bind<std::string>(1, key_tag);
    while (true)
    {
        const auto rc{ database.StepNext() };
        if (rc == SQLiteWrapper::StepDone())
        {
            break;
        }
        if (rc != SQLiteWrapper::StepRow())
        {
            throw std::runtime_error("Step failed: " + database.ErrorMessage());
        }

        auto atom_object{ std::make_unique<AtomObject>() };
        atom_object->SetSerialID(database.GetColumn<int>(0));
        atom_object->SetSequenceID(database.GetColumn<int>(1));
        atom_object->SetComponentID(database.GetColumn<std::string>(2));
        atom_object->SetAtomID(database.GetColumn<std::string>(3));
        atom_object->SetChainID(database.GetColumn<std::string>(4));
        atom_object->SetIndicator(database.GetColumn<std::string>(5));
        atom_object->SetOccupancy(static_cast<float>(database.GetColumn<double>(6)));
        atom_object->SetTemperature(static_cast<float>(database.GetColumn<double>(7)));
        atom_object->SetElement(static_cast<Element>(database.GetColumn<int>(8)));
        atom_object->SetStructure(static_cast<Structure>(database.GetColumn<int>(9)));
        atom_object->SetSpecialAtomFlag(static_cast<bool>(database.GetColumn<int>(10)));
        atom_object->SetPosition(
            static_cast<float>(database.GetColumn<double>(11)),
            static_cast<float>(database.GetColumn<double>(12)),
            static_cast<float>(database.GetColumn<double>(13)));
        atom_object->SetComponentKey(database.GetColumn<ComponentKey>(14));
        atom_object->SetAtomKey(database.GetColumn<AtomKey>(15));
        atom_object->SetSelectedFlag(false);
        atom_object_list.emplace_back(std::move(atom_object));
    }
    return atom_object_list;
}

std::vector<std::unique_ptr<BondObject>> LoadBondObjectList(
    SQLiteWrapper & database,
    const std::string & key_tag,
    const ModelObject & model_obj)
{
    std::vector<std::unique_ptr<BondObject>> bond_object_list;

    database.Prepare(std::string(kSelectModelBondSql));
    SQLiteWrapper::StatementGuard guard(database);
    database.Bind<std::string>(1, key_tag);
    while (true)
    {
        const auto rc{ database.StepNext() };
        if (rc == SQLiteWrapper::StepDone())
        {
            break;
        }
        if (rc != SQLiteWrapper::StepRow())
        {
            throw std::runtime_error("Step failed: " + database.ErrorMessage());
        }

        auto atom_object_1{ model_obj.GetAtomPtr(database.GetColumn<int>(0)) };
        auto atom_object_2{ model_obj.GetAtomPtr(database.GetColumn<int>(1)) };
        auto bond_object{ std::make_unique<BondObject>(atom_object_1, atom_object_2) };
        bond_object->SetBondKey(database.GetColumn<BondKey>(2));
        bond_object->SetBondType(static_cast<BondType>(database.GetColumn<int>(3)));
        bond_object->SetBondOrder(static_cast<BondOrder>(database.GetColumn<int>(4)));
        bond_object->SetSpecialBondFlag(static_cast<bool>(database.GetColumn<int>(5)));
        bond_object->SetSelectedFlag(false);
        bond_object_list.emplace_back(std::move(bond_object));
    }
    return bond_object_list;
}

} // namespace

void SaveStructure(SQLiteWrapper & database, const ModelObject & model_obj, const std::string & key_tag)
{
    SaveModelObjectRow(database, model_obj, key_tag);
    SaveChainMap(database, model_obj, key_tag);
    SaveChemicalComponentEntryList(database, model_obj, key_tag);
    SaveComponentAtomEntryList(database, model_obj, key_tag);
    SaveComponentBondEntryList(database, model_obj, key_tag);
    SaveAtomObjectList(database, model_obj, key_tag);
    SaveBondObjectList(database, model_obj, key_tag);
}

void LoadStructure(SQLiteWrapper & database, ModelObject & model_obj, const std::string & key_tag)
{
    LoadChemicalComponentEntryList(database, model_obj, key_tag);
    LoadComponentAtomEntryList(database, model_obj, key_tag);
    LoadComponentBondEntryList(database, model_obj, key_tag);
    model_obj.SetAtomList(LoadAtomObjectList(database, key_tag));
    model_obj.SetBondList(LoadBondObjectList(database, key_tag, model_obj));
    LoadModelObjectRow(database, model_obj, key_tag);
    LoadChainMap(database, model_obj, key_tag);
}

} // namespace rhbm_gem::model_io
