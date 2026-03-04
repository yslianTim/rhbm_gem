#include "ModelSchemaSql.hpp"
#include "SQLiteStatementBatch.hpp"

#include "AtomObject.hpp"
#include "BondObject.hpp"
#include "ChemicalComponentEntry.hpp"
#include "LocalPotentialEntry.hpp"
#include "ModelObject.hpp"
#include "ModelObjectDAOv2.hpp"
#include "SQLiteWrapper.hpp"

#include <stdexcept>
#include <string>
#include <unordered_map>

namespace rhbm_gem {

void ModelObjectDAOv2::SaveModelObjectRow(const ModelObject * model_obj, const std::string & key_tag)
{
    model_io::SQLiteStatementBatch batch{ *m_database, std::string(model_io::kUpsertModelObjectSql) };
    batch.Execute([&](SQLiteWrapper & database)
    {
        database.Bind<std::string>(1, key_tag);
        database.Bind<int>(2, static_cast<int>(model_obj->GetNumberOfAtom()));
        database.Bind<std::string>(3, model_obj->GetPdbID());
        database.Bind<std::string>(4, model_obj->GetEmdID());
        database.Bind<double>(5, model_obj->GetResolution());
        database.Bind<std::string>(6, model_obj->GetResolutionMethod());
    });
}

void ModelObjectDAOv2::SaveChainMap(const ModelObject * model_obj, const std::string & key_tag)
{
    model_io::SQLiteStatementBatch batch{ *m_database, std::string(model_io::kInsertModelChainMapSql) };
    for (const auto & [entity_id, chain_id_list] : model_obj->GetChainIDListMap())
    {
        for (size_t ordinal = 0; ordinal < chain_id_list.size(); ordinal++)
        {
            batch.Execute([&](SQLiteWrapper & database)
            {
                database.Bind<std::string>(1, key_tag);
                database.Bind<std::string>(2, entity_id);
                database.Bind<int>(3, static_cast<int>(ordinal));
                database.Bind<std::string>(4, chain_id_list.at(ordinal));
            });
        }
    }
}

void ModelObjectDAOv2::SaveChemicalComponentEntryList(const ModelObject * model_obj, const std::string & key_tag)
{
    model_io::SQLiteStatementBatch batch{ *m_database, std::string(model_io::kInsertModelComponentSql) };
    for (const auto & [component_key, component_entry] : model_obj->GetChemicalComponentEntryMap())
    {
        batch.Execute([&](SQLiteWrapper & database)
        {
            database.Bind<std::string>(1, key_tag);
            database.Bind<ComponentKey>(2, component_key);
            database.Bind<std::string>(3, component_entry->GetComponentId());
            database.Bind<std::string>(4, component_entry->GetComponentName());
            database.Bind<std::string>(5, component_entry->GetComponentType());
            database.Bind<std::string>(6, component_entry->GetComponentFormula());
            database.Bind<double>(7, static_cast<double>(component_entry->GetComponentMolecularWeight()));
            database.Bind<int>(8, static_cast<int>(component_entry->IsStandardMonomer()));
        });
    }
}

void ModelObjectDAOv2::SaveComponentAtomEntryList(const ModelObject * model_obj, const std::string & key_tag)
{
    model_io::SQLiteStatementBatch batch{ *m_database, std::string(model_io::kInsertModelComponentAtomSql) };
    for (const auto & [component_key, component_entry] : model_obj->GetChemicalComponentEntryMap())
    {
        for (const auto & [atom_key, atom_entry] : component_entry->GetComponentAtomEntryMap())
        {
            batch.Execute([&](SQLiteWrapper & database)
            {
                database.Bind<std::string>(1, key_tag);
                database.Bind<ComponentKey>(2, component_key);
                database.Bind<AtomKey>(3, atom_key);
                database.Bind<std::string>(4, atom_entry.atom_id);
                database.Bind<int>(5, static_cast<int>(atom_entry.element_type));
                database.Bind<int>(6, static_cast<int>(atom_entry.aromatic_atom_flag));
                database.Bind<int>(7, static_cast<int>(atom_entry.stereo_config));
            });
        }
    }
}

void ModelObjectDAOv2::SaveComponentBondEntryList(const ModelObject * model_obj, const std::string & key_tag)
{
    model_io::SQLiteStatementBatch batch{ *m_database, std::string(model_io::kInsertModelComponentBondSql) };
    for (const auto & [component_key, component_entry] : model_obj->GetChemicalComponentEntryMap())
    {
        for (const auto & [bond_key, bond_entry] : component_entry->GetComponentBondEntryMap())
        {
            batch.Execute([&](SQLiteWrapper & database)
            {
                database.Bind<std::string>(1, key_tag);
                database.Bind<ComponentKey>(2, component_key);
                database.Bind<BondKey>(3, bond_key);
                database.Bind<std::string>(4, bond_entry.bond_id);
                database.Bind<int>(5, static_cast<int>(bond_entry.bond_type));
                database.Bind<int>(6, static_cast<int>(bond_entry.bond_order));
                database.Bind<int>(7, static_cast<int>(bond_entry.aromatic_atom_flag));
                database.Bind<int>(8, static_cast<int>(bond_entry.stereo_config));
            });
        }
    }
}

void ModelObjectDAOv2::SaveAtomObjectList(const ModelObject * model_obj, const std::string & key_tag)
{
    model_io::SQLiteStatementBatch batch{ *m_database, std::string(model_io::kInsertModelAtomSql) };
    for (const auto & atom_object : model_obj->GetAtomList())
    {
        batch.Execute([&](SQLiteWrapper & database)
        {
            database.Bind<std::string>(1, key_tag);
            database.Bind<int>(2, atom_object->GetSerialID());
            database.Bind<int>(3, atom_object->GetSequenceID());
            database.Bind<std::string>(4, atom_object->GetComponentID());
            database.Bind<std::string>(5, atom_object->GetAtomID());
            database.Bind<std::string>(6, atom_object->GetChainID());
            database.Bind<std::string>(7, atom_object->GetIndicator());
            database.Bind<double>(8, static_cast<double>(atom_object->GetOccupancy()));
            database.Bind<double>(9, static_cast<double>(atom_object->GetTemperature()));
            database.Bind<int>(10, static_cast<int>(atom_object->GetElement()));
            database.Bind<int>(11, static_cast<int>(atom_object->GetStructure()));
            database.Bind<int>(12, static_cast<int>(atom_object->GetSpecialAtomFlag()));
            database.Bind<double>(13, static_cast<double>(atom_object->GetPosition().at(0)));
            database.Bind<double>(14, static_cast<double>(atom_object->GetPosition().at(1)));
            database.Bind<double>(15, static_cast<double>(atom_object->GetPosition().at(2)));
            database.Bind<int>(16, static_cast<int>(atom_object->GetComponentKey()));
            database.Bind<int>(17, static_cast<int>(atom_object->GetAtomKey()));
        });
    }
}

void ModelObjectDAOv2::SaveBondObjectList(const ModelObject * model_obj, const std::string & key_tag)
{
    model_io::SQLiteStatementBatch batch{ *m_database, std::string(model_io::kInsertModelBondSql) };
    for (const auto & bond_object : model_obj->GetBondList())
    {
        batch.Execute([&](SQLiteWrapper & database)
        {
            database.Bind<std::string>(1, key_tag);
            database.Bind<int>(2, bond_object->GetAtomSerialID1());
            database.Bind<int>(3, bond_object->GetAtomSerialID2());
            database.Bind<int>(4, bond_object->GetBondKey());
            database.Bind<int>(5, static_cast<int>(bond_object->GetBondType()));
            database.Bind<int>(6, static_cast<int>(bond_object->GetBondOrder()));
            database.Bind<int>(7, static_cast<int>(bond_object->GetSpecialBondFlag()));
        });
    }
}

void ModelObjectDAOv2::LoadModelObjectRow(ModelObject * model_obj, const std::string & key_tag)
{
    m_database->Prepare(std::string(model_io::kSelectModelObjectSql));
    SQLiteWrapper::StatementGuard guard(*m_database);
    m_database->Bind<std::string>(1, key_tag);
    const auto rc{ m_database->StepNext() };
    if (rc == SQLiteWrapper::StepDone())
    {
        throw std::runtime_error("Cannot find the row with key_tag = " + key_tag);
    }
    if (rc != SQLiteWrapper::StepRow())
    {
        throw std::runtime_error("Step failed: " + m_database->ErrorMessage());
    }

    const auto atom_size{ m_database->GetColumn<int>(1) };
    model_obj->SetKeyTag(m_database->GetColumn<std::string>(0));
    model_obj->SetPdbID(m_database->GetColumn<std::string>(2));
    model_obj->SetEmdID(m_database->GetColumn<std::string>(3));
    model_obj->SetResolution(m_database->GetColumn<double>(4));
    model_obj->SetResolutionMethod(m_database->GetColumn<std::string>(5));
    if (atom_size != static_cast<int>(model_obj->GetNumberOfAtom()))
    {
        throw std::runtime_error(
            "The number of atoms in the model object does not match the database record.");
    }
}

void ModelObjectDAOv2::LoadChainMap(ModelObject * model_obj, const std::string & key_tag)
{
    std::unordered_map<std::string, std::vector<std::string>> chain_map;
    m_database->Prepare(std::string(model_io::kSelectModelChainMapSql));
    SQLiteWrapper::StatementGuard guard(*m_database);
    m_database->Bind<std::string>(1, key_tag);
    while (true)
    {
        const auto rc{ m_database->StepNext() };
        if (rc == SQLiteWrapper::StepDone())
        {
            break;
        }
        if (rc != SQLiteWrapper::StepRow())
        {
            throw std::runtime_error("Step failed: " + m_database->ErrorMessage());
        }
        auto & chain_id_list{ chain_map[m_database->GetColumn<std::string>(0)] };
        const auto ordinal{ m_database->GetColumn<int>(1) };
        if (static_cast<int>(chain_id_list.size()) <= ordinal)
        {
            chain_id_list.resize(static_cast<size_t>(ordinal + 1));
        }
        chain_id_list.at(static_cast<size_t>(ordinal)) = m_database->GetColumn<std::string>(2);
    }
    model_obj->SetChainIDListMap(chain_map);
}

void ModelObjectDAOv2::LoadChemicalComponentEntryList(ModelObject * model_obj, const std::string & key_tag)
{
    m_database->Prepare(std::string(model_io::kSelectModelComponentSql));
    SQLiteWrapper::StatementGuard guard(*m_database);
    m_database->Bind<std::string>(1, key_tag);
    while (true)
    {
        const auto rc{ m_database->StepNext() };
        if (rc == SQLiteWrapper::StepDone())
        {
            break;
        }
        if (rc != SQLiteWrapper::StepRow())
        {
            throw std::runtime_error("Step failed: " + m_database->ErrorMessage());
        }

        auto component_entry{ std::make_unique<ChemicalComponentEntry>() };
        const auto component_key{ m_database->GetColumn<ComponentKey>(0) };
        const auto component_id{ m_database->GetColumn<std::string>(1) };
        component_entry->SetComponentId(component_id);
        component_entry->SetComponentName(m_database->GetColumn<std::string>(2));
        component_entry->SetComponentType(m_database->GetColumn<std::string>(3));
        component_entry->SetComponentFormula(m_database->GetColumn<std::string>(4));
        component_entry->SetComponentMolecularWeight(
            static_cast<float>(m_database->GetColumn<double>(5)));
        component_entry->SetStandardMonomerFlag(static_cast<bool>(m_database->GetColumn<int>(6)));
        model_obj->AddChemicalComponentEntry(component_key, std::move(component_entry));
        model_obj->GetComponentKeySystemPtr()->RegisterComponent(component_id, component_key);
    }
}

void ModelObjectDAOv2::LoadComponentAtomEntryList(ModelObject * model_obj, const std::string & key_tag)
{
    m_database->Prepare(std::string(model_io::kSelectModelComponentAtomSql));
    SQLiteWrapper::StatementGuard guard(*m_database);
    m_database->Bind<std::string>(1, key_tag);
    while (true)
    {
        const auto rc{ m_database->StepNext() };
        if (rc == SQLiteWrapper::StepDone())
        {
            break;
        }
        if (rc != SQLiteWrapper::StepRow())
        {
            throw std::runtime_error("Step failed: " + m_database->ErrorMessage());
        }

        const auto component_key{ m_database->GetColumn<ComponentKey>(0) };
        if (model_obj->GetChemicalComponentEntryMap().find(component_key)
            == model_obj->GetChemicalComponentEntryMap().end())
        {
            continue;
        }
        const ComponentAtomEntry atom_entry{
            m_database->GetColumn<std::string>(2),
            static_cast<Element>(m_database->GetColumn<int>(3)),
            static_cast<bool>(m_database->GetColumn<int>(4)),
            static_cast<StereoChemistry>(m_database->GetColumn<int>(5))
        };
        const auto atom_key{ m_database->GetColumn<AtomKey>(1) };
        model_obj->GetChemicalComponentEntryMap().at(component_key)->AddComponentAtomEntry(
            atom_key, atom_entry);
        model_obj->GetAtomKeySystemPtr()->RegisterAtom(atom_entry.atom_id, atom_key);
    }
}

void ModelObjectDAOv2::LoadComponentBondEntryList(ModelObject * model_obj, const std::string & key_tag)
{
    m_database->Prepare(std::string(model_io::kSelectModelComponentBondSql));
    SQLiteWrapper::StatementGuard guard(*m_database);
    m_database->Bind<std::string>(1, key_tag);
    while (true)
    {
        const auto rc{ m_database->StepNext() };
        if (rc == SQLiteWrapper::StepDone())
        {
            break;
        }
        if (rc != SQLiteWrapper::StepRow())
        {
            throw std::runtime_error("Step failed: " + m_database->ErrorMessage());
        }

        const auto component_key{ m_database->GetColumn<ComponentKey>(0) };
        if (model_obj->GetChemicalComponentEntryMap().find(component_key)
            == model_obj->GetChemicalComponentEntryMap().end())
        {
            continue;
        }
        ComponentBondEntry bond_entry;
        const auto bond_key{ m_database->GetColumn<BondKey>(1) };
        bond_entry.bond_id = m_database->GetColumn<std::string>(2);
        bond_entry.bond_type = static_cast<BondType>(m_database->GetColumn<int>(3));
        bond_entry.bond_order = static_cast<BondOrder>(m_database->GetColumn<int>(4));
        bond_entry.aromatic_atom_flag = static_cast<bool>(m_database->GetColumn<int>(5));
        bond_entry.stereo_config = static_cast<StereoChemistry>(m_database->GetColumn<int>(6));
        model_obj->GetChemicalComponentEntryMap().at(component_key)->AddComponentBondEntry(
            bond_key, bond_entry);
        model_obj->GetBondKeySystemPtr()->RegisterBond(bond_entry.bond_id, bond_key);
    }
}

std::vector<std::unique_ptr<AtomObject>> ModelObjectDAOv2::LoadAtomObjectList(const std::string & key_tag)
{
    auto local_entry_map{ LoadAtomLocalPotentialEntryMap(key_tag) };
    std::vector<std::unique_ptr<AtomObject>> atom_object_list;

    m_database->Prepare(std::string(model_io::kSelectModelAtomSql));
    SQLiteWrapper::StatementGuard guard(*m_database);
    m_database->Bind<std::string>(1, key_tag);
    while (true)
    {
        const auto rc{ m_database->StepNext() };
        if (rc == SQLiteWrapper::StepDone())
        {
            break;
        }
        if (rc != SQLiteWrapper::StepRow())
        {
            throw std::runtime_error("Step failed: " + m_database->ErrorMessage());
        }

        auto atom_object{ std::make_unique<AtomObject>() };
        atom_object->SetSerialID(m_database->GetColumn<int>(0));
        atom_object->SetSequenceID(m_database->GetColumn<int>(1));
        atom_object->SetComponentID(m_database->GetColumn<std::string>(2));
        atom_object->SetAtomID(m_database->GetColumn<std::string>(3));
        atom_object->SetChainID(m_database->GetColumn<std::string>(4));
        atom_object->SetIndicator(m_database->GetColumn<std::string>(5));
        atom_object->SetOccupancy(static_cast<float>(m_database->GetColumn<double>(6)));
        atom_object->SetTemperature(static_cast<float>(m_database->GetColumn<double>(7)));
        atom_object->SetElement(static_cast<Element>(m_database->GetColumn<int>(8)));
        atom_object->SetStructure(static_cast<Structure>(m_database->GetColumn<int>(9)));
        atom_object->SetSpecialAtomFlag(static_cast<bool>(m_database->GetColumn<int>(10)));
        atom_object->SetPosition(
            static_cast<float>(m_database->GetColumn<double>(11)),
            static_cast<float>(m_database->GetColumn<double>(12)),
            static_cast<float>(m_database->GetColumn<double>(13)));
        atom_object->SetComponentKey(m_database->GetColumn<ComponentKey>(14));
        atom_object->SetAtomKey(m_database->GetColumn<AtomKey>(15));

        const auto serial_id{ atom_object->GetSerialID() };
        auto iter{ local_entry_map.find(serial_id) };
        if (iter != local_entry_map.end())
        {
            atom_object->SetSelectedFlag(true);
            atom_object->AddLocalPotentialEntry(std::move(iter->second));
        }
        else
        {
            atom_object->SetSelectedFlag(false);
        }
        atom_object_list.emplace_back(std::move(atom_object));
    }
    return atom_object_list;
}

std::vector<std::unique_ptr<BondObject>> ModelObjectDAOv2::LoadBondObjectList(
    const std::string & key_tag, const ModelObject * model_obj)
{
    auto local_entry_map{ LoadBondLocalPotentialEntryMap(key_tag) };
    std::vector<std::unique_ptr<BondObject>> bond_object_list;

    m_database->Prepare(std::string(model_io::kSelectModelBondSql));
    SQLiteWrapper::StatementGuard guard(*m_database);
    m_database->Bind<std::string>(1, key_tag);
    while (true)
    {
        const auto rc{ m_database->StepNext() };
        if (rc == SQLiteWrapper::StepDone())
        {
            break;
        }
        if (rc != SQLiteWrapper::StepRow())
        {
            throw std::runtime_error("Step failed: " + m_database->ErrorMessage());
        }

        auto atom_object_1{ model_obj->GetAtomPtr(m_database->GetColumn<int>(0)) };
        auto atom_object_2{ model_obj->GetAtomPtr(m_database->GetColumn<int>(1)) };
        auto bond_object{ std::make_unique<BondObject>(atom_object_1, atom_object_2) };
        bond_object->SetBondKey(m_database->GetColumn<BondKey>(2));
        bond_object->SetBondType(static_cast<BondType>(m_database->GetColumn<int>(3)));
        bond_object->SetBondOrder(static_cast<BondOrder>(m_database->GetColumn<int>(4)));
        bond_object->SetSpecialBondFlag(static_cast<bool>(m_database->GetColumn<int>(5)));

        const auto serial_id_pair{
            std::make_pair(atom_object_1->GetSerialID(), atom_object_2->GetSerialID()) };
        auto iter{ local_entry_map.find(serial_id_pair) };
        if (iter != local_entry_map.end())
        {
            bond_object->SetSelectedFlag(true);
            bond_object->AddLocalPotentialEntry(std::move(iter->second));
        }
        else
        {
            bond_object->SetSelectedFlag(false);
        }

        bond_object_list.emplace_back(std::move(bond_object));
    }
    return bond_object_list;
}

} // namespace rhbm_gem
