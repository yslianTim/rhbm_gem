#include "ModelObjectDAOv2.hpp"

#include "AtomObject.hpp"
#include "BondObject.hpp"
#include "ChemicalDataHelper.hpp"
#include "DataObjectBase.hpp"
#include "GroupPotentialEntry.hpp"
#include "ModelObject.hpp"
#include "SQLiteWrapper.hpp"

#include "model_io/ModelSchemaSql.hpp"

#include <stdexcept>
#include <string>

namespace rhbm_gem {

ModelObjectDAOv2::ModelObjectDAOv2(SQLiteWrapper * db_manager) :
    m_database{ db_manager }
{
}

ModelObjectDAOv2::~ModelObjectDAOv2()
{
}

void ModelObjectDAOv2::EnsureSchema(SQLiteWrapper & database)
{
    for (const auto create_sql : model_io::kCreateModelTableSqlList)
    {
        database.Execute(std::string(create_sql));
    }
}

void ModelObjectDAOv2::Save(const DataObjectBase * obj, const std::string & key_tag)
{
    auto model_obj{ dynamic_cast<const ModelObject *>(obj) };
    if (!model_obj)
    {
        throw std::runtime_error(
            "ModelObjectDAOv2::Save() failed: object is not a ModelObject instance.");
    }

    EnsureSchema(*m_database);
    for (const auto table_name : model_io::kModelTablesScopedByKey)
    {
        DeleteRowsForKey(std::string(table_name), key_tag);
    }

    SaveModelObjectRow(model_obj, key_tag);
    SaveChainMap(model_obj, key_tag);
    SaveChemicalComponentEntryList(model_obj, key_tag);
    SaveComponentAtomEntryList(model_obj, key_tag);
    SaveComponentBondEntryList(model_obj, key_tag);
    SaveAtomObjectList(model_obj, key_tag);
    SaveBondObjectList(model_obj, key_tag);
    SaveAtomLocalPotentialEntryList(model_obj, key_tag);
    SaveBondLocalPotentialEntryList(model_obj, key_tag);

    for (const auto & [class_key, group_entry] : model_obj->GetAtomGroupPotentialEntryMap())
    {
        if (group_entry == nullptr || group_entry->GetGroupKeySet().empty())
        {
            continue;
        }
        SaveAtomLocalPotentialEntrySubList(model_obj, key_tag, class_key);
        SaveAtomGroupPotentialEntryList(group_entry.get(), key_tag, class_key);
    }
    for (const auto & [class_key, group_entry] : model_obj->GetBondGroupPotentialEntryMap())
    {
        if (group_entry == nullptr || group_entry->GetGroupKeySet().empty())
        {
            continue;
        }
        SaveBondLocalPotentialEntrySubList(model_obj, key_tag, class_key);
        SaveBondGroupPotentialEntryList(group_entry.get(), key_tag, class_key);
    }
}

std::unique_ptr<DataObjectBase> ModelObjectDAOv2::Load(const std::string & key_tag)
{
    auto model_object{ std::make_unique<ModelObject>() };
    LoadChemicalComponentEntryList(model_object.get(), key_tag);
    LoadComponentAtomEntryList(model_object.get(), key_tag);
    LoadComponentBondEntryList(model_object.get(), key_tag);
    model_object->SetAtomList(LoadAtomObjectList(key_tag));
    model_object->SetBondList(LoadBondObjectList(key_tag, model_object.get()));
    LoadModelObjectRow(model_object.get(), key_tag);
    LoadChainMap(model_object.get(), key_tag);

    for (size_t i = 0; i < ChemicalDataHelper::GetGroupAtomClassCount(); i++)
    {
        auto class_key{ ChemicalDataHelper::GetGroupAtomClassKey(i) };
        auto group_entry{ std::make_unique<GroupPotentialEntry>() };
        model_object->AddAtomGroupPotentialEntry(class_key, group_entry);
        LoadAtomGroupPotentialEntryList(model_object.get(), key_tag, class_key);
    }

    for (size_t i = 0; i < ChemicalDataHelper::GetGroupBondClassCount(); i++)
    {
        auto class_key{ ChemicalDataHelper::GetGroupBondClassKey(i) };
        auto group_entry{ std::make_unique<GroupPotentialEntry>() };
        model_object->AddBondGroupPotentialEntry(class_key, group_entry);
        LoadBondGroupPotentialEntryList(model_object.get(), key_tag, class_key);
    }

    return model_object;
}

void ModelObjectDAOv2::DeleteRowsForKey(const std::string & table_name, const std::string & key_tag) const
{
    m_database->Prepare(
        std::string(model_io::kDeleteRowsForKeySqlPrefix)
        + table_name
        + std::string(model_io::kDeleteRowsForKeySqlSuffix));
    SQLiteWrapper::StatementGuard guard(*m_database);
    m_database->Bind<std::string>(1, key_tag);
    m_database->StepOnce();
}

} // namespace rhbm_gem
