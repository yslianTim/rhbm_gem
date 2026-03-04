#pragma once

#include <map>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "DataObjectDAOBase.hpp"

namespace rhbm_gem {

class SQLiteWrapper;
class ModelObject;
class AtomObject;
class BondObject;
class LocalPotentialEntry;
class GroupPotentialEntry;

class ModelObjectDAOv2 : public DataObjectDAOBase
{
    SQLiteWrapper * m_database;

public:
    explicit ModelObjectDAOv2(SQLiteWrapper * db_manager);
    ~ModelObjectDAOv2();

    static void EnsureSchema(SQLiteWrapper & database);

    void Save(const DataObjectBase * obj, const std::string & key_tag) override;
    std::unique_ptr<DataObjectBase> Load(const std::string & key_tag) override;

private:
    void DeleteRowsForKey(const std::string & table_name, const std::string & key_tag) const;

    void SaveModelObjectRow(const ModelObject * model_obj, const std::string & key_tag);
    void SaveChainMap(const ModelObject * model_obj, const std::string & key_tag);
    void SaveChemicalComponentEntryList(const ModelObject * model_obj, const std::string & key_tag);
    void SaveComponentAtomEntryList(const ModelObject * model_obj, const std::string & key_tag);
    void SaveComponentBondEntryList(const ModelObject * model_obj, const std::string & key_tag);
    void SaveAtomObjectList(const ModelObject * model_obj, const std::string & key_tag);
    void SaveBondObjectList(const ModelObject * model_obj, const std::string & key_tag);
    void SaveAtomLocalPotentialEntryList(const ModelObject * model_obj, const std::string & key_tag);
    void SaveBondLocalPotentialEntryList(const ModelObject * model_obj, const std::string & key_tag);
    void SaveAtomLocalPotentialEntrySubList(
        const ModelObject * model_obj, const std::string & key_tag, const std::string & class_key);
    void SaveBondLocalPotentialEntrySubList(
        const ModelObject * model_obj, const std::string & key_tag, const std::string & class_key);
    void SaveAtomGroupPotentialEntryList(
        const GroupPotentialEntry * group_entry,
        const std::string & key_tag,
        const std::string & class_key);
    void SaveBondGroupPotentialEntryList(
        const GroupPotentialEntry * group_entry,
        const std::string & key_tag,
        const std::string & class_key);

    void LoadModelObjectRow(ModelObject * model_obj, const std::string & key_tag);
    void LoadChainMap(ModelObject * model_obj, const std::string & key_tag);
    void LoadChemicalComponentEntryList(ModelObject * model_obj, const std::string & key_tag);
    void LoadComponentAtomEntryList(ModelObject * model_obj, const std::string & key_tag);
    void LoadComponentBondEntryList(ModelObject * model_obj, const std::string & key_tag);
    std::vector<std::unique_ptr<AtomObject>> LoadAtomObjectList(const std::string & key_tag);
    std::vector<std::unique_ptr<BondObject>> LoadBondObjectList(
        const std::string & key_tag, const ModelObject * model_obj);
    std::unordered_map<int, std::unique_ptr<LocalPotentialEntry>>
    LoadAtomLocalPotentialEntryMap(const std::string & key_tag);
    std::map<std::pair<int, int>, std::unique_ptr<LocalPotentialEntry>>
    LoadBondLocalPotentialEntryMap(const std::string & key_tag);
    void LoadAtomLocalPotentialEntrySubList(
        const std::string & key_tag,
        const std::string & class_key,
        std::unordered_map<int, std::unique_ptr<LocalPotentialEntry>> & entry_map);
    void LoadBondLocalPotentialEntrySubList(
        const std::string & key_tag,
        const std::string & class_key,
        std::map<std::pair<int, int>, std::unique_ptr<LocalPotentialEntry>> & entry_map);
    void LoadAtomGroupPotentialEntryList(
        ModelObject * model_obj, const std::string & key_tag, const std::string & class_key);
    void LoadBondGroupPotentialEntryList(
        ModelObject * model_obj, const std::string & key_tag, const std::string & class_key);
};

} // namespace rhbm_gem
