#include "ModelAnalysisPersistence.hpp"

#include "ModelSchemaSql.hpp"
#include "SQLiteStatementBatch.hpp"

#include <rhbm_gem/data/AtomClassifier.hpp>
#include <rhbm_gem/data/AtomObject.hpp>
#include <rhbm_gem/data/BondClassifier.hpp>
#include <rhbm_gem/data/BondObject.hpp>
#include <rhbm_gem/utils/ChemicalDataHelper.hpp>
#include <rhbm_gem/data/GroupPotentialEntry.hpp>
#include <rhbm_gem/data/LocalPotentialEntry.hpp>
#include <rhbm_gem/data/ModelObject.hpp>
#include "internal/SQLiteWrapper.hpp"

#include <map>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>

namespace rhbm_gem::model_io {

namespace {

void SaveAtomLocalPotentialEntryList(
    SQLiteWrapper & database,
    const ModelObject & model_obj,
    const std::string & key_tag)
{
    SQLiteStatementBatch batch{ database, std::string(kInsertModelAtomLocalSql) };
    for (const auto & atom_object : model_obj.GetAtomList())
    {
        auto entry{ atom_object->GetLocalPotentialEntry() };
        if (entry == nullptr) continue;

        batch.Execute([&](SQLiteWrapper & statement_db)
        {
            statement_db.Bind<std::string>(1, key_tag);
            statement_db.Bind<int>(2, atom_object->GetSerialID());
            statement_db.Bind<int>(3, entry->GetDistanceAndMapValueListSize());
            statement_db.Bind<std::vector<std::tuple<float, float>>>(4, entry->GetDistanceAndMapValueList());
            statement_db.Bind<double>(5, entry->GetAmplitudeEstimateOLS());
            statement_db.Bind<double>(6, entry->GetWidthEstimateOLS());
            statement_db.Bind<double>(7, entry->GetAmplitudeEstimateMDPDE());
            statement_db.Bind<double>(8, entry->GetWidthEstimateMDPDE());
            statement_db.Bind<double>(9, entry->GetAlphaR());
        });
    }
}

void SaveBondLocalPotentialEntryList(
    SQLiteWrapper & database,
    const ModelObject & model_obj,
    const std::string & key_tag)
{
    SQLiteStatementBatch batch{ database, std::string(kInsertModelBondLocalSql) };
    for (const auto & bond_object : model_obj.GetBondList())
    {
        auto entry{ bond_object->GetLocalPotentialEntry() };
        if (entry == nullptr) continue;

        batch.Execute([&](SQLiteWrapper & statement_db)
        {
            statement_db.Bind<std::string>(1, key_tag);
            statement_db.Bind<int>(2, bond_object->GetAtomSerialID1());
            statement_db.Bind<int>(3, bond_object->GetAtomSerialID2());
            statement_db.Bind<int>(4, entry->GetDistanceAndMapValueListSize());
            statement_db.Bind<std::vector<std::tuple<float, float>>>(5, entry->GetDistanceAndMapValueList());
            statement_db.Bind<double>(6, entry->GetAmplitudeEstimateOLS());
            statement_db.Bind<double>(7, entry->GetWidthEstimateOLS());
            statement_db.Bind<double>(8, entry->GetAmplitudeEstimateMDPDE());
            statement_db.Bind<double>(9, entry->GetWidthEstimateMDPDE());
            statement_db.Bind<double>(10, entry->GetAlphaR());
        });
    }
}

void SaveAtomLocalPotentialEntrySubList(
    SQLiteWrapper & database,
    const ModelObject & model_obj,
    const std::string & key_tag,
    const std::string & class_key)
{
    SQLiteStatementBatch batch{ database, std::string(kInsertModelAtomPosteriorSql) };
    for (const auto & atom_object : model_obj.GetAtomList())
    {
        auto entry{ atom_object->GetLocalPotentialEntry() };
        if (entry == nullptr) continue;

        batch.Execute([&](SQLiteWrapper & statement_db)
        {
            statement_db.Bind<std::string>(1, key_tag);
            statement_db.Bind<std::string>(2, class_key);
            statement_db.Bind<int>(3, atom_object->GetSerialID());
            statement_db.Bind<double>(4, entry->GetAmplitudeEstimatePosterior(class_key));
            statement_db.Bind<double>(5, entry->GetWidthEstimatePosterior(class_key));
            statement_db.Bind<double>(6, entry->GetAmplitudeVariancePosterior(class_key));
            statement_db.Bind<double>(7, entry->GetWidthVariancePosterior(class_key));
            statement_db.Bind<int>(8, static_cast<int>(entry->GetOutlierTag(class_key)));
            statement_db.Bind<double>(9, entry->GetStatisticalDistance(class_key));
        });
    }
}

void SaveBondLocalPotentialEntrySubList(
    SQLiteWrapper & database,
    const ModelObject & model_obj,
    const std::string & key_tag,
    const std::string & class_key)
{
    SQLiteStatementBatch batch{ database, std::string(kInsertModelBondPosteriorSql) };
    for (const auto & bond_object : model_obj.GetBondList())
    {
        auto entry{ bond_object->GetLocalPotentialEntry() };
        if (entry == nullptr) continue;

        batch.Execute([&](SQLiteWrapper & statement_db)
        {
            statement_db.Bind<std::string>(1, key_tag);
            statement_db.Bind<std::string>(2, class_key);
            statement_db.Bind<int>(3, bond_object->GetAtomSerialID1());
            statement_db.Bind<int>(4, bond_object->GetAtomSerialID2());
            statement_db.Bind<double>(5, entry->GetAmplitudeEstimatePosterior(class_key));
            statement_db.Bind<double>(6, entry->GetWidthEstimatePosterior(class_key));
            statement_db.Bind<double>(7, entry->GetAmplitudeVariancePosterior(class_key));
            statement_db.Bind<double>(8, entry->GetWidthVariancePosterior(class_key));
            statement_db.Bind<int>(9, static_cast<int>(entry->GetOutlierTag(class_key)));
            statement_db.Bind<double>(10, entry->GetStatisticalDistance(class_key));
        });
    }
}

void SaveAtomGroupPotentialEntryList(
    SQLiteWrapper & database,
    const GroupPotentialEntry & group_entry,
    const std::string & key_tag,
    const std::string & class_key)
{
    SQLiteStatementBatch batch{ database, std::string(kInsertModelAtomGroupSql) };
    for (const auto & group_key : group_entry.GetGroupKeySet())
    {
        batch.Execute([&](SQLiteWrapper & statement_db)
        {
            statement_db.Bind<std::string>(1, key_tag);
            statement_db.Bind<std::string>(2, class_key);
            statement_db.Bind<GroupKey>(3, group_key);
            statement_db.Bind<int>(4, group_entry.GetAtomObjectPtrListSize(group_key));
            statement_db.Bind<double>(5, std::get<0>(group_entry.GetGausEstimateMean(group_key)));
            statement_db.Bind<double>(6, std::get<1>(group_entry.GetGausEstimateMean(group_key)));
            statement_db.Bind<double>(7, std::get<0>(group_entry.GetGausEstimateMDPDE(group_key)));
            statement_db.Bind<double>(8, std::get<1>(group_entry.GetGausEstimateMDPDE(group_key)));
            statement_db.Bind<double>(9, std::get<0>(group_entry.GetGausEstimatePrior(group_key)));
            statement_db.Bind<double>(10, std::get<1>(group_entry.GetGausEstimatePrior(group_key)));
            statement_db.Bind<double>(11, std::get<0>(group_entry.GetGausVariancePrior(group_key)));
            statement_db.Bind<double>(12, std::get<1>(group_entry.GetGausVariancePrior(group_key)));
            statement_db.Bind<double>(13, group_entry.GetAlphaG(group_key));
        });
    }
}

void SaveBondGroupPotentialEntryList(
    SQLiteWrapper & database,
    const GroupPotentialEntry & group_entry,
    const std::string & key_tag,
    const std::string & class_key)
{
    SQLiteStatementBatch batch{ database, std::string(kInsertModelBondGroupSql) };
    for (const auto & group_key : group_entry.GetGroupKeySet())
    {
        batch.Execute([&](SQLiteWrapper & statement_db)
        {
            statement_db.Bind<std::string>(1, key_tag);
            statement_db.Bind<std::string>(2, class_key);
            statement_db.Bind<GroupKey>(3, group_key);
            statement_db.Bind<int>(4, group_entry.GetBondObjectPtrListSize(group_key));
            statement_db.Bind<double>(5, std::get<0>(group_entry.GetGausEstimateMean(group_key)));
            statement_db.Bind<double>(6, std::get<1>(group_entry.GetGausEstimateMean(group_key)));
            statement_db.Bind<double>(7, std::get<0>(group_entry.GetGausEstimateMDPDE(group_key)));
            statement_db.Bind<double>(8, std::get<1>(group_entry.GetGausEstimateMDPDE(group_key)));
            statement_db.Bind<double>(9, std::get<0>(group_entry.GetGausEstimatePrior(group_key)));
            statement_db.Bind<double>(10, std::get<1>(group_entry.GetGausEstimatePrior(group_key)));
            statement_db.Bind<double>(11, std::get<0>(group_entry.GetGausVariancePrior(group_key)));
            statement_db.Bind<double>(12, std::get<1>(group_entry.GetGausVariancePrior(group_key)));
            statement_db.Bind<double>(13, group_entry.GetAlphaG(group_key));
        });
    }
}

std::unordered_map<int, std::unique_ptr<LocalPotentialEntry>> LoadAtomLocalPotentialEntryMap(
    SQLiteWrapper & database,
    const std::string & key_tag);

std::map<std::pair<int, int>, std::unique_ptr<LocalPotentialEntry>> LoadBondLocalPotentialEntryMap(
    SQLiteWrapper & database,
    const std::string & key_tag);

void LoadAtomLocalPotentialEntrySubList(
    SQLiteWrapper & database,
    const std::string & key_tag,
    const std::string & class_key,
    std::unordered_map<int, std::unique_ptr<LocalPotentialEntry>> & entry_map)
{
    database.Prepare(std::string(kSelectModelAtomPosteriorSql));
    SQLiteWrapper::StatementGuard guard(database);
    database.Bind<std::string>(1, key_tag);
    database.Bind<std::string>(2, class_key);
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

        const auto serial_id{ database.GetColumn<int>(0) };
        auto iter{ entry_map.find(serial_id) };
        if (iter == entry_map.end()) continue;
        auto & entry{ iter->second };
        entry->AddGausEstimatePosterior(
            class_key, database.GetColumn<double>(1), database.GetColumn<double>(2));
        entry->AddGausVariancePosterior(
            class_key, database.GetColumn<double>(3), database.GetColumn<double>(4));
        entry->AddOutlierTag(class_key, static_cast<bool>(database.GetColumn<int>(5)));
        entry->AddStatisticalDistance(class_key, database.GetColumn<double>(6));
    }
}

void LoadBondLocalPotentialEntrySubList(
    SQLiteWrapper & database,
    const std::string & key_tag,
    const std::string & class_key,
    std::map<std::pair<int, int>, std::unique_ptr<LocalPotentialEntry>> & entry_map)
{
    database.Prepare(std::string(kSelectModelBondPosteriorSql));
    SQLiteWrapper::StatementGuard guard(database);
    database.Bind<std::string>(1, key_tag);
    database.Bind<std::string>(2, class_key);
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

        const auto atom_pair{
            std::make_pair(database.GetColumn<int>(0), database.GetColumn<int>(1)) };
        auto iter{ entry_map.find(atom_pair) };
        if (iter == entry_map.end()) continue;
        auto & entry{ iter->second };
        entry->AddGausEstimatePosterior(
            class_key, database.GetColumn<double>(2), database.GetColumn<double>(3));
        entry->AddGausVariancePosterior(
            class_key, database.GetColumn<double>(4), database.GetColumn<double>(5));
        entry->AddOutlierTag(class_key, static_cast<bool>(database.GetColumn<int>(6)));
        entry->AddStatisticalDistance(class_key, database.GetColumn<double>(7));
    }
}

std::unordered_map<int, std::unique_ptr<LocalPotentialEntry>> LoadAtomLocalPotentialEntryMap(
    SQLiteWrapper & database,
    const std::string & key_tag)
{
    std::unordered_map<int, std::unique_ptr<LocalPotentialEntry>> entry_map;
    database.Prepare(std::string(kSelectModelAtomLocalSql));
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

        auto entry{ std::make_unique<LocalPotentialEntry>() };
        const auto serial_id{ database.GetColumn<int>(0) };
        entry->AddDistanceAndMapValueList(
            database.GetColumn<std::vector<std::tuple<float, float>>>(2));
        entry->AddGausEstimateOLS(database.GetColumn<double>(3), database.GetColumn<double>(4));
        entry->AddGausEstimateMDPDE(database.GetColumn<double>(5), database.GetColumn<double>(6));
        entry->SetAlphaR(database.GetColumn<double>(7));
        entry_map[serial_id] = std::move(entry);
    }

    for (size_t i = 0; i < ChemicalDataHelper::GetGroupAtomClassCount(); i++)
    {
        LoadAtomLocalPotentialEntrySubList(
            database, key_tag, ChemicalDataHelper::GetGroupAtomClassKey(i), entry_map);
    }
    return entry_map;
}

std::map<std::pair<int, int>, std::unique_ptr<LocalPotentialEntry>> LoadBondLocalPotentialEntryMap(
    SQLiteWrapper & database,
    const std::string & key_tag)
{
    std::map<std::pair<int, int>, std::unique_ptr<LocalPotentialEntry>> entry_map;
    database.Prepare(std::string(kSelectModelBondLocalSql));
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

        auto entry{ std::make_unique<LocalPotentialEntry>() };
        const auto key{ std::make_pair(database.GetColumn<int>(0), database.GetColumn<int>(1)) };
        entry->AddDistanceAndMapValueList(
            database.GetColumn<std::vector<std::tuple<float, float>>>(3));
        entry->AddGausEstimateOLS(database.GetColumn<double>(4), database.GetColumn<double>(5));
        entry->AddGausEstimateMDPDE(database.GetColumn<double>(6), database.GetColumn<double>(7));
        entry->SetAlphaR(database.GetColumn<double>(8));
        entry_map[key] = std::move(entry);
    }

    for (size_t i = 0; i < ChemicalDataHelper::GetGroupBondClassCount(); i++)
    {
        LoadBondLocalPotentialEntrySubList(
            database, key_tag, ChemicalDataHelper::GetGroupBondClassKey(i), entry_map);
    }
    return entry_map;
}

void LoadAtomGroupPotentialEntryList(
    SQLiteWrapper & database,
    ModelObject & model_obj,
    const std::string & key_tag,
    const std::string & class_key)
{
    auto group_entry{ model_obj.GetAtomGroupPotentialEntry(class_key) };
    database.Prepare(std::string(kSelectModelAtomGroupSql));
    SQLiteWrapper::StatementGuard guard(database);
    database.Bind<std::string>(1, key_tag);
    database.Bind<std::string>(2, class_key);
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

        const auto group_key{ database.GetColumn<GroupKey>(0) };
        group_entry->InsertGroupKey(group_key);
        group_entry->ReserveAtomObjectPtrList(group_key, database.GetColumn<int>(1));
        group_entry->AddGausEstimateMean(group_key, database.GetColumn<double>(2), database.GetColumn<double>(3));
        group_entry->AddGausEstimateMDPDE(group_key, database.GetColumn<double>(4), database.GetColumn<double>(5));
        group_entry->AddGausEstimatePrior(group_key, database.GetColumn<double>(6), database.GetColumn<double>(7));
        group_entry->AddGausVariancePrior(group_key, database.GetColumn<double>(8), database.GetColumn<double>(9));
        group_entry->AddAlphaG(group_key, database.GetColumn<double>(10));
    }

    for (auto & atom : model_obj.GetSelectedAtomList())
    {
        group_entry->AddAtomObjectPtr(AtomClassifier::GetGroupKeyInClass(atom, class_key), atom);
    }
}

void LoadBondGroupPotentialEntryList(
    SQLiteWrapper & database,
    ModelObject & model_obj,
    const std::string & key_tag,
    const std::string & class_key)
{
    auto group_entry{ model_obj.GetBondGroupPotentialEntry(class_key) };
    database.Prepare(std::string(kSelectModelBondGroupSql));
    SQLiteWrapper::StatementGuard guard(database);
    database.Bind<std::string>(1, key_tag);
    database.Bind<std::string>(2, class_key);
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

        const auto group_key{ database.GetColumn<GroupKey>(0) };
        group_entry->InsertGroupKey(group_key);
        group_entry->ReserveBondObjectPtrList(group_key, database.GetColumn<int>(1));
        group_entry->AddGausEstimateMean(group_key, database.GetColumn<double>(2), database.GetColumn<double>(3));
        group_entry->AddGausEstimateMDPDE(group_key, database.GetColumn<double>(4), database.GetColumn<double>(5));
        group_entry->AddGausEstimatePrior(group_key, database.GetColumn<double>(6), database.GetColumn<double>(7));
        group_entry->AddGausVariancePrior(group_key, database.GetColumn<double>(8), database.GetColumn<double>(9));
        group_entry->AddAlphaG(group_key, database.GetColumn<double>(10));
    }

    for (auto & bond : model_obj.GetSelectedBondList())
    {
        group_entry->AddBondObjectPtr(BondClassifier::GetGroupKeyInClass(bond, class_key), bond);
    }
}

void ApplyAtomLocalPotentialEntries(
    ModelObject & model_obj,
    std::unordered_map<int, std::unique_ptr<LocalPotentialEntry>> entry_map)
{
    for (const auto & atom_object : model_obj.GetAtomList())
    {
        const auto serial_id{ atom_object->GetSerialID() };
        auto iter{ entry_map.find(serial_id) };
        if (iter != entry_map.end())
        {
            atom_object->SetSelectedFlag(true);
            atom_object->AddLocalPotentialEntry(std::move(iter->second));
        }
        else
        {
            atom_object->SetSelectedFlag(false);
        }
    }
}

void ApplyBondLocalPotentialEntries(
    ModelObject & model_obj,
    std::map<std::pair<int, int>, std::unique_ptr<LocalPotentialEntry>> entry_map)
{
    for (const auto & bond_object : model_obj.GetBondList())
    {
        const auto serial_id_pair{
            std::make_pair(bond_object->GetAtomSerialID1(), bond_object->GetAtomSerialID2()) };
        auto iter{ entry_map.find(serial_id_pair) };
        if (iter != entry_map.end())
        {
            bond_object->SetSelectedFlag(true);
            bond_object->AddLocalPotentialEntry(std::move(iter->second));
        }
        else
        {
            bond_object->SetSelectedFlag(false);
        }
    }
}

} // namespace

void SaveAnalysis(SQLiteWrapper & database, const ModelObject & model_obj, const std::string & key_tag)
{
    SaveAtomLocalPotentialEntryList(database, model_obj, key_tag);
    SaveBondLocalPotentialEntryList(database, model_obj, key_tag);

    for (const auto & [class_key, group_entry] : model_obj.GetAtomGroupPotentialEntryMap())
    {
        if (group_entry == nullptr || group_entry->GetGroupKeySet().empty())
        {
            continue;
        }
        SaveAtomLocalPotentialEntrySubList(database, model_obj, key_tag, class_key);
        SaveAtomGroupPotentialEntryList(database, *group_entry, key_tag, class_key);
    }
    for (const auto & [class_key, group_entry] : model_obj.GetBondGroupPotentialEntryMap())
    {
        if (group_entry == nullptr || group_entry->GetGroupKeySet().empty())
        {
            continue;
        }
        SaveBondLocalPotentialEntrySubList(database, model_obj, key_tag, class_key);
        SaveBondGroupPotentialEntryList(database, *group_entry, key_tag, class_key);
    }
}

void LoadAnalysis(SQLiteWrapper & database, ModelObject & model_obj, const std::string & key_tag)
{
    ApplyAtomLocalPotentialEntries(model_obj, LoadAtomLocalPotentialEntryMap(database, key_tag));
    ApplyBondLocalPotentialEntries(model_obj, LoadBondLocalPotentialEntryMap(database, key_tag));
    model_obj.Update();

    for (size_t i = 0; i < ChemicalDataHelper::GetGroupAtomClassCount(); i++)
    {
        auto class_key{ ChemicalDataHelper::GetGroupAtomClassKey(i) };
        auto group_entry{ std::make_unique<GroupPotentialEntry>() };
        model_obj.AddAtomGroupPotentialEntry(class_key, group_entry);
        LoadAtomGroupPotentialEntryList(database, model_obj, key_tag, class_key);
    }

    for (size_t i = 0; i < ChemicalDataHelper::GetGroupBondClassCount(); i++)
    {
        auto class_key{ ChemicalDataHelper::GetGroupBondClassKey(i) };
        auto group_entry{ std::make_unique<GroupPotentialEntry>() };
        model_obj.AddBondGroupPotentialEntry(class_key, group_entry);
        LoadBondGroupPotentialEntryList(database, model_obj, key_tag, class_key);
    }
}

} // namespace rhbm_gem::model_io
