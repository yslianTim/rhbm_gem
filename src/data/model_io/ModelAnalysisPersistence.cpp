#include "ModelSchemaSql.hpp"
#include "SQLiteStatementBatch.hpp"

#include "AtomClassifier.hpp"
#include "AtomObject.hpp"
#include "BondClassifier.hpp"
#include "BondObject.hpp"
#include "ChemicalDataHelper.hpp"
#include "GroupPotentialEntry.hpp"
#include "LocalPotentialEntry.hpp"
#include "ModelObject.hpp"
#include "ModelObjectDAOv2.hpp"
#include "SQLiteWrapper.hpp"

#include <map>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>

namespace rhbm_gem {

void ModelObjectDAOv2::SaveAtomLocalPotentialEntryList(const ModelObject * model_obj, const std::string & key_tag)
{
    model_io::SQLiteStatementBatch batch{ *m_database, std::string(model_io::kInsertModelAtomLocalSql) };
    for (const auto & atom_object : model_obj->GetAtomList())
    {
        auto entry{ atom_object->GetLocalPotentialEntry() };
        if (entry == nullptr) continue;

        batch.Execute([&](SQLiteWrapper & database)
        {
            database.Bind<std::string>(1, key_tag);
            database.Bind<int>(2, atom_object->GetSerialID());
            database.Bind<int>(3, entry->GetDistanceAndMapValueListSize());
            database.Bind<std::vector<std::tuple<float, float>>>(4, entry->GetDistanceAndMapValueList());
            database.Bind<double>(5, entry->GetAmplitudeEstimateOLS());
            database.Bind<double>(6, entry->GetWidthEstimateOLS());
            database.Bind<double>(7, entry->GetAmplitudeEstimateMDPDE());
            database.Bind<double>(8, entry->GetWidthEstimateMDPDE());
            database.Bind<double>(9, entry->GetAlphaR());
        });
    }
}

void ModelObjectDAOv2::SaveBondLocalPotentialEntryList(const ModelObject * model_obj, const std::string & key_tag)
{
    model_io::SQLiteStatementBatch batch{ *m_database, std::string(model_io::kInsertModelBondLocalSql) };
    for (const auto & bond_object : model_obj->GetBondList())
    {
        auto entry{ bond_object->GetLocalPotentialEntry() };
        if (entry == nullptr) continue;

        batch.Execute([&](SQLiteWrapper & database)
        {
            database.Bind<std::string>(1, key_tag);
            database.Bind<int>(2, bond_object->GetAtomSerialID1());
            database.Bind<int>(3, bond_object->GetAtomSerialID2());
            database.Bind<int>(4, entry->GetDistanceAndMapValueListSize());
            database.Bind<std::vector<std::tuple<float, float>>>(5, entry->GetDistanceAndMapValueList());
            database.Bind<double>(6, entry->GetAmplitudeEstimateOLS());
            database.Bind<double>(7, entry->GetWidthEstimateOLS());
            database.Bind<double>(8, entry->GetAmplitudeEstimateMDPDE());
            database.Bind<double>(9, entry->GetWidthEstimateMDPDE());
            database.Bind<double>(10, entry->GetAlphaR());
        });
    }
}

void ModelObjectDAOv2::SaveAtomLocalPotentialEntrySubList(
    const ModelObject * model_obj, const std::string & key_tag, const std::string & class_key)
{
    model_io::SQLiteStatementBatch batch{ *m_database, std::string(model_io::kInsertModelAtomPosteriorSql) };
    for (const auto & atom_object : model_obj->GetAtomList())
    {
        auto entry{ atom_object->GetLocalPotentialEntry() };
        if (entry == nullptr) continue;

        batch.Execute([&](SQLiteWrapper & database)
        {
            database.Bind<std::string>(1, key_tag);
            database.Bind<std::string>(2, class_key);
            database.Bind<int>(3, atom_object->GetSerialID());
            database.Bind<double>(4, entry->GetAmplitudeEstimatePosterior(class_key));
            database.Bind<double>(5, entry->GetWidthEstimatePosterior(class_key));
            database.Bind<double>(6, entry->GetAmplitudeVariancePosterior(class_key));
            database.Bind<double>(7, entry->GetWidthVariancePosterior(class_key));
            database.Bind<int>(8, static_cast<int>(entry->GetOutlierTag(class_key)));
            database.Bind<double>(9, entry->GetStatisticalDistance(class_key));
        });
    }
}

void ModelObjectDAOv2::SaveBondLocalPotentialEntrySubList(
    const ModelObject * model_obj, const std::string & key_tag, const std::string & class_key)
{
    model_io::SQLiteStatementBatch batch{ *m_database, std::string(model_io::kInsertModelBondPosteriorSql) };
    for (const auto & bond_object : model_obj->GetBondList())
    {
        auto entry{ bond_object->GetLocalPotentialEntry() };
        if (entry == nullptr) continue;

        batch.Execute([&](SQLiteWrapper & database)
        {
            database.Bind<std::string>(1, key_tag);
            database.Bind<std::string>(2, class_key);
            database.Bind<int>(3, bond_object->GetAtomSerialID1());
            database.Bind<int>(4, bond_object->GetAtomSerialID2());
            database.Bind<double>(5, entry->GetAmplitudeEstimatePosterior(class_key));
            database.Bind<double>(6, entry->GetWidthEstimatePosterior(class_key));
            database.Bind<double>(7, entry->GetAmplitudeVariancePosterior(class_key));
            database.Bind<double>(8, entry->GetWidthVariancePosterior(class_key));
            database.Bind<int>(9, static_cast<int>(entry->GetOutlierTag(class_key)));
            database.Bind<double>(10, entry->GetStatisticalDistance(class_key));
        });
    }
}

void ModelObjectDAOv2::SaveAtomGroupPotentialEntryList(
    const GroupPotentialEntry * group_entry, const std::string & key_tag, const std::string & class_key)
{
    model_io::SQLiteStatementBatch batch{ *m_database, std::string(model_io::kInsertModelAtomGroupSql) };
    for (const auto & group_key : group_entry->GetGroupKeySet())
    {
        batch.Execute([&](SQLiteWrapper & database)
        {
            database.Bind<std::string>(1, key_tag);
            database.Bind<std::string>(2, class_key);
            database.Bind<GroupKey>(3, group_key);
            database.Bind<int>(4, group_entry->GetAtomObjectPtrListSize(group_key));
            database.Bind<double>(5, std::get<0>(group_entry->GetGausEstimateMean(group_key)));
            database.Bind<double>(6, std::get<1>(group_entry->GetGausEstimateMean(group_key)));
            database.Bind<double>(7, std::get<0>(group_entry->GetGausEstimateMDPDE(group_key)));
            database.Bind<double>(8, std::get<1>(group_entry->GetGausEstimateMDPDE(group_key)));
            database.Bind<double>(9, std::get<0>(group_entry->GetGausEstimatePrior(group_key)));
            database.Bind<double>(10, std::get<1>(group_entry->GetGausEstimatePrior(group_key)));
            database.Bind<double>(11, std::get<0>(group_entry->GetGausVariancePrior(group_key)));
            database.Bind<double>(12, std::get<1>(group_entry->GetGausVariancePrior(group_key)));
            database.Bind<double>(13, group_entry->GetAlphaG(group_key));
        });
    }
}

void ModelObjectDAOv2::SaveBondGroupPotentialEntryList(
    const GroupPotentialEntry * group_entry, const std::string & key_tag, const std::string & class_key)
{
    model_io::SQLiteStatementBatch batch{ *m_database, std::string(model_io::kInsertModelBondGroupSql) };
    for (const auto & group_key : group_entry->GetGroupKeySet())
    {
        batch.Execute([&](SQLiteWrapper & database)
        {
            database.Bind<std::string>(1, key_tag);
            database.Bind<std::string>(2, class_key);
            database.Bind<GroupKey>(3, group_key);
            database.Bind<int>(4, group_entry->GetBondObjectPtrListSize(group_key));
            database.Bind<double>(5, std::get<0>(group_entry->GetGausEstimateMean(group_key)));
            database.Bind<double>(6, std::get<1>(group_entry->GetGausEstimateMean(group_key)));
            database.Bind<double>(7, std::get<0>(group_entry->GetGausEstimateMDPDE(group_key)));
            database.Bind<double>(8, std::get<1>(group_entry->GetGausEstimateMDPDE(group_key)));
            database.Bind<double>(9, std::get<0>(group_entry->GetGausEstimatePrior(group_key)));
            database.Bind<double>(10, std::get<1>(group_entry->GetGausEstimatePrior(group_key)));
            database.Bind<double>(11, std::get<0>(group_entry->GetGausVariancePrior(group_key)));
            database.Bind<double>(12, std::get<1>(group_entry->GetGausVariancePrior(group_key)));
            database.Bind<double>(13, group_entry->GetAlphaG(group_key));
        });
    }
}

std::unordered_map<int, std::unique_ptr<LocalPotentialEntry>>
ModelObjectDAOv2::LoadAtomLocalPotentialEntryMap(const std::string & key_tag)
{
    std::unordered_map<int, std::unique_ptr<LocalPotentialEntry>> entry_map;
    m_database->Prepare(std::string(model_io::kSelectModelAtomLocalSql));
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

        auto entry{ std::make_unique<LocalPotentialEntry>() };
        const auto serial_id{ m_database->GetColumn<int>(0) };
        entry->AddDistanceAndMapValueList(
            m_database->GetColumn<std::vector<std::tuple<float, float>>>(2));
        entry->AddGausEstimateOLS(m_database->GetColumn<double>(3), m_database->GetColumn<double>(4));
        entry->AddGausEstimateMDPDE(m_database->GetColumn<double>(5), m_database->GetColumn<double>(6));
        entry->SetAlphaR(m_database->GetColumn<double>(7));
        entry_map[serial_id] = std::move(entry);
    }

    for (size_t i = 0; i < ChemicalDataHelper::GetGroupAtomClassCount(); i++)
    {
        LoadAtomLocalPotentialEntrySubList(
            key_tag, ChemicalDataHelper::GetGroupAtomClassKey(i), entry_map);
    }
    return entry_map;
}

std::map<std::pair<int, int>, std::unique_ptr<LocalPotentialEntry>>
ModelObjectDAOv2::LoadBondLocalPotentialEntryMap(const std::string & key_tag)
{
    std::map<std::pair<int, int>, std::unique_ptr<LocalPotentialEntry>> entry_map;
    m_database->Prepare(std::string(model_io::kSelectModelBondLocalSql));
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

        auto entry{ std::make_unique<LocalPotentialEntry>() };
        const auto key{ std::make_pair(m_database->GetColumn<int>(0), m_database->GetColumn<int>(1)) };
        entry->AddDistanceAndMapValueList(
            m_database->GetColumn<std::vector<std::tuple<float, float>>>(3));
        entry->AddGausEstimateOLS(m_database->GetColumn<double>(4), m_database->GetColumn<double>(5));
        entry->AddGausEstimateMDPDE(m_database->GetColumn<double>(6), m_database->GetColumn<double>(7));
        entry->SetAlphaR(m_database->GetColumn<double>(8));
        entry_map[key] = std::move(entry);
    }

    for (size_t i = 0; i < ChemicalDataHelper::GetGroupBondClassCount(); i++)
    {
        LoadBondLocalPotentialEntrySubList(
            key_tag, ChemicalDataHelper::GetGroupBondClassKey(i), entry_map);
    }
    return entry_map;
}

void ModelObjectDAOv2::LoadAtomLocalPotentialEntrySubList(
    const std::string & key_tag,
    const std::string & class_key,
    std::unordered_map<int, std::unique_ptr<LocalPotentialEntry>> & entry_map)
{
    m_database->Prepare(std::string(model_io::kSelectModelAtomPosteriorSql));
    SQLiteWrapper::StatementGuard guard(*m_database);
    m_database->Bind<std::string>(1, key_tag);
    m_database->Bind<std::string>(2, class_key);
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

        const auto serial_id{ m_database->GetColumn<int>(0) };
        auto iter{ entry_map.find(serial_id) };
        if (iter == entry_map.end()) continue;
        auto & entry{ iter->second };
        entry->AddGausEstimatePosterior(
            class_key, m_database->GetColumn<double>(1), m_database->GetColumn<double>(2));
        entry->AddGausVariancePosterior(
            class_key, m_database->GetColumn<double>(3), m_database->GetColumn<double>(4));
        entry->AddOutlierTag(class_key, static_cast<bool>(m_database->GetColumn<int>(5)));
        entry->AddStatisticalDistance(class_key, m_database->GetColumn<double>(6));
    }
}

void ModelObjectDAOv2::LoadBondLocalPotentialEntrySubList(
    const std::string & key_tag,
    const std::string & class_key,
    std::map<std::pair<int, int>, std::unique_ptr<LocalPotentialEntry>> & entry_map)
{
    m_database->Prepare(std::string(model_io::kSelectModelBondPosteriorSql));
    SQLiteWrapper::StatementGuard guard(*m_database);
    m_database->Bind<std::string>(1, key_tag);
    m_database->Bind<std::string>(2, class_key);
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

        const auto atom_pair{
            std::make_pair(m_database->GetColumn<int>(0), m_database->GetColumn<int>(1)) };
        auto iter{ entry_map.find(atom_pair) };
        if (iter == entry_map.end()) continue;
        auto & entry{ iter->second };
        entry->AddGausEstimatePosterior(
            class_key, m_database->GetColumn<double>(2), m_database->GetColumn<double>(3));
        entry->AddGausVariancePosterior(
            class_key, m_database->GetColumn<double>(4), m_database->GetColumn<double>(5));
        entry->AddOutlierTag(class_key, static_cast<bool>(m_database->GetColumn<int>(6)));
        entry->AddStatisticalDistance(class_key, m_database->GetColumn<double>(7));
    }
}

void ModelObjectDAOv2::LoadAtomGroupPotentialEntryList(
    ModelObject * model_obj, const std::string & key_tag, const std::string & class_key)
{
    auto group_entry{ model_obj->GetAtomGroupPotentialEntry(class_key) };
    m_database->Prepare(std::string(model_io::kSelectModelAtomGroupSql));
    SQLiteWrapper::StatementGuard guard(*m_database);
    m_database->Bind<std::string>(1, key_tag);
    m_database->Bind<std::string>(2, class_key);
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

        const auto group_key{ m_database->GetColumn<GroupKey>(0) };
        group_entry->InsertGroupKey(group_key);
        group_entry->ReserveAtomObjectPtrList(group_key, m_database->GetColumn<int>(1));
        group_entry->AddGausEstimateMean(group_key, m_database->GetColumn<double>(2), m_database->GetColumn<double>(3));
        group_entry->AddGausEstimateMDPDE(group_key, m_database->GetColumn<double>(4), m_database->GetColumn<double>(5));
        group_entry->AddGausEstimatePrior(group_key, m_database->GetColumn<double>(6), m_database->GetColumn<double>(7));
        group_entry->AddGausVariancePrior(group_key, m_database->GetColumn<double>(8), m_database->GetColumn<double>(9));
        group_entry->AddAlphaG(group_key, m_database->GetColumn<double>(10));
    }

    for (auto & atom : model_obj->GetSelectedAtomList())
    {
        group_entry->AddAtomObjectPtr(AtomClassifier::GetGroupKeyInClass(atom, class_key), atom);
    }
}

void ModelObjectDAOv2::LoadBondGroupPotentialEntryList(
    ModelObject * model_obj, const std::string & key_tag, const std::string & class_key)
{
    auto group_entry{ model_obj->GetBondGroupPotentialEntry(class_key) };
    m_database->Prepare(std::string(model_io::kSelectModelBondGroupSql));
    SQLiteWrapper::StatementGuard guard(*m_database);
    m_database->Bind<std::string>(1, key_tag);
    m_database->Bind<std::string>(2, class_key);
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

        const auto group_key{ m_database->GetColumn<GroupKey>(0) };
        group_entry->InsertGroupKey(group_key);
        group_entry->ReserveBondObjectPtrList(group_key, m_database->GetColumn<int>(1));
        group_entry->AddGausEstimateMean(group_key, m_database->GetColumn<double>(2), m_database->GetColumn<double>(3));
        group_entry->AddGausEstimateMDPDE(group_key, m_database->GetColumn<double>(4), m_database->GetColumn<double>(5));
        group_entry->AddGausEstimatePrior(group_key, m_database->GetColumn<double>(6), m_database->GetColumn<double>(7));
        group_entry->AddGausVariancePrior(group_key, m_database->GetColumn<double>(8), m_database->GetColumn<double>(9));
        group_entry->AddAlphaG(group_key, m_database->GetColumn<double>(10));
    }

    for (auto & bond : model_obj->GetSelectedBondList())
    {
        group_entry->AddBondObjectPtr(BondClassifier::GetGroupKeyInClass(bond, class_key), bond);
    }
}

} // namespace rhbm_gem
