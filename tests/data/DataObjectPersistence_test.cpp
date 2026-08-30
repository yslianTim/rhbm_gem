#include <gtest/gtest.h>

#include <cstddef>
#include <filesystem>
#include <stdexcept>
#include <string>

#include "data/detail/ModelAnalysisData.hpp"
#include <rhbm_gem/data/io/DataRepository.hpp>
#include <rhbm_gem/data/io/ModelMapFileIO.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <rhbm_gem/data/object/ModelAnalysisView.hpp>
#include <rhbm_gem/utils/domain/ChemicalDataHelper.hpp>
#include "io/sqlite/SQLitePersistence.hpp"
#include "support/CommandTestHelpers.hpp"
#include "support/DataObjectTestSupport.hpp"

namespace rg = rhbm_gem;

namespace {

enum class PersistedObjectType
{
    Model,
    Map,
};

enum class RequestedLoadType
{
    Model,
    Map,
};

struct TypedMismatchCase
{
    const char * name;
    PersistedObjectType persisted_type;
    RequestedLoadType requested_type;
};

class SQLitePersistenceTypedMismatchTest
    : public testing::TestWithParam<TypedMismatchCase> {};

TEST_P(SQLitePersistenceTypedMismatchTest, ThrowsWhenCatalogRowTypeDoesNotMatchRequestedLoad)
{
    const auto & test_case{ GetParam() };
    const command_test::ScopedTempDir temp_dir{
        std::string("sqlite_persistence_typed_mismatch_") + test_case.name };
    const auto database_path{ temp_dir.path() / "typed_mismatch.sqlite" };

    rg::SQLitePersistence persistence{ database_path };
    if (test_case.persisted_type == PersistedObjectType::Map)
    {
        auto map{ data_test::MakeMapObject() };
        persistence.SaveMap(map, "shared_key");
    }
    else
    {
        const auto model_path{ command_test::TestDataPath("test_model.cif") };
        auto model{ rg::ReadModel(model_path) };
        persistence.SaveModel(*model, "shared_key");
    }

    if (test_case.requested_type == RequestedLoadType::Model)
    {
        EXPECT_THROW((void)persistence.LoadModel("shared_key"), std::runtime_error);
    }
    else
    {
        EXPECT_THROW((void)persistence.LoadMap("shared_key"), std::runtime_error);
    }
}

INSTANTIATE_TEST_SUITE_P(
    TypedCatalogMismatch,
    SQLitePersistenceTypedMismatchTest,
    testing::Values(
        TypedMismatchCase{
            "ModelLoadAgainstMapRow",
            PersistedObjectType::Map,
            RequestedLoadType::Model },
        TypedMismatchCase{
            "MapLoadAgainstModelRow",
            PersistedObjectType::Model,
            RequestedLoadType::Map }),
    [](const testing::TestParamInfo<TypedMismatchCase> & info)
    {
        return std::string(info.param.name);
    });

} // namespace

TEST(DataObjectPersistenceTest, SaveAndLoadModelWithoutRegistryState)
{
    const command_test::ScopedTempDir temp_dir{ "data_repository_schema_roundtrip" };
    const auto database_path{ temp_dir.path() / "repository.sqlite" };
    const auto model_path{ command_test::TestDataPath("test_model.cif") };

    rg::DataRepository repository{ database_path };
    auto model{ rg::ReadModel(model_path) };
    model->SetKeyTag("model");
    model->SetPdbID("MODEL_REPOSITORY");
    repository.SaveModel(*model, "model");

    auto loaded_model{ repository.LoadModel("model") };
    ASSERT_NE(loaded_model, nullptr);
    EXPECT_EQ(loaded_model->GetPdbID(), "MODEL_REPOSITORY");
    EXPECT_EQ(loaded_model->GetNumberOfAtom(), model->GetNumberOfAtom());
}

TEST(DataObjectPersistenceTest, StandardQScoresAndReferenceParametersRoundTripThroughRepository)
{
    const command_test::ScopedTempDir temp_dir{ "standard_qscore_roundtrip" };
    const auto database_path{ temp_dir.path() / "repository.sqlite" };

    rg::DataRepository repository{ database_path };
    auto model{ data_test::MakeModelWithBond() };
    ASSERT_EQ(model->GetAtomList().size(), 2u);
    model->SetStandardAverageQScore(0.625);
    model->SetReferenceHeight(1.25);
    model->SetReferenceOffset(-0.5);
    model->GetAtomList().at(0)->SetStandardQScore(0.5);
    model->GetAtomList().at(1)->SetStandardQScore(0.75);

    repository.SaveModel(*model, "model");
    auto loaded_model{ repository.LoadModel("model") };

    ASSERT_NE(loaded_model, nullptr);
    ASSERT_EQ(loaded_model->GetAtomList().size(), 2u);
    EXPECT_DOUBLE_EQ(loaded_model->GetStandardAverageQScore(), 0.625);
    EXPECT_DOUBLE_EQ(loaded_model->GetReferenceHeight(), 1.25);
    EXPECT_DOUBLE_EQ(loaded_model->GetReferenceOffset(), -0.5);
    EXPECT_DOUBLE_EQ(
        loaded_model->GetAtomList().at(0)->GetStandardQScore(),
        0.5);
    EXPECT_DOUBLE_EQ(
        loaded_model->GetAtomList().at(1)->GetStandardQScore(),
        0.75);
}

TEST(DataObjectPersistenceTest, SQLiteTypedSaveAndLoadMapRoundTrip)
{
    const command_test::ScopedTempDir temp_dir{ "sqlite_persistence_typed_map" };
    const auto database_path{ temp_dir.path() / "typed_map.sqlite" };

    rg::SQLitePersistence persistence{ database_path };
    auto map{ data_test::MakeMapObject() };
    map.SetKeyTag("map");

    persistence.SaveMap(map, "map");

    auto loaded_map{ persistence.LoadMap("map") };
    ASSERT_NE(loaded_map, nullptr);
    EXPECT_EQ(loaded_map->GetGridSize(), map.GetGridSize());
    EXPECT_FLOAT_EQ(loaded_map->GetMapValue(0), map.GetMapValue(0));
}

TEST(DataObjectPersistenceTest, FinalV2CatalogDatabaseRemainsLoadable)
{
    const command_test::ScopedTempDir temp_dir{ "data_schema_final_v2_loadable" };
    const auto database_path{ temp_dir.path() / "final_v2.sqlite" };
    const auto model_path{ command_test::TestDataPath("test_model.cif") };

    {
        rg::DataRepository repository{ database_path };
        auto model{ rg::ReadModel(model_path) };
        model->SetKeyTag("model");
        repository.SaveModel(*model, "model");
        data_test::SaveTinyMapThroughRepository(repository, "map", 3.0f);
    }

    data_test::ConvertLocalGaussianColumnsToLegacyFinal(database_path);
    data_test::ConvertAtomGroupGaussianColumnsToLegacyFinal(database_path);
    data_test::ConvertSamplingEntryColumnsToLegacyRawOnly(database_path);
    data_test::SetUserVersion(database_path, 5);

    rg::DataRepository repository{ database_path };
    EXPECT_NE(repository.LoadModel("model"), nullptr);
    EXPECT_NE(repository.LoadMap("map"), nullptr);
    EXPECT_EQ(data_test::GetUserVersion(database_path), 12);
}

TEST(DataObjectPersistenceTest, RawSamplingSelectionRoundTripPreservesUnselectedSamples)
{
    const command_test::ScopedTempDir temp_dir{ "data_schema_sampling_selection_roundtrip" };
    const auto database_path{ temp_dir.path() / "sampling_selection.sqlite" };

    {
        rg::DataRepository repository{ database_path };
        auto model{ data_test::MakeModelWithBond() };
        model->SetKeyTag("model");
        auto editor{ model->EditAnalysis().EnsureAtomLocalPotential(*model->GetAtomList().at(0)) };
        editor.SetRawSamplingEntries({
            LocalPotentialSample{ 6.0f, SamplingPoint{ 0.1f, { 0.0f, 0.0f, 0.0f }, true } },
            LocalPotentialSample{ 4.0f, SamplingPoint{ 0.2f, { 0.0f, 0.0f, 0.0f }, false } }
        });

        repository.SaveModel(*model, "model");
    }

    rg::DataRepository repository{ database_path };
    auto loaded_model{ repository.LoadModel("model") };
    ASSERT_NE(loaded_model, nullptr);
    const auto raw_all_entries{
        rg::AtomLocalPotentialView::RequireFor(
            *loaded_model->GetAtomList().at(0)).GetRawSamplingEntries(false)
    };
    const auto raw_selected_entries{
        rg::AtomLocalPotentialView::RequireFor(
            *loaded_model->GetAtomList().at(0)).GetRawSamplingEntries()
    };

    ASSERT_EQ(raw_all_entries.size(), 2u);
    EXPECT_TRUE(raw_all_entries.at(0).point.is_selected);
    EXPECT_FALSE(raw_all_entries.at(1).point.is_selected);
    ASSERT_EQ(raw_selected_entries.size(), 1u);
    EXPECT_FLOAT_EQ(raw_selected_entries.at(0).response, 6.0f);
}

TEST(DataObjectPersistenceTest, RawAndPeelingSamplingEntriesRoundTripPreservesSelection)
{
    const command_test::ScopedTempDir temp_dir{ "data_schema_raw_peeling_sampling_roundtrip" };
    const auto database_path{ temp_dir.path() / "raw_peeling_sampling.sqlite" };

    {
        rg::DataRepository repository{ database_path };
        auto model{ data_test::MakeModelWithBond() };
        model->SetKeyTag("model");
        auto editor{ model->EditAnalysis().EnsureAtomLocalPotential(*model->GetAtomList().at(0)) };
        editor.SetRawSamplingEntries({
            LocalPotentialSample{ 6.0f, SamplingPoint{ 0.1f, { 0.0f, 0.0f, 0.0f }, true } },
            LocalPotentialSample{ 4.0f, SamplingPoint{ 0.2f, { 0.0f, 0.0f, 0.0f }, false } }
        });
        editor.SetPeelingSamplingEntries({
            LocalPotentialSample{ 3.0f, SamplingPoint{ 0.1f, { 0.0f, 0.0f, 0.0f }, true } },
            LocalPotentialSample{ 5.0f, SamplingPoint{ 0.2f, { 0.0f, 0.0f, 0.0f }, false } }
        });
        editor.SetNeighborCountForPeeling(7);

        repository.SaveModel(*model, "model");
    }

    rg::DataRepository repository{ database_path };
    auto loaded_model{ repository.LoadModel("model") };
    ASSERT_NE(loaded_model, nullptr);
    const auto view{ rg::AtomLocalPotentialView::RequireFor(*loaded_model->GetAtomList().at(0)) };
    const auto raw_entries{ view.GetRawSamplingEntries(false) };
    const auto peeling_entries{ view.GetPeelingSamplingEntries() };
    const auto selected_peeling_entries{ view.GetPeelingSamplingEntries(true) };

    ASSERT_EQ(raw_entries.size(), 2u);
    EXPECT_FLOAT_EQ(raw_entries.at(0).response, 6.0f);
    EXPECT_FLOAT_EQ(raw_entries.at(1).response, 4.0f);
    EXPECT_TRUE(raw_entries.at(0).point.is_selected);
    EXPECT_FALSE(raw_entries.at(1).point.is_selected);
    ASSERT_EQ(peeling_entries.size(), 2u);
    EXPECT_FLOAT_EQ(peeling_entries.at(0).response, 3.0f);
    EXPECT_FLOAT_EQ(peeling_entries.at(1).response, 5.0f);
    EXPECT_TRUE(peeling_entries.at(0).point.is_selected);
    EXPECT_FALSE(peeling_entries.at(1).point.is_selected);
    ASSERT_EQ(selected_peeling_entries.size(), 1u);
    EXPECT_FLOAT_EQ(selected_peeling_entries.at(0).response, 3.0f);
    EXPECT_EQ(view.GetNeighborCountForPeeling(), 7);
}

TEST(DataObjectPersistenceTest, VersionEightSamplingColumnsMigrateToRawAndPeelingWithoutDataLoss)
{
    const command_test::ScopedTempDir temp_dir{ "data_schema_v8_raw_peeling_migration" };
    const auto database_path{ temp_dir.path() / "v8_raw_peeling_migration.sqlite" };

    {
        rg::DataRepository repository{ database_path };
        auto model{ data_test::MakeModelWithBond() };
        model->SetKeyTag("model");
        auto editor{
            model->EditAnalysis().EnsureAtomLocalPotential(
                *model->GetAtomList().at(0))
        };
        editor.SetRawSamplingEntries({
            LocalPotentialSample{
                8.0f,
                SamplingPoint{ 0.1f, { 1.0f, 2.0f, 3.0f }, true } },
            LocalPotentialSample{
                6.0f,
                SamplingPoint{ 0.2f, { 4.0f, 5.0f, 6.0f }, false } }
        });
        editor.SetPeelingSamplingEntries({
            LocalPotentialSample{
                3.5f,
                SamplingPoint{ 0.1f, { 1.0f, 2.0f, 3.0f }, false } },
            LocalPotentialSample{
                1.5f,
                SamplingPoint{ 0.2f, { 4.0f, 5.0f, 6.0f }, true } }
        });
        repository.SaveModel(*model, "model");
    }

    data_test::ConvertLocalGaussianColumnsToLegacyFinal(database_path);
    data_test::ConvertAtomGroupGaussianColumnsToLegacyFinal(database_path);
    data_test::ConvertSamplingEntryColumnsToLegacyRawAndUpdated(database_path);
    data_test::SetUserVersion(database_path, 8);
    ASSERT_TRUE(data_test::HasColumn(
        database_path, "model_atom_local_potential", "sampling_size"));
    ASSERT_TRUE(data_test::HasColumn(
        database_path, "model_atom_local_potential", "updated_sampling_size"));

    rg::DataRepository repository{ database_path };
    auto loaded_model{ repository.LoadModel("model") };
    ASSERT_NE(loaded_model, nullptr);

    EXPECT_EQ(data_test::GetUserVersion(database_path), 12);
    EXPECT_TRUE(data_test::HasColumn(
        database_path, "model_atom_local_potential", "raw_sampling_size"));
    EXPECT_TRUE(data_test::HasColumn(
        database_path,
        "model_atom_local_potential",
        "raw_distance_and_map_value_list"));
    EXPECT_TRUE(data_test::HasColumn(
        database_path, "model_atom_local_potential", "peeling_sampling_size"));
    EXPECT_TRUE(data_test::HasColumn(
        database_path,
        "model_atom_local_potential",
        "peeling_distance_and_map_value_list"));
    EXPECT_FALSE(data_test::HasColumn(
        database_path, "model_atom_local_potential", "sampling_size"));
    EXPECT_FALSE(data_test::HasColumn(
        database_path,
        "model_atom_local_potential",
        "distance_and_map_value_list"));
    EXPECT_FALSE(data_test::HasColumn(
        database_path, "model_atom_local_potential", "updated_sampling_size"));
    EXPECT_FALSE(data_test::HasColumn(
        database_path,
        "model_atom_local_potential",
        "updated_distance_and_map_value_list"));

    const auto view{
        rg::AtomLocalPotentialView::RequireFor(
            *loaded_model->GetAtomList().at(0))
    };
    const auto raw_entries{ view.GetRawSamplingEntries(false) };
    const auto peeling_entries{ view.GetPeelingSamplingEntries(false) };
    ASSERT_EQ(raw_entries.size(), 2u);
    EXPECT_FLOAT_EQ(raw_entries.at(0).response, 8.0f);
    EXPECT_FLOAT_EQ(raw_entries.at(1).response, 6.0f);
    EXPECT_TRUE(raw_entries.at(0).point.is_selected);
    EXPECT_FALSE(raw_entries.at(1).point.is_selected);
    ASSERT_EQ(peeling_entries.size(), 2u);
    EXPECT_FLOAT_EQ(peeling_entries.at(0).response, 3.5f);
    EXPECT_FLOAT_EQ(peeling_entries.at(1).response, 1.5f);
    EXPECT_FALSE(peeling_entries.at(0).point.is_selected);
    EXPECT_TRUE(peeling_entries.at(1).point.is_selected);
}

TEST(DataObjectPersistenceTest, GaussianOffsetRoundTripPreservesAnalysisResults)
{
    const command_test::ScopedTempDir temp_dir{ "data_schema_gaussian_offset_roundtrip" };
    const auto database_path{ temp_dir.path() / "gaussian_offset_roundtrip.sqlite" };

    int annotated_serial_id{ 0 };
    GroupKey group_key{};
    {
        rg::DataRepository repository{ database_path };
        auto model{ data_test::MakeModelWithBond() };
        model->SetKeyTag("model");
        model->SelectAllAtoms();
        model->ApplySymmetrySelection(false);

        auto analysis{ model->EditAnalysis() };
        auto local_editor{ analysis.EnsureAtomLocalPotential(*model->GetAtomList().at(0)) };
        rg::LocalGaussianResult local_result;
        local_result.alpha_r = 0.5;
        local_result.ols = rg::GaussianModel3DWithUncertainty{
            rg::GaussianModel3D{ 1.0, 0.6, 0.11 },
            rg::GaussianModel3DUncertainty{}
        };
        local_result.mdpde = rg::GaussianModel3DWithUncertainty{
            rg::GaussianModel3D{ 1.5, 0.8, 0.22 },
            rg::GaussianModel3DUncertainty{}
        };
        auto first_result{ local_result };
        first_result.alpha_r = 0.1;
        first_result.ols = rg::GaussianModel3DWithUncertainty{
            rg::GaussianModel3D{ 1.0, 0.6, 1.11 },
            rg::GaussianModel3DUncertainty{}
        };
        first_result.mdpde = rg::GaussianModel3DWithUncertainty{
            rg::GaussianModel3D{ 1.5, 0.8, 1.22 },
            rg::GaussianModel3DUncertainty{}
        };
        auto second_result{ local_result };
        second_result.alpha_r = 0.2;
        second_result.ols = rg::GaussianModel3DWithUncertainty{
            rg::GaussianModel3D{ 1.0, 0.6, 2.11 },
            rg::GaussianModel3DUncertainty{}
        };
        second_result.mdpde = rg::GaussianModel3DWithUncertainty{
            rg::GaussianModel3D{ 1.5, 0.8, 2.22 },
            rg::GaussianModel3DUncertainty{}
        };
        local_editor.SetGaussianResult(
            rg::FittingStage::First,
            first_result);
        local_editor.SetGaussianResult(
            rg::FittingStage::Second,
            second_result);
        local_editor.SetGaussianResult(
            rg::FittingStage::Third,
            local_result);

        analysis.RebuildAtomGroupsFromSelection();
        const auto analysis_view{ model->GetAnalysisView() };
        const auto group_keys{ analysis_view.CollectAtomGroupKeys(
            rg::FittingStage::Third) };
        ASSERT_FALSE(group_keys.empty());
        group_key = group_keys.front();
        const auto & atom_list{ analysis_view.GetAtomObjectList(
            rg::FittingStage::Third, group_key) };
        ASSERT_FALSE(atom_list.empty());
        annotated_serial_id = atom_list.front()->GetSerialID();

        rg::GroupGaussianResult group_result;
        group_result.alpha_g = 0.25;
        group_result.mean = rg::GaussianModel3D{ 1.1, 1.2, 0.33 };
        group_result.mdpde = rg::GaussianModel3D{ 1.2, 1.1, 0.44 };
        group_result.prior = rg::GaussianModel3DWithUncertainty{
            rg::GaussianModel3D{ 1.3, 1.0, 0.55 },
            rg::GaussianModel3DUncertainty{ 0.1, 0.2, 0.03 }
        };
        group_result.member_results.reserve(atom_list.size());
        for (std::size_t i = 0; i < atom_list.size(); i++)
        {
            rg::GroupGaussianMemberResult member_result;
            member_result.posterior = rg::GaussianModel3DWithUncertainty{
                rg::GaussianModel3D{ 0.4 + 0.1 * static_cast<double>(i), 0.9, 0.66 },
                rg::GaussianModel3DUncertainty{ 0.01, 0.02, 0.04 }
            };
            member_result.is_outlier = (i == 0);
            member_result.statistical_distance = 1.5 + static_cast<double>(i);
            group_result.member_results.emplace_back(member_result);
        }
        auto first_group_result{ group_result };
        first_group_result.alpha_g = 0.1;
        first_group_result.mean = rg::GaussianModel3D{ 1.1, 1.2, 1.33 };
        first_group_result.mdpde = rg::GaussianModel3D{ 1.2, 1.1, 1.44 };
        first_group_result.prior = rg::GaussianModel3DWithUncertainty{
            rg::GaussianModel3D{ 1.3, 1.0, 1.55 },
            rg::GaussianModel3DUncertainty{ 0.1, 0.2, 0.13 }
        };
        analysis.ApplyAtomGroupGaussianResult(
            rg::FittingStage::First,
            group_key,
            first_group_result);
        auto second_group_result{ group_result };
        second_group_result.alpha_g = 0.2;
        second_group_result.mean = rg::GaussianModel3D{ 1.1, 1.2, 2.33 };
        second_group_result.mdpde = rg::GaussianModel3D{ 1.2, 1.1, 2.44 };
        second_group_result.prior = rg::GaussianModel3DWithUncertainty{
            rg::GaussianModel3D{ 1.3, 1.0, 2.55 },
            rg::GaussianModel3DUncertainty{ 0.1, 0.2, 0.23 }
        };
        analysis.ApplyAtomGroupGaussianResult(
            rg::FittingStage::Second,
            group_key,
            second_group_result);
        analysis.ApplyAtomGroupGaussianResult(
            rg::FittingStage::Third,
            group_key,
            group_result);

        repository.SaveModel(*model, "model");
    }

    rg::DataRepository repository{ database_path };
    auto loaded_model{ repository.LoadModel("model") };
    ASSERT_NE(loaded_model, nullptr);

    const auto local_view{
        rg::AtomLocalPotentialView::RequireFor(*loaded_model->GetAtomList().at(0))
    };
    EXPECT_DOUBLE_EQ(
        local_view.GetGaussianResult(rg::FittingStage::First)
            .mdpde.GetModel().GetOffset(),
        1.22);
    EXPECT_DOUBLE_EQ(
        local_view.GetGaussianResult(rg::FittingStage::Second)
            .mdpde.GetModel().GetOffset(),
        2.22);
    EXPECT_DOUBLE_EQ(
        local_view.GetGaussianResult(rg::FittingStage::Third)
            .ols.GetModel().GetOffset(),
        0.11);
    EXPECT_DOUBLE_EQ(
        local_view.GetGaussianResult(rg::FittingStage::Third)
            .mdpde.GetModel().GetOffset(),
        0.22);

    const auto loaded_analysis_view{ loaded_model->GetAnalysisView() };
    EXPECT_DOUBLE_EQ(
        loaded_analysis_view.GetAtomGroupMean(
            rg::FittingStage::First, group_key).GetOffset(),
        1.33);
    EXPECT_DOUBLE_EQ(
        loaded_analysis_view.GetAtomGroupMDPDE(
            rg::FittingStage::First, group_key).GetOffset(),
        1.44);
    EXPECT_DOUBLE_EQ(
        loaded_analysis_view.GetAtomGroupPrior(
            rg::FittingStage::First, group_key).GetOffset(),
        1.55);
    EXPECT_DOUBLE_EQ(
        loaded_analysis_view.GetAtomAlphaG(
            rg::FittingStage::First, group_key),
        0.1);
    EXPECT_DOUBLE_EQ(
        loaded_analysis_view.GetAtomGroupMean(
            rg::FittingStage::Second, group_key).GetOffset(),
        2.33);
    EXPECT_DOUBLE_EQ(
        loaded_analysis_view.GetAtomGroupMDPDE(
            rg::FittingStage::Second, group_key).GetOffset(),
        2.44);
    EXPECT_DOUBLE_EQ(
        loaded_analysis_view.GetAtomGroupPrior(
            rg::FittingStage::Second, group_key).GetOffset(),
        2.55);
    EXPECT_DOUBLE_EQ(
        loaded_analysis_view.GetAtomAlphaG(
            rg::FittingStage::Second, group_key),
        0.2);
    EXPECT_DOUBLE_EQ(
        loaded_analysis_view.GetAtomGroupMean(
            rg::FittingStage::Third, group_key).GetOffset(),
        0.33);
    EXPECT_DOUBLE_EQ(
        loaded_analysis_view.GetAtomGroupMDPDE(
            rg::FittingStage::Third, group_key).GetOffset(),
        0.44);
    EXPECT_DOUBLE_EQ(
        loaded_analysis_view.GetAtomGroupPrior(
            rg::FittingStage::Third, group_key).GetOffset(),
        0.55);
    EXPECT_DOUBLE_EQ(
        loaded_analysis_view.GetAtomGroupPriorWithUncertainty(
            rg::FittingStage::Third, group_key)
            .GetStandardDeviationModel()
            .GetOffset(),
        0.03);

    rg::AtomObject * annotated_atom{ nullptr };
    for (const auto & atom : loaded_model->GetAtomList())
    {
        if (atom->GetSerialID() == annotated_serial_id)
        {
            annotated_atom = atom.get();
            break;
        }
    }
    ASSERT_NE(annotated_atom, nullptr);
    const auto & gaussian_result{
        rg::AtomLocalPotentialView::RequireFor(*annotated_atom)
            .GetGaussianResult(rg::FittingStage::Third)
    };
    ASSERT_TRUE(gaussian_result.posterior.has_value());
    EXPECT_DOUBLE_EQ(gaussian_result.posterior->GetModel().GetOffset(), 0.66);
    EXPECT_DOUBLE_EQ(
        gaussian_result.posterior->GetStandardDeviationModel().GetOffset(),
        0.04);
}

TEST(DataObjectPersistenceTest, LegacyV2SamplingBlobLoadsAsSelectedAndMigratesVersion)
{
    const command_test::ScopedTempDir temp_dir{ "data_schema_legacy_sampling_blob" };
    const auto database_path{ temp_dir.path() / "legacy_sampling.sqlite" };

    {
        rg::DataRepository repository{ database_path };
        auto model{ data_test::MakeModelWithBond() };
        model->SetKeyTag("model");
        repository.SaveModel(*model, "model");
    }

    data_test::ConvertLocalGaussianColumnsToLegacyFinal(database_path);
    data_test::ConvertAtomGroupGaussianColumnsToLegacyFinal(database_path);
    data_test::ConvertSamplingEntryColumnsToLegacyRawOnly(database_path);

    {
        rg::SQLiteWrapper database{ database_path };
        database.Prepare(
            "INSERT OR REPLACE INTO model_atom_local_potential ("
            "key_tag, serial_id, sampling_size, distance_and_map_value_list, "
            "amplitude_estimate_ols, width_estimate_ols, intercept_estimate_ols, "
            "amplitude_estimate_mdpde, width_estimate_mdpde, intercept_estimate_mdpde, alpha_r"
            ") VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?);");
        rg::SQLiteWrapper::StatementGuard guard(database);
        const std::vector<float> legacy_samples{
            0.1f, 6.0f,
            0.2f, 4.0f
        };
        database.Bind<std::string>(1, "model");
        database.Bind<int>(2, 1);
        database.Bind<int>(3, 2);
        database.Bind<std::vector<float>>(4, legacy_samples);
        database.Bind<double>(5, 0.0);
        database.Bind<double>(6, 0.0);
        database.Bind<double>(7, 0.0);
        database.Bind<double>(8, 0.0);
        database.Bind<double>(9, 0.0);
        database.Bind<double>(10, 0.0);
        database.Bind<double>(11, 0.0);
        database.StepOnce();
    }
    data_test::SetUserVersion(database_path, 5);

    rg::DataRepository repository{ database_path };
    auto loaded_model{ repository.LoadModel("model") };
    ASSERT_NE(loaded_model, nullptr);
    const auto raw_entries{
        rg::AtomLocalPotentialView::RequireFor(
            *loaded_model->GetAtomList().at(0)).GetRawSamplingEntries(false)
    };

    ASSERT_EQ(raw_entries.size(), 2u);
    EXPECT_TRUE(raw_entries.at(0).point.is_selected);
    EXPECT_TRUE(raw_entries.at(1).point.is_selected);
    const auto peeling_entries{
        rg::AtomLocalPotentialView::RequireFor(*loaded_model->GetAtomList().at(0))
            .GetPeelingSamplingEntries(false)
    };

    EXPECT_TRUE(peeling_entries.empty());
    EXPECT_EQ(data_test::GetUserVersion(database_path), 12);
}

TEST(DataObjectPersistenceTest, DatabaseRoundTripPreservesChainMetadataAndSymmetryFiltering)
{
    const command_test::ScopedTempDir temp_dir{ "data_schema_roundtrip" };
    const auto database_path{ temp_dir.path() / "roundtrip.sqlite" };
    const auto model_path{ command_test::TestDataPath("test_model_keyvalue_entity.cif") };

    auto original_model{ data_test::LoadFixtureModel(model_path) };
    const auto original_chain_map{ original_model->GetChainIDListMap() };
    original_model->SelectAllAtoms();
    original_model->ApplySymmetrySelection(false);
    const auto original_selected_count{ original_model->GetSelectedAtomCount() };

    rg::DataRepository repository{ database_path };
    auto stored_model{ rg::ReadModel(model_path) };
    stored_model->SetKeyTag("model");
    repository.SaveModel(*stored_model, "model");

    auto loaded_model{ repository.LoadModel("model") };
    EXPECT_EQ(loaded_model->GetChainIDListMap(), original_chain_map);
    EXPECT_GT(data_test::CountRows(database_path, "model_chain_map", "model"), 0);

    loaded_model->SelectAllAtoms();
    loaded_model->ApplySymmetrySelection(false);
    EXPECT_EQ(loaded_model->GetSelectedAtomCount(), original_selected_count);
}

TEST(DataObjectPersistenceTest, DistinctUnsanitizedKeysDoNotCollideInV2Schema)
{
    const command_test::ScopedTempDir temp_dir{ "data_schema_key_collision" };
    const auto database_path{ temp_dir.path() / "collision.sqlite" };
    const auto model_path{ command_test::TestDataPath("test_model.cif") };

    rg::DataRepository repository{ database_path };
    auto model_a{ rg::ReadModel(model_path) };
    auto model_b{ rg::ReadModel(model_path) };
    model_a->SetKeyTag("mem_a");
    model_b->SetKeyTag("mem_b");
    model_a->SetPdbID("MODEL_A");
    model_b->SetPdbID("MODEL_B");

    repository.SaveModel(*model_a, "a-b");
    repository.SaveModel(*model_b, "a_b");

    auto loaded_model_a{ repository.LoadModel("a-b") };
    auto loaded_model_b{ repository.LoadModel("a_b") };
    EXPECT_EQ(loaded_model_a->GetPdbID(), "MODEL_A");
    EXPECT_EQ(loaded_model_b->GetPdbID(), "MODEL_B");
    EXPECT_EQ(data_test::CountRows(database_path, "model_object"), 2);
}

TEST(DataObjectPersistenceTest, SaveRenamedKeyDoesNotRenameInMemoryObject)
{
    const command_test::ScopedTempDir temp_dir{ "data_schema_rename_semantics" };
    const auto database_path{ temp_dir.path() / "rename.sqlite" };
    const auto model_path{ command_test::TestDataPath("test_model.cif") };

    rg::DataRepository repository{ database_path };
    auto model{ rg::ReadModel(model_path) };
    model->SetKeyTag("model");

    repository.SaveModel(*model, "saved_model");

    EXPECT_EQ(model->GetKeyTag(), "model");
    auto loaded_model{ repository.LoadModel("saved_model") };
    ASSERT_NE(loaded_model, nullptr);
    EXPECT_EQ(loaded_model->GetNumberOfAtom(), model->GetNumberOfAtom());
}

TEST(DataObjectPersistenceTest, LoadModelRestoresSelectionFromPersistedLocalEntriesInBulk)
{
    const command_test::ScopedTempDir temp_dir{ "data_schema_bulk_selection_restore" };
    const auto database_path{ temp_dir.path() / "bulk_selection.sqlite" };

    rg::DataRepository repository{ database_path };
    auto model{ data_test::MakeModelWithBond() };
    model->SetKeyTag("model");
    model->SetPdbID("BULK_SELECTION");
    model->SelectAllAtoms(false);
    model->SelectAllBonds(false);

    auto & atoms{ model->GetAtomList() };
    ASSERT_EQ(atoms.size(), 2);

    rg::ModelAnalysisData::Of(*model).EnsureAtomLocalEntry(*atoms.at(0));

    repository.SaveModel(*model, "model");

    auto loaded_model{ repository.LoadModel("model") };
    ASSERT_NE(loaded_model, nullptr);
    ASSERT_EQ(loaded_model->GetSelectedAtomCount(), 1);
    ASSERT_EQ(loaded_model->GetSelectedBondCount(), 0);
    EXPECT_EQ(loaded_model->GetSelectedAtoms().front()->GetSerialID(), atoms.at(0)->GetSerialID());
    EXPECT_NE(
        rg::ModelAnalysisData::Of(*loaded_model).FindAtomLocalEntry(
            *loaded_model->GetSelectedAtoms().front()),
        nullptr);
}

TEST(DataObjectPersistenceTest, LoadModelThrowsWhenDatabaseKeyIsMissing)
{
    const command_test::ScopedTempDir temp_dir{ "data_schema_missing_database_key" };
    const auto database_path{ temp_dir.path() / "missing.sqlite" };

    rg::DataRepository repository{ database_path };
    EXPECT_THROW((void)repository.LoadModel("missing"), std::runtime_error);
}
