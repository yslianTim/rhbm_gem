#include <gtest/gtest.h>

#include <stdexcept>
#include <string>

#include <rhbm_gem/data/io/DataRepository.hpp>
#include <rhbm_gem/data/io/ModelMapFileIO.hpp>
#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <rhbm_gem/data/object/ModelAnalysisView.hpp>
#include "support/CommandTestHelpers.hpp"
#include "support/DataObjectTestSupport.hpp"

namespace rg = rhbm_gem;

TEST(DataObjectPersistenceTest, SaveAndLoadModelWithoutRegistryState)
{
    const command_test::ScopedTempDir temp_dir{ "data_repository_schema_roundtrip" };
    const auto database_path{ temp_dir.path() / "repository.sqlite" };
    const auto model_path{ command_test::TestDataPath("test_model.cif") };

    rg::DataRepository repository{ database_path };
    auto model{ rg::ReadModel(model_path) };
    model->SetKeyTag("memory_model");
    model->SetPdbID("MODEL_REPOSITORY");
    repository.SaveModel(*model, "stored_model");

    auto loaded_model{ repository.LoadModel("stored_model") };
    ASSERT_NE(loaded_model, nullptr);
    EXPECT_EQ(loaded_model->GetKeyTag(), "stored_model");
    EXPECT_EQ(loaded_model->GetPdbID(), "MODEL_REPOSITORY");
    EXPECT_EQ(loaded_model->GetNumberOfAtom(), model->GetNumberOfAtom());
    EXPECT_EQ(loaded_model->GetNumberOfBond(), model->GetNumberOfBond());
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
    EXPECT_DOUBLE_EQ(loaded_model->GetAtomList().at(0)->GetStandardQScore(), 0.5);
    EXPECT_DOUBLE_EQ(loaded_model->GetAtomList().at(1)->GetStandardQScore(), 0.75);
}

TEST(DataObjectPersistenceTest, RawAndPeelingSamplingEntriesRoundTripPreservesSelection)
{
    const command_test::ScopedTempDir temp_dir{ "data_schema_raw_peeling_sampling_roundtrip" };
    const auto database_path{ temp_dir.path() / "raw_peeling_sampling.sqlite" };

    rg::DataRepository repository{ database_path };
    auto model{ data_test::MakeModelWithBond() };
    auto * atom{ model->GetAtomList().at(0).get() };
    auto editor{ model->EditAnalysis() };
    editor.SetAtomLocalRawSamplingEntries(*atom, {
        LocalPotentialSample{
            6.0f, SamplingPoint{ 0.1f, { 0.0f, 0.0f, 0.0f }, true } },
        LocalPotentialSample{
            4.0f, SamplingPoint{ 0.2f, { 0.0f, 0.0f, 0.0f }, false } }
    });
    editor.SetAtomLocalPeelingSamplingEntries(*atom, {
        LocalPotentialSample{
            3.0f, SamplingPoint{ 0.1f, { 0.0f, 0.0f, 0.0f }, true } },
        LocalPotentialSample{
            5.0f, SamplingPoint{ 0.2f, { 0.0f, 0.0f, 0.0f }, false } }
    });
    editor.SetAtomLocalNeighborCountForPeeling(*atom, 7);
    repository.SaveModel(*model, "model");

    auto loaded_model{ repository.LoadModel("model") };
    ASSERT_NE(loaded_model, nullptr);
    const auto view{
        rg::AtomLocalPotentialView::For(*loaded_model->GetAtomList().at(0)) };
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
    EXPECT_EQ(view.GetNeighborCountForPeeling(), 7);
}

TEST(DataObjectPersistenceTest, GaussianOffsetRoundTripPreservesAnalysisResults)
{
    const command_test::ScopedTempDir temp_dir{ "data_schema_gaussian_roundtrip" };
    const auto database_path{ temp_dir.path() / "gaussian.sqlite" };

    rg::DataRepository repository{ database_path };
    auto model{ data_test::MakeModelWithBond() };
    model->SelectAllAtoms();
    auto editor{ model->EditAnalysis() };
    editor.InitializeFromSelection();
    auto * atom{ model->GetAtomList().at(0).get() };

    rg::LocalGaussianResult local_result;
    local_result.alpha_r = 0.5;
    local_result.ols = rg::GaussianModel3DWithUncertainty{
        rg::GaussianModel3D{ 1.0, 0.6, 0.11 },
        rg::GaussianModel3DUncertainty{} };
    local_result.mdpde = rg::GaussianModel3DWithUncertainty{
        rg::GaussianModel3D{ 1.5, 0.7, -0.22 },
        rg::GaussianModel3DUncertainty{} };
    local_result.posterior = rg::GaussianModel3DWithUncertainty{
        rg::GaussianModel3D{ 1.7, 0.8, -0.33 },
        rg::GaussianModel3DUncertainty{ 0.1, 0.2, 0.3 } };
    local_result.is_outlier = true;
    local_result.statistical_distance = 2.5;
    const auto view{ model->GetAnalysisView() };
    const auto group_key{
        view.CollectAtomGroupKeys(rg::FittingStage::Third).front() };
    rg::GroupGaussianResult group_result;
    group_result.alpha_g = 0.25;
    group_result.mean = rg::GaussianModel3D{ 2.0, 0.9, 0.12 };
    group_result.mdpde = rg::GaussianModel3D{ 2.1, 1.0, -0.23 };
    group_result.prior = rg::GaussianModel3DWithUncertainty{
        rg::GaussianModel3D{ 2.2, 1.1, -0.34 },
        rg::GaussianModel3DUncertainty{ 0.4, 0.5, 0.6 } };
    group_result.member_results.resize(
        view.GetAtomObjectList(rg::FittingStage::Third, group_key).size());
    editor.ApplyAtomGroupGaussianResult(
        rg::FittingStage::Third, group_key, group_result);
    editor.SetAtomLocalGaussianResult(
        rg::FittingStage::Third, *atom, local_result);

    repository.SaveModel(*model, "model");
    auto loaded_model{ repository.LoadModel("model") };
    ASSERT_NE(loaded_model, nullptr);

    const auto loaded_local{
        rg::AtomLocalPotentialView::For(*loaded_model->FindAtomPtr(atom->GetSerialID())) };
    const auto loaded_local_result{
        loaded_local.GetGaussianResult(rg::FittingStage::Third) };
    EXPECT_DOUBLE_EQ(loaded_local_result.alpha_r, 0.5);
    EXPECT_DOUBLE_EQ(loaded_local_result.ols.GetModel().GetOffset(), 0.11);
    EXPECT_DOUBLE_EQ(loaded_local_result.mdpde.GetModel().GetOffset(), -0.22);
    ASSERT_TRUE(loaded_local_result.posterior.has_value());
    EXPECT_DOUBLE_EQ(
        loaded_local_result.posterior->GetModel().GetOffset(), -0.33);
    EXPECT_TRUE(loaded_local_result.is_outlier);
    EXPECT_DOUBLE_EQ(loaded_local_result.statistical_distance, 2.5);

    const auto loaded_view{ loaded_model->GetAnalysisView() };
    EXPECT_DOUBLE_EQ(
        loaded_view.GetAtomGroupMean(
            rg::FittingStage::Third, group_key).GetOffset(),
        0.12);
    EXPECT_DOUBLE_EQ(
        loaded_view.GetAtomGroupPriorWithUncertainty(
            rg::FittingStage::Third, group_key)
            .GetStandardDeviationModel().GetOffset(),
        0.6);
}

TEST(DataObjectPersistenceTest, DatabaseRoundTripPreservesChainMetadataAndSymmetryFiltering)
{
    const command_test::ScopedTempDir temp_dir{ "data_schema_chain_roundtrip" };
    const auto database_path{ temp_dir.path() / "roundtrip.sqlite" };
    const auto model_path{ command_test::TestDataPath("test_model_keyvalue_entity.cif") };

    auto original_model{ data_test::LoadFixtureModel(model_path) };
    const auto original_chain_map{ original_model->GetChainIDListMap() };
    original_model->SelectAllAtoms();
    original_model->ApplySymmetrySelection(false);
    const auto original_selected_count{ original_model->GetSelectedAtomCount() };

    rg::DataRepository repository{ database_path };
    repository.SaveModel(*original_model, "model");
    auto loaded_model{ repository.LoadModel("model") };

    EXPECT_EQ(loaded_model->GetChainIDListMap(), original_chain_map);
    EXPECT_GT(data_test::CountRows(database_path, "model_chain_map", "model"), 0);
    loaded_model->SelectAllAtoms();
    loaded_model->ApplySymmetrySelection(false);
    EXPECT_EQ(loaded_model->GetSelectedAtomCount(), original_selected_count);
}

TEST(DataObjectPersistenceTest, DistinctUnsanitizedKeysDoNotCollideInV13Schema)
{
    const command_test::ScopedTempDir temp_dir{ "data_schema_key_collision" };
    const auto database_path{ temp_dir.path() / "collision.sqlite" };
    const auto model_path{ command_test::TestDataPath("test_model.cif") };

    rg::DataRepository repository{ database_path };
    auto model_a{ rg::ReadModel(model_path) };
    auto model_b{ rg::ReadModel(model_path) };
    model_a->SetPdbID("MODEL_A");
    model_b->SetPdbID("MODEL_B");

    repository.SaveModel(*model_a, "a-b");
    repository.SaveModel(*model_b, "a_b");

    EXPECT_EQ(repository.LoadModel("a-b")->GetPdbID(), "MODEL_A");
    EXPECT_EQ(repository.LoadModel("a_b")->GetPdbID(), "MODEL_B");
    EXPECT_EQ(data_test::CountRows(database_path, "model_object"), 2);
}

TEST(DataObjectPersistenceTest, SaveRenamedKeyDoesNotRenameInMemoryObject)
{
    const command_test::ScopedTempDir temp_dir{ "data_schema_rename_semantics" };
    const auto database_path{ temp_dir.path() / "rename.sqlite" };
    auto model{ data_test::MakeModelWithBond() };
    model->SetKeyTag("memory_model");

    rg::DataRepository repository{ database_path };
    repository.SaveModel(*model, "saved_model");

    EXPECT_EQ(model->GetKeyTag(), "memory_model");
    EXPECT_EQ(repository.LoadModel("saved_model")->GetKeyTag(), "saved_model");
}

TEST(DataObjectPersistenceTest, LoadModelRestoresSelectionFromPersistedColumns)
{
    const command_test::ScopedTempDir temp_dir{ "data_schema_selection_restore" };
    const auto database_path{ temp_dir.path() / "selection.sqlite" };

    auto model{ data_test::MakeModelWithBond() };
    auto & atoms{ model->GetAtomList() };
    ASSERT_EQ(atoms.size(), 2u);
    model->SelectAllAtoms(false);
    model->SetAtomSelected(atoms.at(0)->GetSerialID(), true);
    model->SelectAllBonds(false);
    model->SetBondSelected(
        atoms.at(0)->GetSerialID(), atoms.at(1)->GetSerialID(), true);
    model->EditAnalysis().SetAtomLocalAlphaR(
        rg::FittingStage::Third, *atoms.at(1), 0.9);

    rg::DataRepository repository{ database_path };
    repository.SaveModel(*model, "model");
    auto loaded_model{ repository.LoadModel("model") };

    ASSERT_NE(loaded_model, nullptr);
    ASSERT_EQ(loaded_model->GetSelectedAtomCount(), 1u);
    EXPECT_EQ(
        loaded_model->GetSelectedAtoms().front()->GetSerialID(),
        atoms.at(0)->GetSerialID());
    EXPECT_EQ(loaded_model->GetSelectedBondCount(), 1u);
    EXPECT_TRUE(rg::AtomLocalPotentialView::For(
        *loaded_model->FindAtomPtr(atoms.at(1)->GetSerialID())).IsAvailable());
}

TEST(DataObjectPersistenceTest, LoadModelThrowsWhenDatabaseKeyIsMissing)
{
    const command_test::ScopedTempDir temp_dir{ "data_schema_missing_database_key" };
    const auto database_path{ temp_dir.path() / "missing.sqlite" };

    rg::DataRepository repository{ database_path };
    EXPECT_THROW((void)repository.LoadModel("missing"), std::runtime_error);
}
