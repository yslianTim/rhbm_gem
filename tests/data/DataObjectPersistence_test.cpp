#include <gtest/gtest.h>

#include <stdexcept>
#include <string>
#include <fstream>
#include <boost/json.hpp>
#include <rhbm_gem/data/io/JointAnalysisFileIO.hpp>
#include "data/io/detail/JointResultJson.hpp"

#include <rhbm_gem/data/io/DataRepository.hpp>
#include <rhbm_gem/data/io/ModelMapFileIO.hpp>
#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>
#include <rhbm_gem/data/object/ChemicalComponentEntry.hpp>
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
            6.0, SamplingPoint{ 0.1, { 0.0, 0.0, 0.0 }, true } },
        LocalPotentialSample{
            4.0, SamplingPoint{ 0.2, { 0.0, 0.0, 0.0 }, false } }
    });
    editor.SetAtomLocalPeelingSamplingEntries(*atom, {
        LocalPotentialSample{
            3.0, SamplingPoint{ 0.1, { 0.0, 0.0, 0.0 }, true } },
        LocalPotentialSample{
            5.0, SamplingPoint{ 0.2, { 0.0, 0.0, 0.0 }, false } }
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
    EXPECT_DOUBLE_EQ(raw_entries.at(0).response, 6.0);
    EXPECT_DOUBLE_EQ(raw_entries.at(1).response, 4.0);
    EXPECT_TRUE(raw_entries.at(0).point.is_selected);
    EXPECT_FALSE(raw_entries.at(1).point.is_selected);
    ASSERT_EQ(peeling_entries.size(), 2u);
    EXPECT_DOUBLE_EQ(peeling_entries.at(0).response, 3.0);
    EXPECT_DOUBLE_EQ(peeling_entries.at(1).response, 5.0);
    EXPECT_TRUE(peeling_entries.at(0).point.is_selected);
    EXPECT_FALSE(peeling_entries.at(1).point.is_selected);
    ASSERT_EQ(selected_peeling_entries.size(), 1u);
    EXPECT_EQ(view.GetNeighborCountForPeeling(), 7);
}

TEST(DataObjectPersistenceTest, DoublePrecisionDomainValuesRoundTripWithoutNarrowing)
{
    const command_test::ScopedTempDir temp_dir{ "data_schema_double_precision_roundtrip" };
    const auto database_path{ temp_dir.path() / "double_precision.sqlite" };
    constexpr ComponentKey component_key{ 30 };

    rg::ModelObjectParts parts;
    parts.component_key_system->RegisterComponent("DBL", component_key);
    auto component{ std::make_unique<rg::ChemicalComponentEntry>() };
    component->SetComponentId("DBL");
    component->SetComponentName("DOUBLE PRECISION COMPONENT");
    component->SetComponentType("non-polymer");
    component->SetComponentFormula("C1");
    component->SetComponentMolecularWeight(123.456789012345);
    component->SetStandardMonomerFlag(false);
    parts.chemical_component_entry_map.emplace(
        component_key, std::move(component));

    auto atom{ std::make_unique<rg::AtomObject>() };
    atom->SetSerialID(1);
    atom->SetComponentID("DBL");
    atom->SetComponentKey(component_key);
    atom->SetAtomID("C1");
    atom->SetPosition(
        1.0000000000000002,
        16777217.125,
        -0.123456789012345);
    atom->SetOccupancy(0.123456789012345);
    atom->SetTemperature(12.3456789012345);
    parts.atom_list.emplace_back(std::move(atom));

    auto model{ std::make_unique<rg::ModelObject>(
        rg::AssembleModelObject(std::move(parts))) };
    auto * stored_atom{ model->GetAtomList().front().get() };
    const LocalPotentialSampleList raw_samples{
        { 16777217.125,
            SamplingPoint{ 1.0000000000000002, { 0.0, 0.0, 0.0 }, true } },
        { -0.123456789012345,
            SamplingPoint{ 0.3333333333333333, { 0.0, 0.0, 0.0 }, false } }
    };
    const LocalPotentialSampleList peeling_samples{
        { 3.141592653589793,
            SamplingPoint{ 0.987654321098765, { 0.0, 0.0, 0.0 }, true } }
    };
    auto editor{ model->EditAnalysis() };
    editor.SetAtomLocalRawSamplingEntries(*stored_atom, raw_samples);
    editor.SetAtomLocalPeelingSamplingEntries(*stored_atom, peeling_samples);

    rg::DataRepository repository{ database_path };
    repository.SaveModel(*model, "model");
    const auto loaded_model{ repository.LoadModel("model") };
    ASSERT_NE(loaded_model, nullptr);
    ASSERT_EQ(loaded_model->GetAtomList().size(), 1u);

    const auto & loaded_atom{ *loaded_model->GetAtomList().front() };
    EXPECT_DOUBLE_EQ(loaded_atom.GetPosition().at(0), 1.0000000000000002);
    EXPECT_DOUBLE_EQ(loaded_atom.GetPosition().at(1), 16777217.125);
    EXPECT_DOUBLE_EQ(loaded_atom.GetPosition().at(2), -0.123456789012345);
    EXPECT_DOUBLE_EQ(loaded_atom.GetOccupancy(), 0.123456789012345);
    EXPECT_DOUBLE_EQ(loaded_atom.GetTemperature(), 12.3456789012345);

    const auto component_keys{ loaded_model->GetComponentKeyList() };
    ASSERT_EQ(component_keys.size(), 1u);
    const auto * loaded_component{
        loaded_model->FindChemicalComponentEntry(component_keys.front()) };
    ASSERT_NE(loaded_component, nullptr);
    EXPECT_DOUBLE_EQ(
        loaded_component->GetComponentMolecularWeight(),
        123.456789012345);

    const auto view{ rg::AtomLocalPotentialView::For(loaded_atom) };
    const auto loaded_raw{ view.GetRawSamplingEntries(false) };
    const auto loaded_peeling{ view.GetPeelingSamplingEntries(false) };
    ASSERT_EQ(loaded_raw.size(), raw_samples.size());
    ASSERT_EQ(loaded_peeling.size(), peeling_samples.size());
    EXPECT_DOUBLE_EQ(loaded_raw.at(0).point.distance, 1.0000000000000002);
    EXPECT_DOUBLE_EQ(loaded_raw.at(0).response, 16777217.125);
    EXPECT_TRUE(loaded_raw.at(0).point.is_selected);
    EXPECT_DOUBLE_EQ(loaded_raw.at(1).point.distance, 0.3333333333333333);
    EXPECT_DOUBLE_EQ(loaded_raw.at(1).response, -0.123456789012345);
    EXPECT_FALSE(loaded_raw.at(1).point.is_selected);
    EXPECT_DOUBLE_EQ(
        loaded_peeling.at(0).point.distance, 0.987654321098765);
    EXPECT_DOUBLE_EQ(loaded_peeling.at(0).response, 3.141592653589793);

    rg::SQLiteWrapper database{ database_path };
    database.Prepare(
        "SELECT length(raw_distance_and_map_value_list), "
        "length(peeling_distance_and_map_value_list) "
        "FROM model_atom_local_potential "
        "WHERE key_tag = ? AND serial_id = ?;");
    rg::SQLiteWrapper::StatementGuard guard(database);
    database.Bind<std::string>(1, "model");
    database.Bind<int>(2, 1);
    ASSERT_EQ(database.StepNext(), rg::SQLiteWrapper::StepRow());
    EXPECT_EQ(
        database.GetColumn<int>(0),
        static_cast<int>(raw_samples.size() * 3 * sizeof(double)));
    EXPECT_EQ(
        database.GetColumn<int>(1),
        static_cast<int>(peeling_samples.size() * 3 * sizeof(double)));
}

TEST(DataObjectPersistenceTest, InvalidV14SamplingBlobLengthIsRejected)
{
    const command_test::ScopedTempDir temp_dir{ "data_schema_invalid_double_blob" };
    const auto database_path{ temp_dir.path() / "invalid_double_blob.sqlite" };
    auto model{ data_test::MakeModelWithBond() };
    auto * atom{ model->GetAtomList().front().get() };
    model->EditAnalysis().SetAtomLocalRawSamplingEntries(
        *atom,
        { LocalPotentialSample{ 2.0, SamplingPoint{ 0.25 } } });

    {
        rg::DataRepository repository{ database_path };
        repository.SaveModel(*model, "model");
    }
    data_test::ExecuteSql(
        database_path,
        "UPDATE model_atom_local_potential "
        "SET raw_distance_and_map_value_list = X'000102' "
        "WHERE key_tag = 'model' AND serial_id = 1;");

    rg::DataRepository repository{ database_path };
    EXPECT_THROW((void)repository.LoadModel("model"), std::runtime_error);
    EXPECT_EQ(data_test::GetUserVersion(database_path), 18);
    EXPECT_EQ(
        data_test::CountRows(
            database_path, "model_atom_local_potential", "model"),
        1);
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
    rg::GroupGaussianMemberResult member_result;
    member_result.posterior = rg::GaussianModel3DWithUncertainty{
        rg::GaussianModel3D{ 1.7, 0.8, -0.33 },
        rg::GaussianModel3DUncertainty{ 0.1, 0.2, 0.3 } };
    member_result.is_outlier = true;
    member_result.statistical_distance = 2.5;
    const auto view{ model->GetAnalysisView() };
    const auto group_key{ view.CollectAtomGroupKeys().front() };
    atom = view.GetAtomObjectList(group_key).front();
    rg::GroupGaussianResult group_result;
    group_result.alpha_g = 0.25;
    group_result.mean = rg::GaussianModel3D{ 2.0, 0.9, 0.12 };
    group_result.mdpde = rg::GaussianModel3D{ 2.1, 1.0, -0.23 };
    group_result.prior = rg::GaussianModel3DWithUncertainty{
        rg::GaussianModel3D{ 2.2, 1.1, -0.34 },
        rg::GaussianModel3DUncertainty{ 0.4, 0.5, 0.6 } };
    group_result.member_results.resize(
        view.GetAtomObjectList(group_key).size());
    group_result.member_results.front() = member_result;
    editor.SetAtomLocalGaussianResult(
        rg::FittingStage::Second, *atom, local_result);
    editor.ApplyAtomGroupGaussianResult(group_key, group_result);

    rg::LocalGaussianResult first_result{ local_result };
    first_result.alpha_r = 0.1;
    first_result.mdpde = rg::GaussianModel3DWithUncertainty{
        rg::GaussianModel3D{ 4.0, 0.9, 0.4 }, rg::GaussianModel3DUncertainty{} };
    editor.SetAtomLocalGaussianResult(rg::FittingStage::First, *atom, first_result);
    editor.SetAtomLocalNeighborCountForPeeling(*atom, 7);
    repository.SaveModel(*model, "model");
    auto loaded_model{ repository.LoadModel("model") };
    ASSERT_NE(loaded_model, nullptr);

    const auto loaded_local{
        rg::AtomLocalPotentialView::For(*loaded_model->FindAtomPtr(atom->GetSerialID())) };
    EXPECT_DOUBLE_EQ(loaded_local.GetAlphaR(rg::FittingStage::First), 0.1);
    EXPECT_DOUBLE_EQ(loaded_local.GetAlphaR(rg::FittingStage::Second), 0.5);
    EXPECT_EQ(loaded_local.GetNeighborCountForPeeling(), 7);
    EXPECT_DOUBLE_EQ(loaded_local.GetEstimateMDPDE(rg::FittingStage::First).GetOffset(), 0.4);
    const auto loaded_local_result{
        loaded_local.GetGaussianResult(rg::FittingStage::Second) };
    EXPECT_DOUBLE_EQ(loaded_local_result.alpha_r, 0.5);
    EXPECT_DOUBLE_EQ(loaded_local_result.ols.GetModel().GetOffset(), 0.11);
    EXPECT_DOUBLE_EQ(loaded_local_result.mdpde.GetModel().GetOffset(), -0.22);
    const auto & loaded_member_result{ loaded_local.GetGroupMemberResult() };
    ASSERT_TRUE(loaded_member_result.has_value());
    EXPECT_DOUBLE_EQ(
        loaded_member_result->posterior.GetModel().GetOffset(), -0.33);
    EXPECT_TRUE(loaded_member_result->is_outlier);
    EXPECT_DOUBLE_EQ(loaded_member_result->statistical_distance, 2.5);

    const auto loaded_view{ loaded_model->GetAnalysisView() };
    EXPECT_EQ(loaded_view.CollectAtomGroupKeys().size(), view.CollectAtomGroupKeys().size());
    EXPECT_EQ(loaded_view.GetAtomObjectList(group_key).size(), view.GetAtomObjectList(group_key).size());
    EXPECT_DOUBLE_EQ(loaded_view.GetAtomAlphaG(group_key), 0.25);
    EXPECT_DOUBLE_EQ(loaded_view.GetAtomGroupMDPDE(group_key).GetOffset(), -0.23);
    EXPECT_DOUBLE_EQ(loaded_view.GetAtomGroupPrior(group_key).GetOffset(), -0.34);
    for (int parameter = 0; parameter < 3; ++parameter)
    {
        EXPECT_DOUBLE_EQ(loaded_view.GetAtomGroupMean(group_key).GetModelParameter(parameter),
            group_result.mean.GetModelParameter(parameter));
        EXPECT_DOUBLE_EQ(loaded_view.GetAtomGroupMDPDE(group_key).GetModelParameter(parameter),
            group_result.mdpde.GetModelParameter(parameter));
        EXPECT_DOUBLE_EQ(loaded_view.GetAtomGroupPriorWithUncertainty(group_key)
            .GetModelStandardDeviation(parameter), group_result.prior.GetModelStandardDeviation(parameter));
        EXPECT_DOUBLE_EQ(loaded_member_result->posterior.GetModelStandardDeviation(parameter),
            member_result.posterior.GetModelStandardDeviation(parameter));
    }
    EXPECT_DOUBLE_EQ(
        loaded_view.GetAtomGroupMean(group_key).GetOffset(),
        0.12);
    EXPECT_DOUBLE_EQ(
        loaded_view.GetAtomGroupPriorWithUncertainty(group_key)
            .GetStandardDeviationModel().GetOffset(),
        0.6);

    editor.InitializeLocalFittingSeedModels();
    repository.SaveModel(*model, "model");
    EXPECT_EQ(data_test::CountRows(database_path, "model_atom_posterior", "model"), 0);
    auto reset_model{ repository.LoadModel("model") };
    EXPECT_FALSE(rg::AtomLocalPotentialView::For(*reset_model->FindAtomPtr(atom->GetSerialID()))
        .GetGroupMemberResult().has_value());
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

TEST(DataObjectPersistenceTest, DistinctUnsanitizedKeysDoNotCollideInV14Schema)
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
        rg::FittingStage::Second, *atoms.at(1), 0.9);

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

TEST(DataObjectPersistenceTest, GroupWriteFailureRollsBackAnalysisReplacement)
{
    const command_test::ScopedTempDir temp_dir{ "data_group_rollback" };
    const auto database_path{ temp_dir.path() / "group.sqlite" };
    rg::DataRepository repository{ database_path };
    auto model{ data_test::MakeModelWithBond() };
    model->SelectAllAtoms();
    auto editor{ model->EditAnalysis() };
    editor.InitializeFromSelection();
    const auto view{ model->GetAnalysisView() };
    const auto group_key{ view.CollectAtomGroupKeys().front() };
    rg::GroupGaussianResult result;
    result.alpha_g = 0.2;
    result.prior = rg::GaussianModel3DWithUncertainty{
        rg::GaussianModel3D{ 2.0, 0.8, 0.3 }, rg::GaussianModel3DUncertainty{} };
    result.member_results.assign(view.GetAtomObjectList(group_key).size(),
        rg::GroupGaussianMemberResult{ result.prior, true, 2.5 });
    editor.ApplyAtomGroupGaussianResult(group_key, result);
    repository.SaveModel(*model, "model");

    result.alpha_g = 0.9;
    result.member_results.front().statistical_distance = 8.0;
    editor.ApplyAtomGroupGaussianResult(group_key, result);
    data_test::ExecuteSql(database_path,
        "CREATE TRIGGER reject_group_insert BEFORE INSERT ON model_atom_group_potential "
        "BEGIN SELECT RAISE(ABORT, 'test group write failure'); END;");
    EXPECT_THROW(repository.SaveModel(*model, "model"), std::runtime_error);
    data_test::ExecuteSql(database_path, "DROP TRIGGER reject_group_insert;");

    auto loaded{ repository.LoadModel("model") };
    const auto loaded_view{ loaded->GetAnalysisView() };
    EXPECT_DOUBLE_EQ(loaded_view.GetAtomAlphaG(group_key), 0.2);
    const auto & member{ rg::AtomLocalPotentialView::For(
        *loaded_view.GetAtomObjectList(group_key).front()).GetGroupMemberResult() };
    ASSERT_TRUE(member.has_value());
    EXPECT_DOUBLE_EQ(member->statistical_distance, 2.5);
}

namespace {
rg::JointAnalysisResult SavedJointExample()
{
    rg::JointAnalysisResult result;
    result.metadata.model_sha256=std::string(64,'a'); result.metadata.map_sha256=std::string(64,'b');
    result.metadata.software=rg::JointSoftwareProvenance{"saved-version",std::string(64,'c'),std::string(64,'d'),std::string(64,'e')};
    result.metadata.map_normalization=rg::JointMapNormalization{true,true,2.5};
    result.atom_ids={"1","2"}; result.row_ids={"voxel-1","constant"}; result.available_row_mask={true,true};
    result.initialization={true,"",{.5,.6},{}};
    result.selection_domain=rg::JointSelectionDomain{}; result.selection_domain->target_indices={0};
    rg::JointAnalysisComponent component;
    component.id="component-0"; component.atoms={0,1}; component.rows={0};
    component.stop_reason="budget"; component.search_completed=false;
    component.state=rg::JointState{{2.123456789012345,-.3,4,.2},{.5,.6},{-.7,-.5},{1e-14,-2e-14},.012345678901234567};
    component.runtime_convergence=rg::JointCheckStatus::Failed;
    component.evidence.push_back({"local-correction",rg::JointCheckStatus::Failed,rg::JointEvidenceScope::ComponentLocal,.1,1e-10,"too-large"});
    result.components.push_back(component); result.assembled_state=component.state;
    result.objective=.11234567890123456; result.observation_scale=9;
    result.runtime_convergence=rg::JointCheckStatus::Failed;
    return result;
}
}

TEST(DataObjectPersistenceTest, JointResultsRoundTripCopyClearAndReplaceAtomically)
{
    const command_test::ScopedTempDir dir{"joint_persistence"};
    const auto path=dir.path()/"results.sqlite";
    rg::DataRepository repository{path}; auto model=data_test::MakeModelWithBond();
    const auto record=SavedJointExample(); model->EditAnalysis().SetJointResult(record);
    model->EditAnalysis().ClearTransientFitStates();
    rg::ModelObject copied(*model);
    ASSERT_TRUE(copied.GetAnalysisView().GetJointResult());
    EXPECT_EQ(rg::joint_result_io::Encode(*copied.GetAnalysisView().GetJointResult()),rg::joint_result_io::Encode(record));
    copied.EditAnalysis().Clear(); EXPECT_FALSE(copied.GetAnalysisView().GetJointResult());
    repository.SaveModel(*model,"model"); auto loaded=repository.LoadModel("model");
    ASSERT_TRUE(loaded->GetAnalysisView().GetJointResult());
    EXPECT_EQ(rg::joint_result_io::Encode(*loaded->GetAnalysisView().GetJointResult()),rg::joint_result_io::Encode(record));
    auto invalid=record; invalid.available_row_mask.clear(); model->EditAnalysis().SetJointResult(invalid);
    model->SetPdbID("should-roll-back");
    EXPECT_THROW(repository.SaveModel(*model,"model"),std::invalid_argument);
    loaded=repository.LoadModel("model"); EXPECT_NE(loaded->GetPdbID(),"should-roll-back");
    EXPECT_EQ(rg::joint_result_io::Encode(*loaded->GetAnalysisView().GetJointResult()),rg::joint_result_io::Encode(record));
    model->EditAnalysis().ClearJointResult(); repository.SaveModel(*model,"model");
    EXPECT_FALSE(repository.LoadModel("model")->GetAnalysisView().GetJointResult());
    EXPECT_EQ(data_test::CountRows(path,"model_joint_result"),0);
    model->EditAnalysis().SetJointResult(record); repository.SaveModel(*model,"model");
    data_test::ExecuteSql(path,"DELETE FROM model_object WHERE key_tag='model';");
    EXPECT_EQ(data_test::CountRows(path,"model_joint_result"),0);
}

TEST(DataObjectPersistenceTest, JointContributorSubsetRoundTripsAndRejectsForeignOrHydrogenIds)
{
    const command_test::ScopedTempDir dir{"joint_subset"};
    const auto path=dir.path()/"results.sqlite";
    std::vector<std::unique_ptr<rg::AtomObject>> atoms;
    for(int id=1;id<=4;++id)
    {
        auto atom=std::make_unique<rg::AtomObject>(); atom->SetSerialID(id);
        atom->SetElement(id==4 ? Element::HYDROGEN : Element::CARBON);
        atom->SetPosition(id==3 ? 12. : static_cast<double>(id),0,0);
        atoms.push_back(std::move(atom));
    }
    rg::ModelObject model(std::move(atoms));
    model.SelectAtoms([](const auto & atom){return atom.GetSerialID()==1;});
    const auto record=SavedJointExample(); model.EditAnalysis().SetJointResult(record);
    rg::DataRepository repository{path};
    ASSERT_NO_THROW(repository.SaveModel(model,"subset"));
    auto loaded=repository.LoadModel("subset");
    ASSERT_EQ(loaded->GetAtomList().size(),4u);
    ASSERT_TRUE(loaded->GetAnalysisView().GetJointResult());
    EXPECT_EQ(rg::joint_result_io::Encode(*loaded->GetAnalysisView().GetJointResult()),rg::joint_result_io::Encode(record));
    for(const std::string id:{"999","4","1"})
    {
        auto invalid=record; invalid.atom_ids[1]=id;
        model.EditAnalysis().SetJointResult(invalid);
        EXPECT_THROW(repository.SaveModel(model,"subset"),std::invalid_argument);
        EXPECT_EQ(rg::joint_result_io::Encode(*repository.LoadModel("subset")->GetAnalysisView().GetJointResult()),rg::joint_result_io::Encode(record));
    }
    data_test::ExecuteSql(path,"UPDATE model_joint_result SET result_json=replace(result_json,'\"atom_ids\":[\"1\",\"2\"]','\"atom_ids\":[\"1\",\"4\"]');");
    EXPECT_THROW(repository.LoadModel("subset"),std::invalid_argument);
}

TEST(DataObjectPersistenceTest, JointCodecRejectsMalformedMappingsAndPreservesMissingValues)
{
    namespace io=rg::joint_result_io;
    auto record=SavedJointExample();
    auto invalid=record; invalid.components[0].atoms[0]=99;
    EXPECT_THROW(io::Encode(invalid),std::invalid_argument);
    auto json=boost::json::parse(io::Encode(record)).as_object();
    json["schema_version"]=99; EXPECT_THROW(io::Decode(boost::json::serialize(json)),std::invalid_argument);
    json["schema_version"]=3; json.erase("runtime_convergence");
    EXPECT_THROW(io::Decode(boost::json::serialize(json)),std::exception);
    record.components[0].state.reset(); record.assembled_state.reset(); record.objective.reset();
    record.runtime_convergence=rg::JointCheckStatus::Unavailable;
    record.available_row_mask={false,true};
    const auto decoded=io::Decode(io::Encode(record));
    EXPECT_FALSE(decoded.objective); EXPECT_FALSE(decoded.components[0].state);
    EXPECT_EQ(decoded.available_row_mask,record.available_row_mask);
    EXPECT_EQ(decoded.regular_certificate,rg::JointCheckStatus::NotRun);
    const command_test::ScopedTempDir dir{"joint_export"};
    rg::WriteJointAnalysisResult(decoded,dir.path()/"result.json",dir.path()/"atoms.csv");
    std::ifstream csv(dir.path()/"atoms.csv"); std::string header,line;
    std::getline(csv,header); std::getline(csv,line);
    EXPECT_NE(line.find("\"1\",\"component-0\",,,,0,0,"),std::string::npos);
    EXPECT_THROW(rg::WriteJointAnalysisResult(decoded,dir.path()/"missing"/"result.json",dir.path()/"atoms.csv"),std::runtime_error);
}

TEST(DataObjectPersistenceTest, JointMetadataRejectsInvalidValuesAndOldDocumentsWithoutWrites)
{
    namespace io=rg::joint_result_io;
    const command_test::ScopedTempDir dir{"joint_metadata"};
    const auto path=dir.path()/"results.sqlite";
    rg::DataRepository repository{path}; auto model=data_test::MakeModelWithBond();
    const auto record=SavedJointExample(); model->EditAnalysis().SetJointResult(record);
    repository.SaveModel(*model,"model");
    for(int variant=0;variant<6;++variant)
    {
        auto invalid=record;
        if(variant==0) invalid.metadata.map_sha256="not-a-hash";
        if(variant==1) invalid.metadata.map_normalization->divisor=0;
        if(variant==2) invalid.metadata.map_normalization->requested=false;
        if(variant==3) invalid.metadata.simulation=true;
        if(variant==4) invalid.metadata.map_normalization->applied=false;
        if(variant==5) invalid.metadata.software->build_sha256="";
        model->EditAnalysis().SetJointResult(invalid);
        EXPECT_THROW(repository.SaveModel(*model,"model"),std::invalid_argument);
        EXPECT_EQ(io::Encode(*repository.LoadModel("model")->GetAnalysisView().GetJointResult()),io::Encode(record));
    }
    auto malformed=boost::json::parse(io::Encode(record)).as_object();
    malformed.at("metadata").as_object().at("map_normalization").as_object()["divisor"]=0;
    EXPECT_THROW(io::Decode(boost::json::serialize(malformed)),std::invalid_argument);
    malformed=boost::json::parse(io::Encode(record)).as_object();
    malformed.at("metadata").as_object().erase("map_sha256");
    EXPECT_THROW(io::Decode(boost::json::serialize(malformed)),std::exception);
    auto old=boost::json::parse(io::Encode(record)).as_object(); old["schema_version"]=2;
    const auto old_text=boost::json::serialize(old);
    EXPECT_THROW(io::Decode(old_text),std::invalid_argument);
    data_test::ExecuteSql(path,"UPDATE model_joint_result SET result_json='"+old_text+"' WHERE key_tag='model';");
    std::ifstream before_file(path,std::ios::binary);
    const std::string before((std::istreambuf_iterator<char>(before_file)),{});
    EXPECT_THROW(repository.LoadModel("model"),std::invalid_argument);
    std::ifstream after_file(path,std::ios::binary);
    EXPECT_EQ(std::string((std::istreambuf_iterator<char>(after_file)),{}),before);
    auto wrong_units=boost::json::parse(io::Encode(record)).as_object();
    wrong_units.at("metadata").as_object().at("units").as_object()["B"]="angstrom^2";
    EXPECT_THROW(io::Decode(boost::json::serialize(wrong_units)),std::invalid_argument);
}

TEST(DataObjectPersistenceTest, JointSelectionMetadataValidationAndCsvRoles)
{
    namespace io=rg::joint_result_io;
    const command_test::ScopedTempDir dir{"joint_selection"};
    const auto path=dir.path()/"results.sqlite";
    rg::DataRepository repository{path}; auto model=data_test::MakeModelWithBond();
    const auto saved=SavedJointExample(); model->EditAnalysis().SetJointResult(saved); repository.SaveModel(*model,"model");
    for(int variant=0;variant<7;++variant)
    {
        auto bad=saved;
        if(variant==0) bad.selection_domain->target_indices={0,0};
        if(variant==1) bad.selection_domain->target_indices={1,0};
        if(variant==2) bad.selection_domain->target_indices={2};
        if(variant==3) bad.selection_domain->target_indices.clear();
        if(variant==4) bad.selection_domain->support_radius=3;
        if(variant==5) bad.selection_domain->contract="other";
        if(variant==6) bad.initialization.data_scope="unknown";
        model->EditAnalysis().SetJointResult(bad);
        EXPECT_THROW(repository.SaveModel(*model,"model"),std::invalid_argument);
        EXPECT_EQ(io::Encode(*repository.LoadModel("model")->GetAnalysisView().GetJointResult()),io::Encode(saved));
    }
    rg::WriteJointAnalysisResult(saved,dir.path()/"out.json",dir.path()/"out.csv");
    std::ifstream csv(dir.path()/"out.csv"); std::string line;
    std::getline(csv,line); EXPECT_TRUE(line.ends_with(",SelectionRole"));
    std::getline(csv,line); EXPECT_TRUE(line.ends_with(",target"));
    std::getline(csv,line); EXPECT_TRUE(line.ends_with(",halo"));
    auto raw=saved; raw.selection_domain.reset();
    rg::WriteJointAnalysisResult(raw,dir.path()/"raw.json",dir.path()/"raw.csv");
    std::ifstream raw_csv(dir.path()/"raw.csv"); std::getline(raw_csv,line); std::getline(raw_csv,line);
    EXPECT_TRUE(line.ends_with(",not-recorded"));
}

#include "support/JointPartialSelection.hpp"
#include "data/detail/LocalPotentialEntry.hpp"
#include "data/detail/ModelAnalysisData.hpp"
#include <rhbm_gem/core/GaussianEstimator.hpp>
TEST(DataObjectPersistenceTest, JointNeutralRoundTripPreservesSourcesGeometryAndUnavailableStates)
{
    const command_test::ScopedTempDir dir{"joint_neutral_roundtrip"};
    auto f=joint_partial_test::Make("partial");
    rg::core::FitOptions options; options.estimator=rg::core::PotentialEstimator::JOINT_COMPONENTS; options.quiet_mode=true;
    rg::core::RunPotentialFittingWorkflow(*f.map,*f.model,options);
    rg::DataRepository repository(dir.path()/"result.sqlite");
    repository.SaveModel(*f.model,"joint");
    auto loaded=repository.LoadModel("joint");
    for(int id:{1,2})
    {
        const auto before=rg::AtomLocalPotentialView::For(*f.model->FindAtomPtr(id));
        const auto after=rg::AtomLocalPotentialView::For(*loaded->FindAtomPtr(id));
        const auto & a=before.GetStageEstimate(rg::FittingStage::Second);
        const auto & b=after.GetStageEstimate(rg::FittingStage::Second);
        ASSERT_TRUE(b.point); EXPECT_EQ(a.point->ToVector(),b.point->ToVector());
        EXPECT_EQ(a.source.run_id,b.source.run_id); EXPECT_EQ(a.source.role,b.source.role);
        EXPECT_EQ(a.uncertainty.status,b.uncertainty.status); EXPECT_EQ(a.uncertainty.reason,b.uncertainty.reason);
        EXPECT_EQ(a.uncertainty.covariance.has_value(),b.uncertainty.covariance.has_value());
        if(a.uncertainty.covariance) EXPECT_EQ(*a.uncertainty.covariance,*b.uncertainty.covariance);
        ASSERT_TRUE(after.GetPostFitPeeling()); EXPECT_EQ(before.GetPostFitPeeling()->mode,after.GetPostFitPeeling()->mode);
        const auto raw=before.GetRawSamplingEntries(false), saved=after.GetRawSamplingEntries(false);
        ASSERT_EQ(raw.size(),saved.size()); EXPECT_TRUE(after.HasSampleGeometry());
        for(std::size_t i=0;i<raw.size();++i) {EXPECT_EQ(raw[i].point.position,saved[i].point.position); EXPECT_EQ(raw[i].response,saved[i].response); EXPECT_EQ(before.GetPostFitPeeling()->samples[i].response,after.GetPostFitPeeling()->samples[i].response); EXPECT_EQ(before.GetPostFitPeeling()->samples[i].reason,after.GetPostFitPeeling()->samples[i].reason);}
        EXPECT_EQ(before.GetGroupEvidence()->status,after.GetGroupEvidence()->status);
        EXPECT_THROW(after.GetEstimateMDPDE(rg::FittingStage::Second),std::runtime_error);
    }
    {
        rg::SQLiteWrapper db(dir.path()/"result.sqlite");
        db.Prepare("SELECT amplitude_estimate_mdpde_2nd IS NULL FROM model_atom_local_potential WHERE key_tag='joint';");
        rg::SQLiteWrapper::StatementGuard guard(db); ASSERT_EQ(db.StepNext(),rg::SQLiteWrapper::StepRow()); EXPECT_EQ(db.GetColumn<int>(0),1);
    }
    auto * target=f.model->FindAtomPtr(1);
    auto changed=rg::AtomLocalPotentialView::For(*target).GetStageEstimate(rg::FittingStage::Second);
    const auto original=*changed.point; changed.point=original.WithAmplitude(original.GetAmplitude()+1);
    const auto snapshot = *f.model->GetAnalysisView().GetJointResult();
    f.model->EditAnalysis().SetAtomStageEstimate(rg::FittingStage::Second,*target,changed);
    EXPECT_FALSE(f.model->GetAnalysisView().GetJointResult());
    f.model->EditAnalysis().SetJointResult(snapshot); // An explicitly reattached stale snapshot must still fail validation.
    EXPECT_THROW(repository.SaveModel(*f.model,"joint"),std::invalid_argument);
    const auto preserved=repository.LoadModel("joint");
    EXPECT_EQ(rg::AtomLocalPotentialView::For(*preserved->FindAtomPtr(1)).GetFinalModel(rg::FittingStage::Second).ToVector(),original.ToVector());
}

TEST(DataObjectPersistenceTest, V17ReadsRemainUnmodifiedAndFirstSaveUpgradeRollsBackAtomically)
{
    const command_test::ScopedTempDir dir{"v17_transactional_upgrade"}; const auto path=dir.path()/"legacy.sqlite";
    auto model=data_test::MakeModelWithBond(); model->SelectAllAtoms(); model->EditAnalysis().InitializeFromSelection();
    model->EditAnalysis().SetAtomLocalRawSamplingEntries(*model->FindAtomPtr(1),{{2.0,{0.3,{1,2,3},true}}});
    {rg::DataRepository repository(path); repository.SaveModel(*model,"old");}
    data_test::ExecuteSql(path,"DROP TABLE model_stage_result;"); data_test::ExecuteSql(path,"PRAGMA user_version=17;");
    const auto read_bytes=[&]() { std::ifstream file(path,std::ios::binary); return std::string(std::istreambuf_iterator<char>(file),{}); };
    const auto before=read_bytes();
    rg::DataRepository repository(path); auto legacy=repository.LoadModel("old");
    EXPECT_EQ(read_bytes(),before);
    EXPECT_EQ(data_test::GetUserVersion(path),17); EXPECT_FALSE(data_test::HasTable(path,"model_stage_result"));
    EXPECT_FALSE(rg::AtomLocalPotentialView::For(*legacy->FindAtomPtr(1)).HasSampleGeometry());
    auto & entry=rg::ModelAnalysisData::Of(*model).EnsureAtomLocalEntry(*model->FindAtomPtr(1));
    auto stage=entry.StageEstimate(rg::FittingStage::Second); stage.source.atom_id="wrong-identity"; entry.SetStageEstimate(rg::FittingStage::Second,stage);
    EXPECT_THROW(repository.SaveModel(*model,"new"),std::invalid_argument);
    EXPECT_EQ(data_test::GetUserVersion(path),17); EXPECT_FALSE(data_test::HasTable(path,"model_stage_result"));
    EXPECT_EQ(data_test::CountRows(path,"model_object","old"),1);
    repository.SaveModel(*legacy,"new");
    EXPECT_EQ(data_test::GetUserVersion(path),18); EXPECT_TRUE(data_test::HasTable(path,"model_stage_result"));
    EXPECT_EQ(data_test::CountRows(path,"model_object","old"),1);
    EXPECT_FALSE(rg::AtomLocalPotentialView::For(*repository.LoadModel("new")->FindAtomPtr(1)).HasSampleGeometry());
}

TEST(DataObjectPersistenceTest, LegacyJointSnapshotMapsPointsWithoutInventingDerivedEvidence)
{
    const command_test::ScopedTempDir dir{"legacy_joint_adaptation"}; const auto path=dir.path()/"legacy.sqlite";
    auto model=data_test::MakeModelWithBond(); const auto snapshot=SavedJointExample();
    model->EditAnalysis().SetJointResult(snapshot);
    { rg::DataRepository repository(path); repository.SaveModel(*model,"joint"); }
    data_test::ExecuteSql(path,"DROP TABLE model_stage_result;"); data_test::ExecuteSql(path,"PRAGMA user_version=17;");
    rg::DataRepository repository(path); const auto loaded=repository.LoadModel("joint");
    const auto target=rg::AtomLocalPotentialView::For(*loaded->FindAtomPtr(1));
    EXPECT_DOUBLE_EQ(target.GetFinalModel(rg::FittingStage::Second).GetAmplitude(),snapshot.components[0].state->ac[0]);
    EXPECT_EQ(target.GetStageEstimate(rg::FittingStage::Second).uncertainty.status,rg::EvidenceStatus::NotRun);
    EXPECT_FALSE(target.GetPostFitPeeling()); EXPECT_FALSE(target.GetGroupMemberResult()); EXPECT_FALSE(target.HasSampleGeometry());
    EXPECT_FALSE(target.GetLocalFittingPeelingRatio(1,2)); EXPECT_EQ(data_test::GetUserVersion(path),17);
}
