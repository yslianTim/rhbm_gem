#include <gtest/gtest.h>

#include <stdexcept>
#include <string>
#include <fstream>
#include <limits>
#include <map>
#include <boost/json.hpp>
#include <rhbm_gem/data/io/JointAnalysisFileIO.hpp>
#include "data/io/detail/JointResultJson.hpp"
#include "data/io/detail/StageResultJson.hpp"
#include "core/detail/FirstStageInitialization.hpp"
#include "core/detail/PotentialFittingWorkflow.hpp"

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

    EXPECT_FALSE(data_test::HasTable(database_path,"model_atom_local_potential"));
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
    EXPECT_FALSE(data_test::HasTable(database_path, "model_atom_posterior"));
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
        "CREATE TRIGGER reject_group_insert BEFORE INSERT ON model_stage_result "
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
void UseHistoricalGlobalRankKeys(boost::json::object & json)
{
    if(json.contains("legacy_global_rank_evidence"))
    {
        json["ranks"]=json.at("legacy_global_rank_evidence");
        json.erase("legacy_global_rank_evidence");
    }
    for(auto & item:json.at("components").as_array())
    {
        auto & component=item.as_object();
        if(component.contains("legacy_global_rank_evidence"))
        {
            component["ranks"]=component.at("legacy_global_rank_evidence");
            component.erase("legacy_global_rank_evidence");
        }
    }
}
rg::JointAnalysisResult CurrentJointExample(bool with_halo=true)
{
    rg::JointAnalysisResult result;
    result.parameterization_contract=std::string(rg::JointPureHaloSupportExclusionParameterizationContract);
    result.metadata.model_sha256=std::string(64,'a');
    result.metadata.map_sha256=std::string(64,'b');
    result.metadata.software=rg::JointSoftwareProvenance{"saved-version",std::string(64,'c'),
        std::string(64,'d'),std::string(64,'e')};
    result.metadata.map_normalization=rg::JointMapNormalization{true,true,2.5};
    result.metadata.solver=rg::JointSolverProvenance{};
    auto & solver=*result.metadata.solver;
    solver.search_method="fixed-neighbor";
    solver.fixed_neighbor_core_atoms=12;
    solver.contract_version=std::string(rg::JointSolverProvenanceContractVersion);
    solver.sparse_backend="EIGEN";
    solver.fixed_neighbor_policy_version=std::string(rg::FixedNeighborPolicyContractVersion);
    solver.fixed_neighbor_maximum_sweeps=200;
    solver.fixed_neighbor_order="forward";
    solver.fixed_neighbor_local_search=std::string(rg::FixedNeighborBoundedLocalSearchName);
    result.atom_ids=with_halo ? std::vector<std::string>{"1","2"} : std::vector<std::string>{"1"};
    result.row_ids={"voxel-1","constant"};
    result.fitting_row_ids={"voxel-1"};
    result.fitting_to_original_contributors={0};
    result.fitting_to_original_rows={0};
    result.available_row_mask={true,false};
    result.selection_domain=rg::JointSelectionDomain{};
    result.selection_domain->target_indices={0};
    result.selection_domain->contributor_policy=std::string(rg::PureHaloSupportExclusionContributorPolicy);
    result.selection_domain->pure_halo_support_exclusion=rg::JointPureHaloSupportExclusionProvenance{
        2,1,1,with_halo ? 1u:0u,with_halo ? 1u:0u,0};
    result.initialization.valid=true;
    result.initialization.reason="valid-widths";
    result.initialization.b={0.5};
    result.objective=0.012345678901234567;
    result.observation_scale=9;
    result.runtime_convergence=rg::JointCheckStatus::Failed;
    result.target_runtime_convergence=rg::JointCheckStatus::Failed;
    rg::JointAnalysisComponent component;
    component.id="component-0";
    component.atoms={0};
    component.rows={0};
    component.stop_reason="budget";
    component.search_sweeps=120;
    component.first_order_stationarity_sweep=107;
    component.confirmed_stationarity_sweep=108;
    component.endpoint_certification_sweep=108;
    component.state=rg::JointState{{2.123456789012345,-0.3},{0.5},{-0.7},{1e-14},0.012345678901234567};
    component.runtime_convergence=rg::JointCheckStatus::Failed;
    component.target_runtime_convergence=rg::JointCheckStatus::Failed;
    component.evidence.push_back({"active-set-lm-subproblems",rg::JointCheckStatus::Passed,
        rg::JointEvidenceScope::ComponentLocal,1.0,{},"solves=1"});
    result.components.push_back(component);
    result.assembled_state=component.state;
    result.evidence.push_back({"global-width-calibration",rg::JointCheckStatus::Passed,
        rg::JointEvidenceScope::AssembledGlobal,0.6,2.4,"qualified-first-stage-widths"});
    return result;
}
rg::JointAnalysisResult SavedJointExample() { return CurrentJointExample(); }
}

TEST(DataObjectPersistenceTest, JointCurrentSchemaV19RoundTripsUnqualifiedResult)
{
    namespace io=rg::joint_result_io;
    const auto result=CurrentJointExample();
    const auto encoded=io::Encode(result);
    const auto object=boost::json::parse(encoded).as_object();
    EXPECT_EQ(boost::json::value_to<int>(object.at("schema_version")),19);
    for(const auto key:{"fixed_halos","layout","legacy_global_rank_evidence","ranks"})
        EXPECT_FALSE(object.contains(key));
    EXPECT_FALSE(object.at("costs").as_object().contains("search_reference_seconds"));
    const auto & component=object.at("components").as_array().front().as_object();
    for(const auto key:{"reference_evaluations","legacy_global_rank_evidence","legacy_candidate_evidence","layout"})
        EXPECT_FALSE(component.contains(key));
    EXPECT_FALSE(component.at("state").as_object().contains("nuisance_amplitudes"));
    const auto decoded=io::Decode(encoded);
    EXPECT_EQ(decoded.parameterization_contract,rg::JointPureHaloSupportExclusionParameterizationContract);
    EXPECT_EQ(decoded.runtime_convergence,rg::JointCheckStatus::Failed);
    ASSERT_TRUE(decoded.components.front().state);
    EXPECT_DOUBLE_EQ(decoded.components.front().state->ac.front(),2.123456789012345);
    ASSERT_TRUE(decoded.metadata.solver);
    EXPECT_EQ(decoded.metadata.solver->fixed_neighbor_local_search,rg::FixedNeighborBoundedLocalSearchName);
}

TEST(DataObjectPersistenceTest, JointCurrentSchemaV19RoundTripsQualifiedTargetWithoutHalo)
{
    namespace io=rg::joint_result_io;
    auto result=CurrentJointExample(false);
    result.search_completed=true;
    result.runtime_convergence=rg::JointCheckStatus::Passed;
    result.target_runtime_convergence=rg::JointCheckStatus::Passed;
    result.components.front().search_completed=true;
    result.components.front().runtime_convergence=rg::JointCheckStatus::Passed;
    result.components.front().target_runtime_convergence=rg::JointCheckStatus::Passed;
    ASSERT_TRUE(rg::IsQualifiedJointAnalysisResult(result));
    const auto decoded=io::Decode(io::Encode(result));
    EXPECT_TRUE(rg::IsQualifiedJointAnalysisResult(decoded));
    EXPECT_FALSE(decoded.halo_median_reconstruction);
}

TEST(DataObjectPersistenceTest, JointCurrentSchemaV19RoundTripsGeometryLocalHaloMedian)
{
    namespace io=rg::joint_result_io;
    auto result=CurrentJointExample();
    result.search_completed=true;
    result.runtime_convergence=rg::JointCheckStatus::Passed;
    result.target_runtime_convergence=rg::JointCheckStatus::Passed;
    result.components.front().search_completed=true;
    result.components.front().runtime_convergence=rg::JointCheckStatus::Passed;
    result.components.front().target_runtime_convergence=rg::JointCheckStatus::Passed;
    rg::JointHaloMedianSnapshot snapshot;
    snapshot.joint_run_id="deterministic-run";
    snapshot.global_median={2.123456789012345,0.5,-0.3};
    snapshot.valid_target_count=1;
    rg::JointHaloMedianReconstruction halo;
    halo.halo_contributor_index=1;
    halo.atom_id="2";
    halo.amplitude=2.123456789012345;
    halo.width=0.5;
    halo.charge=-0.3;
    halo.source=rg::HaloMedianSource::GeometryLocalTargetMedian;
    halo.donor_count=1;
    halo.joint_run_id=snapshot.joint_run_id;
    halo.donor_atom_ids={"1"};
    snapshot.halos.push_back(halo);
    result.halo_median_reconstruction=snapshot;
    const auto decoded=io::Decode(io::Encode(result));
    ASSERT_TRUE(decoded.halo_median_reconstruction);
    ASSERT_EQ(decoded.halo_median_reconstruction->halos.size(),1u);
    EXPECT_EQ(decoded.halo_median_reconstruction->halos.front().source,
        rg::HaloMedianSource::GeometryLocalTargetMedian);
    EXPECT_EQ(decoded.halo_median_reconstruction->halos.front().donor_atom_ids,
        std::vector<std::string>{"1"});
}

TEST(DataObjectPersistenceTest, JointCurrentSchemaV19PreservesUnavailableAndNullableDiagnostics)
{
    namespace io=rg::joint_result_io;
    auto unavailable=CurrentJointExample();
    unavailable.components.clear();
    unavailable.initialization.valid=false;
    unavailable.initialization.b.clear();
    unavailable.initialization.reason="width-calibration-unavailable";
    unavailable.assembled_state.reset();
    unavailable.objective.reset();
    unavailable.runtime_convergence=rg::JointCheckStatus::Unavailable;
    unavailable.target_runtime_convergence=rg::JointCheckStatus::NotRun;
    unavailable.evidence.push_back({"global-width-calibration",rg::JointCheckStatus::Unavailable,
        rg::JointEvidenceScope::AssembledGlobal,{}, {},"calibration-unavailable:no-qualified-donors"});
    const auto decoded=io::Decode(io::Encode(unavailable));
    EXPECT_FALSE(decoded.initialization.valid);
    EXPECT_EQ(decoded.runtime_convergence,rg::JointCheckStatus::Unavailable);
    EXPECT_FALSE(decoded.objective);

    auto halo_unavailable=CurrentJointExample();
    halo_unavailable.evidence.push_back({std::string(rg::JointHaloMedianReconstructionEvidenceName),
        rg::JointCheckStatus::Unavailable,rg::JointEvidenceScope::AssembledGlobal,{}, {},
        "halo-median-no-valid-targets"});
    const auto halo_round_trip=io::Decode(io::Encode(halo_unavailable));
    EXPECT_TRUE(rg::HasUnavailableJointHaloMedianReconstruction(halo_round_trip));

    auto nullable=CurrentJointExample();
    nullable.components.front().state.reset();
    nullable.components.front().candidate_diagnostic=rg::JointCandidateDiagnostic{};
    nullable.assembled_state.reset();
    nullable.objective.reset();
    const auto encoded=io::Encode(nullable);
    const auto object=boost::json::parse(encoded).as_object();
    EXPECT_TRUE(object.at("objective").is_null());
    EXPECT_TRUE(object.at("components").as_array().front().as_object().at("candidate_diagnostic")
        .as_object().at("candidate_objective").is_null());
    const auto round_trip=io::Decode(encoded);
    EXPECT_TRUE(round_trip.components.front().candidate_diagnostic);
    EXPECT_FALSE(round_trip.components.front().candidate_diagnostic->candidate_objective);
}

TEST(DataObjectPersistenceTest, LocalStabilityBoundedProbeDiagnosticsRoundTrip)
{
    namespace io=rg::joint_result_io;
    auto record=SavedJointExample();
    rg::JointParameterStability stability;
    rg::JointBlockParameterStability block; block.block_id="component-0/block-1";
    block.status="Sensitive"; block.reason="bounded-local-directional-perturbation-assessed";
    block.baseline_numerical_validity_status="Valid";
    block.baseline_conditional_optimality_status="PracticalOnly";
    block.baseline_eligibility_status="Eligible";
    block.probe_selected=true;
    rg::JointStabilityProbe probe; probe.direction="common-mode";
    probe.source="normalized-all-width-direction"; probe.status="OneSidedAssessed";
    probe.unit_direction={1.}; probe.minus_seed_status="Assessed";
    probe.minus_requested_displacement_inf=.001; probe.minus_effective_displacement_inf=.001;
    probe.plus_seed_status="ZeroDisplacement";
    probe.plus_seed_reason="upper-bound-projection-produced-zero-displacement";
    probe.plus_requested_displacement_inf=.001; probe.plus_effective_displacement_inf=0.;
    block.probes.push_back(probe); stability.blocks.push_back(block);
    record.components.front().parameter_stability=stability;
    const auto decoded=io::Decode(io::Encode(record));
    ASSERT_TRUE(decoded.components.front().parameter_stability);
    const auto & saved=decoded.components.front().parameter_stability->blocks.front();
    EXPECT_EQ(saved.baseline_conditional_optimality_status,"PracticalOnly");
    EXPECT_EQ(saved.baseline_eligibility_status,"Eligible");
    EXPECT_TRUE(saved.probe_selected);
    ASSERT_EQ(saved.probes.size(),1u);
    EXPECT_EQ(saved.probes.front().status,"OneSidedAssessed");
    EXPECT_EQ(saved.probes.front().plus_seed_status,"ZeroDisplacement");
    EXPECT_DOUBLE_EQ(*saved.probes.front().plus_effective_displacement_inf,0.);
}

TEST(DataObjectPersistenceTest, JointResultsRoundTripCopyClearAndReplaceAtomically)
{
    const command_test::ScopedTempDir dir{"joint_persistence"};
    const auto path=dir.path()/"results.sqlite";
    rg::DataRepository repository{path}; auto model=data_test::MakeModelWithBond();
    const auto record=SavedJointExample(); model->EditAnalysis().ApplyJointResult(record,"imported-joint");
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
    model->EditAnalysis().ApplyJointResult(record,"imported-joint"); repository.SaveModel(*model,"model");
    data_test::ExecuteSql(path,"DELETE FROM model_object WHERE key_tag='model';");
    EXPECT_EQ(data_test::CountRows(path,"model_joint_result"),0);
}

TEST(DataObjectPersistenceTest, JointSolverProvenanceRoundTripsAndMissingFieldRemainsReadable)
{
    namespace io=rg::joint_result_io;
    const auto fixed=SavedJointExample();
    const auto decoded=io::Decode(io::Encode(fixed));
    ASSERT_TRUE(decoded.metadata.solver);
    EXPECT_EQ(decoded.metadata.solver->search_method,"fixed-neighbor");
    EXPECT_EQ(decoded.metadata.solver->fixed_neighbor_core_atoms,12u);
    EXPECT_EQ(decoded.metadata.solver->fixed_neighbor_local_search,"active-set-lm");

    auto old=boost::json::parse(io::Encode(fixed)).as_object();
    old.at("metadata").as_object().erase("solver");
    EXPECT_ANY_THROW(io::Decode(boost::json::serialize(old)));
}

TEST(DataObjectPersistenceTest, FixedNeighborFailureEvidenceAndSweepCountsRoundTrip)
{
    namespace io=rg::joint_result_io;
    auto result=SavedJointExample();
    result.components.front().evidence.push_back({"endpoint-trust",rg::JointCheckStatus::Failed,
        rg::JointEvidenceScope::ComponentLocal,{},{},"reference-disagreement"});
    const auto decoded=io::Decode(io::Encode(result));
    ASSERT_EQ(decoded.components.size(),1u);
    const auto & component=decoded.components.front();
    EXPECT_EQ(component.search_sweeps,120u);
    EXPECT_EQ(component.first_order_stationarity_sweep,107u);
    EXPECT_EQ(component.confirmed_stationarity_sweep,108u);
    EXPECT_EQ(component.endpoint_certification_sweep,108u);
    ASSERT_FALSE(component.evidence.empty());
    EXPECT_EQ(component.evidence.back().reason,"reference-disagreement");

    auto missing_sweeps=boost::json::parse(io::Encode(result)).as_object();
    auto & saved_component=missing_sweeps.at("components").as_array().front().as_object();
    for(const auto key:{"search_sweeps","first_order_stationarity_sweep",
        "confirmed_stationarity_sweep","endpoint_certification_sweep"}) saved_component.erase(key);
    EXPECT_ANY_THROW(io::Decode(boost::json::serialize(missing_sweeps)));
}

TEST(DataObjectPersistenceTest, CurrentFixedNeighborSolverProvenanceRoundTrips)
{
    namespace io=rg::joint_result_io;
    auto fixed=SavedJointExample();
    rg::JointSolverProvenance fixed_provenance{"fixed-neighbor",12u,std::nullopt};
    fixed_provenance.contract_version=rg::JointSolverProvenanceContractVersion;
    fixed_provenance.sparse_backend="EIGEN";
    fixed_provenance.fixed_neighbor_policy_version=rg::FixedNeighborPolicyContractVersion;
    fixed_provenance.fixed_neighbor_maximum_sweeps=200;
    fixed_provenance.fixed_neighbor_order="forward";
    fixed_provenance.fixed_neighbor_local_search="active-set-lm";
    fixed.metadata.solver=fixed_provenance;
    const auto fixed_json=io::Encode(fixed);
    const auto fixed_object=boost::json::parse(fixed_json).as_object();
    auto downgraded_schema=fixed_object;
    UseHistoricalGlobalRankKeys(downgraded_schema);
    downgraded_schema["schema_version"]=7;
    EXPECT_ANY_THROW(io::Decode(boost::json::serialize(downgraded_schema)));
    EXPECT_EQ(fixed_object.at("metadata").as_object().at("solver").as_object()
        .at("fixed_neighbor_local_search").as_string(),"active-set-lm");
    const auto fixed_decoded=io::Decode(fixed_json);
    ASSERT_TRUE(fixed_decoded.metadata.solver);
    EXPECT_EQ(fixed_decoded.metadata.solver->contract_version,rg::JointSolverProvenanceContractVersion);
    EXPECT_EQ(fixed_decoded.metadata.solver->sparse_backend,"EIGEN");
    EXPECT_FALSE(fixed_decoded.metadata.solver->fixed_neighbor_local_work);
    EXPECT_EQ(fixed_decoded.metadata.solver->fixed_neighbor_maximum_sweeps,200u);
    EXPECT_EQ(fixed_decoded.metadata.solver->fixed_neighbor_order,"forward");
    EXPECT_EQ(fixed_decoded.metadata.solver->fixed_neighbor_local_search,"active-set-lm");
    EXPECT_FALSE(fixed_decoded.metadata.solver->preconditioner);
    EXPECT_FALSE(fixed_decoded.metadata.solver->operator_rank_mode);
    EXPECT_FALSE(fixed_decoded.metadata.solver->operator_rank_backend);
    EXPECT_FALSE(fixed_decoded.metadata.solver->operator_pcg_iterations);
    EXPECT_FALSE(fixed_decoded.metadata.solver->operator_damping_trials);
    EXPECT_FALSE(fixed_decoded.metadata.solver->operator_rank_budget_seconds);
    EXPECT_FALSE(fixed_decoded.metadata.solver->operator_rank_budget_entries);
    EXPECT_FALSE(fixed_decoded.metadata.solver->operator_rank_budget_workspace_bytes);
    EXPECT_FALSE(fixed_decoded.metadata.solver->schwarz_core_atoms);
    EXPECT_FALSE(fixed_decoded.metadata.solver->schwarz_overlap_hops);
    EXPECT_FALSE(fixed_decoded.metadata.solver->schwarz_max_block_atoms);
    EXPECT_FALSE(fixed_decoded.metadata.solver->schwarz_storage_bytes);
    EXPECT_FALSE(fixed_decoded.metadata.solver->schwarz_scratch_bytes);

    auto bounded=fixed;
    bounded.metadata.solver->fixed_neighbor_local_search="active-set-lm";
    const auto bounded_decoded=io::Decode(io::Encode(bounded));
    ASSERT_TRUE(bounded_decoded.metadata.solver);
    EXPECT_EQ(bounded_decoded.metadata.solver->fixed_neighbor_local_search,"active-set-lm");

    auto retired_search=boost::json::parse(fixed_json).as_object();
    retired_search.at("metadata").as_object().at("solver").as_object()["fixed_neighbor_local_search"]=
        "legacy-compact";
    EXPECT_THROW(io::Decode(boost::json::serialize(retired_search)),std::invalid_argument);
}

TEST(DataObjectPersistenceTest, MalformedCurrentSolverProvenanceIsRejected)
{
    namespace io=rg::joint_result_io;
    auto fixed=SavedJointExample();
    rg::JointSolverProvenance provenance{"fixed-neighbor",12u,std::nullopt};
    provenance.contract_version=rg::JointSolverProvenanceContractVersion;
    provenance.sparse_backend="EIGEN";
    provenance.fixed_neighbor_policy_version=rg::FixedNeighborPolicyContractVersion;
    provenance.fixed_neighbor_maximum_sweeps=200;
    provenance.fixed_neighbor_order="forward";
    provenance.fixed_neighbor_local_search="active-set-lm";
    fixed.metadata.solver=provenance;

    auto malformed=boost::json::parse(io::Encode(fixed)).as_object();
    malformed.at("metadata").as_object().at("solver").as_object().erase("sparse_backend");
    EXPECT_THROW(io::Decode(boost::json::serialize(malformed)),std::invalid_argument);

    malformed=boost::json::parse(io::Encode(fixed)).as_object();
    malformed.at("metadata").as_object().at("solver").as_object()["sparse_backend"]="unsupported";
    EXPECT_THROW(io::Decode(boost::json::serialize(malformed)),std::invalid_argument);

    malformed=boost::json::parse(io::Encode(fixed)).as_object();
    malformed.at("metadata").as_object().at("solver").as_object()["search_method"]="unsupported";
    EXPECT_THROW(io::Decode(boost::json::serialize(malformed)),std::invalid_argument);

    malformed=boost::json::parse(io::Encode(fixed)).as_object();
    malformed.at("metadata").as_object().at("solver").as_object()["preconditioner"]="Identity";
    EXPECT_THROW(io::Decode(boost::json::serialize(malformed)),std::invalid_argument);

    malformed=boost::json::parse(io::Encode(fixed)).as_object();
    malformed.at("metadata").as_object().at("solver").as_object()["fixed_neighbor_order"]="sideways";
    EXPECT_THROW(io::Decode(boost::json::serialize(malformed)),std::invalid_argument);
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
    const auto record=SavedJointExample(); model.EditAnalysis().ApplyJointResult(record,"imported-joint");
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
    UseHistoricalGlobalRankKeys(json);
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

TEST(DataObjectPersistenceTest, UncertifiedCandidateDiagnosticsRoundTripWithoutEnteringScientificExports)
{
    namespace io=rg::joint_result_io;
    const command_test::ScopedTempDir dir{"joint_candidate_diagnostic"};
    auto record=SavedJointExample();
    auto & component=record.components.front();
    component.state.reset(); record.assembled_state.reset(); record.objective.reset();
    component.runtime_convergence=rg::JointCheckStatus::Unavailable;
    record.runtime_convergence=rg::JointCheckStatus::Unavailable;
    component.target_evidence.reset(); record.target_evidence.reset();
    rg::JointCandidateDiagnostic candidate;
    candidate.candidate_available=true;
    candidate.candidate_eta={-.7}; candidate.candidate_ac={999.125,-.3};
    candidate.candidate_objective=.0123; candidate.candidate_ac_kkt=.1;
    candidate.candidate_width_gradient_inf_norm=.2; candidate.candidate_feasible=true;
    candidate.convergence_criteria_id="practical-v1";
    candidate.global_ac_kkt_threshold=1e-7; candidate.global_width_gradient_threshold=1e-9;
    candidate.confirmed_log_width_change_threshold=1e-6; candidate.candidate_log_width_change_inf_norm=5e-7;
    candidate.certification_status=rg::JointCheckStatus::Unavailable;
    candidate.endpoint_verification_status=rg::JointCheckStatus::Passed;
    candidate.numerical_qualification_status=rg::JointCheckStatus::Unavailable;
    candidate.verification_reason="endpoint-verified-but-unqualified";
    candidate.resource_limit_reason="resource-limit:qualification:budget-exceeded";
    component.candidate_diagnostic=candidate;

    const auto json=boost::json::parse(io::Encode(record));
    EXPECT_EQ(json.at("schema_version").as_int64(),19);
    const auto & saved_candidate=json.at("components").as_array()[0].at("candidate_diagnostic");
    EXPECT_EQ(saved_candidate.at("candidate_status").as_string(),"uncertified");
    EXPECT_EQ(saved_candidate.at("convergence_criteria_id").as_string(),"practical-v1");
    EXPECT_DOUBLE_EQ(saved_candidate.at("global_ac_kkt_threshold").as_double(),1e-7);
    EXPECT_DOUBLE_EQ(saved_candidate.at("global_width_gradient_threshold").as_double(),1e-9);
    EXPECT_DOUBLE_EQ(saved_candidate.at("confirmed_log_width_change_threshold").as_double(),1e-6);
    EXPECT_DOUBLE_EQ(saved_candidate.at("candidate_log_width_change_inf_norm").as_double(),5e-7);
    const auto decoded=io::Decode(boost::json::serialize(json));
    ASSERT_TRUE(decoded.components.front().candidate_diagnostic);
    EXPECT_EQ(decoded.components.front().id,component.id);
    EXPECT_EQ(decoded.components.front().search_sweeps,component.search_sweeps);
    EXPECT_EQ(decoded.components.front().stop_reason,component.stop_reason);
    EXPECT_FALSE(decoded.components.front().state);
    EXPECT_EQ(decoded.components.front().runtime_convergence,rg::JointCheckStatus::Unavailable);
    EXPECT_EQ(decoded.metadata.model_sha256,record.metadata.model_sha256);
    ASSERT_TRUE(decoded.metadata.software);
    EXPECT_EQ(decoded.metadata.software->configuration_sha256,record.metadata.software->configuration_sha256);
    EXPECT_EQ(decoded.components.front().candidate_diagnostic->candidate_ac,candidate.candidate_ac);
    EXPECT_EQ(decoded.components.front().candidate_diagnostic->certification_status,rg::JointCheckStatus::Unavailable);
    EXPECT_EQ(decoded.components.front().candidate_diagnostic->endpoint_verification_status,rg::JointCheckStatus::Passed);
    EXPECT_EQ(decoded.components.front().candidate_diagnostic->numerical_qualification_status,rg::JointCheckStatus::Unavailable);
    EXPECT_EQ(decoded.components.front().candidate_diagnostic->convergence_criteria_id,"practical-v1");
    ASSERT_TRUE(decoded.components.front().candidate_diagnostic->global_ac_kkt_threshold);
    EXPECT_DOUBLE_EQ(*decoded.components.front().candidate_diagnostic->global_ac_kkt_threshold,1e-7);
    ASSERT_TRUE(decoded.components.front().candidate_diagnostic->candidate_log_width_change_inf_norm);
    EXPECT_DOUBLE_EQ(*decoded.components.front().candidate_diagnostic->candidate_log_width_change_inf_norm,5e-7);
    EXPECT_EQ(decoded.components.front().candidate_diagnostic->verification_reason,
        "endpoint-verified-but-unqualified");
    EXPECT_FALSE(decoded.components.front().target_evidence);

    auto missing_current_criteria=json.as_object();
    missing_current_criteria.at("components").as_array()[0].as_object()
        .at("candidate_diagnostic").as_object().erase("convergence_criteria_id");
    EXPECT_ANY_THROW(io::Decode(boost::json::serialize(missing_current_criteria)));

    rg::WriteJointAnalysisResult(decoded,dir.path()/"result.json",dir.path()/"atoms.csv");
    std::ifstream csv(dir.path()/"atoms.csv"); std::string line;
    while(std::getline(csv,line)) EXPECT_EQ(line.find("999.125"),std::string::npos);
}

TEST(DataObjectPersistenceTest, CandidateDiagnosticsRejectInvalidOrUntrustedCheckpointParameters)
{
    namespace io=rg::joint_result_io;
    auto record=SavedJointExample();
    auto & component=record.components.front();
    component.candidate_diagnostic=rg::JointCandidateDiagnostic{};
    component.candidate_diagnostic->candidate_available=true;
    component.candidate_diagnostic->candidate_eta={-.7};
    component.candidate_diagnostic->candidate_ac={1,2};
    component.candidate_diagnostic->candidate_objective=.1;
    EXPECT_THROW(io::Encode(record),std::invalid_argument);

    component.state.reset(); record.assembled_state.reset(); record.objective.reset();
    component.candidate_diagnostic->candidate_eta={std::numeric_limits<double>::quiet_NaN(),-.5};
    component.candidate_diagnostic->candidate_ac={1,2,3,4};
    EXPECT_THROW(io::Encode(record),std::invalid_argument);

    component.candidate_diagnostic->candidate_available=false;
    component.candidate_diagnostic->candidate_eta.clear();
    component.candidate_diagnostic->candidate_ac.clear();
    component.candidate_diagnostic->candidate_objective.reset();
    EXPECT_NO_THROW(io::Encode(record));
}

TEST(DataObjectPersistenceTest, JointMetadataRejectsInvalidValuesAndOldDocumentsWithoutWrites)
{
    namespace io=rg::joint_result_io;
    const command_test::ScopedTempDir dir{"joint_metadata"};
    const auto path=dir.path()/"results.sqlite";
    rg::DataRepository repository{path}; auto model=data_test::MakeModelWithBond();
    const auto record=SavedJointExample(); model->EditAnalysis().ApplyJointResult(record,"imported-joint");
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
    const auto saved=SavedJointExample(); model->EditAnalysis().ApplyJointResult(saved,"imported-joint"); repository.SaveModel(*model,"model");
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
    std::getline(csv,line); EXPECT_TRUE(line.ends_with(",SelectionRole,Parameterization,ContributionGroupID,TargetRuntimeConvergence,PointQuality"));
    std::getline(csv,line); EXPECT_TRUE(line.ends_with(",target,JointEstimated,\"\",failed,not-run"));
    std::getline(csv,line); EXPECT_TRUE(line.ends_with(",halo,HaloReconstructionUnavailable,\"\",not-run,NotApplicable"));
    auto raw=saved; raw.selection_domain.reset();
    EXPECT_THROW(rg::WriteJointAnalysisResult(raw,dir.path()/"raw.json",dir.path()/"raw.csv"),std::invalid_argument);
}

#include "support/JointPartialSelection.hpp"
#include "data/detail/LocalPotentialEntry.hpp"
#include "data/detail/ModelAnalysisData.hpp"
TEST(DataObjectPersistenceTest, JointNeutralRoundTripPreservesCurrentJointAndStageResults)
{
    const command_test::ScopedTempDir dir{"joint_neutral_roundtrip"};
    auto fixture=joint_partial_test::Make("partial");
    rg::core::FitOptions options;
    options.estimator=rg::core::PotentialEstimator::JOINT_COMPONENTS;
    options.quiet_mode=true;
    rg::core::RunPotentialFittingWorkflow(*fixture.map,*fixture.model,options);
    rg::DataRepository repository(dir.path()/"result.sqlite");
    repository.SaveModel(*fixture.model,"joint");
    auto loaded=repository.LoadModel("joint");
    ASSERT_TRUE(fixture.model->GetAnalysisView().GetJointResult());
    ASSERT_TRUE(loaded->GetAnalysisView().GetJointResult());
    EXPECT_EQ(rg::joint_result_io::Encode(*loaded->GetAnalysisView().GetJointResult()),
        rg::joint_result_io::Encode(*fixture.model->GetAnalysisView().GetJointResult()));
    EXPECT_EQ(rg::stage_result_io::Encode(*loaded),rg::stage_result_io::Encode(*fixture.model));

    loaded->EditAnalysis().ClearJointResult();
    repository.SaveModel(*loaded,"no-snapshot");
    const auto restored=repository.LoadModel("no-snapshot");
    EXPECT_EQ(rg::stage_result_io::Encode(*restored),rg::stage_result_io::Encode(*loaded));
    EXPECT_FALSE(restored->GetAnalysisView().GetJointResult());
    EXPECT_FALSE(data_test::HasTable(dir.path()/"result.sqlite","model_atom_local_potential"));
}

TEST(DataObjectPersistenceTest, VersionTwoJointEvidenceIsInvalidatedDuringNeutralMigration)
{
    auto model = data_test::MakeModelWithBond();
    model->SelectAllAtoms();
    auto editor = model->EditAnalysis();
    editor.InitializeFromSelection();
    std::map<int, rg::LocalStageEstimate> stages;
    for (const auto & atom : model->GetAtomList())
    {
        rg::LocalStageEstimate stage;
        stage.source = {rg::EstimateMethod::JointComponents, std::to_string(atom->GetSerialID()),
            "legacy-component", "legacy-v2-run", rg::FittingRole::Target};
        stage.point = rg::GaussianModel3D(2.0, 0.5, 0.1);
        stage.convergence = rg::JointCheckStatus::Passed;
        stages.emplace(atom->GetSerialID(), std::move(stage));
    }
    editor.ApplySecondStageEstimates(std::move(stages));
    auto old_document = boost::json::parse(rg::stage_result_io::Encode(*model)).as_object();
    old_document["version"] = 2;
    for (auto & atom : old_document["atoms"].as_array())
        atom.as_object()["evidence"] = boost::json::object{
            {"atom", "legacy"}, {"component", "legacy-component"}, {"source", "legacy-v2-run"},
            {"status", 1}, {"reason", ""}, {"method", "local-conditional-iid-ls"},
            {"correlation", "block-diagonal-by-atom"}, {"estimate", boost::json::array{1.0, 0.5}},
            {"covariance", boost::json::array{1.0, 0.0, 0.0, 1.0}}};

    rg::stage_result_io::Decode(*model, boost::json::serialize(old_document), 3, {});
    for (const auto & atom : model->GetAtomList())
    {
        const auto view = rg::AtomLocalPotentialView::For(*atom);
        EXPECT_FALSE(view.GetPostJointLocalResult());
        EXPECT_FALSE(view.GetGroupMemberResult());
    }
    for (const auto key : model->GetAnalysisView().CollectAtomGroupKeys())
    {
        const auto summary = model->GetAnalysisView().GetGroupParameterSummary(key);
        if (!summary) continue;
        EXPECT_EQ(summary->status, rg::EvidenceStatus::Unavailable);
        EXPECT_EQ(summary->reason, "post-joint-local-results-required");
        EXPECT_FALSE(summary->inference);
    }
}

TEST(DataObjectPersistenceTest, CanonicalNativeDiagnosticsAndMissingStatesRoundTripExactly)
{
    const command_test::ScopedTempDir dir{"canonical_diagnostics"}; const auto path=dir.path()/"model.sqlite";
    auto model=data_test::MakeModelWithBond(); auto editor=model->EditAnalysis();
    rg::LocalGaussianResult result;
    result.alpha_r=.23;
    result.ols={rg::GaussianModel3D{2.123456789012345,.4,-.2},rg::GaussianModel3DUncertainty{.12,.13,.14}};
    result.mdpde={rg::GaussianModel3D{3.123456789012345,.5,-.3},rg::GaussianModel3DUncertainty{.21,.22,.23}};
    result.fit_result=rg::RHBMBetaEstimateResult{};
    auto & fit=*result.fit_result; fit.status=rg::RHBMEstimationStatus::MAX_ITERATIONS_REACHED;
    fit.sigma_square=.123456789012345; fit.diagnostics={17,.0001,.0002};
    fit.refinement=rg::RHBMEndpointRefinementDiagnostics{false,"rejected",12,25,"converged",.3,.2,.1,std::array<double,3>{.01,.02,.03},std::numeric_limits<double>::infinity(),false};
    for(auto stage:{rg::FittingStage::First,rg::FittingStage::Second})
        editor.SetAtomLocalGaussianResult(stage,*model->FindAtomPtr(1),result);
    editor.ClearTransientFitStates();
    rg::LocalStageEstimate missing; missing.source.method=rg::EstimateMethod::JointComponents;
    missing.source.run_id="no-snapshot"; missing.reason="initialization-failed"; missing.source.role=rg::FittingRole::Halo;
    editor.SetAtomStageEstimate(rg::FittingStage::Second,*model->FindAtomPtr(2),missing);
    model->SelectAllAtoms(); editor.RebuildAtomGroupsFromSelection();
    const auto group=model->GetAnalysisView().CollectAtomGroupKeys().front();
    editor.SetAtomGroupAlphaG(group,.87);
    rg::GroupParameterSummary summary; summary.status=rg::EvidenceStatus::Unavailable;
    summary.reason="insufficient-evidence"; summary.point_count=1; summary.excluded_count=1;
    summary.descriptive_mean=rg::GaussianModel3D{3.,.5,.1};
    editor.ApplyAtomGroupParameterSummary(group,summary);
    const auto expected=rg::stage_result_io::Encode(*model);
    rg::DataRepository repo(path); repo.SaveModel(*model,"mixed"); const auto loaded=repo.LoadModel("mixed");
    EXPECT_EQ(rg::stage_result_io::Encode(*loaded),expected);
    const auto actual=rg::AtomLocalPotentialView::For(*loaded->FindAtomPtr(1)).GetGaussianResult(rg::FittingStage::Second);
    EXPECT_DOUBLE_EQ(loaded->GetAnalysisView().GetAtomAlphaG(group),.87);
    EXPECT_FALSE(actual.fit_result); ASSERT_TRUE(actual.diagnostics); EXPECT_EQ(actual.diagnostics->iterations.iterations,17);
    EXPECT_EQ(actual.diagnostics->refinement->relative_coordinate_difference,fit.refinement->relative_coordinate_difference);
    const auto document=boost::json::parse(expected).as_object(); EXPECT_EQ(document.at("version").as_int64(),4);
    EXPECT_FALSE(loaded->GetAnalysisView().GetJointResult());
}
