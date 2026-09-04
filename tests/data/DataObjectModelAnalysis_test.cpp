#include <gtest/gtest.h>

#include <array>
#include <limits>
#include <memory>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

#include "data/detail/LocalPotentialEntry.hpp"
#include "data/detail/ModelAnalysisData.hpp"
#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <rhbm_gem/data/object/ModelAnalysisView.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/domain/AtomSelector.hpp>
#include <rhbm_gem/utils/domain/ChemicalDataHelper.hpp>
#include "support/DataObjectTestSupport.hpp"

namespace rg = rhbm_gem;

TEST(DataObjectModelAnalysisTest, StandardQScoresDefaultToZeroAndCanBeSet)
{
    auto model{ data_test::MakeModelWithBond() };
    ASSERT_FALSE(model->GetAtomList().empty());

    EXPECT_DOUBLE_EQ(model->GetStandardAverageQScore(), 0.0);
    EXPECT_DOUBLE_EQ(model->GetAtomList().front()->GetStandardQScore(), 0.0);

    model->SetStandardAverageQScore(0.75);
    model->GetAtomList().front()->SetStandardQScore(0.5);

    EXPECT_DOUBLE_EQ(model->GetStandardAverageQScore(), 0.75);
    EXPECT_DOUBLE_EQ(model->GetAtomList().front()->GetStandardQScore(), 0.5);
}

TEST(DataObjectModelAnalysisTest, ReferenceGaussianParametersDefaultToZeroAndCanBeSet)
{
    auto model{ data_test::MakeModelWithBond() };

    EXPECT_DOUBLE_EQ(model->GetReferenceHeight(), 0.0);
    EXPECT_DOUBLE_EQ(model->GetReferenceOffset(), 0.0);

    model->SetReferenceHeight(1.25);
    model->SetReferenceOffset(-0.5);

    EXPECT_DOUBLE_EQ(model->GetReferenceHeight(), 1.25);
    EXPECT_DOUBLE_EQ(model->GetReferenceOffset(), -0.5);
}

TEST(DataObjectModelAnalysisTest, SelectedModelEntriesCanBeInitializedForTypedWorkflows)
{
    auto model{ data_test::MakeModelWithBond() };
    model->SelectAllAtoms(false);
    model->SelectAllBonds(false);
    ASSERT_EQ(model->GetSelectedAtomCount(), 0);

    model->SelectAllAtoms();
    model->ApplySymmetrySelection(false);
    for (auto * atom : model->GetSelectedAtoms())
    {
        rg::ModelAnalysisData::Of(*model).EnsureAtomLocalEntry(*atom);
    }

    EXPECT_EQ(model->GetSelectedAtomCount(), model->GetNumberOfAtom());
    for (const auto * atom : model->GetSelectedAtoms())
    {
        ASSERT_NE(atom, nullptr);
        EXPECT_TRUE(rg::AtomLocalPotentialView::For(*atom).IsAvailable());
    }
}

TEST(DataObjectModelAnalysisTest, AtomSelectorCanDriveModelSelectionState)
{
    auto model{ data_test::MakeModelWithBond() };
    auto & atom_list{ model->GetAtomList() };
    ASSERT_EQ(atom_list.size(), 2);
    atom_list.at(0)->SetChainID("A");
    atom_list.at(1)->SetChainID("B");

    ::AtomSelector selector;
    selector.PickChainID("A");
    model->SelectAtoms(
        [&selector](const rg::AtomObject & atom)
        {
            return selector.GetSelectionFlag(
                atom.GetChainID(),
                atom.GetResidue(),
                atom.GetElement());
        });

    ASSERT_EQ(model->GetSelectedAtomCount(), 1);
    EXPECT_EQ(model->GetSelectedAtoms().front(), atom_list.at(0).get());
}

TEST(DataObjectModelAnalysisTest, ModelSelectionAndLocalEntriesRemainDirectlyQueryable)
{
    auto model{ data_test::MakeModelWithBond() };
    auto & atoms{ model->GetAtomList() };
    ASSERT_EQ(atoms.size(), 2);
    model->SetAtomSelected(1, true);
    model->SetAtomSelected(2, false);
    rg::ModelAnalysisData::Of(*model).EnsureAtomLocalEntry(*atoms[0]);

    const auto & selected_only_atoms{ model->GetSelectedAtoms() };
    ASSERT_EQ(selected_only_atoms.size(), 1);
    EXPECT_EQ(selected_only_atoms.front(), atoms[0].get());

    std::vector<rg::AtomObject *> require_entry_atoms;
    for (auto & atom : model->GetAtomList())
    {
        if (rg::AtomLocalPotentialView::For(*atom).IsAvailable())
        {
            require_entry_atoms.emplace_back(atom.get());
        }
    }
    ASSERT_EQ(require_entry_atoms.size(), 1);
    EXPECT_EQ(require_entry_atoms.front(), atoms[0].get());
}

TEST(DataObjectModelAnalysisTest, SelectedAtomListCanBeQueriedByResidueId)
{
    std::vector<std::unique_ptr<rg::AtomObject>> atom_list;
    atom_list.reserve(4);

    auto make_atom =
        [](int serial_id, int sequence_id, const std::string & chain_id)
        {
            auto atom{ std::make_unique<rg::AtomObject>() };
            atom->SetSerialID(serial_id);
            atom->SetSequenceID(sequence_id);
            atom->SetChainID(chain_id);
            return atom;
        };
    atom_list.emplace_back(make_atom(1, 10, "A"));
    atom_list.emplace_back(make_atom(2, 10, "A"));
    atom_list.emplace_back(make_atom(3, 20, "A"));
    atom_list.emplace_back(make_atom(4, 10, "B"));

    rg::ModelObject model(std::move(atom_list));
    const auto & atoms{ model.GetAtomList() };
    model.SetAtomSelected(1, true);
    model.SetAtomSelected(2, false);
    model.SetAtomSelected(3, true);
    model.SetAtomSelected(4, true);

    const auto & residue_10_atoms{ model.GetSelectedAtomList(10) };
    ASSERT_EQ(residue_10_atoms.size(), 2u);
    EXPECT_EQ(residue_10_atoms.at(0), atoms.at(0).get());
    EXPECT_EQ(residue_10_atoms.at(1), atoms.at(3).get());

    const auto & residue_20_atoms{ model.GetSelectedAtomList(20) };
    ASSERT_EQ(residue_20_atoms.size(), 1u);
    EXPECT_EQ(residue_20_atoms.front(), atoms.at(2).get());
    EXPECT_TRUE(model.GetSelectedAtomList(999).empty());

    model.SetAtomSelected(1, false);
    model.SetAtomSelected(4, false);
    EXPECT_TRUE(model.GetSelectedAtomList(10).empty());

    model.SelectAllAtoms(false);
    EXPECT_TRUE(model.GetSelectedAtomList(20).empty());
}

TEST(DataObjectModelAnalysisTest, SequenceIDListReportsSortedUniqueModelResidues)
{
    std::vector<std::unique_ptr<rg::AtomObject>> atom_list;
    atom_list.reserve(4);

    auto make_atom =
        [](int serial_id, int sequence_id)
        {
            auto atom{ std::make_unique<rg::AtomObject>() };
            atom->SetSerialID(serial_id);
            atom->SetSequenceID(sequence_id);
            return atom;
        };
    atom_list.emplace_back(make_atom(1, 20));
    atom_list.emplace_back(make_atom(2, 10));
    atom_list.emplace_back(make_atom(3, 10));
    atom_list.emplace_back(make_atom(4, 30));

    rg::ModelObject model(std::move(atom_list));
    model.SelectAllAtoms(false);

    EXPECT_EQ(model.GetSequenceIDList(), (std::vector<int>{ 10, 20, 30 }));

    model.GetAtomList().at(0)->SetSequenceID(5);
    EXPECT_EQ(model.GetSequenceIDList(), (std::vector<int>{ 5, 10, 30 }));
}

TEST(DataObjectModelAnalysisTest, SequenceIDListIsEmptyForEmptyModel)
{
    const rg::ModelObject model;

    EXPECT_TRUE(model.GetSequenceIDList().empty());
}

TEST(DataObjectModelAnalysisTest, ModelObjectCanApplyElementSelectionAsExclusion)
{
    auto model{ data_test::MakeModelWithBond() };
    auto & atoms{ model->GetAtomList() };
    atoms.at(0)->SetElement(Element::CARBON);
    atoms.at(1)->SetElement(Element::HYDROGEN);
    model->SelectAllAtoms();

    model->ApplyElementSelection(Element::HYDROGEN, true);

    ASSERT_EQ(model->GetSelectedAtomCount(), 1u);
    EXPECT_EQ(model->GetSelectedAtoms().front(), atoms.at(0).get());
}

TEST(DataObjectModelAnalysisTest, ModelObjectElementSelectionExclusionDoesNotWidenSelection)
{
    auto model{ data_test::MakeModelWithBond() };
    auto & atoms{ model->GetAtomList() };
    atoms.at(0)->SetElement(Element::CARBON);
    atoms.at(1)->SetElement(Element::HYDROGEN);
    model->SelectAllAtoms(false);
    model->SetAtomSelected(atoms.at(1)->GetSerialID(), true);

    model->ApplyElementSelection(Element::HYDROGEN, true);

    EXPECT_EQ(model->GetSelectedAtomCount(), 0u);
}

TEST(DataObjectModelAnalysisTest, ModelObjectElementSelectionNoOpKeepsSelection)
{
    auto model{ data_test::MakeModelWithBond() };
    auto & atoms{ model->GetAtomList() };
    atoms.at(0)->SetElement(Element::CARBON);
    atoms.at(1)->SetElement(Element::HYDROGEN);
    model->SelectAllAtoms(false);
    model->SetAtomSelected(atoms.at(0)->GetSerialID(), true);

    model->ApplyElementSelection(Element::HYDROGEN, false);

    ASSERT_EQ(model->GetSelectedAtomCount(), 1u);
    EXPECT_EQ(model->GetSelectedAtoms().front(), atoms.at(0).get());
}

TEST(DataObjectModelAnalysisTest, ModelObjectCanApplySpotSelectionAsExclusion)
{
    auto model{ data_test::MakeModelWithBond() };
    auto & atoms{ model->GetAtomList() };
    atoms.at(0)->SetSpot(Spot::CA);
    atoms.at(1)->SetSpot(Spot::H);
    model->SelectAllAtoms();

    model->ApplySpotSelection(Spot::H, true);

    ASSERT_EQ(model->GetSelectedAtomCount(), 1u);
    EXPECT_EQ(model->GetSelectedAtoms().front(), atoms.at(0).get());
}

TEST(DataObjectModelAnalysisTest, ModelObjectSpotSelectionExclusionDoesNotWidenSelection)
{
    auto model{ data_test::MakeModelWithBond() };
    auto & atoms{ model->GetAtomList() };
    atoms.at(0)->SetSpot(Spot::CA);
    atoms.at(1)->SetSpot(Spot::H);
    model->SelectAllAtoms(false);
    model->SetAtomSelected(atoms.at(1)->GetSerialID(), true);

    model->ApplySpotSelection(Spot::H, true);

    EXPECT_EQ(model->GetSelectedAtomCount(), 0u);
}

TEST(DataObjectModelAnalysisTest, ModelObjectSpotSelectionNoOpKeepsSelection)
{
    auto model{ data_test::MakeModelWithBond() };
    auto & atoms{ model->GetAtomList() };
    atoms.at(0)->SetSpot(Spot::CA);
    atoms.at(1)->SetSpot(Spot::H);
    model->SelectAllAtoms(false);
    model->SetAtomSelected(atoms.at(0)->GetSerialID(), true);

    model->ApplySpotSelection(Spot::H, false);

    ASSERT_EQ(model->GetSelectedAtomCount(), 1u);
    EXPECT_EQ(model->GetSelectedAtoms().front(), atoms.at(0).get());
}

TEST(DataObjectModelAnalysisTest, ModelObjectCanApplyBackboneSelectionAsExclusion)
{
    std::vector<std::unique_ptr<rg::AtomObject>> atom_list;
    const std::vector<Spot> spots{
        Spot::C, Spot::CA, Spot::N, Spot::O, Spot::H, Spot::HA, Spot::CB
    };
    for (size_t i = 0; i < spots.size(); ++i)
    {
        auto atom{ std::make_unique<rg::AtomObject>() };
        atom->SetSerialID(static_cast<int>(i + 1));
        atom->SetSpot(spots.at(i));
        atom_list.emplace_back(std::move(atom));
    }
    rg::ModelObject model{ std::move(atom_list) };
    model.SelectAllAtoms();

    model.ApplyBackboneSelection(true);

    ASSERT_EQ(model.GetSelectedAtomCount(), 6u);
    for (auto * atom : model.GetSelectedAtoms())
    {
        EXPECT_NE(atom->GetSpot(), Spot::CB);
    }
}

TEST(DataObjectModelAnalysisTest, ModelObjectBackboneSelectionDoesNotWidenSelection)
{
    auto model{ data_test::MakeModelWithBond() };
    auto & atoms{ model->GetAtomList() };
    atoms.at(0)->SetSpot(Spot::CA);
    atoms.at(1)->SetSpot(Spot::CB);
    model->SelectAllAtoms(false);
    model->SetAtomSelected(atoms.at(1)->GetSerialID(), true);

    model->ApplyBackboneSelection(true);

    EXPECT_EQ(model->GetSelectedAtomCount(), 0u);
}

TEST(DataObjectModelAnalysisTest, ModelObjectBackboneSelectionNoOpKeepsSelection)
{
    auto model{ data_test::MakeModelWithBond() };
    auto & atoms{ model->GetAtomList() };
    atoms.at(0)->SetSpot(Spot::CA);
    atoms.at(1)->SetSpot(Spot::CB);
    model->SelectAllAtoms(false);
    model->SetAtomSelected(atoms.at(1)->GetSerialID(), true);

    model->ApplyBackboneSelection(false);

    ASSERT_EQ(model->GetSelectedAtomCount(), 1u);
    EXPECT_EQ(model->GetSelectedAtoms().front(), atoms.at(1).get());
}

TEST(DataObjectModelAnalysisTest, ModelObjectCanApplyComponentIDSelectionAsExclusion)
{
    auto model{ data_test::MakeModelWithBond() };
    auto & atoms{ model->GetAtomList() };
    atoms.at(0)->SetComponentID("ALA");
    atoms.at(1)->SetComponentID("HOH");
    model->SelectAllAtoms();

    model->ApplyComponentIDSelection("HOH", true);

    ASSERT_EQ(model->GetSelectedAtomCount(), 1u);
    EXPECT_EQ(model->GetSelectedAtoms().front(), atoms.at(0).get());
}

TEST(DataObjectModelAnalysisTest, ModelObjectComponentIDSelectionExclusionDoesNotWidenSelection)
{
    auto model{ data_test::MakeModelWithBond() };
    auto & atoms{ model->GetAtomList() };
    atoms.at(0)->SetComponentID("ALA");
    atoms.at(1)->SetComponentID("HOH");
    model->SelectAllAtoms(false);
    model->SetAtomSelected(atoms.at(1)->GetSerialID(), true);

    model->ApplyComponentIDSelection("HOH", true);

    EXPECT_EQ(model->GetSelectedAtomCount(), 0u);
}

TEST(DataObjectModelAnalysisTest, ModelObjectComponentIDSelectionNoOpKeepsSelection)
{
    auto model{ data_test::MakeModelWithBond() };
    auto & atoms{ model->GetAtomList() };
    atoms.at(0)->SetComponentID("ALA");
    atoms.at(1)->SetComponentID("HOH");
    model->SelectAllAtoms(false);
    model->SetAtomSelected(atoms.at(0)->GetSerialID(), true);

    model->ApplyComponentIDSelection("HOH", false);

    ASSERT_EQ(model->GetSelectedAtomCount(), 1u);
    EXPECT_EQ(model->GetSelectedAtoms().front(), atoms.at(0).get());
}

TEST(DataObjectModelAnalysisTest, ModelObjectAppliesSimulationMetadata)
{
    auto model{ data_test::MakeModelWithBond() };

    model->ApplySimulationMetadata(1.5);

    EXPECT_EQ(model->GetEmdID(), "Simulation");
    EXPECT_DOUBLE_EQ(model->GetResolution(), 1.5);
    EXPECT_EQ(model->GetResolutionMethod(), "Blurring Width");
}

TEST(DataObjectModelAnalysisTest, ModelObjectInitializesLocalPotentialAnalysis)
{
    auto model{ data_test::MakeModelWithBond() };
    auto * first_atom{ model->GetAtomList().at(0).get() };
    auto * second_atom{ model->GetAtomList().at(1).get() };
    auto & analysis_data{ rg::ModelAnalysisData::Of(*model) };
    model->SelectAllAtoms(false);
    model->SetAtomSelected(first_atom->GetSerialID(), true);
    analysis_data.EnsureAtomLocalEntry(*second_atom).SetAlphaR(
        rg::FittingStage::Third,
        0.7);
    analysis_data.AtomGroupEntry().AddMember(999, *second_atom);

    model->EditAnalysis().InitializeFromSelection();

    size_t member_count{ 0 };
    for (const auto group_key : analysis_data.AtomGroupEntry().CollectGroupKeys())
    {
        member_count += analysis_data.AtomGroupEntry().GetMemberCount(group_key);
        for (const auto * atom : analysis_data.AtomGroupEntry().GetMembers(group_key))
        {
            EXPECT_EQ(atom, first_atom);
        }
    }
    EXPECT_EQ(member_count, 1U);
    EXPECT_DOUBLE_EQ(
        0.0,
        rg::AtomLocalPotentialView::For(*first_atom).GetAlphaR(
            rg::FittingStage::Third));
    EXPECT_EQ(analysis_data.FindAtomLocalEntry(*second_atom), nullptr);

    const auto analysis_view{ model->GetAnalysisView() };
    for (const auto group_key : analysis_view.CollectAtomGroupKeys())
    {
        EXPECT_DOUBLE_EQ(
            0.0,
            analysis_view.GetAtomAlphaG(group_key));
    }
}

TEST(DataObjectModelAnalysisTest, ModelAnalysisEditorCanClearTransientFitStatesWithoutDroppingEntries)
{
    auto model{ data_test::MakeModelWithBond() };
    auto * atom{ model->GetAtomList().at(0).get() };
    auto & analysis_data{ rg::ModelAnalysisData::Of(*model) };
    auto analysis{ model->EditAnalysis() };

    auto & atom_entry{ analysis_data.EnsureAtomLocalEntry(*atom) };
    analysis_data.AtomGroupEntry().AddMember(999, *atom);
    rg::LocalGaussianResult atom_result;
    atom_result.alpha_r = 0.2;
    atom_result.fit_result = rg::RHBMBetaEstimateResult{};
    atom_entry.SetGaussianResult(rg::FittingStage::Third, atom_result);
    atom_entry.SetGroupMemberResult(rg::GroupGaussianMemberResult{ atom_result.mdpde, true, 2.5 });

    ASSERT_NE(analysis_data.FindAtomLocalEntry(*atom), nullptr);
    ASSERT_FALSE(analysis_data.AtomGroupEntry().CollectGroupKeys().empty());

    analysis.ClearTransientFitStates();

    const auto * cleared_atom_entry{ analysis_data.FindAtomLocalEntry(*atom) };
    ASSERT_NE(cleared_atom_entry, nullptr);
    EXPECT_FALSE(analysis_data.AtomGroupEntry().CollectGroupKeys().empty());
    EXPECT_DOUBLE_EQ(
        0.2,
        cleared_atom_entry->GaussianResult(rg::FittingStage::Third).alpha_r);
    ASSERT_TRUE(cleared_atom_entry->GroupMemberResult().has_value());
    EXPECT_TRUE(cleared_atom_entry->GroupMemberResult()->is_outlier);
    EXPECT_DOUBLE_EQ(cleared_atom_entry->GroupMemberResult()->statistical_distance, 2.5);
    EXPECT_FALSE(cleared_atom_entry->GaussianResult(
        rg::FittingStage::Third).fit_result.has_value());
}

TEST(DataObjectModelAnalysisTest, ModelAnalysisDataClearDropsEntriesAndFitStates)
{
    auto model{ data_test::MakeModelWithBond() };
    auto * atom{ model->GetAtomList().at(0).get() };
    auto & analysis_data{ rg::ModelAnalysisData::Of(*model) };

    auto & atom_entry{ analysis_data.EnsureAtomLocalEntry(*atom) };
    analysis_data.AtomGroupEntry().AddMember(999, *atom);
    atom_entry.SetAlphaR(rg::FittingStage::Third, 0.2);
    atom_entry.SetGroupMemberResult(rg::GroupGaussianMemberResult{});

    analysis_data.Clear();

    EXPECT_EQ(analysis_data.FindAtomLocalEntry(*atom), nullptr);
    EXPECT_TRUE(analysis_data.AtomGroupEntry().CollectGroupKeys().empty());
    EXPECT_THROW(rg::AtomLocalPotentialView::For(*atom).GetGroupMemberResult(), std::runtime_error);
    EXPECT_FALSE(analysis_data.EnsureAtomLocalEntry(*atom).GroupMemberResult().has_value());
}

TEST(DataObjectModelAnalysisTest, LocalPotentialEntryClearTransientFitStateKeepsGaussianResult)
{
    rg::LocalPotentialEntry entry;
    rg::LocalGaussianResult result;
    result.alpha_r = 0.4;
    result.mdpde = rg::GaussianModel3DWithUncertainty{
        rg::GaussianModel3D{ 2.0, 0.7 },
        rg::GaussianModel3DUncertainty{}
    };
    result.fit_result = rg::RHBMBetaEstimateResult{};
    entry.SetGaussianResult(rg::FittingStage::Third, result);
    entry.SetPeelingSamplingEntries({
        LocalPotentialSample{ 3.0, SamplingPoint{ 0.1 } }
    });
    ASSERT_TRUE(entry.GaussianResult(
        rg::FittingStage::Third).fit_result.has_value());
    ASSERT_FALSE(entry.PeelingSamplingEntries().empty());

    entry.ClearTransientFitState(rg::FittingStage::Third);

    EXPECT_DOUBLE_EQ(
        0.4,
        entry.GaussianResult(rg::FittingStage::Third).alpha_r);
    EXPECT_DOUBLE_EQ(
        2.0,
        entry.GaussianResult(rg::FittingStage::Third)
            .mdpde.GetModel().GetAmplitude());
    EXPECT_DOUBLE_EQ(
        0.7,
        entry.GaussianResult(rg::FittingStage::Third)
            .mdpde.GetModel().GetWidth());
    EXPECT_FALSE(entry.GaussianResult(
        rg::FittingStage::Third).fit_result.has_value());
    ASSERT_EQ(entry.PeelingSamplingEntries().size(), 1u);
    EXPECT_DOUBLE_EQ(entry.PeelingSamplingEntries().front().response, 3.0);
}

TEST(DataObjectModelAnalysisTest, LocalPotentialEntryStoresRawAndPeelingSamplingEntries)
{
    rg::LocalPotentialEntry entry;
    EXPECT_EQ(entry.NeighborCountForPeeling(), 0);
    entry.SetRawSamplingEntries({
        LocalPotentialSample{ 6.0, SamplingPoint{ 0.1 } },
        LocalPotentialSample{ 4.0, SamplingPoint{ 0.2 } }
    });
    entry.SetPeelingSamplingEntries({
        LocalPotentialSample{ 3.0, SamplingPoint{ 0.1 } }
    });
    entry.SetNeighborCountForPeeling(7);

    EXPECT_EQ(entry.RawSamplingEntryCount(), 2);
    EXPECT_EQ(entry.PeelingSamplingEntryCount(), 1);
    ASSERT_EQ(entry.RawSamplingEntries().size(), 2u);
    ASSERT_EQ(entry.PeelingSamplingEntries().size(), 1u);
    EXPECT_DOUBLE_EQ(entry.RawSamplingEntries().front().response, 6.0);
    EXPECT_DOUBLE_EQ(entry.PeelingSamplingEntries().front().response, 3.0);
    EXPECT_EQ(entry.NeighborCountForPeeling(), 7);
}

TEST(DataObjectModelAnalysisTest, LocalPotentialEntryStoresGaussianResult)
{
    rg::LocalPotentialEntry entry;
    rg::LocalGaussianResult result;
    result.alpha_r = 0.5;
    result.ols = rg::GaussianModel3DWithUncertainty{
        rg::GaussianModel3D{ 1.0, 0.6 },
        rg::GaussianModel3DUncertainty{}
    };
    result.mdpde = rg::GaussianModel3DWithUncertainty{
        rg::GaussianModel3D{ 1.5, 0.8 },
        rg::GaussianModel3DUncertainty{}
    };

    entry.SetGaussianResult(rg::FittingStage::Third, result);

    EXPECT_DOUBLE_EQ(
        0.5,
        entry.GaussianResult(rg::FittingStage::Third).alpha_r);
    EXPECT_DOUBLE_EQ(
        1.0,
        entry.GaussianResult(rg::FittingStage::Third)
            .ols.GetModel().GetAmplitude());
    EXPECT_DOUBLE_EQ(
        1.5,
        entry.GaussianResult(rg::FittingStage::Third)
            .mdpde.GetModel().GetAmplitude());
}

TEST(DataObjectModelAnalysisTest, LocalPotentialEntryKeepsGaussianStagesIndependent)
{
    rg::LocalPotentialEntry entry;
    for (const auto [stage, amplitude] : {
             std::pair{ rg::FittingStage::First, 1.0 },
             std::pair{ rg::FittingStage::Second, 2.0 },
             std::pair{ rg::FittingStage::Third, 3.0 } })
    {
        rg::LocalGaussianResult result;
        result.alpha_r = amplitude / 10.0;
        result.mdpde = rg::GaussianModel3DWithUncertainty{
            rg::GaussianModel3D{ amplitude, 0.5, -amplitude },
            rg::GaussianModel3DUncertainty{}
        };
        entry.SetGaussianResult(stage, result);
    }

    EXPECT_DOUBLE_EQ(
        entry.GaussianResult(rg::FittingStage::First)
            .mdpde.GetModel().GetAmplitude(),
        1.0);
    EXPECT_DOUBLE_EQ(
        entry.GaussianResult(rg::FittingStage::Second)
            .mdpde.GetModel().GetAmplitude(),
        2.0);
    EXPECT_DOUBLE_EQ(
        entry.GaussianResult(rg::FittingStage::Third)
            .mdpde.GetModel().GetAmplitude(),
        3.0);
    entry.SetAlphaR(rg::FittingStage::Second, 0.75);
    EXPECT_DOUBLE_EQ(
        entry.GaussianResult(rg::FittingStage::First).alpha_r,
        0.1);
    EXPECT_DOUBLE_EQ(
        entry.GaussianResult(rg::FittingStage::Second).alpha_r,
        0.75);
    EXPECT_DOUBLE_EQ(
        entry.GaussianResult(rg::FittingStage::Third).alpha_r,
        0.3);
}

TEST(DataObjectModelAnalysisTest, GroupPotentialEntryKeepsSingleResultPerGroup)
{
    auto model{ data_test::MakeModelWithBond() };
    auto * first_atom{ model->GetAtomList().at(0).get() };
    auto * second_atom{ model->GetAtomList().at(1).get() };
    rg::AtomGroupPotentialEntry entry;
    constexpr GroupKey group_key{ 42 };

    EXPECT_FALSE(entry.HasGroup(group_key));
    EXPECT_THROW(entry.GetPrior(group_key), std::runtime_error);
    entry.AddMember(group_key, *first_atom);
    rg::GroupGaussianResult result;
    result.alpha_g = 0.1;
    result.mean = rg::GaussianModel3D{ 1.0, 0.5, 0.1 };
    result.mdpde = rg::GaussianModel3D{ 1.1, 0.6, 0.2 };
    result.prior = rg::GaussianModel3DWithUncertainty{
        rg::GaussianModel3D{ 1.2, 0.7, 0.3 },
        rg::GaussianModel3DUncertainty{ 0.01, 0.02, 0.03 }
    };
    entry.SetGaussianResult(group_key, result);
    EXPECT_DOUBLE_EQ(entry.GetPrior(group_key).GetAmplitude(), 1.2);
    entry.AddMember(group_key, *second_atom);
    entry.SetAlphaG(group_key, 0.8);
    EXPECT_DOUBLE_EQ(entry.GetAlphaG(group_key), 0.8);

    result.mean = rg::GaussianModel3D{ 2.0, 0.6, 0.2 };
    result.mdpde = rg::GaussianModel3D{ 2.1, 0.7, 0.3 };
    result.prior = rg::GaussianModel3DWithUncertainty{
        rg::GaussianModel3D{ 2.2, 0.8, 0.4 },
        rg::GaussianModel3DUncertainty{ 0.04, 0.05, 0.06 }
    };
    result.alpha_g = 0.3;
    entry.SetGaussianResult(group_key, result);

    EXPECT_EQ(entry.GroupCount(), 1u);
    EXPECT_EQ(entry.CollectGroupKeys(), std::vector<GroupKey>{ group_key });
    ASSERT_EQ(entry.GetMemberCount(group_key), 2u);
    EXPECT_EQ(entry.GetMembers(group_key).at(0), first_atom);
    EXPECT_EQ(entry.GetMembers(group_key).at(1), second_atom);
    EXPECT_DOUBLE_EQ(entry.GetMean(group_key).GetAmplitude(), 2.0);
    EXPECT_DOUBLE_EQ(entry.GetMDPDE(group_key).GetAmplitude(), 2.1);
    EXPECT_DOUBLE_EQ(entry.GetPrior(group_key).GetAmplitude(), 2.2);
    EXPECT_DOUBLE_EQ(entry.GetPriorStandardDeviation(group_key).GetOffset(), 0.06);
    EXPECT_DOUBLE_EQ(entry.GetAlphaG(group_key), 0.3);
}

TEST(DataObjectModelAnalysisTest, ModelCopyPreservesGroupResultAndRebindsMembers)
{
    auto model{ data_test::MakeModelWithBond() };
    model->SelectAllAtoms();
    model->ApplySymmetrySelection(false);
    auto analysis{ model->EditAnalysis() };
    analysis.RebuildAtomGroupsFromSelection();
    const auto source_view{ model->GetAnalysisView() };
    const auto group_key{ source_view.CollectAtomGroupKeys().front() };
    const auto & source_members{ source_view.GetAtomObjectList(group_key) };

    rg::GroupGaussianResult result;
    result.alpha_g = 0.3;
    result.mean = rg::GaussianModel3D{ 3.0, 0.5, 0.1 };
    result.mdpde = rg::GaussianModel3D{ 3.1, 0.6, 0.2 };
    result.prior = rg::GaussianModel3DWithUncertainty{
        rg::GaussianModel3D{ 3.2, 0.7, 0.3 },
        rg::GaussianModel3DUncertainty{ 0.1, 0.2, 0.3 }
    };
    for (const auto * atom : source_members)
    {
        result.member_results.emplace_back(rg::GroupGaussianMemberResult{
            result.prior, true, static_cast<double>(atom->GetSerialID())
        });
    }
    analysis.ApplyAtomGroupGaussianResult(group_key, result);

    rg::ModelObject copied_model{ *model };
    const auto copied_view{ copied_model.GetAnalysisView() };
    EXPECT_EQ(copied_view.CollectAtomGroupKeys().size(), source_view.CollectAtomGroupKeys().size());
    EXPECT_DOUBLE_EQ(copied_view.GetAtomGroupMean(group_key).GetAmplitude(), 3.0);
    EXPECT_DOUBLE_EQ(copied_view.GetAtomGroupMDPDE(group_key).GetAmplitude(), 3.1);
    EXPECT_DOUBLE_EQ(copied_view.GetAtomGroupPrior(group_key).GetAmplitude(), 3.2);
    EXPECT_DOUBLE_EQ(copied_view.GetAtomAlphaG(group_key), 0.3);
    EXPECT_DOUBLE_EQ(copied_view.GetAtomGroupPriorWithUncertainty(group_key)
        .GetStandardDeviationModel().GetOffset(), 0.3);
    const auto & copied_members{ copied_view.GetAtomObjectList(group_key) };
    ASSERT_EQ(copied_members.size(), source_members.size());
    for (std::size_t i = 0; i < source_members.size(); ++i)
    {
        EXPECT_NE(copied_members.at(i), source_members.at(i));
        EXPECT_EQ(copied_members.at(i)->GetSerialID(), source_members.at(i)->GetSerialID());
        const auto & member{ rg::AtomLocalPotentialView::For(*copied_members.at(i)).GetGroupMemberResult() };
        ASSERT_TRUE(member.has_value());
        EXPECT_TRUE(member->is_outlier);
        EXPECT_DOUBLE_EQ(member->statistical_distance, result.member_results.at(i).statistical_distance);
        EXPECT_DOUBLE_EQ(member->posterior.GetModel().GetOffset(), 0.3);
    }

    result.alpha_g = 0.8;
    result.member_results.front().is_outlier = false;
    copied_model.EditAnalysis().ApplyAtomGroupGaussianResult(group_key, result);
    EXPECT_DOUBLE_EQ(source_view.GetAtomAlphaG(group_key), 0.3);
    EXPECT_TRUE(rg::AtomLocalPotentialView::For(*source_members.front())
        .GetGroupMemberResult()->is_outlier);
}

TEST(DataObjectModelAnalysisTest, ModelAnalysisEditorCanSetAlphaRAndCreateEntry)
{
    auto model{ data_test::MakeModelWithBond() };
    auto * atom{ model->GetAtomList().at(0).get() };
    auto analysis{ model->EditAnalysis() };

    analysis.SetAtomLocalAlphaR(rg::FittingStage::Third, *atom, 0.37);

    EXPECT_DOUBLE_EQ(
        0.37,
        rg::AtomLocalPotentialView::For(*atom).GetAlphaR(
            rg::FittingStage::Third));
}

TEST(DataObjectModelAnalysisTest, ModelAnalysisEditorCanSetPeelingNeighborCount)
{
    auto model{ data_test::MakeModelWithBond() };
    auto * atom{ model->GetAtomList().at(0).get() };

    model->EditAnalysis().SetAtomLocalNeighborCountForPeeling(*atom, 5);

    EXPECT_EQ(
        rg::AtomLocalPotentialView::For(*atom)
            .GetNeighborCountForPeeling(),
        5);
}

TEST(DataObjectModelAnalysisTest, ModelAnalysisEditorSetsPeelingNeighborCountAndCreatesEntry)
{
    auto model{ data_test::MakeModelWithBond() };
    auto * atom{ model->GetAtomList().at(0).get() };
    auto analysis{ model->EditAnalysis() };

    EXPECT_FALSE(rg::AtomLocalPotentialView::For(*atom).IsAvailable());

    analysis.SetAtomLocalNeighborCountForPeeling(*atom, 5);

    const auto view{ rg::AtomLocalPotentialView::For(*atom) };
    EXPECT_EQ(view.GetNeighborCountForPeeling(), 5);
}

TEST(DataObjectModelAnalysisTest, ModelAnalysisEditorAppliesAtomLocalSecondStageResult)
{
    auto model{ data_test::MakeModelWithBond() };
    auto * atom{ model->GetAtomList().at(0).get() };
    auto analysis{ model->EditAnalysis() };

    const LocalPotentialSampleList raw_sampling_entries{
        LocalPotentialSample{ 8.0, SamplingPoint{ 0.2 } }
    };
    const LocalPotentialSampleList peeling_sampling_entries{
        LocalPotentialSample{ 3.0, SamplingPoint{ 0.4 } }
    };
    analysis.SetAtomLocalRawSamplingEntries(*atom, raw_sampling_entries);

    rg::LocalGaussianResult first_result;
    first_result.mdpde = rg::GaussianModel3DWithUncertainty{
        rg::GaussianModel3D{ 1.0, 0.6 },
        rg::GaussianModel3DUncertainty{}
    };
    analysis.SetAtomLocalGaussianResult(
        rg::FittingStage::First, *atom, first_result);

    rg::LocalGaussianResult third_result;
    third_result.mdpde = rg::GaussianModel3DWithUncertainty{
        rg::GaussianModel3D{ 3.0, 0.8 },
        rg::GaussianModel3DUncertainty{}
    };
    analysis.SetAtomLocalGaussianResult(
        rg::FittingStage::Third, *atom, third_result);

    rg::LocalGaussianResult second_result;
    second_result.alpha_r = 0.42;
    second_result.mdpde = rg::GaussianModel3DWithUncertainty{
        rg::GaussianModel3D{ 2.0, 0.7, -0.1 },
        rg::GaussianModel3DUncertainty{}
    };
    analysis.ApplyAtomLocalSecondStageResult(
        *atom,
        second_result,
        peeling_sampling_entries);

    const auto view{ rg::AtomLocalPotentialView::For(*atom) };
    EXPECT_EQ(view.GetRawSamplingEntries(false).size(), 1u);
    EXPECT_DOUBLE_EQ(view.GetRawSamplingEntries(false).front().response, 8.0);
    EXPECT_EQ(view.GetPeelingSamplingEntries(false).size(), 1u);
    EXPECT_DOUBLE_EQ(view.GetPeelingSamplingEntries(false).front().response, 3.0);
    EXPECT_DOUBLE_EQ(
        view.GetGaussianResult(rg::FittingStage::Second).alpha_r,
        0.42);
    EXPECT_DOUBLE_EQ(
        view.GetEstimateMDPDE(rg::FittingStage::Second).GetAmplitude(),
        2.0);
    EXPECT_DOUBLE_EQ(
        view.GetEstimateMDPDE(rg::FittingStage::First).GetAmplitude(),
        1.0);
    EXPECT_DOUBLE_EQ(
        view.GetEstimateMDPDE(rg::FittingStage::Third).GetAmplitude(),
        3.0);
}

TEST(DataObjectModelAnalysisTest, ModelAnalysisEditorSetGaussianResultUpdatesViewEstimates)
{
    auto model{ data_test::MakeModelWithBond() };
    auto * atom{ model->GetAtomList().at(0).get() };
    auto analysis{ model->EditAnalysis() };

    LocalPotentialSampleList sampling_entries{
        LocalPotentialSample{ 6.0, SamplingPoint{ 0.0 } },
        LocalPotentialSample{ 4.0, SamplingPoint{ 0.5 } },
        LocalPotentialSample{ 2.0, SamplingPoint{ 0.9 } }
    };
    analysis.SetAtomLocalRawSamplingEntries(*atom, sampling_entries);

    rg::LocalGaussianResult gaussian_result;
    gaussian_result.alpha_r = 0.6;
    gaussian_result.ols = rg::GaussianModel3DWithUncertainty{
        rg::GaussianModel3D{ 0.0, 0.0 },
        rg::GaussianModel3DUncertainty{}
    };
    gaussian_result.mdpde = rg::GaussianModel3DWithUncertainty{
        rg::GaussianModel3D{ 0.0, 0.0 },
        rg::GaussianModel3DUncertainty{}
    };

    analysis.SetAtomLocalGaussianResult(
        rg::FittingStage::Third, *atom, gaussian_result);

    EXPECT_DOUBLE_EQ(
        0.6,
        rg::AtomLocalPotentialView::For(*atom).GetGaussianResult(
            rg::FittingStage::Third).alpha_r);
    EXPECT_DOUBLE_EQ(
        0.0,
        rg::AtomLocalPotentialView::For(*atom).GetEstimateOLS(
            rg::FittingStage::Third).GetWidth());
    EXPECT_DOUBLE_EQ(
        0.0,
        rg::AtomLocalPotentialView::For(*atom).GetEstimateMDPDE(
            rg::FittingStage::Third).GetWidth());
}

TEST(DataObjectModelAnalysisTest, AtomLocalPotentialViewCanApplyRawSamplingSelection)
{
    auto model{ data_test::MakeModelWithBond() };
    auto * atom{ model->GetAtomList().at(0).get() };
    auto analysis{ model->EditAnalysis() };

    analysis.SetAtomLocalRawSamplingEntries(*atom, {
        LocalPotentialSample{ 6.0, SamplingPoint{ 0.0, { 0.0, 0.0, 0.0 }, true } },
        LocalPotentialSample{ 4.0, SamplingPoint{ 0.5, { 0.0, 0.0, 0.0 }, false } },
        LocalPotentialSample{ 2.0, SamplingPoint{ 0.9, { 0.0, 0.0, 0.0 }, true } }
    });

    const auto view{ rg::AtomLocalPotentialView::For(*atom) };
    const auto raw_selected_entries{ view.GetRawSamplingEntries() };
    const auto raw_all_entries{ view.GetRawSamplingEntries(false) };

    ASSERT_EQ(raw_selected_entries.size(), 2u);
    EXPECT_DOUBLE_EQ(raw_selected_entries.at(0).response, 6.0);
    EXPECT_DOUBLE_EQ(raw_selected_entries.at(1).response, 2.0);
    ASSERT_EQ(raw_all_entries.size(), 3u);
    EXPECT_FALSE(raw_all_entries.at(1).point.is_selected);
}

TEST(DataObjectModelAnalysisTest, AtomLocalPotentialViewCanApplyPeelingSamplingSelection)
{
    auto model{ data_test::MakeModelWithBond() };
    auto * atom{ model->GetAtomList().at(0).get() };
    auto analysis{ model->EditAnalysis() };

    analysis.SetAtomLocalRawSamplingEntries(*atom, {
        LocalPotentialSample{ 6.0, SamplingPoint{ 0.0, { 0.0, 0.0, 0.0 }, true } },
        LocalPotentialSample{ 4.0, SamplingPoint{ 0.5, { 0.0, 0.0, 0.0 }, false } },
        LocalPotentialSample{ 2.0, SamplingPoint{ 0.9, { 0.0, 0.0, 0.0 }, true } }
    });
    analysis.SetAtomLocalPeelingSamplingEntries(*atom, {
        LocalPotentialSample{ 3.0, SamplingPoint{ 0.0, { 0.0, 0.0, 0.0 }, true } },
        LocalPotentialSample{ 5.0, SamplingPoint{ 0.5, { 0.0, 0.0, 0.0 }, false } },
        LocalPotentialSample{ 7.0, SamplingPoint{ 0.9, { 0.0, 0.0, 0.0 }, true } }
    });

    const auto view{ rg::AtomLocalPotentialView::For(*atom) };
    const auto raw_selected_entries{ view.GetRawSamplingEntries() };
    const auto peeling_selected_entries{ view.GetPeelingSamplingEntries(true) };
    const auto peeling_all_entries{ view.GetPeelingSamplingEntries() };

    ASSERT_EQ(raw_selected_entries.size(), 2u);
    EXPECT_DOUBLE_EQ(raw_selected_entries.at(0).response, 6.0);
    EXPECT_DOUBLE_EQ(raw_selected_entries.at(1).response, 2.0);
    ASSERT_EQ(peeling_selected_entries.size(), 2u);
    EXPECT_DOUBLE_EQ(peeling_selected_entries.at(0).response, 3.0);
    EXPECT_DOUBLE_EQ(peeling_selected_entries.at(1).response, 7.0);
    ASSERT_EQ(peeling_all_entries.size(), 3u);
    EXPECT_DOUBLE_EQ(peeling_all_entries.at(1).response, 5.0);
    EXPECT_FALSE(peeling_all_entries.at(1).point.is_selected);
}

TEST(
    DataObjectModelAnalysisTest,
    AtomLocalPotentialViewComputesPeelingRatioInInclusiveDistanceRange)
{
    auto model{ data_test::MakeModelWithBond() };
    auto * atom{ model->GetAtomList().at(0).get() };
    auto analysis{ model->EditAnalysis() };
    const auto non_finite_distance{ std::numeric_limits<double>::quiet_NaN() };

    analysis.SetAtomLocalRawSamplingEntries(*atom, {
        LocalPotentialSample{ 10.0, SamplingPoint{ 0.5 } },
        LocalPotentialSample{ 8.0, SamplingPoint{ 1.0 } },
        LocalPotentialSample{
            4.0,
            SamplingPoint{ 1.5, { 0.0, 0.0, 0.0 }, false } },
        LocalPotentialSample{ 8.0, SamplingPoint{ 2.0 } },
        LocalPotentialSample{ 100.0, SamplingPoint{ 2.5 } },
        LocalPotentialSample{ 100.0, SamplingPoint{ non_finite_distance } }
    });
    analysis.SetAtomLocalPeelingSamplingEntries(*atom, {
        LocalPotentialSample{ 5.0, SamplingPoint{ 0.5 } },
        LocalPotentialSample{ 2.0, SamplingPoint{ 1.0 } },
        LocalPotentialSample{
            1.0,
            SamplingPoint{ 1.5, { 0.0, 0.0, 0.0 }, false } },
        LocalPotentialSample{ 3.0, SamplingPoint{ 2.0 } },
        LocalPotentialSample{ 50.0, SamplingPoint{ 2.5 } },
        LocalPotentialSample{ 50.0, SamplingPoint{ non_finite_distance } }
    });

    const auto view{ rg::AtomLocalPotentialView::For(*atom) };
    const auto ratio{ view.GetLocalFittingPeelingRatio(true, 1.0, 2.0) };
    ASSERT_TRUE(ratio.has_value());
    EXPECT_DOUBLE_EQ(*ratio, 0.7);

    const auto boundary_ratio{
        view.GetLocalFittingPeelingRatio(true, 1.0, 1.0)
    };
    ASSERT_TRUE(boundary_ratio.has_value());
    EXPECT_DOUBLE_EQ(*boundary_ratio, 0.75);
}

TEST(
    DataObjectModelAnalysisTest,
    AtomLocalPotentialViewReturnsNulloptWhenPeelingRatioCannotBeComputed)
{
    auto model{ data_test::MakeModelWithBond() };
    auto * atom{ model->GetAtomList().at(0).get() };
    auto analysis{ model->EditAnalysis() };
    const auto view{ rg::AtomLocalPotentialView::For(*atom) };

    analysis.SetAtomLocalRawSamplingEntries(*atom, {
        LocalPotentialSample{ 4.0, SamplingPoint{ 0.5 } }
    });
    analysis.SetAtomLocalPeelingSamplingEntries(*atom, {
        LocalPotentialSample{ 2.0, SamplingPoint{ 0.5 } }
    });
    EXPECT_FALSE(view.GetLocalFittingPeelingRatio(false, 0.0, 1.0).has_value());

    analysis.SetAtomLocalRawSamplingEntries(*atom, {
        LocalPotentialSample{ 4.0, SamplingPoint{ 1.5 } }
    });
    EXPECT_FALSE(view.GetLocalFittingPeelingRatio(true, 0.0, 1.0).has_value());

    analysis.SetAtomLocalRawSamplingEntries(*atom, {
        LocalPotentialSample{ 4.0, SamplingPoint{ 0.5 } }
    });
    analysis.SetAtomLocalPeelingSamplingEntries(*atom, {
        LocalPotentialSample{ 2.0, SamplingPoint{ 1.5 } }
    });
    EXPECT_FALSE(view.GetLocalFittingPeelingRatio(true, 0.0, 1.0).has_value());

    analysis.SetAtomLocalRawSamplingEntries(*atom, {
        LocalPotentialSample{ 1.0, SamplingPoint{ 0.5 } },
        LocalPotentialSample{ -1.0, SamplingPoint{ 0.75 } }
    });
    analysis.SetAtomLocalPeelingSamplingEntries(*atom, {
        LocalPotentialSample{ 1.0, SamplingPoint{ 0.5 } }
    });
    EXPECT_FALSE(view.GetLocalFittingPeelingRatio(true, 0.0, 1.0).has_value());

    analysis.SetAtomLocalRawSamplingEntries(*atom, {
        LocalPotentialSample{
            std::numeric_limits<double>::infinity(),
            SamplingPoint{ 0.5 } }
    });
    EXPECT_FALSE(view.GetLocalFittingPeelingRatio(true, 0.0, 1.0).has_value());

    analysis.SetAtomLocalRawSamplingEntries(*atom, {
        LocalPotentialSample{ 1.0, SamplingPoint{ 0.5 } }
    });
    analysis.SetAtomLocalPeelingSamplingEntries(*atom, {
        LocalPotentialSample{
            std::numeric_limits<double>::infinity(),
            SamplingPoint{ 0.5 } }
    });
    EXPECT_FALSE(view.GetLocalFittingPeelingRatio(true, 0.0, 1.0).has_value());

    analysis.SetAtomLocalRawSamplingEntries(*atom, {
        LocalPotentialSample{
            std::numeric_limits<double>::min(),
            SamplingPoint{ 0.5 } }
    });
    analysis.SetAtomLocalPeelingSamplingEntries(*atom, {
        LocalPotentialSample{
            std::numeric_limits<double>::max(),
            SamplingPoint{ 0.5 } }
    });
    EXPECT_FALSE(view.GetLocalFittingPeelingRatio(true, 0.0, 1.0).has_value());
}

TEST(
    DataObjectModelAnalysisTest,
    AtomLocalPotentialViewRejectsInvalidPeelingRatioDistanceRange)
{
    auto model{ data_test::MakeModelWithBond() };
    auto * atom{ model->GetAtomList().at(0).get() };
    const auto view{ rg::AtomLocalPotentialView::For(*atom) };

    EXPECT_THROW(
        view.GetLocalFittingPeelingRatio(false, -0.1, 1.0),
        std::invalid_argument);
    EXPECT_THROW(
        view.GetLocalFittingPeelingRatio(false, 2.0, 1.0),
        std::invalid_argument);
    EXPECT_THROW(
        view.GetLocalFittingPeelingRatio(
            false,
            std::numeric_limits<double>::quiet_NaN(),
            1.0),
        std::invalid_argument);
    EXPECT_THROW(
        view.GetLocalFittingPeelingRatio(
            false,
            0.0,
            std::numeric_limits<double>::infinity()),
        std::invalid_argument);
}

TEST(DataObjectModelAnalysisTest, AtomLocalPotentialViewGetsSamplingEntriesByFittingStage)
{
    auto model{ data_test::MakeModelWithBond() };
    auto * atom{ model->GetAtomList().at(0).get() };
    auto analysis{ model->EditAnalysis() };

    analysis.SetAtomLocalRawSamplingEntries(*atom, {
        LocalPotentialSample{ 6.0, SamplingPoint{ 0.0, { 0.0, 0.0, 0.0 }, true } },
        LocalPotentialSample{ 4.0, SamplingPoint{ 0.5, { 0.0, 0.0, 0.0 }, false } }
    });
    analysis.SetAtomLocalPeelingSamplingEntries(*atom, {
        LocalPotentialSample{ 3.0, SamplingPoint{ 0.0, { 0.0, 0.0, 0.0 }, true } },
        LocalPotentialSample{ 5.0, SamplingPoint{ 0.5, { 0.0, 0.0, 0.0 }, false } }
    });

    const auto view{ rg::AtomLocalPotentialView::For(*atom) };
    const auto first_entries{ view.GetSamplingEntries(rg::FittingStage::First) };
    const auto second_entries{ view.GetSamplingEntries(rg::FittingStage::Second) };
    const auto third_entries{ view.GetSamplingEntries(rg::FittingStage::Third) };

    ASSERT_EQ(first_entries.size(), 1u);
    EXPECT_DOUBLE_EQ(first_entries.front().response, 6.0);
    ASSERT_EQ(second_entries.size(), 2u);
    EXPECT_DOUBLE_EQ(second_entries.at(0).response, 3.0);
    EXPECT_DOUBLE_EQ(second_entries.at(1).response, 5.0);
    ASSERT_EQ(third_entries.size(), 2u);
    EXPECT_DOUBLE_EQ(third_entries.at(0).response, 3.0);
    EXPECT_DOUBLE_EQ(third_entries.at(1).response, 5.0);
    EXPECT_THROW(
        view.GetSamplingEntries(static_cast<rg::FittingStage>(3)),
        std::invalid_argument);
}

TEST(DataObjectModelAnalysisTest, EmptyPeelingSamplingEntriesDoNotFallBackToRaw)
{
    auto model{ data_test::MakeModelWithBond() };
    auto * atom{ model->GetAtomList().at(0).get() };
    model->EditAnalysis().SetAtomLocalRawSamplingEntries(*atom, {
        LocalPotentialSample{ 6.0, SamplingPoint{ 0.0 } }
    });

    const auto view{ rg::AtomLocalPotentialView::For(*atom) };
    EXPECT_EQ(view.GetRawSamplingEntries(false).size(), 1u);
    EXPECT_TRUE(view.GetPeelingSamplingEntries(false).empty());
}

TEST(DataObjectModelAnalysisTest, ModelAnalysisEditorAppliesAtomGroupGaussianResultToStatisticsAndPosterior)
{
    auto model{ data_test::MakeModelWithBond() };
    model->SelectAllAtoms();
    model->ApplySymmetrySelection(false);
    auto analysis{ model->EditAnalysis() };
    analysis.RebuildAtomGroupsFromSelection();
    const auto analysis_view{ model->GetAnalysisView() };
    const auto group_keys{ analysis_view.CollectAtomGroupKeys() };
    ASSERT_FALSE(group_keys.empty());

    const auto group_key{ group_keys.front() };
    const auto & atom_list{ analysis_view.GetAtomObjectList(group_key) };
    ASSERT_FALSE(atom_list.empty());

    rg::LocalGaussianResult first_atom_local_result;
    first_atom_local_result.mdpde = rg::GaussianModel3DWithUncertainty{
        rg::GaussianModel3D{ 9.0, 0.5 },
        rg::GaussianModel3DUncertainty{}
    };
    analysis.SetAtomLocalGaussianResult(
        rg::FittingStage::Third,
        *atom_list.front(),
        first_atom_local_result);

    constexpr double alpha_g{ 0.25 };
    rg::GroupGaussianResult result;
    result.alpha_g = alpha_g;
    result.mean = rg::GaussianModel3D{ 1.1, 1.2 };
    result.mdpde = rg::GaussianModel3D{ 1.2, 1.1 };
    result.prior = rg::GaussianModel3DWithUncertainty{
        rg::GaussianModel3D{ 1.3, 1.0 },
        rg::GaussianModel3DUncertainty{ 0.1, 0.2 }
    };
    result.member_results.reserve(atom_list.size());
    for (std::size_t i = 0; i < atom_list.size(); i++)
    {
        rg::GroupGaussianMemberResult member_result;
        member_result.posterior = rg::GaussianModel3DWithUncertainty{
            rg::GaussianModel3D{ 0.4 + 0.1 * static_cast<double>(i), 0.9 },
            rg::GaussianModel3DUncertainty{ 0.01, 0.02 }
        };
        member_result.is_outlier = (i == 0);
        member_result.statistical_distance = 1.5 + static_cast<double>(i);
        result.member_results.emplace_back(member_result);
    }

    analysis.ApplyAtomGroupGaussianResult(group_key,
        result);

    EXPECT_NEAR(result.mean.GetAmplitude(), analysis_view.GetAtomGroupMean(group_key).GetAmplitude(), 1e-12);
    EXPECT_NEAR(result.mean.GetWidth(), analysis_view.GetAtomGroupMean(group_key).GetWidth(), 1e-12);
    EXPECT_NEAR(result.mdpde.GetAmplitude(), analysis_view.GetAtomGroupMDPDE(group_key).GetAmplitude(), 1e-12);
    EXPECT_NEAR(result.mdpde.GetWidth(), analysis_view.GetAtomGroupMDPDE(group_key).GetWidth(), 1e-12);
    EXPECT_NEAR(result.prior.GetModel().GetAmplitude(), analysis_view.GetAtomGroupPrior(group_key).GetAmplitude(), 1e-12);
    EXPECT_NEAR(result.prior.GetModel().GetWidth(), analysis_view.GetAtomGroupPrior(group_key).GetWidth(), 1e-12);
    EXPECT_NEAR(result.prior.GetStandardDeviationModel().GetAmplitude(),
        analysis_view.GetAtomGroupPriorWithUncertainty(group_key)
            .GetStandardDeviationModel().GetAmplitude(), 1e-12);
    EXPECT_NEAR(result.prior.GetStandardDeviationModel().GetWidth(),
        analysis_view.GetAtomGroupPriorWithUncertainty(group_key)
            .GetStandardDeviationModel().GetWidth(), 1e-12);
    EXPECT_DOUBLE_EQ(alpha_g, analysis_view.GetAtomAlphaG(group_key));

    const auto & gaussian_result{
        rg::AtomLocalPotentialView::For(*atom_list.front())
            .GetGaussianResult(rg::FittingStage::Third)
    };
    EXPECT_NEAR(
        first_atom_local_result.mdpde.GetModel().GetAmplitude(),
        gaussian_result.mdpde.GetModel().GetAmplitude(),
        1e-12);
    EXPECT_NEAR(
        first_atom_local_result.mdpde.GetModel().GetWidth(),
        gaussian_result.mdpde.GetModel().GetWidth(),
        1e-12);
    const auto & member_result{
        rg::AtomLocalPotentialView::For(*atom_list.front()).GetGroupMemberResult()
    };
    ASSERT_TRUE(member_result.has_value());
    const auto expected_gaussian{ result.member_results.front().posterior };
    EXPECT_NEAR(
        expected_gaussian.GetModel().GetAmplitude(),
        member_result->posterior.GetModel().GetAmplitude(),
        1e-12);
    EXPECT_NEAR(
        expected_gaussian.GetModel().GetWidth(),
        member_result->posterior.GetModel().GetWidth(),
        1e-12);
    EXPECT_NEAR(
        expected_gaussian.GetStandardDeviationModel().GetAmplitude(),
        member_result->posterior.GetStandardDeviationModel().GetAmplitude(),
        1e-12);
    EXPECT_NEAR(
        expected_gaussian.GetStandardDeviationModel().GetWidth(),
        member_result->posterior.GetStandardDeviationModel().GetWidth(),
        1e-12);
    EXPECT_TRUE(member_result->is_outlier);
    EXPECT_DOUBLE_EQ(1.5, member_result->statistical_distance);
}

TEST(DataObjectModelAnalysisTest, ModelAnalysisEditorRejectsAtomGroupGaussianResultWithMismatchedMemberCount)
{
    auto model{ data_test::MakeModelWithBond() };
    model->SelectAllAtoms();
    model->ApplySymmetrySelection(false);
    auto analysis{ model->EditAnalysis() };
    analysis.RebuildAtomGroupsFromSelection();
    const auto analysis_view{ model->GetAnalysisView() };
    const auto group_key{ analysis_view.CollectAtomGroupKeys().front() };

    rg::GroupGaussianResult result;
    result.alpha_g = 0.4;
    result.member_results.assign(analysis_view.GetAtomObjectList(group_key).size(),
        rg::GroupGaussianMemberResult{ rg::GaussianModel3DWithUncertainty{}, true, 2.5 });
    analysis.ApplyAtomGroupGaussianResult(group_key, result);
    result.alpha_g = 0.9;
    result.member_results.clear();

    EXPECT_THROW(
        analysis.ApplyAtomGroupGaussianResult(group_key,
            result),
        std::invalid_argument);
    EXPECT_DOUBLE_EQ(analysis_view.GetAtomAlphaG(group_key), 0.4);
    for (const auto * atom : analysis_view.GetAtomObjectList(group_key))
    {
        const auto & member{ rg::AtomLocalPotentialView::For(*atom).GetGroupMemberResult() };
        ASSERT_TRUE(member.has_value());
        EXPECT_TRUE(member->is_outlier);
        EXPECT_DOUBLE_EQ(member->statistical_distance, 2.5);
    }
    EXPECT_THROW(analysis.ApplyAtomGroupGaussianResult(999, result), std::runtime_error);
}

TEST(DataObjectModelAnalysisTest, ModelAnalysisEditorRebuildsAtomGroupsFromSelection)
{
    auto model{ data_test::MakeModelWithBond() };
    auto * first_atom{ model->GetAtomList().at(0).get() };
    auto * second_atom{ model->GetAtomList().at(1).get() };
    auto & analysis_data{ rg::ModelAnalysisData::Of(*model) };
    auto analysis{ model->EditAnalysis() };

    model->SelectAllAtoms(false);
    model->SetAtomSelected(first_atom->GetSerialID(), true);
    analysis_data.AtomGroupEntry().AddMember(999, *second_atom);

    analysis.RebuildAtomGroupsFromSelection();

    const auto & component_group_entry{ analysis_data.AtomGroupEntry() };

    size_t member_count{ 0 };
    for (const auto group_key : component_group_entry.CollectGroupKeys())
    {
        member_count += component_group_entry.GetMemberCount(group_key);
        for (const auto * atom : component_group_entry.GetMembers(group_key))
        {
            ASSERT_NE(atom, nullptr);
            EXPECT_EQ(atom, first_atom);
        }
    }
    EXPECT_EQ(member_count, 1U);

    model->SetAtomSelected(first_atom->GetSerialID(), false);
    model->SetAtomSelected(second_atom->GetSerialID(), true);
    analysis.RebuildAtomGroupsFromSelection();

    member_count = 0;
    for (const auto group_key : component_group_entry.CollectGroupKeys())
    {
        member_count += component_group_entry.GetMemberCount(group_key);
        for (const auto * atom : component_group_entry.GetMembers(group_key))
        {
            ASSERT_NE(atom, nullptr);
            EXPECT_EQ(atom, second_atom);
        }
    }
    EXPECT_EQ(member_count, 1U);
}

TEST(DataObjectModelAnalysisTest, ModelAnalysisEditorInitializesLocalAlphaForSelectedAtoms)
{
    auto model{ data_test::MakeModelWithBond() };
    auto * first_atom{ model->GetAtomList().at(0).get() };
    auto * second_atom{ model->GetAtomList().at(1).get() };
    auto analysis{ model->EditAnalysis() };

    model->SelectAllAtoms(false);
    model->SetAtomSelected(first_atom->GetSerialID(), true);
    model->SetAtomSelected(second_atom->GetSerialID(), false);
    analysis.SetAtomLocalAlphaR(rg::FittingStage::Third, *second_atom, 0.9);

    analysis.InitializeLocalAlpha(rg::FittingStage::Third, 0.4);

    EXPECT_DOUBLE_EQ(
        0.4,
        rg::AtomLocalPotentialView::For(*first_atom).GetAlphaR(
            rg::FittingStage::Third));
    EXPECT_DOUBLE_EQ(
        0.9,
        rg::AtomLocalPotentialView::For(*second_atom).GetAlphaR(
        rg::FittingStage::Third));
}

TEST(DataObjectModelAnalysisTest, ModelAnalysisEditorEnsuresSelectedAtomLocalPotentials)
{
    auto model{ data_test::MakeModelWithBond() };
    auto * selected_atom{ model->GetAtomList().at(0).get() };
    auto * unselected_atom{ model->GetAtomList().at(1).get() };
    model->SelectAllAtoms(false);
    model->SetAtomSelected(selected_atom->GetSerialID(), true);
    auto analysis{ model->EditAnalysis() };

    EXPECT_FALSE(rg::AtomLocalPotentialView::For(*selected_atom).IsAvailable());
    EXPECT_FALSE(rg::AtomLocalPotentialView::For(*unselected_atom).IsAvailable());

    analysis.EnsureSelectedAtomLocalPotentials();

    EXPECT_TRUE(rg::AtomLocalPotentialView::For(*selected_atom).IsAvailable());
    EXPECT_FALSE(rg::AtomLocalPotentialView::For(*unselected_atom).IsAvailable());

    analysis.SetAtomLocalAlphaR(rg::FittingStage::Third, *selected_atom, 0.42);
    analysis.EnsureSelectedAtomLocalPotentials();
    EXPECT_DOUBLE_EQ(
        0.42,
        rg::AtomLocalPotentialView::For(*selected_atom).GetAlphaR(
            rg::FittingStage::Third));
}

TEST(DataObjectModelAnalysisTest, ModelAnalysisEditorEnsuresAtomGroupLocalPotentials)
{
    auto model{ data_test::MakeModelWithBond() };
    auto * first_atom{ model->GetAtomList().at(0).get() };
    auto * second_atom{ model->GetAtomList().at(1).get() };
    auto & analysis_data{ rg::ModelAnalysisData::Of(*model) };
    auto analysis{ model->EditAnalysis() };
    analysis_data.AtomGroupEntry().AddMember(101,
        *first_atom);
    analysis_data.AtomGroupEntry().AddMember(202,
        *second_atom);

    analysis.EnsureAtomGroupLocalPotentials(101);

    EXPECT_TRUE(rg::AtomLocalPotentialView::For(*first_atom).IsAvailable());
    EXPECT_FALSE(rg::AtomLocalPotentialView::For(*second_atom).IsAvailable());
    EXPECT_THROW(
        analysis.EnsureAtomGroupLocalPotentials(999),
        std::runtime_error);
}

TEST(DataObjectModelAnalysisTest, ModelAnalysisEditorSetsAtomGroupAlphaR)
{
    auto model{ data_test::MakeModelWithBond() };
    auto * first_atom{ model->GetAtomList().at(0).get() };
    auto * second_atom{ model->GetAtomList().at(1).get() };
    auto & analysis_data{ rg::ModelAnalysisData::Of(*model) };
    auto analysis{ model->EditAnalysis() };
    analysis_data.AtomGroupEntry().AddMember(101,
        *first_atom);
    analysis_data.AtomGroupEntry().AddMember(202,
        *second_atom);
    analysis.SetAtomLocalAlphaR(rg::FittingStage::Third, *first_atom, 0.8);

    analysis.SetAtomGroupAlphaR(rg::FittingStage::First, 101, 0.42);

    const auto first_view{ rg::AtomLocalPotentialView::For(*first_atom) };
    EXPECT_DOUBLE_EQ(0.42, first_view.GetAlphaR(rg::FittingStage::First));
    EXPECT_DOUBLE_EQ(0.8, first_view.GetAlphaR(rg::FittingStage::Third));
    EXPECT_FALSE(rg::AtomLocalPotentialView::For(*second_atom).IsAvailable());
}

TEST(DataObjectModelAnalysisTest, ModelAnalysisEditorAppliesAtomLocalGaussianResult)
{
    auto model{ data_test::MakeModelWithBond() };
    auto * target_atom{ model->GetAtomList().at(0).get() };
    auto * missing_atom{ model->GetAtomList().at(1).get() };
    auto analysis{ model->EditAnalysis() };

    rg::LocalGaussianResult result;
    result.alpha_r = 0.42;
    result.ols = rg::GaussianModel3DWithUncertainty{
        rg::GaussianModel3D{ 1.0, 0.6, -0.1 },
        rg::GaussianModel3DUncertainty{ 0.1, 0.2, 0.3 }
    };
    result.mdpde = rg::GaussianModel3DWithUncertainty{
        rg::GaussianModel3D{ 1.5, 0.8, -0.2 },
        rg::GaussianModel3DUncertainty{ 0.4, 0.5, 0.6 }
    };

    analysis.SetAtomLocalGaussianResult(
        rg::FittingStage::Second,
        *target_atom,
        result);

    const auto & saved_result{
        rg::AtomLocalPotentialView::For(*target_atom)
            .GetGaussianResult(rg::FittingStage::Second)
    };
    EXPECT_DOUBLE_EQ(0.42, saved_result.alpha_r);
    EXPECT_DOUBLE_EQ(1.0, saved_result.ols.GetModel().GetAmplitude());
    EXPECT_DOUBLE_EQ(0.8, saved_result.mdpde.GetModel().GetWidth());
    EXPECT_FALSE(rg::AtomLocalPotentialView::For(*target_atom)
        .GetGroupMemberResult().has_value());
    EXPECT_FALSE(rg::AtomLocalPotentialView::For(*missing_atom).IsAvailable());
    analysis.SetAtomLocalGaussianResult(
        rg::FittingStage::First,
        *missing_atom,
        result);
    EXPECT_TRUE(rg::AtomLocalPotentialView::For(*missing_atom).IsAvailable());
}

TEST(DataObjectModelAnalysisTest, ModelAnalysisEditorInitializesLocalFittingSeedModels)
{
    auto model{ data_test::MakeModelWithBond() };
    auto * selected_atom{ model->GetAtomList().at(0).get() };
    auto * unselected_atom{ model->GetAtomList().at(1).get() };
    model->SelectAllAtoms(false);
    model->SetAtomSelected(selected_atom->GetSerialID(), true);
    auto analysis{ model->EditAnalysis() };

    rg::LocalGaussianResult selected_result;
    selected_result.alpha_r = 0.42;
    selected_result.ols = rg::GaussianModel3DWithUncertainty{
        rg::GaussianModel3D{ 2.0, 0.7, -0.3 },
        rg::GaussianModel3DUncertainty{ 0.1, 0.2, 0.3 }
    };
    selected_result.mdpde = selected_result.ols;
    selected_result.fit_result = rg::RHBMBetaEstimateResult{};
    analysis.SetAtomLocalRawSamplingEntries(
        *selected_atom,
        { LocalPotentialSample{ 5.0, SamplingPoint{ 0.4 } } });
    analysis.SetAtomLocalPeelingSamplingEntries(
        *selected_atom,
        { LocalPotentialSample{ 3.0, SamplingPoint{ 0.6 } } });
    analysis.SetAtomLocalNeighborCountForPeeling(*selected_atom, 7);
    analysis.SetAtomLocalGaussianResult(
        rg::FittingStage::First, *selected_atom, selected_result);

    rg::LocalGaussianResult unselected_result;
    unselected_result.alpha_r = 0.9;
    unselected_result.mdpde = rg::GaussianModel3DWithUncertainty{
        rg::GaussianModel3D{ 4.0, 0.8, 0.2 },
        rg::GaussianModel3DUncertainty{}
    };
    analysis.SetAtomLocalGaussianResult(
        rg::FittingStage::Third, *unselected_atom, unselected_result);

    auto & analysis_data{ rg::ModelAnalysisData::Of(*model) };
    const rg::GroupGaussianMemberResult group_result{ selected_result.mdpde, true, 2.5 };
    analysis_data.EnsureAtomLocalEntry(*selected_atom).SetGroupMemberResult(group_result);
    analysis_data.EnsureAtomLocalEntry(*unselected_atom).SetGroupMemberResult(group_result);
    analysis.InitializeLocalFittingSeedModels();

    const auto selected_view{ rg::AtomLocalPotentialView::For(*selected_atom) };
    for (const auto stage : {
             rg::FittingStage::First,
             rg::FittingStage::Second,
             rg::FittingStage::Third })
    {
        const auto result{ selected_view.GetGaussianResult(stage) };
        EXPECT_DOUBLE_EQ(0.42, result.alpha_r);
        EXPECT_DOUBLE_EQ(0.0, result.ols.GetModel().GetAmplitude());
        EXPECT_DOUBLE_EQ(1.0, result.ols.GetModel().GetWidth());
        EXPECT_DOUBLE_EQ(0.0, result.ols.GetModel().GetOffset());
        EXPECT_DOUBLE_EQ(0.0, result.mdpde.GetModel().GetAmplitude());
        EXPECT_DOUBLE_EQ(1.0, result.mdpde.GetModel().GetWidth());
        EXPECT_DOUBLE_EQ(0.0, result.mdpde.GetModel().GetOffset());
        EXPECT_FALSE(result.fit_result.has_value());
    }
    EXPECT_EQ(1u, selected_view.GetRawSamplingEntries(false).size());
    EXPECT_EQ(1u, selected_view.GetPeelingSamplingEntries(false).size());
    EXPECT_EQ(7, selected_view.GetNeighborCountForPeeling());
    EXPECT_FALSE(selected_view.GetGroupMemberResult().has_value());

    const auto unselected_view{ rg::AtomLocalPotentialView::For(*unselected_atom) };
    ASSERT_TRUE(unselected_view.GetGroupMemberResult().has_value());
    EXPECT_TRUE(unselected_view.GetGroupMemberResult()->is_outlier);
    EXPECT_DOUBLE_EQ(unselected_view.GetGroupMemberResult()->statistical_distance, 2.5);
    EXPECT_DOUBLE_EQ(
        0.9,
        unselected_view.GetGaussianResult(rg::FittingStage::Third).alpha_r);
    EXPECT_DOUBLE_EQ(
        4.0,
        unselected_view.GetGaussianResult(rg::FittingStage::Third)
            .mdpde.GetModel().GetAmplitude());
}

TEST(DataObjectModelAnalysisTest, ModelAnalysisEditorInitializesMissingLocalFittingSeedEntry)
{
    auto model{ data_test::MakeModelWithBond() };
    auto * atom{ model->GetAtomList().at(0).get() };
    model->SelectAllAtoms(false);
    model->SetAtomSelected(atom->GetSerialID(), true);

    const auto view_before{ rg::AtomLocalPotentialView::For(*atom) };
    EXPECT_FALSE(view_before.IsAvailable());

    model->EditAnalysis().InitializeLocalFittingSeedModels();

    const auto view_after{ rg::AtomLocalPotentialView::For(*atom) };
    EXPECT_DOUBLE_EQ(
        1.0,
        view_after.GetGaussianResult(rg::FittingStage::Third)
            .ols.GetModel().GetWidth());
}

TEST(DataObjectModelAnalysisTest, ModelAnalysisEditorInitializesGroupAlphaForExistingGroups)
{
    auto model{ data_test::MakeModelWithBond() };
    model->SelectAllAtoms();
    model->ApplySymmetrySelection(false);
    auto analysis{ model->EditAnalysis() };

    EXPECT_NO_THROW(analysis.InitializeGroupAlpha(0.3));

    analysis.RebuildAtomGroupsFromSelection();
    analysis.InitializeGroupAlpha(0.3);
    const auto analysis_view{ model->GetAnalysisView() };
    const auto group_key_list{ analysis_view.CollectAtomGroupKeys() };
    for (const auto group_key : group_key_list)
    {
        EXPECT_DOUBLE_EQ(0.3, analysis_view.GetAtomAlphaG(group_key));
    }
    ASSERT_FALSE(group_key_list.empty());
    analysis.SetAtomGroupAlphaG(group_key_list.front(),
        0.8);
    EXPECT_DOUBLE_EQ(
        analysis_view.GetAtomAlphaG(group_key_list.front()),
        0.8);
}

TEST(DataObjectModelAnalysisTest, ModelAnalysisEditorCopiesLocalFittingStageResult)
{
    auto model{ data_test::MakeModelWithBond() };
    auto * selected_atom{ model->GetAtomList().at(0).get() };
    auto * unselected_atom{ model->GetAtomList().at(1).get() };
    model->SelectAllAtoms(false);
    model->SetAtomSelected(selected_atom->GetSerialID(), true);
    auto analysis{ model->EditAnalysis() };
    analysis.RebuildAtomGroupsFromSelection();

    rg::LocalGaussianResult source_result;
    source_result.alpha_r = 0.7;
    source_result.ols = rg::GaussianModel3DWithUncertainty{
        rg::GaussianModel3D{ 1.0, 0.6, -0.1 },
        rg::GaussianModel3DUncertainty{ 0.1, 0.2, 0.3 }
    };
    source_result.mdpde = rg::GaussianModel3DWithUncertainty{
        rg::GaussianModel3D{ 1.5, 0.8, -0.2 },
        rg::GaussianModel3DUncertainty{ 0.4, 0.5, 0.6 }
    };
    rg::GroupGaussianMemberResult source_group_result;
    source_group_result.posterior = rg::GaussianModel3DWithUncertainty{
        rg::GaussianModel3D{ 1.8, 0.9, -0.3 },
        rg::GaussianModel3DUncertainty{ 0.7, 0.8, 0.9 }
    };
    source_group_result.is_outlier = true;
    source_group_result.statistical_distance = 2.5;
    rg::ModelAnalysisData::Of(*model).EnsureAtomLocalEntry(*selected_atom)
        .SetGroupMemberResult(source_group_result);
    source_result.fit_result = rg::RHBMBetaEstimateResult{};
    analysis.SetAtomLocalGaussianResult(
        rg::FittingStage::First,
        *selected_atom,
        source_result);

    rg::LocalGaussianResult unselected_result;
    unselected_result.alpha_r = 0.9;
    analysis.SetAtomLocalGaussianResult(
        rg::FittingStage::Second,
        *unselected_atom,
        unselected_result);

    analysis.CopyLocalFittingStageResult(
        rg::FittingStage::First,
        rg::FittingStage::Second);

    const auto copied_result{
        rg::AtomLocalPotentialView::For(*selected_atom)
            .GetGaussianResult(rg::FittingStage::Second)
    };
    EXPECT_DOUBLE_EQ(copied_result.alpha_r, source_result.alpha_r);
    EXPECT_DOUBLE_EQ(
        copied_result.ols.GetModel().GetOffset(),
        source_result.ols.GetModel().GetOffset());
    EXPECT_DOUBLE_EQ(
        copied_result.ols.GetStandardDeviationModel().GetWidth(),
        source_result.ols.GetStandardDeviationModel().GetWidth());
    EXPECT_DOUBLE_EQ(
        copied_result.mdpde.GetModel().GetAmplitude(),
        source_result.mdpde.GetModel().GetAmplitude());
    const auto & group_result{
        rg::AtomLocalPotentialView::For(*selected_atom).GetGroupMemberResult()
    };
    ASSERT_TRUE(group_result.has_value());
    EXPECT_DOUBLE_EQ(group_result->posterior.GetModel().GetWidth(),
        source_group_result.posterior.GetModel().GetWidth());
    EXPECT_TRUE(group_result->is_outlier);
    EXPECT_DOUBLE_EQ(group_result->statistical_distance, source_group_result.statistical_distance);
    EXPECT_TRUE(copied_result.fit_result.has_value());

    EXPECT_DOUBLE_EQ(
        rg::AtomLocalPotentialView::For(*unselected_atom)
            .GetGaussianResult(rg::FittingStage::Second)
            .alpha_r,
        unselected_result.alpha_r);

    analysis.SetAtomLocalAlphaR(rg::FittingStage::Second, *selected_atom, 0.95);
    EXPECT_DOUBLE_EQ(
        rg::AtomLocalPotentialView::For(*selected_atom)
            .GetGaussianResult(rg::FittingStage::First)
            .alpha_r,
        source_result.alpha_r);
}

TEST(DataObjectModelAnalysisTest, CollectAtomGroupKeysReturnsRebuiltGroupKeySet)
{
    auto model{ data_test::MakeModelWithBond() };
    auto & analysis_data{ rg::ModelAnalysisData::Of(*model) };
    auto analysis{ model->EditAnalysis() };

    model->SelectAllAtoms();
    analysis.RebuildAtomGroupsFromSelection();

    const auto & group_entry{ analysis_data.AtomGroupEntry() };
    const auto group_keys{ group_entry.CollectGroupKeys() };
    EXPECT_EQ(
        group_keys.size(),
        group_entry.GroupCount());
    for (const auto & group_key : group_keys)
    {
        EXPECT_TRUE(group_entry.HasGroup(group_key));
    }
}

TEST(DataObjectModelAnalysisTest, RebuildAtomGroupsUsesComponentAtomGroupKeys)
{
    auto model{ data_test::MakeModelWithBond() };
    auto analysis{ model->EditAnalysis() };

    model->SelectAllAtoms();
    analysis.RebuildAtomGroupsFromSelection();

    const auto analysis_view{ model->GetAnalysisView() };
    for (const auto group_key : analysis_view.CollectAtomGroupKeys())
    {
        EXPECT_TRUE(analysis_view.HasAtomGroup(group_key));
    }
}

TEST(DataObjectModelAnalysisTest, AtomGroupKeyCollectionCoversSingleGroupEntry)
{
    auto model{ data_test::MakeModelWithBond() };
    model->SelectAllAtoms();
    model->ApplySymmetrySelection(false);
    auto analysis{ model->EditAnalysis() };
    analysis.RebuildAtomGroupsFromSelection();
    const auto & analysis_data{ rg::ModelAnalysisData::Of(*model) };
    const auto analysis_view{ model->GetAnalysisView() };

    EXPECT_EQ(
        analysis_view.CollectAtomGroupKeys().size(),
        analysis_data.AtomGroupEntry().GroupCount());
}

TEST(DataObjectModelAnalysisTest, GroupKeyCollectionsAreSafeBeforeRebuildAndEmpty)
{
    auto model{ data_test::MakeModelWithBond() };
    model->SelectAllAtoms(false);

    const auto analysis_view{ model->GetAnalysisView() };
    EXPECT_TRUE(analysis_view.CollectAtomGroupKeys().empty());
}

TEST(DataObjectModelAnalysisTest, ModelAtomsExposeStableSerialAndPositionInputsForTypedWorkflows)
{
    auto model{ data_test::MakeModelWithBond() };
    std::vector<rg::AtomObject *> atom_list;
    std::unordered_map<int, double> atom_charge_map;
    std::array<double, 3> range_minimum{
        std::numeric_limits<double>::max(),
        std::numeric_limits<double>::max(),
        std::numeric_limits<double>::max() };
    std::array<double, 3> range_maximum{
        std::numeric_limits<double>::lowest(),
        std::numeric_limits<double>::lowest(),
        std::numeric_limits<double>::lowest() };
    for (auto & atom : model->GetAtomList())
    {
        atom_list.emplace_back(atom.get());
        atom_charge_map.emplace(atom->GetSerialID(), 0.0);
        const auto & atom_position{ atom->GetPositionRef() };
        range_minimum[0] = std::min(range_minimum[0], atom_position[0]);
        range_minimum[1] = std::min(range_minimum[1], atom_position[1]);
        range_minimum[2] = std::min(range_minimum[2], atom_position[2]);
        range_maximum[0] = std::max(range_maximum[0], atom_position[0]);
        range_maximum[1] = std::max(range_maximum[1], atom_position[1]);
        range_maximum[2] = std::max(range_maximum[2], atom_position[2]);
    }

    ASSERT_FALSE(atom_list.empty());
    EXPECT_EQ(atom_list.size(), 2);
    ASSERT_EQ(atom_charge_map.size(), 2);
    EXPECT_DOUBLE_EQ(atom_charge_map.at(1), 0.0);
    EXPECT_DOUBLE_EQ(atom_charge_map.at(2), 0.0);
    EXPECT_DOUBLE_EQ(range_minimum[0], 0.0);
    EXPECT_DOUBLE_EQ(range_maximum[0], 1.0);
}

TEST(DataObjectModelAnalysisTest, SelectedAtomsAndBondsRemainQueryableForContextAssembly)
{
    auto model{ data_test::MakeModelWithBond() };
    auto & atoms{ model->GetAtomList() };
    auto & bonds{ model->GetBondList() };
    ASSERT_EQ(atoms.size(), 2);
    ASSERT_EQ(bonds.size(), 1);
    model->SelectAllAtoms();
    model->SelectAllBonds();

    std::unordered_map<int, rg::AtomObject *> atom_map;
    std::unordered_map<int, std::vector<rg::BondObject *>> bond_map;
    for (auto * atom : model->GetSelectedAtoms())
    {
        atom_map.emplace(atom->GetSerialID(), atom);
    }
    for (auto * bond : model->GetSelectedBonds())
    {
        bond_map[bond->GetAtomSerialID1()].emplace_back(bond);
        bond_map[bond->GetAtomSerialID2()].emplace_back(bond);
    }

    ASSERT_EQ(atom_map.size(), 2);
    ASSERT_EQ(bond_map.size(), 2);
    EXPECT_EQ(atom_map.at(1), atoms[0].get());
    EXPECT_EQ(atom_map.at(2), atoms[1].get());
    EXPECT_EQ(bond_map.at(1).size(), 1);
    EXPECT_EQ(bond_map.at(2).size(), 1);
}
