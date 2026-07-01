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
#include <rhbm_gem/data/object/ModelAnalysisView.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/domain/AtomSelector.hpp>
#include <rhbm_gem/utils/domain/ChemicalDataHelper.hpp>
#include "support/DataObjectTestSupport.hpp"

namespace rg = rhbm_gem;

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

TEST(DataObjectModelAnalysisTest, AtomCountingSummaryReportsSelectedElementCounts)
{
    auto model{ data_test::MakeModelWithBond() };
    auto & atoms{ model->GetAtomList() };
    atoms.at(0)->SetElement(Element::CARBON);
    atoms.at(1)->SetElement(Element::HYDROGEN);
    model->SelectAllAtoms();

    const auto summary{ model->GetAnalysisView().GetAtomCountingSummary() };

    EXPECT_NE(summary.find("Number of selected atom = 2"), std::string::npos);
    EXPECT_NE(summary.find(" - Element type: H include 1 atoms."), std::string::npos);
    EXPECT_NE(summary.find(" - Element type: C include 1 atoms."), std::string::npos);
}

TEST(DataObjectModelAnalysisTest, AtomCountingSummaryOmitsUnselectedElements)
{
    auto model{ data_test::MakeModelWithBond() };
    auto & atoms{ model->GetAtomList() };
    atoms.at(0)->SetElement(Element::CARBON);
    atoms.at(1)->SetElement(Element::HYDROGEN);
    model->SelectAllAtoms(false);
    model->SetAtomSelected(atoms.at(0)->GetSerialID(), true);

    const auto summary{ model->GetAnalysisView().GetAtomCountingSummary() };

    EXPECT_NE(summary.find("Number of selected atom = 1"), std::string::npos);
    EXPECT_NE(summary.find(" - Element type: C include 1 atoms."), std::string::npos);
    EXPECT_EQ(summary.find(" - Element type: H include"), std::string::npos);
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
    analysis_data.EnsureAtomLocalEntry(*second_atom).SetAlphaR(0.7);
    analysis_data.AtomGroupEntry().AddMember(999, *second_atom);

    model->LocalPotentialInitialization();

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
    EXPECT_DOUBLE_EQ(0.0, rg::AtomLocalPotentialView::RequireFor(*first_atom).GetAlphaR());
    EXPECT_EQ(analysis_data.FindAtomLocalEntry(*second_atom), nullptr);

    const auto analysis_view{ model->GetAnalysisView() };
    for (const auto group_key : analysis_view.CollectAtomGroupKeys())
    {
        EXPECT_DOUBLE_EQ(0.0, analysis_view.GetAtomAlphaG(group_key));
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
    atom_entry.SetGaussianResult(atom_result);

    ASSERT_NE(analysis_data.FindAtomLocalEntry(*atom), nullptr);
    ASSERT_FALSE(analysis_data.AtomGroupEntry().CollectGroupKeys().empty());

    analysis.ClearTransientFitStates();

    const auto * cleared_atom_entry{ analysis_data.FindAtomLocalEntry(*atom) };
    ASSERT_NE(cleared_atom_entry, nullptr);
    EXPECT_FALSE(analysis_data.AtomGroupEntry().CollectGroupKeys().empty());
    EXPECT_DOUBLE_EQ(0.2, cleared_atom_entry->GaussianResult().alpha_r);
    EXPECT_FALSE(cleared_atom_entry->GaussianResult().fit_result.has_value());
}

TEST(DataObjectModelAnalysisTest, ModelAnalysisDataClearDropsEntriesAndFitStates)
{
    auto model{ data_test::MakeModelWithBond() };
    auto * atom{ model->GetAtomList().at(0).get() };
    auto & analysis_data{ rg::ModelAnalysisData::Of(*model) };

    auto & atom_entry{ analysis_data.EnsureAtomLocalEntry(*atom) };
    analysis_data.AtomGroupEntry().AddMember(999, *atom);
    atom_entry.SetAlphaR(0.2);

    analysis_data.Clear();

    EXPECT_EQ(analysis_data.FindAtomLocalEntry(*atom), nullptr);
    EXPECT_TRUE(analysis_data.AtomGroupEntry().CollectGroupKeys().empty());
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
    entry.SetGaussianResult(result);
    ASSERT_TRUE(entry.GaussianResult().fit_result.has_value());

    entry.ClearTransientFitState();

    EXPECT_DOUBLE_EQ(0.4, entry.GaussianResult().alpha_r);
    EXPECT_DOUBLE_EQ(2.0, entry.GaussianResult().mdpde.GetModel().GetAmplitude());
    EXPECT_DOUBLE_EQ(0.7, entry.GaussianResult().mdpde.GetModel().GetWidth());
    EXPECT_FALSE(entry.GaussianResult().fit_result.has_value());
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

    entry.SetGaussianResult(result);

    EXPECT_DOUBLE_EQ(0.5, entry.GaussianResult().alpha_r);
    EXPECT_DOUBLE_EQ(1.0, entry.GaussianResult().ols.GetModel().GetAmplitude());
    EXPECT_DOUBLE_EQ(1.5, entry.GaussianResult().mdpde.GetModel().GetAmplitude());
}

TEST(DataObjectModelAnalysisTest, AtomLocalPotentialEditorCanSetAlphaR)
{
    auto model{ data_test::MakeModelWithBond() };
    auto * atom{ model->GetAtomList().at(0).get() };
    auto analysis{ model->EditAnalysis() };
    auto editor{ analysis.EnsureAtomLocalPotential(*atom) };

    editor.SetAlphaR(0.37);

    EXPECT_DOUBLE_EQ(0.37, rg::AtomLocalPotentialView::RequireFor(*atom).GetAlphaR());
}

TEST(DataObjectModelAnalysisTest, AtomLocalPotentialEditorSetGaussianResultUpdatesViewEstimates)
{
    auto model{ data_test::MakeModelWithBond() };
    auto * atom{ model->GetAtomList().at(0).get() };
    auto analysis{ model->EditAnalysis() };
    auto editor{ analysis.EnsureAtomLocalPotential(*atom) };

    LocalPotentialSampleList sampling_entries{
        LocalPotentialSample{ 6.0f, SamplingPoint{ 0.0f } },
        LocalPotentialSample{ 4.0f, SamplingPoint{ 0.5f } },
        LocalPotentialSample{ 2.0f, SamplingPoint{ 0.9f } }
    };
    editor.SetSamplingEntries(sampling_entries);

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

    editor.SetGaussianResult(gaussian_result);

    EXPECT_DOUBLE_EQ(0.6, rg::AtomLocalPotentialView::RequireFor(*atom).GetGaussianResult().alpha_r);
    EXPECT_DOUBLE_EQ(0.0, rg::AtomLocalPotentialView::RequireFor(*atom).GetEstimateOLS().GetWidth());
    EXPECT_DOUBLE_EQ(0.0, rg::AtomLocalPotentialView::RequireFor(*atom).GetEstimateMDPDE().GetWidth());
}

TEST(DataObjectModelAnalysisTest, AtomLocalPotentialViewCanApplySamplingSelection)
{
    auto model{ data_test::MakeModelWithBond() };
    auto * atom{ model->GetAtomList().at(0).get() };
    auto analysis{ model->EditAnalysis() };
    auto editor{ analysis.EnsureAtomLocalPotential(*atom) };

    editor.SetSamplingEntries({
        LocalPotentialSample{ 6.0f, SamplingPoint{ 0.0f, { 0.0f, 0.0f, 0.0f }, true } },
        LocalPotentialSample{ 4.0f, SamplingPoint{ 0.5f, { 0.0f, 0.0f, 0.0f }, false } },
        LocalPotentialSample{ 2.0f, SamplingPoint{ 0.9f, { 0.0f, 0.0f, 0.0f }, true } }
    });

    const auto view{ rg::AtomLocalPotentialView::RequireFor(*atom) };
    const auto selected_entries{ view.GetSamplingEntries() };
    const auto all_entries{ view.GetSamplingEntries(false) };

    ASSERT_EQ(selected_entries.size(), 2u);
    EXPECT_FLOAT_EQ(selected_entries.at(0).response, 6.0f);
    EXPECT_FLOAT_EQ(selected_entries.at(1).response, 2.0f);
    ASSERT_EQ(all_entries.size(), 3u);
    EXPECT_FALSE(all_entries.at(1).point.is_selected);
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
    analysis.EnsureAtomLocalPotential(*atom_list.front()).SetGaussianResult(first_atom_local_result);

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
        rg::LocalGaussianResult member_result;
        member_result.mdpde = rg::GaussianModel3DWithUncertainty{
            rg::GaussianModel3D{ 0.4 + 0.1 * static_cast<double>(i), 0.9 },
            rg::GaussianModel3DUncertainty{ 0.01, 0.02 }
        };
        member_result.is_outlier = (i == 0);
        member_result.statistical_distance = 1.5 + static_cast<double>(i);
        result.member_results.emplace_back(member_result);
    }

    analysis.ApplyAtomGroupGaussianResult(group_key, result);

    EXPECT_NEAR(result.mean.GetAmplitude(), analysis_view.GetAtomGroupMean(group_key).GetAmplitude(), 1e-12);
    EXPECT_NEAR(result.mean.GetWidth(), analysis_view.GetAtomGroupMean(group_key).GetWidth(), 1e-12);
    EXPECT_NEAR(result.mdpde.GetAmplitude(), analysis_view.GetAtomGroupMDPDE(group_key).GetAmplitude(), 1e-12);
    EXPECT_NEAR(result.mdpde.GetWidth(), analysis_view.GetAtomGroupMDPDE(group_key).GetWidth(), 1e-12);
    EXPECT_NEAR(result.prior.GetModel().GetAmplitude(), analysis_view.GetAtomGroupPrior(group_key).GetAmplitude(), 1e-12);
    EXPECT_NEAR(result.prior.GetModel().GetWidth(), analysis_view.GetAtomGroupPrior(group_key).GetWidth(), 1e-12);
    EXPECT_NEAR(result.prior.GetStandardDeviationModel().GetAmplitude(), analysis_view.GetAtomGroupPriorWithUncertainty(group_key).GetStandardDeviationModel().GetAmplitude(), 1e-12);
    EXPECT_NEAR(result.prior.GetStandardDeviationModel().GetWidth(), analysis_view.GetAtomGroupPriorWithUncertainty(group_key).GetStandardDeviationModel().GetWidth(), 1e-12);
    EXPECT_DOUBLE_EQ(alpha_g, analysis_view.GetAtomAlphaG(group_key));

    const auto & gaussian_result{
        rg::AtomLocalPotentialView::RequireFor(*atom_list.front()).GetGaussianResult()
    };
    EXPECT_NEAR(
        first_atom_local_result.mdpde.GetModel().GetAmplitude(),
        gaussian_result.mdpde.GetModel().GetAmplitude(),
        1e-12);
    EXPECT_NEAR(
        first_atom_local_result.mdpde.GetModel().GetWidth(),
        gaussian_result.mdpde.GetModel().GetWidth(),
        1e-12);
    ASSERT_TRUE(gaussian_result.posterior.has_value());
    const auto expected_gaussian{ result.member_results.front().mdpde };
    EXPECT_NEAR(
        expected_gaussian.GetModel().GetAmplitude(),
        gaussian_result.posterior->GetModel().GetAmplitude(),
        1e-12);
    EXPECT_NEAR(
        expected_gaussian.GetModel().GetWidth(),
        gaussian_result.posterior->GetModel().GetWidth(),
        1e-12);
    EXPECT_NEAR(
        expected_gaussian.GetStandardDeviationModel().GetAmplitude(),
        gaussian_result.posterior->GetStandardDeviationModel().GetAmplitude(),
        1e-12);
    EXPECT_NEAR(
        expected_gaussian.GetStandardDeviationModel().GetWidth(),
        gaussian_result.posterior->GetStandardDeviationModel().GetWidth(),
        1e-12);
    EXPECT_TRUE(gaussian_result.is_outlier);
    EXPECT_DOUBLE_EQ(1.5, gaussian_result.statistical_distance);
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

    EXPECT_THROW(
        analysis.ApplyAtomGroupGaussianResult(group_key, result),
        std::invalid_argument);
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
    analysis.EnsureAtomLocalPotential(*second_atom).SetAlphaR(0.9);

    analysis.InitializeLocalAlpha(0.4);

    EXPECT_DOUBLE_EQ(0.4, rg::AtomLocalPotentialView::RequireFor(*first_atom).GetAlphaR());
    EXPECT_DOUBLE_EQ(0.9, rg::AtomLocalPotentialView::RequireFor(*second_atom).GetAlphaR());
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
    for (const auto group_key : analysis_view.CollectAtomGroupKeys())
    {
        EXPECT_DOUBLE_EQ(0.3, analysis_view.GetAtomAlphaG(group_key));
    }
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
    EXPECT_EQ(group_keys.size(), group_entry.GroupCount());
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

TEST(DataObjectModelAnalysisTest, AtomGroupingSummaryReportsAtomGroupCount)
{
    auto model{ data_test::MakeModelWithBond() };
    model->SelectAllAtoms();
    model->ApplySymmetrySelection(false);
    auto analysis{ model->EditAnalysis() };
    analysis.RebuildAtomGroupsFromSelection();
    const auto analysis_view{ model->GetAnalysisView() };

    const auto summary{ analysis_view.GetAtomGroupingSummary() };

    const std::string expected_line{
        "Atomic model includes "
        + std::to_string(analysis_view.CollectAtomGroupKeys().size()) + " atom groups."
    };
    EXPECT_NE(summary.find(expected_line), std::string::npos);
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
    std::array<float, 3> range_minimum{
        std::numeric_limits<float>::max(),
        std::numeric_limits<float>::max(),
        std::numeric_limits<float>::max() };
    std::array<float, 3> range_maximum{
        std::numeric_limits<float>::lowest(),
        std::numeric_limits<float>::lowest(),
        std::numeric_limits<float>::lowest() };
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
    EXPECT_FLOAT_EQ(range_minimum[0], 0.0f);
    EXPECT_FLOAT_EQ(range_maximum[0], 1.0f);
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
