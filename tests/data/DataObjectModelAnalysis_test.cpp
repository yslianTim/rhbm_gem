#include <gtest/gtest.h>

#include <array>
#include <limits>
#include <memory>
#include <tuple>
#include <unordered_map>
#include <vector>

#include "data/detail/LocalPotentialEntry.hpp"
#include "data/detail/ModelAnalysisData.hpp"
#include <rhbm_gem/data/object/ModelAnalysisView.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/domain/AtomSelector.hpp>
#include <rhbm_gem/utils/domain/ChemicalDataHelper.hpp>
#include <rhbm_gem/utils/hrl/LinearizationService.hpp>
#include "support/DataObjectTestSupport.hpp"

namespace rg = rhbm_gem;
namespace ls = rhbm_gem::linearization_service;

namespace
{

const ls::LinearizationSpec & AtomGroupDecodeSpec()
{
    static const auto spec{ ls::LinearizationSpec::AtomGroupDecode() };
    return spec;
}

const ls::LinearizationSpec & BondGroupDecodeSpec()
{
    static const auto spec{ ls::LinearizationSpec::BondGroupDecode() };
    return spec;
}

} // namespace

TEST(DataObjectModelAnalysisTest, SelectedModelEntriesCanBeInitializedForTypedWorkflows)
{
    auto model{ data_test::MakeModelWithBond() };
    model->SelectAllAtoms(false);
    model->SelectAllBonds(false);
    ASSERT_EQ(model->GetSelectedAtomCount(), 0);
    ASSERT_EQ(model->GetSelectedBondCount(), 0);

    model->SelectAllAtoms();
    model->SelectAllBonds();
    model->ApplySymmetrySelection(false);
    for (auto * atom : model->GetSelectedAtoms())
    {
        rg::ModelAnalysisData::Of(*model).EnsureAtomLocalEntry(*atom);
    }
    for (auto * bond : model->GetSelectedBonds())
    {
        rg::ModelAnalysisData::Of(*model).EnsureBondLocalEntry(*bond);
    }

    EXPECT_EQ(model->GetSelectedAtomCount(), model->GetNumberOfAtom());
    EXPECT_EQ(model->GetSelectedBondCount(), model->GetNumberOfBond());
    for (const auto * atom : model->GetSelectedAtoms())
    {
        ASSERT_NE(atom, nullptr);
        EXPECT_NE(rg::ModelAnalysisData::FindLocalEntry(*atom), nullptr);
    }
    for (const auto * bond : model->GetSelectedBonds())
    {
        ASSERT_NE(bond, nullptr);
        EXPECT_NE(rg::ModelAnalysisData::FindLocalEntry(*bond), nullptr);
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
        if (rg::ModelAnalysisData::FindLocalEntry(*atom) != nullptr)
        {
            require_entry_atoms.emplace_back(atom.get());
        }
    }
    ASSERT_EQ(require_entry_atoms.size(), 1);
    EXPECT_EQ(require_entry_atoms.front(), atoms[0].get());
}

TEST(DataObjectModelAnalysisTest, ModelAnalysisDataCanClearTransientFitStatesWithoutDroppingEntries)
{
    auto model{ data_test::MakeModelWithBond() };
    auto * atom{ model->GetAtomList().at(0).get() };
    auto * bond{ model->GetBondList().at(0).get() };
    auto & analysis_data{ rg::ModelAnalysisData::Of(*model) };

    auto & atom_entry{ analysis_data.EnsureAtomLocalEntry(*atom) };
    auto & bond_entry{ analysis_data.EnsureBondLocalEntry(*bond) };
    analysis_data.EnsureAtomGroupEntry("atom_class");
    analysis_data.EnsureBondGroupEntry("bond_class");
    atom_entry.SetDataset({});
    atom_entry.SetFitResult({});
    bond_entry.SetDataset({});
    bond_entry.SetFitResult({});

    ASSERT_NE(analysis_data.FindAtomLocalEntry(*atom), nullptr);
    ASSERT_NE(analysis_data.FindBondLocalEntry(*bond), nullptr);
    ASSERT_NE(analysis_data.FindAtomGroupEntry("atom_class"), nullptr);
    ASSERT_NE(analysis_data.FindBondGroupEntry("bond_class"), nullptr);
    ASSERT_TRUE(atom_entry.HasDataset());
    ASSERT_TRUE(atom_entry.HasFitResult());
    ASSERT_TRUE(bond_entry.HasDataset());
    ASSERT_TRUE(bond_entry.HasFitResult());

    analysis_data.ClearTransientFitStates();

    const auto * cleared_atom_entry{ analysis_data.FindAtomLocalEntry(*atom) };
    const auto * cleared_bond_entry{ analysis_data.FindBondLocalEntry(*bond) };
    ASSERT_NE(cleared_atom_entry, nullptr);
    ASSERT_NE(cleared_bond_entry, nullptr);
    EXPECT_NE(analysis_data.FindAtomGroupEntry("atom_class"), nullptr);
    EXPECT_NE(analysis_data.FindBondGroupEntry("bond_class"), nullptr);
    EXPECT_FALSE(cleared_atom_entry->HasDataset());
    EXPECT_FALSE(cleared_atom_entry->HasFitResult());
    EXPECT_FALSE(cleared_bond_entry->HasDataset());
    EXPECT_FALSE(cleared_bond_entry->HasFitResult());
}

TEST(DataObjectModelAnalysisTest, ModelAnalysisDataClearDropsEntriesAndFitStates)
{
    auto model{ data_test::MakeModelWithBond() };
    auto * atom{ model->GetAtomList().at(0).get() };
    auto * bond{ model->GetBondList().at(0).get() };
    auto & analysis_data{ rg::ModelAnalysisData::Of(*model) };

    auto & atom_entry{ analysis_data.EnsureAtomLocalEntry(*atom) };
    auto & bond_entry{ analysis_data.EnsureBondLocalEntry(*bond) };
    analysis_data.EnsureAtomGroupEntry("atom_class");
    analysis_data.EnsureBondGroupEntry("bond_class");
    atom_entry.SetDataset({});
    atom_entry.SetFitResult({});
    bond_entry.SetDataset({});
    bond_entry.SetFitResult({});

    analysis_data.Clear();

    EXPECT_EQ(analysis_data.FindAtomLocalEntry(*atom), nullptr);
    EXPECT_EQ(analysis_data.FindBondLocalEntry(*bond), nullptr);
    EXPECT_EQ(analysis_data.FindAtomGroupEntry("atom_class"), nullptr);
    EXPECT_EQ(analysis_data.FindBondGroupEntry("bond_class"), nullptr);
}

TEST(DataObjectModelAnalysisTest, LocalPotentialEntryCanClearTransientFitPayload)
{
    rg::LocalPotentialEntry entry;
    entry.SetDataset({});
    entry.SetFitResult({});

    ASSERT_TRUE(entry.HasDataset());
    ASSERT_TRUE(entry.HasFitResult());

    entry.ClearTransientFitState();

    EXPECT_FALSE(entry.HasDataset());
    EXPECT_FALSE(entry.HasFitResult());
}

TEST(DataObjectModelAnalysisTest, LocalPotentialEntryPreservesFitResultStatus)
{
    rg::LocalPotentialEntry entry;
    RHBMBetaEstimateResult fit_result;
    fit_result.status = RHBMEstimationStatus::NUMERICAL_FALLBACK;

    entry.SetFitResult(fit_result);

    ASSERT_TRUE(entry.HasFitResult());
    EXPECT_EQ(RHBMEstimationStatus::NUMERICAL_FALLBACK, entry.GetFitResult().status);
}

TEST(DataObjectModelAnalysisTest, MutableLocalPotentialViewSetFitResultUpdatesFitAndEstimates)
{
    auto model{ data_test::MakeModelWithBond() };
    auto * atom{ model->GetAtomList().at(0).get() };
    auto analysis{ model->EditAnalysis() };
    auto entry{ analysis.EnsureAtomLocalPotential(*atom) };

    LocalPotentialSampleList sampling_entries{
        LocalPotentialSample{ 0.0f, 6.0f, 1.0f },
        LocalPotentialSample{ 0.5f, 4.0f, 1.0f },
        LocalPotentialSample{ 0.9f, 2.0f, 1.0f }
    };
    entry.SetSamplingEntries(sampling_entries);

    RHBMBetaEstimateResult fit_result;
    fit_result.status = RHBMEstimationStatus::SUCCESS;
    fit_result.beta_ols = Eigen::VectorXd::Zero(3);
    fit_result.beta_mdpde = Eigen::VectorXd::Zero(3);

    entry.SetFitResult(fit_result);

    ASSERT_TRUE(entry.HasFitResult());
    EXPECT_EQ(RHBMEstimationStatus::SUCCESS, entry.GetFitResult().status);
    EXPECT_TRUE(entry.GetFitResult().beta_ols.isApprox(Eigen::VectorXd::Zero(3), 1e-12));
    EXPECT_TRUE(entry.GetFitResult().beta_mdpde.isApprox(Eigen::VectorXd::Zero(3), 1e-12));
    EXPECT_DOUBLE_EQ(0.0, rg::LocalPotentialView::RequireFor(*atom).GetEstimateOLS().width);
    EXPECT_DOUBLE_EQ(0.0, rg::LocalPotentialView::RequireFor(*atom).GetEstimateMDPDE().width);
}

TEST(DataObjectModelAnalysisTest, ModelAnalysisEditorAppliesAtomGroupEstimateResultToStatisticsAndAnnotations)
{
    auto model{ data_test::MakeModelWithBond() };
    model->SelectAllAtoms();
    model->ApplySymmetrySelection(false);
    auto analysis{ model->EditAnalysis() };
    analysis.RebuildAtomGroupsFromSelection();
    const auto analysis_view{ model->GetAnalysisView() };
    const auto & class_key{ ChemicalDataHelper::GetSimpleAtomClassKey() };
    const auto group_keys{ analysis_view.CollectAtomGroupKeys(class_key) };
    ASSERT_FALSE(group_keys.empty());

    const auto group_key{ group_keys.front() };
    const auto & atom_list{ analysis_view.GetAtomObjectList(group_key, class_key) };
    ASSERT_FALSE(atom_list.empty());

    RHBMGroupEstimationResult result;
    result.mu_mean = (Eigen::Vector2d() << 0.1, 1.2).finished();
    result.mu_mdpde = (Eigen::Vector2d() << 0.2, 1.1).finished();
    result.mu_prior = (Eigen::Vector2d() << 0.3, 1.0).finished();
    result.capital_lambda = Eigen::Matrix2d::Identity();
    result.beta_posterior_matrix = Eigen::MatrixXd::Zero(2, static_cast<Eigen::Index>(atom_list.size()));
    result.capital_sigma_posterior_list.assign(atom_list.size(), Eigen::Matrix2d::Identity());
    result.outlier_flag_array =
        Eigen::Array<bool, Eigen::Dynamic, 1>::Constant(static_cast<Eigen::Index>(atom_list.size()), false);
    result.statistical_distance_array =
        Eigen::ArrayXd::Zero(static_cast<Eigen::Index>(atom_list.size()));
    for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(atom_list.size()); i++)
    {
        result.beta_posterior_matrix.col(i) =
            (Eigen::Vector2d() << 0.4 + 0.1 * static_cast<double>(i), 0.9).finished();
        result.outlier_flag_array(i) = (i == 0);
        result.statistical_distance_array(i) = 1.5 + static_cast<double>(i);
    }

    constexpr double alpha_g{ 0.25 };
    analysis.ApplyAtomGroupEstimateResult(group_key, class_key, result, alpha_g);

    const auto expected_mean{ ls::DecodeGroupEstimate(AtomGroupDecodeSpec(), result.mu_mean) };
    const auto expected_mdpde{ ls::DecodeGroupEstimate(AtomGroupDecodeSpec(), result.mu_mdpde) };
    const auto expected_prior{
        ls::DecodePosteriorEstimate(AtomGroupDecodeSpec(), result.mu_prior, result.capital_lambda)
    };
    EXPECT_NEAR(expected_mean.amplitude, analysis_view.GetAtomGroupMean(group_key, class_key).amplitude, 1e-12);
    EXPECT_NEAR(expected_mean.width, analysis_view.GetAtomGroupMean(group_key, class_key).width, 1e-12);
    EXPECT_NEAR(expected_mdpde.amplitude, analysis_view.GetAtomGroupMDPDE(group_key, class_key).amplitude, 1e-12);
    EXPECT_NEAR(expected_mdpde.width, analysis_view.GetAtomGroupMDPDE(group_key, class_key).width, 1e-12);
    EXPECT_NEAR(expected_prior.estimate.amplitude, analysis_view.GetAtomGroupPrior(group_key, class_key).amplitude, 1e-12);
    EXPECT_NEAR(expected_prior.estimate.width, analysis_view.GetAtomGroupPrior(group_key, class_key).width, 1e-12);
    EXPECT_NEAR(expected_prior.variance.amplitude, analysis_view.GetAtomGroupPriorPosterior(group_key, class_key).variance.amplitude, 1e-12);
    EXPECT_NEAR(expected_prior.variance.width, analysis_view.GetAtomGroupPriorPosterior(group_key, class_key).variance.width, 1e-12);
    EXPECT_DOUBLE_EQ(alpha_g, analysis_view.GetAtomAlphaG(group_key, class_key));

    const auto annotation{ rg::LocalPotentialView::RequireFor(*atom_list.front()).FindAnnotation(class_key) };
    ASSERT_TRUE(annotation.has_value());
    const auto expected_posterior{
        ls::DecodePosteriorEstimate(AtomGroupDecodeSpec(),
            result.beta_posterior_matrix.col(0),
            result.capital_sigma_posterior_list.front())
    };
    EXPECT_NEAR(expected_posterior.estimate.amplitude, annotation->posterior.estimate.amplitude, 1e-12);
    EXPECT_NEAR(expected_posterior.estimate.width, annotation->posterior.estimate.width, 1e-12);
    EXPECT_NEAR(expected_posterior.variance.amplitude, annotation->posterior.variance.amplitude, 1e-12);
    EXPECT_NEAR(expected_posterior.variance.width, annotation->posterior.variance.width, 1e-12);
    EXPECT_TRUE(annotation->is_outlier);
    EXPECT_DOUBLE_EQ(1.5, annotation->statistical_distance);
}

TEST(DataObjectModelAnalysisTest, ModelAnalysisEditorRejectsAtomGroupEstimateResultWithMismatchedMemberCount)
{
    auto model{ data_test::MakeModelWithBond() };
    model->SelectAllAtoms();
    model->ApplySymmetrySelection(false);
    auto analysis{ model->EditAnalysis() };
    analysis.RebuildAtomGroupsFromSelection();
    const auto analysis_view{ model->GetAnalysisView() };
    const auto & class_key{ ChemicalDataHelper::GetSimpleAtomClassKey() };
    const auto group_key{ analysis_view.CollectAtomGroupKeys(class_key).front() };

    RHBMGroupEstimationResult result;
    result.beta_posterior_matrix = Eigen::MatrixXd::Zero(2, 0);
    result.outlier_flag_array = Eigen::Array<bool, Eigen::Dynamic, 1>::Zero(0);
    result.statistical_distance_array = Eigen::ArrayXd::Zero(0);

    EXPECT_THROW(
        analysis.ApplyAtomGroupEstimateResult(group_key, class_key, result, 0.0),
        std::invalid_argument);
}

TEST(DataObjectModelAnalysisTest, ModelAnalysisEditorAppliesBondGroupEstimateResultToStatisticsAndAnnotations)
{
    auto model{ data_test::MakeModelWithBond() };
    model->SelectAllBonds();
    auto analysis{ model->EditAnalysis() };
    analysis.RebuildBondGroupsFromSelection();
    const auto analysis_view{ model->GetAnalysisView() };
    const auto & class_key{ ChemicalDataHelper::GetSimpleBondClassKey() };
    const auto group_keys{ analysis_view.CollectBondGroupKeys(class_key) };
    ASSERT_FALSE(group_keys.empty());

    const auto group_key{ group_keys.front() };
    const auto & bond_list{ analysis_view.GetBondObjectList(group_key, class_key) };
    ASSERT_FALSE(bond_list.empty());

    RHBMGroupEstimationResult result;
    result.mu_mean = (Eigen::Vector2d() << 0.1, 1.2).finished();
    result.mu_mdpde = (Eigen::Vector2d() << 0.2, 1.1).finished();
    result.mu_prior = (Eigen::Vector2d() << 0.3, 1.0).finished();
    result.capital_lambda = Eigen::Matrix2d::Identity();
    result.beta_posterior_matrix = Eigen::MatrixXd::Zero(2, static_cast<Eigen::Index>(bond_list.size()));
    result.capital_sigma_posterior_list.assign(bond_list.size(), Eigen::Matrix2d::Identity());
    result.outlier_flag_array =
        Eigen::Array<bool, Eigen::Dynamic, 1>::Constant(static_cast<Eigen::Index>(bond_list.size()), false);
    result.statistical_distance_array =
        Eigen::ArrayXd::Zero(static_cast<Eigen::Index>(bond_list.size()));
    for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(bond_list.size()); i++)
    {
        result.beta_posterior_matrix.col(i) =
            (Eigen::Vector2d() << 0.4 + 0.1 * static_cast<double>(i), 0.9).finished();
        result.outlier_flag_array(i) = (i == 0);
        result.statistical_distance_array(i) = 2.5 + static_cast<double>(i);
    }

    constexpr double alpha_g{ 0.5 };
    analysis.ApplyBondGroupEstimateResult(group_key, class_key, result, alpha_g);

    const auto expected_mean{ ls::DecodeGroupEstimate(BondGroupDecodeSpec(), result.mu_mean) };
    const auto expected_prior{
        ls::DecodePosteriorEstimate(BondGroupDecodeSpec(), result.mu_prior, result.capital_lambda)
    };
    EXPECT_NEAR(expected_mean.amplitude, analysis_view.GetBondGroupMean(group_key, class_key).amplitude, 1e-12);
    EXPECT_NEAR(expected_mean.width, analysis_view.GetBondGroupMean(group_key, class_key).width, 1e-12);
    EXPECT_NEAR(expected_prior.estimate.amplitude, analysis_view.GetBondGroupPrior(group_key, class_key).amplitude, 1e-12);
    EXPECT_NEAR(expected_prior.estimate.width, analysis_view.GetBondGroupPrior(group_key, class_key).width, 1e-12);
    EXPECT_NEAR(expected_prior.variance.amplitude, analysis_view.GetBondGroupPriorPosterior(group_key, class_key).variance.amplitude, 1e-12);
    EXPECT_NEAR(expected_prior.variance.width, analysis_view.GetBondGroupPriorPosterior(group_key, class_key).variance.width, 1e-12);
    EXPECT_DOUBLE_EQ(alpha_g, analysis_view.GetBondAlphaG(group_key, class_key));

    const auto annotation{ rg::LocalPotentialView::RequireFor(*bond_list.front()).FindAnnotation(class_key) };
    ASSERT_TRUE(annotation.has_value());
    EXPECT_TRUE(annotation->is_outlier);
    EXPECT_DOUBLE_EQ(2.5, annotation->statistical_distance);
}

TEST(DataObjectModelAnalysisTest, ModelAnalysisEditorRejectsBondGroupEstimateResultWithMismatchedCovarianceSize)
{
    auto model{ data_test::MakeModelWithBond() };
    model->SelectAllBonds();
    auto analysis{ model->EditAnalysis() };
    analysis.RebuildBondGroupsFromSelection();
    const auto analysis_view{ model->GetAnalysisView() };
    const auto & class_key{ ChemicalDataHelper::GetSimpleBondClassKey() };
    const auto group_key{ analysis_view.CollectBondGroupKeys(class_key).front() };
    const auto member_count{ analysis_view.GetBondObjectList(group_key, class_key).size() };

    RHBMGroupEstimationResult result;
    result.beta_posterior_matrix =
        Eigen::MatrixXd::Zero(2, static_cast<Eigen::Index>(member_count));
    result.capital_sigma_posterior_list.assign(member_count == 0 ? 0 : member_count - 1, Eigen::Matrix2d::Identity());
    result.outlier_flag_array =
        Eigen::Array<bool, Eigen::Dynamic, 1>::Constant(static_cast<Eigen::Index>(member_count), false);
    result.statistical_distance_array =
        Eigen::ArrayXd::Zero(static_cast<Eigen::Index>(member_count));

    EXPECT_THROW(
        analysis.ApplyBondGroupEstimateResult(group_key, class_key, result, 0.0),
        std::invalid_argument);
}

TEST(DataObjectModelAnalysisTest, RebuildAtomGroupEntriesFromSelectionTracksOnlySelectedAtoms)
{
    auto model{ data_test::MakeModelWithBond() };
    auto * first_atom{ model->GetAtomList().at(0).get() };
    auto * second_atom{ model->GetAtomList().at(1).get() };
    auto & analysis_data{ rg::ModelAnalysisData::Of(*model) };
    const auto & simple_class_key{ ChemicalDataHelper::GetSimpleAtomClassKey() };

    model->SelectAllAtoms(false);
    model->SetAtomSelected(first_atom->GetSerialID(), true);
    analysis_data.EnsureAtomGroupEntry("stale_atom_class").AddMember(999, *second_atom);

    analysis_data.RebuildAtomGroupEntriesFromSelection(*model);

    EXPECT_EQ(analysis_data.FindAtomGroupEntry("stale_atom_class"), nullptr);
    const auto * simple_group_entry{ analysis_data.FindAtomGroupEntry(simple_class_key) };
    ASSERT_NE(simple_group_entry, nullptr);

    size_t member_count{ 0 };
    for (const auto group_key : simple_group_entry->CollectGroupKeys())
    {
        member_count += simple_group_entry->GetMemberCount(group_key);
        for (const auto * atom : simple_group_entry->GetMembers(group_key))
        {
            ASSERT_NE(atom, nullptr);
            EXPECT_EQ(atom, first_atom);
        }
    }
    EXPECT_EQ(member_count, 1U);

    model->SetAtomSelected(first_atom->GetSerialID(), false);
    model->SetAtomSelected(second_atom->GetSerialID(), true);
    analysis_data.RebuildAtomGroupEntriesFromSelection(*model);

    simple_group_entry = analysis_data.FindAtomGroupEntry(simple_class_key);
    ASSERT_NE(simple_group_entry, nullptr);
    member_count = 0;
    for (const auto group_key : simple_group_entry->CollectGroupKeys())
    {
        member_count += simple_group_entry->GetMemberCount(group_key);
        for (const auto * atom : simple_group_entry->GetMembers(group_key))
        {
            ASSERT_NE(atom, nullptr);
            EXPECT_EQ(atom, second_atom);
        }
    }
    EXPECT_EQ(member_count, 1U);
}

TEST(DataObjectModelAnalysisTest, RebuildBondGroupEntriesFromSelectionTracksOnlySelectedBonds)
{
    auto model{ data_test::MakeModelWithBond() };
    auto * bond{ model->GetBondList().at(0).get() };
    auto & analysis_data{ rg::ModelAnalysisData::Of(*model) };
    const auto & simple_class_key{ ChemicalDataHelper::GetSimpleBondClassKey() };

    model->SelectAllBonds(false);
    analysis_data.EnsureBondGroupEntry("stale_bond_class").AddMember(888, *bond);
    analysis_data.RebuildBondGroupEntriesFromSelection(*model);

    EXPECT_EQ(analysis_data.FindBondGroupEntry("stale_bond_class"), nullptr);
    const auto * simple_group_entry{ analysis_data.FindBondGroupEntry(simple_class_key) };
    ASSERT_NE(simple_group_entry, nullptr);
    size_t member_count{ 0 };
    for (const auto group_key : simple_group_entry->CollectGroupKeys())
    {
        member_count += simple_group_entry->GetMemberCount(group_key);
    }
    EXPECT_EQ(member_count, 0U);

    model->SetBondSelected(
        bond->GetAtomSerialID1(),
        bond->GetAtomSerialID2(),
        true);
    analysis_data.RebuildBondGroupEntriesFromSelection(*model);

    simple_group_entry = analysis_data.FindBondGroupEntry(simple_class_key);
    ASSERT_NE(simple_group_entry, nullptr);
    member_count = 0;
    for (const auto group_key : simple_group_entry->CollectGroupKeys())
    {
        member_count += simple_group_entry->GetMemberCount(group_key);
        for (const auto * group_bond : simple_group_entry->GetMembers(group_key))
        {
            ASSERT_NE(group_bond, nullptr);
            EXPECT_EQ(group_bond, bond);
        }
    }
    EXPECT_EQ(member_count, 1U);
}

TEST(DataObjectModelAnalysisTest, CollectAtomGroupKeysReturnsRebuiltGroupKeySet)
{
    auto model{ data_test::MakeModelWithBond() };
    auto & analysis_data{ rg::ModelAnalysisData::Of(*model) };
    const auto & simple_class_key{ ChemicalDataHelper::GetSimpleAtomClassKey() };

    model->SelectAllAtoms();
    analysis_data.RebuildAtomGroupEntriesFromSelection(*model);

    const auto group_keys{ analysis_data.CollectAtomGroupKeys(simple_class_key) };
    const auto * group_entry{ analysis_data.FindAtomGroupEntry(simple_class_key) };
    ASSERT_NE(group_entry, nullptr);
    EXPECT_EQ(group_keys.size(), group_entry->GroupCount());
    for (const auto & group_key : group_keys)
    {
        EXPECT_TRUE(group_entry->HasGroup(group_key));
    }

    EXPECT_TRUE(analysis_data.CollectAtomGroupKeys("missing_atom_class").empty());
}

TEST(DataObjectModelAnalysisTest, CollectBondGroupKeysReturnsRebuiltGroupKeySet)
{
    auto model{ data_test::MakeModelWithBond() };
    auto & analysis_data{ rg::ModelAnalysisData::Of(*model) };
    const auto & simple_class_key{ ChemicalDataHelper::GetSimpleBondClassKey() };

    model->SelectAllBonds();
    analysis_data.RebuildBondGroupEntriesFromSelection(*model);

    const auto group_keys{ analysis_data.CollectBondGroupKeys(simple_class_key) };
    const auto * group_entry{ analysis_data.FindBondGroupEntry(simple_class_key) };
    ASSERT_NE(group_entry, nullptr);
    EXPECT_EQ(group_keys.size(), group_entry->GroupCount());
    for (const auto & group_key : group_keys)
    {
        EXPECT_TRUE(group_entry->HasGroup(group_key));
    }

    EXPECT_TRUE(analysis_data.CollectBondGroupKeys("missing_bond_class").empty());
}

TEST(DataObjectModelAnalysisTest, AtomGroupingSummaryIncludesAllConfiguredClassesInHelperOrder)
{
    auto model{ data_test::MakeModelWithBond() };
    model->SelectAllAtoms();
    model->ApplySymmetrySelection(false);
    auto analysis{ model->EditAnalysis() };
    analysis.RebuildAtomGroupsFromSelection();

    const auto summary{ model->GetAnalysisView().CollectAtomGroupingSummary() };

    EXPECT_EQ(summary.title, "Atom Grouping Summary:");
    ASSERT_EQ(summary.items.size(), ChemicalDataHelper::GetGroupAtomClassCount());
    for (size_t i = 0; i < summary.items.size(); i++)
    {
        const auto & expected_class_key{ ChemicalDataHelper::GetGroupAtomClassKey(i) };
        EXPECT_EQ(summary.items.at(i).class_key, expected_class_key);
        EXPECT_EQ(
            summary.items.at(i).group_count,
            model->GetAnalysisView().CollectAtomGroupKeys(expected_class_key).size());
    }
}

TEST(DataObjectModelAnalysisTest, BondGroupingSummaryIncludesAllConfiguredClassesInHelperOrder)
{
    auto model{ data_test::MakeModelWithBond() };
    model->SelectAllBonds();
    auto analysis{ model->EditAnalysis() };
    analysis.RebuildBondGroupsFromSelection();

    const auto summary{ model->GetAnalysisView().CollectBondGroupingSummary() };

    EXPECT_EQ(summary.title, "Bond Classification Summary:");
    ASSERT_EQ(summary.items.size(), ChemicalDataHelper::GetGroupBondClassCount());
    for (size_t i = 0; i < summary.items.size(); i++)
    {
        const auto & expected_class_key{ ChemicalDataHelper::GetGroupBondClassKey(i) };
        EXPECT_EQ(summary.items.at(i).class_key, expected_class_key);
        EXPECT_EQ(
            summary.items.at(i).group_count,
            model->GetAnalysisView().CollectBondGroupKeys(expected_class_key).size());
    }
}

TEST(DataObjectModelAnalysisTest, GroupingSummariesAreSafeBeforeRebuildAndReportZeroCounts)
{
    auto model{ data_test::MakeModelWithBond() };
    model->SelectAllAtoms(false);
    model->SelectAllBonds(false);

    const auto analysis_view{ model->GetAnalysisView() };
    const auto atom_summary{ analysis_view.CollectAtomGroupingSummary() };
    const auto bond_summary{ analysis_view.CollectBondGroupingSummary() };

    ASSERT_EQ(atom_summary.items.size(), ChemicalDataHelper::GetGroupAtomClassCount());
    for (const auto & item : atom_summary.items)
    {
        EXPECT_EQ(item.group_count, 0U);
    }

    ASSERT_EQ(bond_summary.items.size(), ChemicalDataHelper::GetGroupBondClassCount());
    for (const auto & item : bond_summary.items)
    {
        EXPECT_EQ(item.group_count, 0U);
    }
}

TEST(DataObjectModelAnalysisTest, GroupingDescriptionsMatchExpectedOutputFormat)
{
    auto model{ data_test::MakeModelWithBond() };
    model->SelectAllAtoms();
    model->SelectAllBonds();
    model->ApplySymmetrySelection(false);
    auto analysis{ model->EditAnalysis() };
    analysis.RebuildAtomGroupsFromSelection();
    analysis.RebuildBondGroupsFromSelection();

    const auto analysis_view{ model->GetAnalysisView() };
    std::string expected_atom_description{ "Atom Grouping Summary:" };
    for (size_t i = 0; i < ChemicalDataHelper::GetGroupAtomClassCount(); i++)
    {
        const auto & class_key{ ChemicalDataHelper::GetGroupAtomClassKey(i) };
        expected_atom_description +=
            "\n - Class type: " + class_key + " include "
            + std::to_string(analysis_view.CollectAtomGroupKeys(class_key).size()) + " groups.";
    }

    std::string expected_bond_description{ "Bond Classification Summary:" };
    for (size_t i = 0; i < ChemicalDataHelper::GetGroupBondClassCount(); i++)
    {
        const auto & class_key{ ChemicalDataHelper::GetGroupBondClassKey(i) };
        expected_bond_description +=
            "\n - Class type: " + class_key + " include "
            + std::to_string(analysis_view.CollectBondGroupKeys(class_key).size()) + " groups.";
    }

    EXPECT_EQ(analysis_view.DescribeAtomGrouping(), expected_atom_description);
    EXPECT_EQ(analysis_view.DescribeBondGrouping(), expected_bond_description);
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
