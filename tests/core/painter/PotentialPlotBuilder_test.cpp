#include <gtest/gtest.h>

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "support/CommandTestHelpers.hpp"
#include "painter/detail/PotentialPlotBuilder.hpp"
#include <rhbm_gem/data/io/ModelMapFileIO.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <rhbm_gem/data/object/ModelAnalysisView.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/domain/ChemicalDataHelper.hpp>
#include <rhbm_gem/utils/domain/KeyPacker.hpp>

#ifdef HAVE_ROOT
#include <TF1.h>
#include <TGraphErrors.h>
#include <TH1.h>
#endif

namespace rg = rhbm_gem;

namespace {

std::shared_ptr<rg::ModelObject> LoadModelFixture(const std::string & fixture_name)
{
    auto model{ rg::ReadModel(command_test::TestDataPath(fixture_name)) };
    model->SetKeyTag("model");
    return std::shared_ptr<rg::ModelObject>{ std::move(model) };
}

void EnsureLocalPotentialEntries(rg::ModelObject & model)
{
    auto analysis{ model.EditAnalysis() };
    for (auto & atom : model.GetAtomList())
    {
        if (!rg::AtomLocalPotentialView::For(*atom).IsAvailable())
        {
            analysis.EnsureAtomLocalPotential(*atom);
        }
    }
}

} // namespace

TEST(PotentialPlotBuilderTest, ConstructorsKeepQueryObjectsReachable)
{
    auto model{ LoadModelFixture("test_model_auth_seq_alnum_struct_conn.cif") };
    ASSERT_NE(model, nullptr);
    EnsureLocalPotentialEntries(*model);

    ASSERT_FALSE(model->GetAtomList().empty());

    auto * atom{ model->GetAtomList().front().get() };

    rg::PotentialPlotBuilder model_builder{ model.get() };
    rg::PotentialPlotBuilder atom_builder{ atom };

    (void)model_builder;
    (void)atom_builder;
    SUCCEED();
}

#ifdef HAVE_ROOT
TEST(PotentialPlotBuilderTest, RepresentativeBuildersProduceRootObjects)
{
    auto model{ LoadModelFixture("test_model_auth_seq_alnum_struct_conn.cif") };
    ASSERT_NE(model, nullptr);
    EnsureLocalPotentialEntries(*model);
    ASSERT_FALSE(model->GetAtomList().empty());

    auto * atom{ model->GetAtomList().front().get() };

    rg::PotentialPlotBuilder model_builder{ model.get() };
    rg::PotentialPlotBuilder atom_builder{ atom };

    auto tomography_graph{ model_builder.CreateAtomXYPositionTomographyGraph() };
    auto distance_graph{ atom_builder.CreateDistanceToMapValueGraph() };
    auto gaus_function{ atom_builder.CreateAtomLocalGausFunctionMDPDE() };

    EXPECT_NE(tomography_graph, nullptr);
    EXPECT_NE(distance_graph, nullptr);
    EXPECT_NE(gaus_function, nullptr);
}

TEST(PotentialPlotBuilderTest, LinearModelDataBuildersUseFiniteSamplingRange)
{
    auto model{ LoadModelFixture("test_model_auth_seq_alnum_struct_conn.cif") };
    ASSERT_NE(model, nullptr);
    ASSERT_FALSE(model->GetAtomList().empty());

    auto * atom{ model->GetAtomList().front().get() };
    auto analysis{ model->EditAnalysis() };
    analysis.EnsureAtomLocalPotential(*atom).SetRawSamplingEntries({
        {4.0f, SamplingPoint{ 0.0f }},
        {3.0f, SamplingPoint{ 0.4f }},
        {2.0f, SamplingPoint{ 0.8f }},
    });

    rg::PotentialPlotBuilder atom_builder{ atom };

    auto basis_hist{ atom_builder.CreateLinearModelDataHistogram(0) };
    auto response_hist{ atom_builder.CreateLinearModelDataHistogram(1) };
    auto graph{ atom_builder.CreateLinearModelDistanceToMapValueGraph() };

    EXPECT_NE(basis_hist, nullptr);
    EXPECT_NE(response_hist, nullptr);
    EXPECT_NE(graph, nullptr);
    EXPECT_EQ(graph->GetN(), 3);
}

TEST(PotentialPlotBuilderTest, GaussianFunctionsPreserveModelOffset)
{
    auto model{ LoadModelFixture("test_model_auth_seq_alnum_struct_conn.cif") };
    ASSERT_NE(model, nullptr);
    model->SelectAllAtoms();
    model->ApplySymmetrySelection(false);
    EnsureLocalPotentialEntries(*model);
    ASSERT_FALSE(model->GetSelectedAtoms().empty());

    auto analysis{ model->EditAnalysis() };
    auto * atom{ model->GetSelectedAtoms().front() };
    const rg::GaussianModel3D local_model{ 1.25, 0.55, -0.15 };
    rg::LocalGaussianResult local_result;
    local_result.mdpde = rg::GaussianModel3DWithUncertainty{
        local_model,
        rg::GaussianModel3DUncertainty{}
    };
    analysis.EnsureAtomLocalPotential(*atom).SetGaussianResult(
        rg::LocalFittingStage::Third,
        local_result);

    rg::PotentialPlotBuilder atom_builder{ atom };
    auto local_function{ atom_builder.CreateAtomLocalGausFunctionMDPDE() };
    ASSERT_NE(local_function, nullptr);
    EXPECT_EQ(local_function->GetNpar(), 3);
    EXPECT_DOUBLE_EQ(local_function->GetParameter(2), local_model.GetOffset());
    EXPECT_NEAR(local_function->Eval(0.75), local_model.ResponseAtDistance(0.75), 1e-12);
    EXPECT_NE(local_function->Eval(0.75), local_model.SignalAtDistance(0.75));

    analysis.RebuildAtomGroupsFromSelection();
    const auto analysis_view{ model->GetAnalysisView() };
    const auto group_keys{ analysis_view.CollectAtomGroupKeys() };
    ASSERT_FALSE(group_keys.empty());
    const auto group_key{ group_keys.front() };
    const auto & atom_list{ analysis_view.GetAtomObjectList(group_key) };
    ASSERT_FALSE(atom_list.empty());

    const rg::GaussianModel3D mean_model{ 1.10, 0.70, 0.21 };
    const rg::GaussianModel3D prior_model{ 1.35, 0.65, 0.32 };
    rg::GroupGaussianResult group_result;
    group_result.mean = mean_model;
    group_result.prior = rg::GaussianModel3DWithUncertainty{
        prior_model,
        rg::GaussianModel3DUncertainty{}
    };
    group_result.member_results.resize(atom_list.size());
    analysis.ApplyAtomGroupGaussianResult(
        rg::LocalFittingStage::Third,
        group_key,
        group_result);

    rg::PotentialPlotBuilder model_builder{ model.get() };
    auto mean_function{ model_builder.CreateAtomGroupGausFunctionMean(group_key) };
    auto prior_function{ model_builder.CreateAtomGroupGausFunctionPrior(group_key) };

    ASSERT_NE(mean_function, nullptr);
    EXPECT_EQ(mean_function->GetNpar(), 3);
    EXPECT_DOUBLE_EQ(mean_function->GetParameter(2), mean_model.GetOffset());
    EXPECT_NEAR(mean_function->Eval(0.75), mean_model.ResponseAtDistance(0.75), 1e-12);
    EXPECT_NE(mean_function->Eval(0.75), mean_model.SignalAtDistance(0.75));

    ASSERT_NE(prior_function, nullptr);
    EXPECT_EQ(prior_function->GetNpar(), 3);
    EXPECT_DOUBLE_EQ(prior_function->GetParameter(2), prior_model.GetOffset());
    EXPECT_NEAR(prior_function->Eval(0.75), prior_model.ResponseAtDistance(0.75), 1e-12);
    EXPECT_NE(prior_function->Eval(0.75), prior_model.SignalAtDistance(0.75));
}

TEST(PotentialPlotBuilderTest, ComponentAtomAveragePriorUsesEqualComponentWeight)
{
    auto model{ LoadModelFixture("test_model_auth_seq_alnum_struct_conn.cif") };
    ASSERT_NE(model, nullptr);
    ASSERT_GE(model->GetAtomList().size(), 2u);
    const auto shared_atom_key{ static_cast<AtomKey>(Spot::CA) };
    model->GetAtomList().at(0)->SetAtomKey(shared_atom_key);
    model->GetAtomList().at(0)->SetComponentKey(1);
    model->GetAtomList().at(1)->SetAtomKey(shared_atom_key);
    model->GetAtomList().at(1)->SetComponentKey(2);
    model->SelectAllAtoms();
    model->ApplySymmetrySelection(false);
    EnsureLocalPotentialEntries(*model);

    auto analysis{ model->EditAnalysis() };
    analysis.RebuildAtomGroupsFromSelection();
    const auto analysis_view{ model->GetAnalysisView() };

    std::map<AtomKey, std::vector<GroupKey>> groups_by_atom_key;
    for (const auto group_key : analysis_view.CollectAtomGroupKeys())
    {
        const auto unpacked_key{ KeyPackerComponentAtomClass::Unpack(group_key) };
        groups_by_atom_key[std::get<1>(unpacked_key)].emplace_back(group_key);
    }

    auto target_iter{ groups_by_atom_key.end() };
    for (auto iter{ groups_by_atom_key.begin() }; iter != groups_by_atom_key.end(); ++iter)
    {
        if (iter->second.size() >= 2)
        {
            target_iter = iter;
            break;
        }
    }
    ASSERT_NE(target_iter, groups_by_atom_key.end());

    const auto atom_key{ target_iter->first };
    const auto first_group_key{ target_iter->second.at(0) };
    const auto second_group_key{ target_iter->second.at(1) };

    rg::GroupGaussianResult first_result;
    first_result.prior = rg::GaussianModel3DWithUncertainty{
        rg::GaussianModel3D{ 1.0, 2.0, 0.1 },
        rg::GaussianModel3DUncertainty{ 0.2, 0.4, 0.06 }
    };
    first_result.member_results.resize(
        analysis_view.GetAtomObjectList(first_group_key).size());
    analysis.ApplyAtomGroupGaussianResult(
        rg::LocalFittingStage::Third,
        first_group_key,
        first_result);

    rg::GroupGaussianResult second_result;
    second_result.prior = rg::GaussianModel3DWithUncertainty{
        rg::GaussianModel3D{ 3.0, 4.0, 0.3 },
        rg::GaussianModel3DUncertainty{ 0.6, 0.8, 0.10 }
    };
    second_result.member_results.resize(
        analysis_view.GetAtomObjectList(second_group_key).size());
    analysis.ApplyAtomGroupGaussianResult(
        rg::LocalFittingStage::Third,
        second_group_key,
        second_result);

    const auto average_prior{
        rg::PotentialPlotBuilder::ComputeComponentAtomAveragePrior(analysis_view, atom_key)
    };
    ASSERT_TRUE(average_prior.has_value());
    EXPECT_DOUBLE_EQ(2.0, average_prior->GetModel().GetAmplitude());
    EXPECT_DOUBLE_EQ(3.0, average_prior->GetModel().GetWidth());
    EXPECT_DOUBLE_EQ(0.2, average_prior->GetModel().GetOffset());
    EXPECT_DOUBLE_EQ(0.4, average_prior->GetStandardDeviationModel().GetAmplitude());
    EXPECT_DOUBLE_EQ(0.6, average_prior->GetStandardDeviationModel().GetWidth());
    EXPECT_DOUBLE_EQ(0.08, average_prior->GetStandardDeviationModel().GetOffset());

    rg::PotentialPlotBuilder model_builder{ model.get() };
    auto prior_function{ model_builder.CreateComponentAtomAverageGausFunctionPrior(atom_key) };
    ASSERT_NE(prior_function, nullptr);
    EXPECT_DOUBLE_EQ(0.2, prior_function->GetParameter(2));
    EXPECT_NEAR(
        prior_function->Eval(0.75),
        average_prior->GetModel().ResponseAtDistance(0.75),
        1e-12);
}

TEST(PotentialPlotBuilderTest, MissingComponentAtomAverageReturnsEmptyOutputs)
{
    auto model{ LoadModelFixture("test_model_auth_seq_alnum_struct_conn.cif") };
    ASSERT_NE(model, nullptr);
    model->SelectAllAtoms();
    model->ApplySymmetrySelection(false);
    EnsureLocalPotentialEntries(*model);
    model->EditAnalysis().RebuildAtomGroupsFromSelection();

    const auto missing_atom_key{ static_cast<AtomKey>(65535) };
    const auto analysis_view{ model->GetAnalysisView() };

    EXPECT_FALSE(
        rg::PotentialPlotBuilder::ComputeComponentAtomAveragePrior(
            analysis_view, missing_atom_key).has_value());

    rg::PotentialPlotBuilder model_builder{ model.get() };
    EXPECT_EQ(model_builder.CreateComponentAtomAverageGausFunctionPrior(missing_atom_key), nullptr);

    auto scatter_graph{
        rg::PotentialPlotBuilder::CreateMapValueScatterGraph(
            missing_atom_key, model.get(), model.get())
    };
    ASSERT_NE(scatter_graph, nullptr);
    EXPECT_EQ(scatter_graph->GetN(), 0);
}
#endif
