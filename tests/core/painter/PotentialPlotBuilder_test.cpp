#include <gtest/gtest.h>

#include <memory>
#include <string>

#include "support/CommandTestHelpers.hpp"
#include "painter/detail/PotentialPlotBuilder.hpp"
#include <rhbm_gem/data/io/ModelMapFileIO.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <rhbm_gem/data/object/ModelAnalysisView.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/domain/ChemicalDataHelper.hpp>

#ifdef HAVE_ROOT
#include <TF1.h>
#include <TGraphErrors.h>
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

TEST(PotentialPlotBuilderTest, GaussianFunctionsPreserveModelIntercept)
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
    analysis.EnsureAtomLocalPotential(*atom).SetGaussianResult(local_result);

    rg::PotentialPlotBuilder atom_builder{ atom };
    auto local_function{ atom_builder.CreateAtomLocalGausFunctionMDPDE() };
    ASSERT_NE(local_function, nullptr);
    EXPECT_EQ(local_function->GetNpar(), 3);
    EXPECT_DOUBLE_EQ(local_function->GetParameter(2), local_model.GetIntercept());
    EXPECT_NEAR(local_function->Eval(0.75), local_model.ResponseAtDistance(0.75), 1e-12);

    analysis.RebuildAtomGroupsFromSelection();
    const auto analysis_view{ model->GetAnalysisView() };
    const auto & class_key{ ChemicalDataHelper::GetSimpleAtomClassKey() };
    const auto group_keys{ analysis_view.CollectAtomGroupKeys(class_key) };
    ASSERT_FALSE(group_keys.empty());
    const auto group_key{ group_keys.front() };
    const auto & atom_list{ analysis_view.GetAtomObjectList(group_key, class_key) };
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
    analysis.ApplyAtomGroupGaussianResult(group_key, class_key, group_result);

    rg::PotentialPlotBuilder model_builder{ model.get() };
    auto mean_function{ model_builder.CreateAtomGroupGausFunctionMean(group_key, class_key) };
    auto prior_function{ model_builder.CreateAtomGroupGausFunctionPrior(group_key, class_key) };

    ASSERT_NE(mean_function, nullptr);
    EXPECT_EQ(mean_function->GetNpar(), 3);
    EXPECT_DOUBLE_EQ(mean_function->GetParameter(2), mean_model.GetIntercept());
    EXPECT_NEAR(mean_function->Eval(0.75), mean_model.ResponseAtDistance(0.75), 1e-12);

    ASSERT_NE(prior_function, nullptr);
    EXPECT_EQ(prior_function->GetNpar(), 3);
    EXPECT_DOUBLE_EQ(prior_function->GetParameter(2), prior_model.GetIntercept());
    EXPECT_NEAR(prior_function->Eval(0.75), prior_model.ResponseAtDistance(0.75), 1e-12);
}
#endif
