#include <gtest/gtest.h>

#include <memory>
#include <string>

#include "support/CommandTestHelpers.hpp"
#include "painter/detail/PotentialPlotBuilder.hpp"
#include <rhbm_gem/data/io/ModelMapFileIO.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <rhbm_gem/data/object/ModelAnalysisView.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>

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
#endif
