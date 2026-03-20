#include <gtest/gtest.h>

#include <memory>
#include <string>

#include "support/CommandTestHelpers.hpp"
#include "PotentialPlotBuilder.hpp"
#include <rhbm_gem/data/io/DataObjectManager.hpp>
#include <rhbm_gem/data/object/BondObject.hpp>
#include <rhbm_gem/data/object/LocalPotentialEntry.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>

#ifdef HAVE_ROOT
#include <TF1.h>
#include <TGraphErrors.h>
#endif

namespace rg = rhbm_gem;

namespace {

std::shared_ptr<rg::ModelObject> LoadModelFixture(const std::string & fixture_name)
{
    rg::DataObjectManager manager{};
    manager.ProcessFile(command_test::TestDataPath(fixture_name), "model");
    return manager.GetTypedDataObject<rg::ModelObject>("model");
}

void EnsureLocalPotentialEntries(rg::ModelObject & model)
{
    for (auto & atom : model.GetAtomList())
    {
        if (atom->GetLocalPotentialEntry() == nullptr)
        {
            atom->AddLocalPotentialEntry(std::make_unique<rg::LocalPotentialEntry>());
        }
    }
    for (auto & bond : model.GetBondList())
    {
        if (bond->GetLocalPotentialEntry() == nullptr)
        {
            bond->AddLocalPotentialEntry(std::make_unique<rg::LocalPotentialEntry>());
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
    ASSERT_FALSE(model->GetBondList().empty());

    auto * atom{ model->GetAtomList().front().get() };
    auto * bond{ model->GetBondList().front().get() };

    rg::PotentialPlotBuilder model_builder{ model.get() };
    rg::PotentialPlotBuilder atom_builder{ atom };
    rg::PotentialPlotBuilder bond_builder{ bond };

    (void)model_builder;
    (void)atom_builder;
    (void)bond_builder;
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
    rg::PotentialPlotBuilder bond_builder{ model->GetBondList().front().get() };

    auto tomography_graph{ model_builder.CreateAtomXYPositionTomographyGraph() };
    auto distance_graph{ atom_builder.CreateDistanceToMapValueGraph() };
    auto bond_distance_graph{ bond_builder.CreateDistanceToMapValueGraph() };
    auto gaus_function{ atom_builder.CreateAtomLocalGausFunctionMDPDE() };

    EXPECT_NE(tomography_graph, nullptr);
    EXPECT_NE(distance_graph, nullptr);
    EXPECT_NE(bond_distance_graph, nullptr);
    EXPECT_NE(gaus_function, nullptr);
}
#endif
