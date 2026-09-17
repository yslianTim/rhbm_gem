#include <gtest/gtest.h>
#include "support/FixedBOracle.hpp"
#include "support/AtomCenteredVoxelUnion.hpp"
#include <rhbm_gem/data/object/MapObject.hpp>
#include <cmath>

namespace {
namespace f=second_stage_test::matched::fixed_b;
namespace ac=second_stage_test::matched::joint_ac;
Eigen::MatrixXd Design()
{
    Eigen::MatrixXd x(8301,4);
    for (Eigen::Index k=0;k<x.rows();++k)
    {
        const double t=static_cast<double>(k)/static_cast<double>(x.rows());
        x.row(k)<<1+t,std::sin(7*t),std::cos(3*t),t*t;
    }
    return x;
}
Eigen::VectorXd Beta(const boost::json::object & fit)
{
    const auto & a=fit.at("primary").at("beta"); Eigen::VectorXd out(4);
    for (std::size_t k=0;k<4;++k) out(static_cast<Eigen::Index>(k))=boost::json::value_to<double>(a.at(k));
    return out;
}
}
TEST(FixedBOracleTest, ExactRecoveryAndNegativeChargeUseMeanOnlyCertificate)
{
    const Eigen::SparseMatrix<double> x=Design().sparseView(0,0);
    const Eigen::VectorXd truth=Eigen::Vector4d(2,-.4,1,.2),y=x*truth;
    const auto spectrum=ac::SparseSpectrum(x,Eigen::VectorXd::Ones(y.size()));
    const auto fit=f::Fit(x,y,spectrum);
    EXPECT_TRUE(fit.at("qualified").as_bool()); EXPECT_LT((Beta(fit)-truth).norm(),1e-10);
    const auto exact=f::Certificate(x,y,truth);
    EXPECT_TRUE(exact.at("kkt_passed").as_bool()); EXPECT_EQ(exact.at("rss"),0.0);
    EXPECT_EQ(exact.at("residual_scale"),0.0);
    const auto zero=f::Fit(x,Eigen::VectorXd::Zero(y.size()),spectrum);
    EXPECT_TRUE(zero.at("qualified").as_bool()); EXPECT_EQ(zero.at("primary").at("residual_scale"),0.0);
}
TEST(FixedBOracleTest, BoundaryConstraintAndIndependentDenseReference)
{
    const auto dense=Design(); const Eigen::SparseMatrix<double> x=dense.sparseView(0,0);
    const Eigen::VectorXd y=dense*Eigen::Vector4d(-.3,.7,1.2,-.2),w=Eigen::VectorXd::Ones(y.size());
    const auto fit=f::Fit(x,y,ac::SparseSpectrum(x,w));
    const auto reference=ac::WeightedSolve(dense,y,w,true);
    ASSERT_TRUE(reference.valid); EXPECT_TRUE(fit.at("qualified").as_bool());
    EXPECT_EQ(Beta(fit)(0),0); EXPECT_LT((Beta(fit)-reference.beta).norm(),1e-10);
    const auto wrong=f::Certificate(x,y,Eigen::Vector4d(0,0,0,0));
    EXPECT_FALSE(wrong.at("kkt_passed").as_bool());
    EXPECT_FALSE(f::Certificate(x,y,Eigen::Vector4d(-1,0,1,0)).at("feasible").as_bool());
}
TEST(FixedBOracleTest, RankDeficiencyCannotQualify)
{
    auto dense=Design(); dense.col(2)=dense.col(0);
    const Eigen::SparseMatrix<double> x=dense.sparseView(0,0);
    const Eigen::VectorXd y=dense*Eigen::Vector4d(1,.2,1,-.3),w=Eigen::VectorXd::Ones(y.size());
    const auto fit=f::Fit(x,y,ac::SparseSpectrum(x,w));
    EXPECT_FALSE(fit.at("qualified").as_bool());
    EXPECT_EQ(fit.at("primary").at("reason"),"rank-deficient");
}
TEST(FixedBOracleTest, EmptyStencilUnionPreservesGeometryAndWidth)
{
    namespace m=second_stage_test::matched;
    rhbm_gem::MapObject generation({11,11,11},{.5,.5,.5},{-2.5,-2.5,-2.5});
    rhbm_gem::MapObject header({11,11,11},{.50001,.50001,.50001},{-2.5001,-2.5001,-2.5001});
    header.ClearMapValueArray();
    std::vector<m::Atom> atoms{{{0,0,0},0,.500364,0},{{.5,0,0},0,.4992,0}};
    const auto grid=m::atom_union::BuildGrid(atoms,{},generation,header);
    const auto first=m::atom_union::BuildDesign(grid,atoms,generation);
    const auto repeated=m::atom_union::BuildGrid(atoms,{},generation,header);
    EXPECT_EQ(grid.voxels.size(),repeated.voxels.size()); EXPECT_TRUE(grid.sample_rows.empty());
    const double saved=atoms[0].width; atoms[0].width=.5;
    const auto second=m::atom_union::BuildDesign(grid,atoms,generation);
    EXPECT_NE(saved,atoms[0].width); EXPECT_GT((first-second).norm(),0);
    bool boundary=false;
    for (std::size_t k=0;k<grid.voxels.size();++k)
    {
        EXPECT_EQ(grid.voxels[k].index,repeated.voxels[k].index);
        EXPECT_EQ(grid.voxels[k].position,generation.GetGridPosition(grid.voxels[k].index));
        if (k) EXPECT_LT(grid.voxels[k-1].index,grid.voxels[k].index);
        boundary |= grid.voxels[k].nearest_distance==2.5;
    }
    EXPECT_TRUE(boundary);
}
