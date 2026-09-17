#include <gtest/gtest.h>
#include "support/AtomCenteredVoxelUnion.hpp"
#include <rhbm_gem/data/object/MapObject.hpp>
#include <cmath>
#include <set>

namespace {
namespace m=second_stage_test::matched;
namespace u=m::atom_union;
namespace ac=m::joint_ac;
Eigen::VectorXd V(const boost::json::value & a)
{
    Eigen::VectorXd v(static_cast<Eigen::Index>(a.as_array().size()));
    for (Eigen::Index k=0;k<v.size();++k) v(k)=boost::json::value_to<double>(a.at(static_cast<std::size_t>(k)));
    return v;
}
}
TEST(AtomCenteredVoxelUnionTest, ExactSphereBoundaryOverlapAndClippedEdges)
{
    rhbm_gem::MapObject map({7,8,9},{.5,.75,1.0},{-1.5,-2.25,-4});
    auto values=std::make_unique<double[]>(map.GetMapValueArraySize());
    for (std::size_t k=0;k<map.GetMapValueArraySize();++k) values[k]=k%2 ? -1 : 0;
    map.SetMapValueArray(std::move(values));
    const std::vector<m::Atom> atoms{{{0,0,0},2,.5,-.2},{{.5,0,0},1,.6,.1},{{-1.5,-2.25,-4},1,.4,0}};
    const auto grid=u::BuildGrid(atoms,{},map,map,1.5); std::size_t p{}; bool boundary{},overlap{};
    for (std::size_t i=0;i<map.GetMapValueArraySize();++i)
    {
        std::size_t count{}; double nearest=100;
        for (const auto & a:atoms)
        {
            const double d=m::SquareDistance(map.GetGridPosition(i),a.position);
            if (d<=2.25) {++count; nearest=std::min(nearest,d); boundary|=d==2.25;}
        }
        if (!count) continue;
        ASSERT_LT(p,grid.voxels.size()); const auto & v=grid.voxels[p++];
        EXPECT_EQ(v.index,i); EXPECT_EQ(v.multiplicity,count); EXPECT_DOUBLE_EQ(v.nearest_distance,std::sqrt(nearest));
        EXPECT_EQ(v.observed,map.GetMapValue(i)); overlap|=count>1;
    }
    EXPECT_EQ(p,grid.voxels.size()); EXPECT_TRUE(overlap); EXPECT_TRUE(boundary);
}
TEST(AtomCenteredVoxelUnionTest, HeaderGeometryAndSparseDesignReplay)
{
    rhbm_gem::MapObject generation({9,9,9},{.2,.3,.4},{-.8,-1.2,-1.6});
    rhbm_gem::MapObject observed({9,9,9},{.2001,.3001,.4001},{-.8001,-1.2001,-1.6001}); observed.ClearMapValueArray();
    const std::vector<m::Atom> atoms{{{0,0,0},2,.5,-.2},{{.4,.3,.4},1,.6,0}};
    const std::vector<m::Stencil> stencils{m::MakeStencil(generation,observed,{0,0,0})};
    const auto grid=u::BuildGrid(atoms,stencils,generation,observed);
    const auto sparse=u::BuildDesign(grid,atoms,generation);
    const auto dense=m::unique_grid::BuildDesign(grid,atoms,2.5);
    EXPECT_EQ(sparse.nonZeros(),(dense.array()!=0).count()); EXPECT_TRUE(Eigen::MatrixXd(sparse).isApprox(dense,0));
    for (const auto & v:grid.voxels) EXPECT_EQ(v.position,generation.GetGridPosition(v.index));
    EXPECT_THROW(u::BuildGrid(atoms,stencils,generation,observed,.01),std::invalid_argument);
}
TEST(AtomCenteredVoxelUnionTest, IndependentTSQRAndConstraintsMatchDenseSVD)
{
    Eigen::MatrixXd x=Eigen::MatrixXd::Zero(8301,6); Eigen::VectorXd y(x.rows()),w(x.rows());
    for (Eigen::Index i=0;i<x.rows();++i)
    {
        const double t=static_cast<double>(i)/static_cast<double>(x.rows());
        x.row(i)<<1,t,std::sin(7*t),std::cos(6*t),std::cos(13*t),t*t;
        if (i%3==0) x(i,2)=0;
        y(i)=-.3*x(i,0)+.7*x(i,1)+1.2*x(i,2)-.2*x(i,3)+.8*x(i,4)+.01*std::sin(3.1*static_cast<double>(i));
        w(i)=i%7 ? .3+t : 0;
    }
    const Eigen::SparseMatrix<double> sparse=x.sparseView(0,0);
    const auto reference=ac::WeightedSolve(x,y,w,true),actual=ac::WeightedSolve(sparse,y,w,true),primary=ac::WeightedSolve(sparse,y,w);
    ASSERT_TRUE(reference.valid); ASSERT_TRUE(actual.valid); ASSERT_TRUE(primary.valid);
    EXPECT_LT((actual.beta-reference.beta).norm(),1e-10); EXPECT_LT((primary.beta-reference.beta).norm(),1e-10);
    EXPECT_EQ(actual.beta(0),0); EXPECT_EQ(actual.rank,reference.rank);
    const auto spectrum=ac::SparseSpectrum(sparse,w); EXPECT_EQ(spectrum.at("rank"),6);
    Eigen::MatrixXd deficient=x; deficient.col(4)=deficient.col(0);
    const Eigen::SparseMatrix<double> bad=deficient.sparseView(0,0);
    EXPECT_FALSE(ac::WeightedSolve(bad,y,w,true).valid);
}
TEST(AtomCenteredVoxelUnionTest, AllAlphaBranchesAndQualificationMatchDense)
{
    Eigen::MatrixXd x(180,4); Eigen::VectorXd y(180);
    for (Eigen::Index i=0;i<x.rows();++i)
    {
        const double t=static_cast<double>(i)/180; x.row(i)<<1+t,std::sin(7*t),std::cos(4*t),t*t;
        y(i)=x.row(i).dot(Eigen::Vector4d(2,-.4,1,.2))+.05*std::sin(3.1*static_cast<double>(i));
    }
    const Eigen::SparseMatrix<double> sparse=x.sparseView(0,0); const Eigen::VectorXd initial=Eigen::Vector4d(2,-.4,1,.2);
    for (double alpha:{0.,.1,.5,1.})
    {
        const auto dense=m::unique_grid::Fit(x,y,initial,alpha),actual=m::unique_grid::Fit(sparse,y,initial,alpha);
        EXPECT_EQ(actual.at("reason"),dense.at("reason")); EXPECT_EQ(actual.at("qualified"),dense.at("qualified"));
        EXPECT_LT((V(actual.at("beta"))-V(dense.at("beta"))).norm(),1e-8);
        for (const auto & b:actual.at("branches").as_array())
        {
            EXPECT_LE(b.at("primary").at("iterations").as_int64(),100);
            EXPECT_LE(b.at("reference").at("iterations").as_int64(),100);
        }
        const auto blocks=m::unique_grid::GlobalBlock(y.size(),alpha);
        const Eigen::VectorXd variance=Eigen::VectorXd::Constant(1,.01);
        const auto a=ac::Evaluate(x,y,initial,variance,blocks),b=ac::Evaluate(sparse,y,initial,variance,blocks);
        EXPECT_NEAR(a.objective,b.objective,1e-12); EXPECT_NEAR(a.stationarity,b.stationarity,1e-12);
    }
    const auto exact=ac::Evaluate(sparse,Eigen::VectorXd(sparse*initial),initial,Eigen::VectorXd::Ones(1),m::unique_grid::GlobalBlock(y.size(),.1));
    EXPECT_EQ(exact.reason,"exact-fit-boundary");
    const auto invalid=ac::Evaluate(sparse,y,initial,Eigen::VectorXd::Constant(1,-1),m::unique_grid::GlobalBlock(y.size(),.1));
    EXPECT_EQ(invalid.reason,"variance-boundary");
    const auto tiny=ac::Evaluate(sparse,y,initial,Eigen::VectorXd::Constant(1,1e-30),m::unique_grid::GlobalBlock(y.size(),1));
    EXPECT_EQ(tiny.reason,"invalid-denominator");
}
