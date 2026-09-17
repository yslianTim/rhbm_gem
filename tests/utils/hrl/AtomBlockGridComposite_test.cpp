#include <gtest/gtest.h>
#include "support/AtomBlockGridComposite.hpp"
#include <rhbm_gem/data/object/MapObject.hpp>
#include <cmath>
#include <numeric>

namespace {
namespace m=second_stage_test::matched;
namespace ac=m::joint_ac;
using Sparse=Eigen::SparseMatrix<double>;
struct Fixture
{
    Eigen::MatrixXd x{128,4}; Eigen::VectorXd y{128},beta{4},v{2}; ac::Blocks blocks{{0,0,{}},{1,.5,{}}};
    Fixture()
    {
        beta<<1,.2,.7,-.1; v<<.02,.03;
        for (Eigen::Index p=0;p<x.rows();++p)
        {
            const double t=static_cast<double>(p)/128.; x.row(p)<<1,t,std::sin(9*t),std::cos(7*t);
            y(p)=x.row(p).dot(beta)+.11*std::sin(1.7*static_cast<double>(p));
            if (p<90) blocks[0].rows.push_back(p);
            if (p>=60) blocks[1].rows.push_back(p);
        }
    }
};
Eigen::VectorXd V(const boost::json::value & a)
{
    Eigen::VectorXd out(static_cast<Eigen::Index>(a.as_array().size()));
    for (Eigen::Index i=0;i<out.size();++i) out(i)=boost::json::value_to<double>(a.at(static_cast<std::size_t>(i)));
    return out;
}
}
TEST(AtomBlockGridCompositeTest, ExplicitExpansionMatchesObjectiveEquationsUpdatesAndMass)
{
    Fixture f; const Sparse x=f.x.sparseView(0,0);
    for (double alpha:{0.,.1,.5,1.})
    {
        f.blocks[1].alpha=alpha;
        const auto e=ac::EvaluateComposite(x,f.y,f.beta,f.v,f.blocks); ASSERT_TRUE(e.valid)<<e.reason;
        Eigen::MatrixXd expanded(158,4); Eigen::VectorXd y(158); ac::Blocks blocks=f.blocks; Eigen::Index cursor{};
        for (auto & block:blocks) for (auto & row:block.rows) {expanded.row(cursor)=f.x.row(row); y(cursor)=f.y(row); row=cursor++;}
        const auto reference=ac::Evaluate(expanded,y,f.beta,f.v,blocks); ASSERT_TRUE(reference.valid);
        EXPECT_NEAR(e.objective,reference.objective,1e-13);
        EXPECT_LT((e.scaled-reference.scaled).norm(),1e-13);
        EXPECT_LT((e.membership_weights-reference.weights).norm(),16*std::numeric_limits<double>::epsilon()*reference.weights.norm());
        EXPECT_LT((e.denominators-reference.denominators).norm(),1e-13);
        Eigen::VectorXd omega=Eigen::VectorXd::Zero(128); cursor=0; double mass{};
        for (std::size_t i=0;i<f.blocks.size();++i) for (auto row:f.blocks[i].rows)
        {omega(row)+=reference.linear_weights(cursor++); mass+=e.block_prefactors(static_cast<Eigen::Index>(i))*e.membership_weights(cursor-1);}
        EXPECT_LT((omega-e.linear_weights).norm(),1e-13); EXPECT_NEAR(mass/e.linear_weights.sum(),1,1e-14);
        const auto unique=ac::WeightedSolve(x,f.y,e.linear_weights,true);
        // Use the unique design's column scaling and original row-count threshold in the expanded reference.
        Eigen::VectorXd scales=f.x.colwise().norm();
        Eigen::MatrixXd z=reference.linear_weights.cwiseSqrt().asDiagonal()*expanded*scales.cwiseInverse().asDiagonal();
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(z,Eigen::ComputeThinU|Eigen::ComputeThinV);
        svd.setThreshold(std::numeric_limits<double>::epsilon()*128.);
        const Eigen::VectorXd b=svd.solve(reference.linear_weights.cwiseSqrt().cwiseProduct(y)).cwiseQuotient(scales);
        ASSERT_TRUE(unique.valid); ASSERT_GT(b(0),0); ASSERT_GT(b(2),0); EXPECT_LT((unique.beta-b).norm(),1e-12);
        const Eigen::VectorXd residual=f.y-f.x*unique.beta, full=y-expanded*unique.beta;
        cursor=0;
        for (std::size_t i=0;i<f.blocks.size();++i)
        {
            double a{},bvar{};
            for (auto row:f.blocks[i].rows)
            {a+=e.membership_weights(cursor)*residual(row)*residual(row); bvar+=reference.weights(cursor)*full(cursor)*full(cursor); ++cursor;}
            EXPECT_NEAR(a/e.denominators(static_cast<Eigen::Index>(i)),bvar/reference.denominators(static_cast<Eigen::Index>(i)),1e-14);
        }
    }
}
TEST(AtomBlockGridCompositeTest, PartitionEntryRejectsOverlapAndCompositeRejectsWithinBlockDuplicates)
{
    Fixture f; const Sparse x=f.x.sparseView(0,0);
    EXPECT_FALSE(ac::Evaluate(x,f.y,f.beta,f.v,f.blocks).valid);
    auto duplicate=f.blocks; duplicate[0].rows.push_back(duplicate[0].rows[0]);
    EXPECT_EQ(ac::EvaluateComposite(x,f.y,f.beta,f.v,duplicate).reason,"invalid-input");
    auto missing=f.blocks; missing[0].rows.erase(missing[0].rows.begin());
    EXPECT_EQ(ac::EvaluateComposite(x,f.y,f.beta,f.v,missing).reason,"invalid-input");
}
TEST(AtomBlockGridCompositeTest, DisjointAndGlobalBlocksKeepLegacyEndpoints)
{
    Fixture f; const Sparse x=f.x.sparseView(0,0);
    f.blocks[0].rows.resize(60);
    for (double alpha:{0.,.1,.5,1.})
    {
        f.blocks[0].alpha=f.blocks[1].alpha=alpha;
        for (const auto & blocks:{f.blocks,m::unique_grid::GlobalBlock(128,alpha)})
        {
            const auto old=ac::Fit(x,f.y,f.beta,blocks,100,100),fresh=ac::FitComposite(x,f.y,f.beta,blocks);
            EXPECT_EQ(old.at("reason"),fresh.at("reason")); EXPECT_EQ(old.at("beta"),fresh.at("beta"));
            EXPECT_EQ(old.at("variances"),fresh.at("variances")); EXPECT_EQ(old.at("stationarity"),fresh.at("stationarity"));
            for (std::size_t k=0;k<2;++k)
            {
                const auto & branch=fresh.at("branches").at(k);
                EXPECT_LE(branch.at("primary").at("iterations").as_int64(),100);
                EXPECT_LE(branch.at("reference").at("iterations").as_int64(),100);
            }
        }
    }
}
TEST(AtomBlockGridCompositeTest, CompositeGradientMatchesObjectiveAndAlphaZeroIsNotOrdinaryLS)
{
    Fixture f; const Sparse x=f.x.sparseView(0,0); auto e=ac::EvaluateComposite(x,f.y,f.beta,f.v,f.blocks);
    const Eigen::VectorXd derivative=-std::exp(e.log_prefactors.maxCoeff())*(x.transpose()*(e.linear_weights.array()*(f.y-x*f.beta).array()).matrix());
    for (Eigen::Index k=0;k<4;++k)
    {
        auto lo=f.beta,hi=f.beta; lo(k)-=1e-6; hi(k)+=1e-6;
        const double delta=(ac::EvaluateComposite(x,f.y,hi,f.v,f.blocks).objective-ac::EvaluateComposite(x,f.y,lo,f.v,f.blocks).objective)/2e-6;
        EXPECT_NEAR(delta,derivative(k),1e-8);
    }
    f.blocks[1].alpha=0; f.v<<.01,1.; e=ac::EvaluateComposite(x,f.y,f.beta,f.v,f.blocks);
    EXPECT_GT(e.linear_weights.maxCoeff()/e.linear_weights.minCoeff(),50);
    const auto composite=ac::WeightedSolve(x,f.y,e.linear_weights),ls=ac::WeightedSolve(x,f.y,Eigen::VectorXd::Ones(128));
    EXPECT_GT((composite.beta-ls.beta).norm(),1e-4);
}
TEST(AtomBlockGridCompositeTest, BoundaryDenominatorRankConstraintsAndExhaustion)
{
    Fixture f; const Sparse x=f.x.sparseView(0,0);
    auto y=f.y; for (auto row:f.blocks[0].rows) y(row)=f.x.row(row).dot(f.beta);
    const auto boundary=ac::EvaluateComposite(x,y,f.beta,f.v,f.blocks);
    EXPECT_EQ(boundary.reason,"exact-fit-boundary"); EXPECT_EQ(boundary.failure_owner,0);
    f.blocks[0].alpha=f.blocks[1].alpha=1; f.v.setConstant(1e-30);
    EXPECT_EQ(ac::EvaluateComposite(x,f.y,f.beta,f.v,f.blocks).reason,"invalid-denominator");
    const auto fit=ac::FitComposite(x,f.y,f.beta,f.blocks,0,0);
    EXPECT_FALSE(fit.at("qualified").as_bool()); EXPECT_EQ(fit.at("reason"),"budget-exhausted");
    auto deficient=f.x; deficient.col(2)=deficient.col(0);
    EXPECT_EQ(ac::FitComposite(Sparse(deficient.sparseView(0,0)),f.y,f.beta,f.blocks).at("reason"),"rank-deficient");
    Eigen::VectorXd negative=f.beta; negative(0)=-1; negative(1)=0; negative(3)=.2;
    y=f.x*negative+(f.y-f.x*f.beta); const auto constrained=ac::FitComposite(x,y,f.beta,m::unique_grid::GlobalBlock(128,0));
    EXPECT_DOUBLE_EQ(V(constrained.at("beta"))(0),0);
}
TEST(AtomBlockGridCompositeTest, MembershipGeometryCoverageAndDiagnosticsAreIndependentOfBasis)
{
    rhbm_gem::MapObject map({8,9,10},{.5,.75,1.},{-1.5,-2.25,-4}); map.ClearMapValueArray();
    const std::vector<m::Atom> atoms{{{0,0,0},0,1e-10,0},{{.5,0,0},0,1e-10,0},{{-1.5,-2.25,-4},0,1e-10,0}};
    auto grid=m::atom_union::BuildGrid(atoms,{},map,map,1.5);
    auto blocks=m::atom_union::BuildBlocks(grid,atoms,map,{0,.1,1},1.5); std::vector<std::size_t> coverage(grid.voxels.size());
    for (const auto & b:blocks)
    {
        EXPECT_TRUE(std::is_sorted(b.rows.begin(),b.rows.end()));
        for (auto p:b.rows) {++coverage[static_cast<std::size_t>(p)]; EXPECT_LE(m::SquareDistance(grid.voxels[static_cast<std::size_t>(p)].position,atoms[b.owner].position),2.25);}
        for (std::size_t p=0;p<grid.voxels.size();++p)
            EXPECT_EQ(std::binary_search(b.rows.begin(),b.rows.end(),static_cast<Eigen::Index>(p)),m::SquareDistance(grid.voxels[p].position,atoms[b.owner].position)<=2.25);
    }
    for (std::size_t p=0;p<coverage.size();++p) EXPECT_EQ(coverage[p],grid.voxels[p].multiplicity);
}
TEST(AtomBlockGridCompositeTest, StagnationAndUnverifiedRefinementCannotQualify)
{
    Eigen::MatrixXd x(24,2); Eigen::VectorXd y(24);
    for (int i=0;i<24;++i) {x(i,0)=1; x(i,1)=std::sin(.71*i); y(i)=1e9+.1*std::sin(2.17*i);}
    auto blocks=m::unique_grid::GlobalBlock(24,0); blocks.push_back(blocks[0]); blocks[1].owner=1;
    auto fit=ac::FitComposite(Sparse(x.sparseView(0,0)),y,Eigen::Vector2d(1e9,0),blocks);
    EXPECT_FALSE(fit.at("qualified").as_bool()); EXPECT_EQ(fit.at("reason"),"stalled");
    Fixture f; const Sparse sparse=f.x.sparseView(0,0);
    fit=ac::FitComposite(sparse,f.y,f.beta,m::unique_grid::GlobalBlock(128,.1),100,0);
    EXPECT_FALSE(fit.at("qualified").as_bool()); EXPECT_EQ(fit.at("reason"),"reference-unverified");
}
