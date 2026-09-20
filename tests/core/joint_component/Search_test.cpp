#include <gtest/gtest.h>
#include "support/JointTestNumerics.hpp"
#include "support/InstrumentedLM.hpp"
#include "core/detail/joint_component/TiledQR.hpp"
#include "support/CommandTestHelpers.hpp"
#include "core/command/detail/SimulationGeometry.hpp"
#include <rhbm_gem/utils/math/ElectricPotential.hpp>

namespace {
namespace m=second_stage_test::matched;
namespace p=m::joint_abc;
namespace j=boost::json;
using Vector=Eigen::VectorXd;
struct Sample
{
    m::unique_grid::Grid grid;
    std::vector<m::Atom> atoms{{{-.6,0,0},2,.42,-.3},{{.7,.2,0},1.5,.67,.4}};
    Vector y;
    Sample()
    {
        for(int x=-5;x<=5;++x) for(int z=-2;z<=2;++z) for(int y0=-3;y0<=3;++y0)
        {m::unique_grid::Voxel v; v.position={x*.4,y0*.4,z*.4}; v.index=grid.voxels.size(); grid.voxels.push_back(v);}
        y.resize(static_cast<Eigen::Index>(grid.voxels.size()));
        for(Eigen::Index k=0;k<y.size();++k) y(k)=m::unique_grid::Direct(grid.voxels[static_cast<std::size_t>(k)].position,atoms,2.5);
    }
};
struct RejectOnce
{
    int calls{},rejected{}; bool always{},reject_all{},tiny_step{}; std::string failure;
    bool retry() const {return calls<40;}
    int values() const {return 2;}
    int operator()(const Vector & x,Vector & r) {++calls; r=Vector::Constant(2,tiny_step ? -1 : x(0)-1); return 0;}
    int df(const Vector &,Eigen::MatrixXd & d) {d=Eigen::MatrixXd::Constant(2,1,tiny_step ? 1e100 : 1.); return 0;}
    int linearize(const Vector & x,const Vector & residual,Eigen::MatrixXd & d,Vector & response,Vector & norms)
    {df(x,d); response=residual; norms=d.colwise().blueNorm(); return 0;}
    bool Trial(const Vector &,const Vector &,const Vector &,double,double,double,double,double,bool proposed)
    {if(reject_all || (proposed && (always || rejected==0))) {++rejected; return false;} return true;}
};
}

TEST(JointComponentSearchTest, GuardedSearchOnlyAcceptsReplayedStates)
{
    Sample s; const p::Domain domain(s.grid,s.atoms); const Vector b=Eigen::Vector2d(.48,.59);
    const auto fit=p::Fit(domain,s.y,b); ASSERT_EQ(fit.at("runtime_convergence"),"passed");
    for(const auto & trial:fit.at("trials").as_array()) if(trial.at("accepted").as_bool())
        EXPECT_TRUE(trial.at("trust").at("passed").as_bool());
    EXPECT_LE(j::value_to<int>(fit.at("profile_evaluations")),200);
}

TEST(JointComponentSearchTest, InvalidTrialShrinksWithoutChangingAcceptedState)
{
    RejectOnce functor; p::InstrumentedLM<RejectOnce> lm(functor); Vector x=Vector::Zero(1);
    lm.parameters.factor=.1; auto status=lm.minimizeInit(x);
    while(status==Eigen::LevenbergMarquardtSpace::NotStarted || status==Eigen::LevenbergMarquardtSpace::Running) status=lm.minimizeOneStep(x);
    EXPECT_EQ(functor.rejected,1); EXPECT_NEAR(x(0),1,1e-10);
    RejectOnce fail; fail.always=true; p::InstrumentedLM<RejectOnce> rejecting(fail); x.setZero();
    rejecting.minimizeInit(x); status=rejecting.minimizeOneStep(x);
    EXPECT_EQ(status,Eigen::LevenbergMarquardtSpace::UserAsked); EXPECT_EQ(x(0),0); EXPECT_LE(fail.calls,40);
}

TEST(JointComponentSearchTest, ExplicitSupportArithmeticKeepsStrictBoundary)
{
    namespace sim=rhbm_gem::core::simulation;
    const std::array<double,3> origin{-3.6,-3.6,-3.6},spacing{.3,.3,.3},center{.1,0,0};
    EXPECT_GT(sim::SupportSquare(sim::GridPosition({7,16,7},spacing,origin),center),6.25);
    EXPECT_LE(sim::SupportSquare(sim::GridPosition({19,9,8},spacing,origin),center),6.25);
    EXPECT_EQ(sim::SupportSquare({2.5,0,0},{0,0,0}),6.25);
}

TEST(JointComponentSearchTest, EffectiveWidthsFollowAllModelSemantics)
{
    ElectricPotential potential; potential.SetBlurringWidth(.5); potential.SetModelChoice(0);
    EXPECT_EQ(potential.GetEffectiveWidths(Element::OXYGEN).gaussian,.4);
    EXPECT_EQ(potential.GetEffectiveWidths(Element::NITROGEN).charge,.45);
    EXPECT_EQ(potential.GetEffectiveWidths(Element::CARBON).gaussian,.5);
    potential.SetModelChoice(1); EXPECT_FALSE(potential.GetEffectiveWidths(Element::OXYGEN).gaussian);
    EXPECT_EQ(potential.GetEffectiveWidths(Element::OXYGEN).charge,.5);
    potential.SetBlurringWidth(0); EXPECT_FALSE(potential.GetEffectiveWidths(Element::OXYGEN).charge);
    potential.SetModelChoice(2); EXPECT_FALSE(potential.GetEffectiveWidths(Element::OXYGEN).gaussian);
}

TEST(JointComponentSearchTest, BudgetAndUnrepresentableStepsReturnWithoutAnUpdate)
{
    RejectOnce budget; budget.reject_all=true; p::InstrumentedLM<RejectOnce> lm(budget);
    Vector x=Vector::Zero(1); lm.minimizeInit(x);
    EXPECT_EQ(lm.minimizeOneStep(x),Eigen::LevenbergMarquardtSpace::UserAsked);
    EXPECT_EQ(budget.calls,40); EXPECT_EQ(x(0),0); EXPECT_EQ(lm.iter,1);
    RejectOnce tiny; tiny.tiny_step=true; p::InstrumentedLM<RejectOnce> small(tiny);
    small.useExternalScaling=true; small.diag=Vector::Ones(1); x.setOnes(); small.minimizeInit(x);
    EXPECT_EQ(small.minimizeOneStep(x),Eigen::LevenbergMarquardtSpace::UserAsked);
    EXPECT_EQ(tiny.failure,"unrepresentable-step"); EXPECT_EQ(x(0),1); EXPECT_EQ(tiny.calls,1);
}

namespace {
struct LinearResidual
{
    Eigen::MatrixXd design;
    Vector target,step;
    bool tiled;
    double predicted{},actual{},ratio{};
    std::string failure;
    explicit LinearResidual(bool compact):design(7,2),target(7),tiled(compact)
    {design<<1,2, 3,-1, 2,4, -1,3, 4,2, 2,-2, 0,1; target<<1,3,-2,5,2,-3,4;}
    int values() const {return 7;}
    bool retry() const {return false;}
    int operator()(const Vector & x,Vector & r) {r=design*x-target; return 0;}
    int linearize(const Vector &,const Vector & residual,Eigen::MatrixXd & factor,Vector & response,Vector & norms)
    {
        norms=design.colwise().blueNorm();
        if(!tiled) {factor=design; response=residual; return 0;}
        p::runtime::TiledQR reduced(2,1);
        for(Eigen::Index first=0;first<7;first+=2)
        {const auto count=std::min<Eigen::Index>(2,7-first); reduced.Append(design.middleRows(first,count),residual.segment(first,count));}
        factor=reduced.r; response=reduced.target.col(0); return 0;
    }
    bool Trial(const Vector &,const Vector & s,const Vector &,double,double,double a,double p0,double r,bool)
    {step=s; actual=a; predicted=p0; ratio=r; return true;}
};
}
TEST(JointComponentSearchTest, CompactLinearizationPreservesStepAndFullResidualReduction)
{
    LinearResidual dense(false),tiled(true); p::InstrumentedLM<LinearResidual> a(dense),b(tiled);
    Vector x=Vector::Zero(2),y=x;
    a.parameters.factor=b.parameters.factor=.1;
    a.minimizeInit(x); b.minimizeInit(y);
    EXPECT_EQ(a.fjac.rows(),2); EXPECT_EQ(b.fjac.rows(),2);
    EXPECT_EQ(a.minimizeOneStep(x),b.minimizeOneStep(y));
    EXPECT_LT((x-y).norm(),1e-12); EXPECT_LT((dense.step-tiled.step).norm(),1e-12);
    EXPECT_NEAR(dense.predicted,tiled.predicted,1e-12);
    EXPECT_NEAR(dense.actual,tiled.actual,1e-12); EXPECT_NEAR(dense.ratio,tiled.ratio,1e-12);
    EXPECT_NEAR(b.fnorm,(tiled.design*y-tiled.target).stableNorm(),1e-12);
    EXPECT_GT(b.fnorm,1); // Nonzero orthogonal residual must survive row reduction.
    EXPECT_EQ(b.fjac.rows(),2);
}
