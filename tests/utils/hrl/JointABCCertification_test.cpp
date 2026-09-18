#include <gtest/gtest.h>
#include "support/JointABCCertification.hpp"
#include "support/JointABCPrecision.hpp"
#include "support/InstrumentedLM.hpp"
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
j::object Science(j::object value)
{
    for(const char * key:{"seconds","trust","lm"}) value.erase(key);
    return value;
}
struct RejectOnce
{
    int calls{},rejected{}; bool always{},reject_all{},tiny_step{}; std::string failure;
    bool guarded() const {return true;}
    bool retry() const {return calls<40;}
    int values() const {return 2;}
    int operator()(const Vector & x,Vector & r) {++calls; r=Vector::Constant(2,tiny_step ? -1 : x(0)-1); return 0;}
    int df(const Vector &,Eigen::MatrixXd & d) {d=Eigen::MatrixXd::Constant(2,1,tiny_step ? 1e100 : 1.); return 0;}
    bool Trial(const Vector &,const Vector &,const Vector &,double,double,double,double,double,bool proposed)
    {if(reject_all || (proposed && (always || rejected==0))) {++rejected; return false;} return true;}
};
}

TEST(JointABCCertificationTest, InstrumentedLegacyPreservesEveryTrial)
{
    Sample s; const p::Domain domain(s.grid,s.atoms); const Vector b=Eigen::Vector2d(.48,.59);
    const auto original=p::Fit(domain,s.y,b),legacy=p::Fit(domain,s.y,b,nullptr,"legacy");
    ASSERT_EQ(original.at("primary"),legacy.at("primary"));
    ASSERT_EQ(original.at("trials").as_array().size(),legacy.at("trials").as_array().size());
    for(std::size_t i=0;i<original.at("trials").as_array().size();++i)
        EXPECT_EQ(Science(original.at("trials").at(i).as_object()),Science(legacy.at("trials").at(i).as_object()));
}

TEST(JointABCCertificationTest, GuardedSearchOnlyAcceptsReplayedStates)
{
    Sample s; const p::Domain domain(s.grid,s.atoms); const Vector b=Eigen::Vector2d(.48,.59);
    for(const std::string variant:{"guarded","guarded-log"})
    {
        const auto fit=p::Fit(domain,s.y,b,nullptr,variant); ASSERT_TRUE(fit.at("joint_qualified").as_bool());
        for(const auto & trial:fit.at("trials").as_array()) if(trial.at("accepted").as_bool())
            EXPECT_TRUE(trial.at("trust").at("passed").as_bool());
        EXPECT_LE(j::value_to<int>(fit.at("profile_evaluations")),200);
    }
}

TEST(JointABCCertificationTest, InvalidTrialShrinksWithoutChangingAcceptedState)
{
    RejectOnce functor; p::InstrumentedLM<RejectOnce> lm(functor); Vector x=Vector::Zero(1);
    lm.parameters.factor=.1; auto status=lm.minimizeInit(x);
    while(status==Eigen::LevenbergMarquardtSpace::NotStarted || status==Eigen::LevenbergMarquardtSpace::Running) status=lm.minimizeOneStep(x);
    EXPECT_EQ(functor.rejected,1); EXPECT_NEAR(x(0),1,1e-10);
    RejectOnce fail; fail.always=true; p::InstrumentedLM<RejectOnce> rejecting(fail); x.setZero();
    rejecting.minimizeInit(x); status=rejecting.minimizeOneStep(x);
    EXPECT_EQ(status,Eigen::LevenbergMarquardtSpace::UserAsked); EXPECT_EQ(x(0),0); EXPECT_LE(fail.calls,40);
}

TEST(JointABCCertificationTest, HighPrecisionReferenceIsIndependentAndStable)
{
    Sample s; const p::Domain domain(s.grid,s.atoms); const Vector eta=Eigen::Vector2d(.48,.59).array().log();
    const auto e=p::Evaluate(domain,s.y,eta); ASSERT_TRUE(e.valid);
    const auto result=m::certification::PrecisionAudit(domain,s.y,e,Eigen::MatrixXd::Identity(2,2));
    EXPECT_TRUE(result.at("agreement_passed").as_bool())<<j::serialize(result);
    EXPECT_TRUE(result.at("derivative_passed").as_bool())<<j::serialize(result);
}

TEST(JointABCCertificationTest, ExplicitSupportArithmeticKeepsStrictBoundary)
{
    namespace sim=rhbm_gem::core::simulation;
    const std::array<double,3> origin{-3.6,-3.6,-3.6},spacing{.3,.3,.3},center{.1,0,0};
    EXPECT_GT(sim::SupportSquare(sim::GridPosition({7,16,7},spacing,origin),center),6.25);
    EXPECT_LE(sim::SupportSquare(sim::GridPosition({19,9,8},spacing,origin),center),6.25);
    EXPECT_EQ(sim::SupportSquare({2.5,0,0},{0,0,0}),6.25);
}

TEST(JointABCCertificationTest, EffectiveWidthsFollowAllModelSemantics)
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

TEST(JointABCCertificationTest, BoundaryReductionPreservesFeasibleCompensationScans)
{
    m::unique_grid::Grid grid; std::vector<m::Atom> atoms;
    Vector eta=Vector::Constant(12,std::log(.5)),beta=Vector::Zero(24);
    for(int a=0;a<12;++a)
    {
        atoms.push_back({{5.0*a,0,0},0,.5,.2}); beta(2*a+1)=.2;
        for(int k=0;k<8;++k)
        {m::unique_grid::Voxel v; v.position={5.0*a,.2*k,0}; v.index=grid.voxels.size(); grid.voxels.push_back(v);}
    }
    const p::Domain domain(grid,atoms);
    Vector y(static_cast<Eigen::Index>(grid.voxels.size()));
    for(Eigen::Index k=0;k<y.size();++k) y(k)=m::unique_grid::Direct(grid.voxels[static_cast<std::size_t>(k)].position,atoms,2.5);
    const auto audit=m::certification::BoundaryAudit(domain,y,eta,beta);
    ASSERT_TRUE(audit.at("agreement_passed").as_bool());
    for(const auto & row:audit.at("precision100").at("rows").as_array())
    {
        EXPECT_TRUE(row.at("path_feasible").as_bool());
        EXPECT_GT(j::value_to<int>(row.at("near_zero_rows")),0);
        EXPECT_LT(std::stod(j::value_to<std::string>(row.at("first_order_prediction_norm"))),1e-20);
        ASSERT_TRUE(row.at("profile_valid").as_bool());
        EXPECT_LE(std::stod(j::value_to<std::string>(row.at("profile_objective"))),
            std::stod(j::value_to<std::string>(row.at("path_objective")))+1e-40);
    }
}

TEST(JointABCCertificationTest, BudgetAndUnrepresentableStepsReturnWithoutAnUpdate)
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
