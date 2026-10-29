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
namespace joint=rhbm_gem::core::joint_component;
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
    int calls{},rejected{},retry_limit{40}; bool always{},reject_all{},tiny_step{}; std::string failure;
    bool retry() const {return calls<retry_limit;}
    int values() const {return 2;}
    int operator()(const Vector & x,Vector & r) {++calls; r=Vector::Constant(2,tiny_step ? -1e-100 : x(0)-1); return 0;}
    int df(const Vector &,Eigen::MatrixXd & d) {d=Eigen::MatrixXd::Constant(2,1,tiny_step ? 1e100 : 1.); return 0;}
    int linearize(const Vector & x,const Vector & residual,Eigen::MatrixXd & d,Vector & response,Vector & norms)
    {df(x,d); response=residual; norms=d.colwise().blueNorm(); return 0;}
    joint::ActiveSetStep solve_bounded_step(const Eigen::MatrixXd & j,const Vector & response,
        const Vector & x,const Vector & diagonal,double radius,double damping)
    {return joint::SolveBoundedTrustRegionStep(j,response,diagonal,Vector::Constant(x.size(),2.0)-x,radius,damping);}
    bool Trial(const Vector &,const Vector &,const Vector &,double,double,double,double,double,bool proposed)
    {if(reject_all || (proposed && (always || rejected==0))) {++rejected; return false;} return true;}
};
}

TEST(JointComponentSearchTest, GuardedSearchOnlyAcceptsReplayedStates)
{
    Sample s; const p::Domain domain(s.grid,s.atoms); const Vector b=Eigen::Vector2d(.48,.59);
    const auto fit=p::Fit(domain,s.y,b); ASSERT_EQ(fit.at("runtime_convergence"),"passed");
    auto context=p::MakeContext(s.y,static_cast<Eigen::Index>(s.atoms.size()));
    context.eta_upper_bound=std::log(4.0*b.maxCoeff());
    const auto search=joint::SearchProfile(domain,s.y,b,context);
    ASSERT_FALSE(search.trials.empty());
    for(const auto & trial:search.trials) if(trial.accepted)
    {
        ASSERT_TRUE(trial.trust.has_value());
        EXPECT_TRUE(trial.trust->passed);
    }
    EXPECT_LE(search.evaluations,200);
}

TEST(JointComponentSearchTest, InvalidTrialShrinksWithoutChangingAcceptedState)
{
    RejectOnce functor; p::InstrumentedLM<RejectOnce> lm(functor); Vector x=Vector::Zero(1);
    lm.parameters.factor=.1; auto status=lm.minimizeInit(x);
    while(status==Eigen::LevenbergMarquardtSpace::NotStarted || status==Eigen::LevenbergMarquardtSpace::Running) status=lm.minimizeOneStep(x);
    EXPECT_EQ(functor.rejected,1); EXPECT_NEAR(x(0),1,1e-10);
    RejectOnce fail; fail.always=true; fail.retry_limit=2;
    p::InstrumentedLM<RejectOnce> rejecting(fail); x.setZero();
    rejecting.minimizeInit(x); status=rejecting.minimizeOneStep(x);
    EXPECT_EQ(status,Eigen::LevenbergMarquardtSpace::UserAsked); EXPECT_EQ(x(0),0); EXPECT_EQ(fail.calls,2);
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

TEST(JointComponentSearchTest, BudgetAndTinyBoundedStepsReturnWithoutAnUpdate)
{
    RejectOnce budget; budget.reject_all=true; budget.retry_limit=2; p::InstrumentedLM<RejectOnce> lm(budget);
    Vector x=Vector::Zero(1); lm.minimizeInit(x);
    EXPECT_EQ(lm.minimizeOneStep(x),Eigen::LevenbergMarquardtSpace::UserAsked);
    EXPECT_EQ(budget.calls,2); EXPECT_EQ(x(0),0); EXPECT_EQ(lm.iterations(),1);
    RejectOnce tiny; tiny.tiny_step=true; p::InstrumentedLM<RejectOnce> small(tiny);
    x.setOnes(); small.minimizeInit(x);
    EXPECT_EQ(small.minimizeOneStep(x),Eigen::LevenbergMarquardtSpace::UserAsked);
    EXPECT_EQ(tiny.failure,"bounded-tiny-step-nonstationary"); EXPECT_EQ(x(0),1); EXPECT_EQ(tiny.calls,1);
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

struct BoundedLinearResidual
{
    Eigen::MatrixXd design;
    Vector target,upper_eta;
    Vector step;
    double radius{},damping{},actual{},predicted{},ratio{};
    std::string failure;
    BoundedLinearResidual():design(2,2),target(2),upper_eta(2)
    {design.setIdentity(); target<<20.0,0.0; upper_eta<<10.2,12.0;}
    int values() const {return 2;}
    bool retry() const {return false;}
    int operator()(const Vector & x,Vector & residual) {residual=design*x-target; return 0;}
    int linearize(const Vector &,const Vector & residual,Eigen::MatrixXd & factor,Vector & response,Vector & norms)
    {factor=design; response=residual; norms=design.colwise().blueNorm(); return 0;}
    joint::ActiveSetStep solve_bounded_step(const Eigen::MatrixXd & j,const Vector & response,
        const Vector & x,const Vector & diagonal,double trust_radius,double lambda)
    {return joint::SolveBoundedTrustRegionStep(j,response,diagonal,upper_eta-x,trust_radius,lambda);}
    bool Trial(const Vector &,const Vector & proposed_step,const Vector &,double used_radius,double used_damping,
        double used_actual,double used_predicted,double used_ratio,bool)
    {step=proposed_step; radius=used_radius; damping=used_damping; actual=used_actual;
        predicted=used_predicted; ratio=used_ratio; return true;}
};

}

TEST(JointComponentSearchTest, InstrumentedLmUsesTheConstrainedStepForTrustRatio)
{
    BoundedLinearResidual functor;
    p::InstrumentedLM<BoundedLinearResidual> lm(functor);
    Vector x(2); x<<10.0,10.0;
    lm.parameters.factor=.5; // Reach the active bound with an origin-independent radius.
    ASSERT_EQ(lm.minimizeInit(x),Eigen::LevenbergMarquardtSpace::NotStarted);
    const auto status=lm.minimizeOneStep(x);
    EXPECT_TRUE(status==Eigen::LevenbergMarquardtSpace::Running ||
        status==Eigen::LevenbergMarquardtSpace::RelativeErrorTooSmall ||
        status==Eigen::LevenbergMarquardtSpace::RelativeReductionTooSmall ||
        status==Eigen::LevenbergMarquardtSpace::RelativeErrorAndReductionTooSmall);
    ASSERT_EQ(functor.step.size(),2);
    EXPECT_NEAR(functor.step(0),.2,1e-13);
    EXPECT_LE(functor.step.norm(),functor.radius+1e-13);
    const Vector initial_residual=Eigen::Vector2d(-10.0,10.0);
    const double expected_actual=1.0-(initial_residual+functor.step).squaredNorm()/initial_residual.squaredNorm();
    const double expected_predicted=(functor.step.squaredNorm()+
        2.0*functor.damping*functor.step.squaredNorm())/initial_residual.squaredNorm();
    EXPECT_NEAR(functor.actual,expected_actual,1e-13);
    EXPECT_NEAR(functor.predicted,expected_predicted,1e-13)<<"step="<<functor.step.transpose()
        <<" damping="<<functor.damping<<" radius="<<functor.radius;
    EXPECT_NEAR(functor.ratio,expected_actual/expected_predicted,1e-12)<<"step="<<functor.step.transpose()
        <<" damping="<<functor.damping<<" radius="<<functor.radius;
}

TEST(JointComponentSearchTest, BoundedLmStopsAtAValidBoundaryKktPoint)
{
    BoundedLinearResidual functor;
    p::InstrumentedLM<BoundedLinearResidual> lm(functor);
    Vector x(2); x<<10.0,10.0;
    lm.parameters.factor=.1; lm.parameters.maxfev=200;
    ASSERT_EQ(lm.minimizeInit(x),Eigen::LevenbergMarquardtSpace::NotStarted);
    auto status=Eigen::LevenbergMarquardtSpace::Running;
    for(int iteration=0;iteration<100 && status==Eigen::LevenbergMarquardtSpace::Running;++iteration)
        status=lm.minimizeOneStep(x);
    EXPECT_NE(status,Eigen::LevenbergMarquardtSpace::UserAsked)<<functor.failure;
    EXPECT_NEAR(x(0),10.2,1e-12);
    EXPECT_NEAR(x(1),0.0,1e-8);
    EXPECT_EQ(status,Eigen::LevenbergMarquardtSpace::CosinusTooSmall);
    EXPECT_EQ(functor.failure,"bounded-stationary");
}

TEST(JointComponentSearchTest, BoundedProfileSearchKeepsEveryEvaluatedWidthFeasible)
{
    Sample sample;
    const p::Domain domain(sample.grid,sample.atoms);
    auto context=p::MakeContext(sample.y,static_cast<Eigen::Index>(sample.atoms.size()));
    const Vector initial=Eigen::Vector2d(.48,.59);
    const double eta_max=std::log(.6);
    context.eta_upper_bound=eta_max;
    const auto search=joint::SearchProfile(domain,sample.y,initial,context);
    EXPECT_FALSE(search.trials.empty())<<search.stop_reason;
    bool reached_bound=false;
    for(const auto & trial:search.trials)
    {
        if(trial.endpoint.eta.size()) EXPECT_LE(trial.endpoint.eta.maxCoeff(),eta_max+1e-14);
        if(trial.lm && trial.lm->bounded_subproblem)
        {
            EXPECT_EQ(trial.lm->subproblem_status,"solved");
            EXPECT_LE(trial.lm->subproblem_kkt_inf_norm,1e-9);
            EXPECT_LE(trial.lm->diagonal.cwiseProduct(trial.lm->step).norm(),trial.lm->radius+1e-12);
            reached_bound|=!trial.lm->active_bound_indices.empty();
        }
    }
    EXPECT_TRUE(reached_bound)<<search.stop_reason;
}

TEST(JointComponentSearchTest, ProfileSearchRequiresFiniteWidthUpperBound)
{
    Sample sample;
    const p::Domain domain(sample.grid,sample.atoms);
    const auto context=p::MakeContext(sample.y,static_cast<Eigen::Index>(sample.atoms.size()));
    const auto search=joint::SearchProfile(domain,sample.y,Vector::Constant(2,.5),context);
    EXPECT_TRUE(search.stopped);
    EXPECT_TRUE(search.trials.empty());
    EXPECT_EQ(search.stop_reason,"invalid-width-upper-bound");
}

namespace {
// The same physical log-width t=x+reference is expressed about different origins.
struct ShiftedWidthResidual
{
    double reference{},target{2.1},sensitivity{2.199712482e-4},background{1e-3};
    int evaluations{};
    double initial_radius{};
    std::string failure;
    int values() const {return 2;}
    bool retry() const {return false;}
    int operator()(const Vector & x,Vector & residual)
    {++evaluations; residual.resize(2); residual<<sensitivity*(x(0)+reference-target),background; return 0;}
    int linearize(const Vector &,const Vector & residual,Eigen::MatrixXd & j,Vector & response,Vector & norms)
    {j=Eigen::MatrixXd::Zero(2,1); j(0,0)=sensitivity; response=residual; norms=Vector::Constant(1,sensitivity); return 0;}
    joint::ActiveSetStep solve_bounded_step(const Eigen::MatrixXd & j,const Vector & r,
        const Vector & x,const Vector & diagonal,double radius,double damping)
    {initial_radius=radius; return joint::SolveBoundedTrustRegionStep(j,r,diagonal,
        Vector::Constant(1,std::log(3.0660309206)-reference)-x,radius,damping);}
    bool Trial(const Vector &,const Vector &,const Vector &,double,double,double,double,double,bool) {return true;}
};
}
TEST(JointComponentSearchTest, NearUnitWidthHasOriginIndependentRadiusAndEvaluatesCandidates)
{
    double first_radius=0.0,first_eta=0.0;
    for(const double reference:{0.0,1.0,-4.0})
        for(const double physical_eta:{0.0,std::log(1.000000000626)})
        {
            ShiftedWidthResidual f; f.reference=reference;
            p::InstrumentedLM<ShiftedWidthResidual> lm(f);
            Vector x=Vector::Constant(1,physical_eta-reference); lm.parameters.factor=.1;
            lm.parameters.gtol=1e-12;
            lm.minimizeInit(x); auto status=lm.minimizeOneStep(x);
            EXPECT_GT(f.evaluations,1)<<f.failure;
            EXPECT_TRUE(x.allFinite());
            EXPECT_GT(x(0)+reference,physical_eta);
            EXPECT_NEAR(f.initial_radius,.1*f.sensitivity,1e-18);
            if(first_radius==0.0) {first_radius=f.initial_radius; first_eta=x(0)+reference;}
            EXPECT_NEAR(f.initial_radius,first_radius,1e-18);
            EXPECT_NEAR(x(0)+reference,first_eta,1e-8);
            EXPECT_LT(std::abs(x(0)+reference-f.target),std::abs(physical_eta-f.target));
            for(int i=0;i<30 && status==Eigen::LevenbergMarquardtSpace::Running;++i)
                status=lm.minimizeOneStep(x);
            const double eta=x(0)+reference;
            const double upper=std::log(3.0660309206);
            const double gradient=f.sensitivity*f.sensitivity*(eta-f.target);
            EXPECT_LE(eta,upper+1e-14);
            EXPECT_LE(eta>=upper-1e-14 ? std::max(0.0,gradient) : std::abs(gradient),1e-9);
        }
}

TEST(JointComponentSearchTest, InitialAcAcceptanceDoesNotHideAWidthSubproblemFailure)
{
    Sample sample; const p::Domain domain(sample.grid,sample.atoms);
    auto context=p::MakeContext(sample.y,2); context.eta_upper_bound=std::log(2.2);
    context.force_width_subproblem_failure_for_testing=true;
    const auto search=joint::SearchProfile(domain,sample.y,Eigen::Vector2d(.48,.59),context);
    EXPECT_TRUE(search.initial_accepted);
    EXPECT_TRUE(search.width_failed);
    EXPECT_EQ(search.accepted,0);
    EXPECT_EQ(search.profile_work.bounded_lm.failed_subproblems,1u);
    EXPECT_EQ(search.stop_reason,"bounded-subproblem:injected-width-subproblem-failure");
    EXPECT_TRUE(search.eta.allFinite());
}

TEST(JointComponentSearchTest, SmallWeightedNormStillEvaluatesAMeaningfulPhysicalStep)
{
    ShiftedWidthResidual f; f.sensitivity=1e-16; f.background=1e-20;
    p::InstrumentedLM<ShiftedWidthResidual> lm(f); lm.parameters.factor=.1;
    Vector x=Vector::Zero(1); lm.minimizeInit(x); lm.minimizeOneStep(x);
    EXPECT_GT(f.evaluations,1)<<f.failure;
    EXPECT_TRUE(x.allFinite());
    EXPECT_GT(x(0),.09);
    EXPECT_LE(x(0),.1+1e-14);
    EXPECT_LT(f.sensitivity*std::abs(x(0)),64*std::numeric_limits<double>::epsilon());
}
