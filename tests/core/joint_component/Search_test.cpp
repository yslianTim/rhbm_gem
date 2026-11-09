#include <gtest/gtest.h>
#include "support/JointTestNumerics.hpp"
#include "core/detail/joint_component/BoundedLmDriver.hpp"
#include "support/JointSyntheticWorkload.hpp"
#include "core/detail/joint_component/ActiveSetBoundedLM.hpp"
#include "core/detail/joint_component/Problem.hpp"
#include "core/detail/joint_component/TiledQR.hpp"
#include "support/CommandTestHelpers.hpp"
#include "core/command/detail/SimulationGeometry.hpp"
#include <rhbm_gem/utils/math/ElectricPotential.hpp>

namespace {
namespace m=second_stage_test::matched;
namespace p=m::joint_component;
namespace joint=rhbm_gem::core::joint_component;
using Vector=Eigen::VectorXd;
using joint::BoundedLmDriver;
using joint::BoundedLmStatus;
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
    bool worsen_first_trial{};
    std::vector<double> trial_ratios, solve_radii, solve_eta, solve_response;
    bool CanRetry() const {return calls<retry_limit;}
    int ResidualCount() const {return 2;}
    bool Evaluate(const Vector & x,Vector & r)
    {
        ++calls; r=Vector::Constant(2,tiny_step ? -1e-100 : x(0)-1);
        if(worsen_first_trial && calls==2) r.setConstant(-2);
        return true;
    }
    bool Linearize(const Vector &,const Vector & residual,Eigen::MatrixXd & d,Vector & response,Vector & norms)
    {
        d=Eigen::MatrixXd::Constant(2,1,tiny_step ? 1e100 : 1.);
        response=residual; norms=d.colwise().blueNorm(); return true;
    }
    bool IsStationary(const Vector & eta) const {return !tiny_step && eta(0)==1.;}
    joint::ActiveSetStep SolveBoundedStep(const Eigen::MatrixXd & j,const Vector & response,
        const Vector & x,const Vector & diagonal,double radius,double damping)
    {
        solve_radii.push_back(radius); solve_eta.push_back(x(0)); solve_response.push_back(response(0));
        return joint::SolveBoundedTrustRegionStep(j,response,diagonal,Vector::Constant(x.size(),2.0)-x,radius,damping);
    }
    bool AssessTrial(const Vector &,const Vector &,const Vector &,double,double,double,double,double ratio,bool proposed)
    {
        trial_ratios.push_back(ratio);
        if(reject_all || (proposed && (always || (!worsen_first_trial && rejected==0))))
        {++rejected; return false;}
        return true;
    }
};
}

TEST(JointComponentSearchTest, GuardedSearchOnlyAcceptsReplayedStates)
{
    Sample s; const p::Domain domain(s.grid,s.atoms); const Vector tau=Eigen::Vector2d(.48,.59);
    const auto fit=p::Fit(domain,s.y,tau); ASSERT_EQ(fit.at("runtime_convergence"),"passed");
    auto context=p::MakeContext(s.y,static_cast<Eigen::Index>(s.atoms.size()));
    context.eta_upper_bound=std::log(4.0*tau.maxCoeff());
    const auto search=joint::SearchProfile(domain,s.y,tau,context);
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
    RejectOnce functor; BoundedLmDriver lm(functor,{.initial_radius_factor=.1}); Vector x=Vector::Zero(1);
    auto status=lm.Initialize(x);
    while(status==BoundedLmStatus::Ready || status==BoundedLmStatus::Running) status=lm.Advance(x);
    EXPECT_EQ(functor.rejected,1); EXPECT_NEAR(x(0),1,1e-10);
    RejectOnce fail; fail.always=true; fail.retry_limit=2;
    BoundedLmDriver rejecting(fail); x.setZero();
    rejecting.Initialize(x); status=rejecting.Advance(x);
    EXPECT_EQ(status,BoundedLmStatus::ProblemStopped); EXPECT_EQ(x(0),0); EXPECT_EQ(fail.calls,2);
}

TEST(JointComponentSearchTest, TrustedLowRatioAndUntrustedTrialUseDifferentRetries)
{
    for(const bool trusted_low_ratio:{false,true})
    {
        SCOPED_TRACE(trusted_low_ratio);
        RejectOnce problem; problem.worsen_first_trial=trusted_low_ratio;
        BoundedLmDriver lm(problem,{.initial_radius_factor=.1});
        Vector eta=Vector::Zero(1);
        ASSERT_EQ(lm.Initialize(eta),BoundedLmStatus::Ready);
        EXPECT_EQ(lm.Advance(eta),BoundedLmStatus::Running);
        ASSERT_EQ(problem.solve_radii.size(),2u);
        ASSERT_EQ(problem.trial_ratios.size(),2u);
        EXPECT_EQ(problem.solve_eta,(std::vector<double>{0,0}));
        EXPECT_EQ(problem.solve_response,(std::vector<double>{-1,-1}));
        EXPECT_GT(eta(0),0);
        EXPECT_LT(eta(0),.1);
        EXPECT_EQ(lm.AcceptedIterations(),1);
        EXPECT_EQ(lm.FunctionEvaluations(),3);
        EXPECT_EQ(problem.rejected,trusted_low_ratio ? 0 : 1);
        if(trusted_low_ratio) EXPECT_LT(problem.trial_ratios.front(),0);
        else EXPECT_GE(problem.trial_ratios.front(),1e-4);
        EXPECT_GE(problem.trial_ratios.back(),1e-4);
        EXPECT_NEAR(problem.solve_radii.back()/problem.solve_radii.front(),
            trusted_low_ratio ? .1 : .25,1e-12);
    }
}

TEST(JointComponentSearchTest, ProfileSearchPreservesAcceptedTrajectory)
{
    const rhbm_gem::core::JointProblem problem(second_stage_test::SyntheticJointWorkload("chain",3));
    const auto & data=rhbm_gem::core::JointProblemAccess::Get(problem);
    auto context=data.context;
    context.eta_upper_bound=std::log(.53);
    context.update_budget=3;
    const Vector tau=Eigen::Vector3d(.46,.49,.52);
    const auto search=joint::SearchProfile(data.domain,data.y,tau,context);
    const std::array<Eigen::Vector3d,4> expected_eta{{
        {-.776528789498996,-.713349887877465,-.653926467406664},
        {-.693859784826252,-.692966728106180,-.692837467146539},
        {-.693143774255377,-.693153781751112,-.693147361351597},
        {-.693143834605029,-.693153796340366,-.693147370056698}}};
    const std::array expected_radius{0.,.0290531450193031,.0581062900386062,.000503430554232906};
    const std::array expected_ratio{0.,1.01772153847589,1.00000563441645,1.00000488429624};
    ASSERT_EQ(search.trials.size(),expected_eta.size());
    for(std::size_t i=0;i<search.trials.size();++i)
    {
        SCOPED_TRACE(i);
        const auto & trial=search.trials[i];
        EXPECT_EQ(trial.evaluation,static_cast<int>(i+1));
        EXPECT_TRUE(trial.accepted);
        EXPECT_EQ(trial.accepted_update,static_cast<int>(i));
        ASSERT_TRUE(trial.trust && trial.trust->passed);
        ASSERT_TRUE(trial.lm && trial.lm->proposed_acceptance);
        EXPECT_TRUE(trial.endpoint.eta.isApprox(expected_eta[i],1e-11));
        EXPECT_DOUBLE_EQ(trial.lm->damping,0.);
        EXPECT_NEAR(trial.lm->radius,expected_radius[i],1e-12);
        EXPECT_NEAR(trial.lm->ratio,expected_ratio[i],1e-7);
        if(i>0)
        {
            EXPECT_TRUE(trial.lm->accepted_eta.isApprox(expected_eta[i-1],1e-11));
            EXPECT_TRUE((trial.lm->accepted_eta+trial.lm->step).isApprox(trial.endpoint.eta,1e-14));
        }
    }
    EXPECT_TRUE(search.eta.isApprox(expected_eta.back(),1e-11));
    EXPECT_EQ(search.lm_status,-1);
    EXPECT_EQ(search.stop_reason,"accepted-update-budget");
    EXPECT_TRUE(search.stopped);
    EXPECT_FALSE(search.width_failed);
    EXPECT_EQ(search.evaluations,4);
    EXPECT_EQ(search.derivatives,3);
    EXPECT_EQ(search.accepted,3);
    EXPECT_EQ(search.bounded_lm_work.subproblem_solves,3u);
    EXPECT_EQ(search.bounded_lm_work.damping_qr_factorizations,3u);
    const auto endpoint=joint::EvaluateProfile(data.domain,data.y,search.eta,&context);
    EXPECT_NEAR(endpoint.certificate.objective/(context.scale*context.scale),2.01990781328648e-10,1e-18);
    EXPECT_LE(endpoint.certificate.width_bound_projected_kkt,1e-12);
    const auto replay=joint::CheckIndependentReplay(data.domain,data.y,endpoint,context);
    EXPECT_TRUE(replay.passed);
    EXPECT_TRUE(joint::VerifyJointEndpoint(data.domain,data.y,endpoint,replay,context).numerical_qualified);
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
    RejectOnce budget; budget.reject_all=true; budget.retry_limit=2; BoundedLmDriver lm(budget);
    Vector x=Vector::Zero(1); lm.Initialize(x);
    EXPECT_EQ(lm.Advance(x),BoundedLmStatus::ProblemStopped);
    EXPECT_EQ(budget.calls,2); EXPECT_EQ(x(0),0); EXPECT_EQ(lm.AcceptedIterations(),0);
    RejectOnce tiny; tiny.tiny_step=true; BoundedLmDriver small(tiny);
    x.setOnes(); small.Initialize(x);
    EXPECT_EQ(small.Advance(x),BoundedLmStatus::ProblemStopped);
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
    int ResidualCount() const {return 7;}
    bool CanRetry() const {return false;}
    bool Evaluate(const Vector & x,Vector & r) {r=design*x-target; return true;}
    bool Linearize(const Vector &,const Vector & residual,Eigen::MatrixXd & factor,Vector & response,Vector & norms)
    {
        norms=design.colwise().blueNorm();
        if(!tiled) {factor=design; response=residual; return true;}
        p::runtime::TiledQR reduced(2,1);
        for(Eigen::Index first=0;first<7;first+=2)
        {const auto count=std::min<Eigen::Index>(2,7-first); reduced.Append(design.middleRows(first,count),residual.segment(first,count));}
        factor=reduced.r; response=reduced.target.col(0); return true;
    }
    joint::ActiveSetStep SolveBoundedStep(const Eigen::MatrixXd & j,const Vector & response,
        const Vector & x,const Vector & diagonal,double radius,double damping)
    {return joint::SolveBoundedTrustRegionStep(j,response,diagonal,Vector::Constant(x.size(),2.)-x,radius,damping);}
    bool IsStationary(const Vector &) const {return false;}
    bool AssessTrial(const Vector &,const Vector & s,const Vector &,double,double,double a,double p0,double r,bool)
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
    int ResidualCount() const {return 2;}
    bool CanRetry() const {return false;}
    bool Evaluate(const Vector & x,Vector & residual) {residual=design*x-target; return true;}
    bool Linearize(const Vector &,const Vector & residual,Eigen::MatrixXd & factor,Vector & response,Vector & norms)
    {factor=design; response=residual; norms=design.colwise().blueNorm(); return true;}
    joint::ActiveSetStep SolveBoundedStep(const Eigen::MatrixXd & j,const Vector & response,
        const Vector & x,const Vector & diagonal,double trust_radius,double lambda)
    {return joint::SolveBoundedTrustRegionStep(j,response,diagonal,upper_eta-x,trust_radius,lambda);}
    bool IsStationary(const Vector & eta) const
    {
        Vector gradient=design.transpose()*(design*eta-target);
        for(Eigen::Index i=0;i<eta.size();++i)
            if(eta(i)>=upper_eta(i)) gradient(i)=std::max(0.,gradient(i));
        return gradient.isZero(0.);
    }
    bool AssessTrial(const Vector &,const Vector & proposed_step,const Vector &,double used_radius,double used_damping,
        double used_actual,double used_predicted,double used_ratio,bool)
    {step=proposed_step; radius=used_radius; damping=used_damping; actual=used_actual;
        predicted=used_predicted; ratio=used_ratio; return true;}
};

}

TEST(JointComponentSearchTest, CompactModelUsesFullResidualForActualReduction)
{
    LinearResidual dense(false),compact(true);
    for(auto * problem:{&dense,&compact})
    {
        BoundedLmDriver lm(*problem,{.initial_radius_factor=.1});
        Vector eta=Vector::Zero(2);
        ASSERT_EQ(lm.Initialize(eta),BoundedLmStatus::Ready);
        EXPECT_EQ(lm.Advance(eta),BoundedLmStatus::Running);
        EXPECT_EQ(lm.AcceptedIterations(),1);
        EXPECT_TRUE(eta.isApprox(problem->step,1e-14));
        EXPECT_NEAR(problem->actual,1.-(problem->design*eta-problem->target).squaredNorm()/
            problem->target.squaredNorm(),1e-14);
    }
    EXPECT_TRUE(compact.step.isApprox(dense.step,1e-12));
    EXPECT_NEAR(compact.predicted,dense.predicted,1e-12);
    EXPECT_NEAR(compact.ratio,dense.ratio,1e-12);
}

TEST(JointComponentSearchTest, BoundedLmUsesTheConstrainedStepForTrustRatio)
{
    BoundedLinearResidual functor;
    // Reach the active bound with an origin-independent radius.
    BoundedLmDriver lm(functor,{.initial_radius_factor=.5});
    Vector x(2); x<<10.0,10.0;
    ASSERT_EQ(lm.Initialize(x),BoundedLmStatus::Ready);
    const auto status=lm.Advance(x);
    EXPECT_TRUE(status==BoundedLmStatus::Running ||
        status==BoundedLmStatus::StepTolerance ||
        status==BoundedLmStatus::ReductionTolerance ||
        status==BoundedLmStatus::ReductionAndStepTolerance);
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
    BoundedLmDriver lm(functor,{.initial_radius_factor=.1,.evaluation_budget=200});
    Vector x(2); x<<10.0,10.0;
    ASSERT_EQ(lm.Initialize(x),BoundedLmStatus::Ready);
    auto status=BoundedLmStatus::Running;
    for(int iteration=0;iteration<100 && status==BoundedLmStatus::Running;++iteration)
        status=lm.Advance(x);
    EXPECT_NE(status,BoundedLmStatus::ProblemStopped)<<functor.failure;
    EXPECT_NEAR(x(0),10.2,1e-12);
    EXPECT_NEAR(x(1),0.0,1e-8);
    EXPECT_EQ(status,BoundedLmStatus::Stationary);
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

TEST(JointComponentSearchTest, InfeasibleRoundedTrialKeepsInitialAcceptanceEvidence)
{
    Sample sample;
    for(auto & atom:sample.atoms) atom.width=1.2;
    for(Eigen::Index row=0;row<sample.y.size();++row)
        sample.y(row)=m::unique_grid::Direct(sample.grid.voxels[static_cast<std::size_t>(row)].position,sample.atoms,2.5);
    const p::Domain domain(sample.grid,sample.atoms);
    auto context=p::MakeContext(sample.y,2);
    // Adding the full bound increment rounds to zero from log(.98), above the bound.
    context.eta_upper_bound=-1e-18;
    const auto search=joint::SearchProfile(domain,sample.y,Vector::Constant(2,.98),context);
    ASSERT_TRUE(search.initial_accepted);
    ASSERT_FALSE(search.trials.empty());
    const auto & initial=search.trials.front();
    ASSERT_TRUE(initial.lm && initial.trust && initial.trust->passed);
    EXPECT_TRUE(initial.lm->proposed_acceptance);
    EXPECT_TRUE(initial.lm->step.isZero(0.));
    EXPECT_DOUBLE_EQ(initial.lm->radius,0.);
}

namespace {
// The same physical log-width t=x+reference is expressed about different origins.
struct ShiftedWidthResidual
{
    double reference{},target{2.1},sensitivity{2.199712482e-4},background{1e-3};
    int evaluations{};
    double initial_radius{};
    std::string failure;
    int ResidualCount() const {return 2;}
    bool CanRetry() const {return false;}
    bool Evaluate(const Vector & x,Vector & residual)
    {++evaluations; residual.resize(2); residual<<sensitivity*(x(0)+reference-target),background; return true;}
    bool Linearize(const Vector &,const Vector & residual,Eigen::MatrixXd & j,Vector & response,Vector & norms)
    {j=Eigen::MatrixXd::Zero(2,1); j(0,0)=sensitivity; response=residual; norms=Vector::Constant(1,sensitivity); return true;}
    joint::ActiveSetStep SolveBoundedStep(const Eigen::MatrixXd & j,const Vector & r,
        const Vector & x,const Vector & diagonal,double radius,double damping)
    {initial_radius=radius; return joint::SolveBoundedTrustRegionStep(j,r,diagonal,
        Vector::Constant(1,std::log(3.0660309206)-reference)-x,radius,damping);}
    bool IsStationary(const Vector & eta) const
    {return eta(0)+reference>=std::log(3.0660309206) || eta(0)+reference==target;}
    bool AssessTrial(const Vector &,const Vector &,const Vector &,double,double,double,double,double,bool) {return true;}
};
}
TEST(JointComponentSearchTest, NearUnitWidthHasOriginIndependentRadiusAndEvaluatesCandidates)
{
    double first_radius=0.0,first_eta=0.0;
    for(const double reference:{0.0,1.0,-4.0})
        for(const double physical_eta:{0.0,std::log(1.000000000626)})
        {
            ShiftedWidthResidual f; f.reference=reference;
            BoundedLmDriver lm(f,{.initial_radius_factor=.1,.gtol=1e-12});
            Vector x=Vector::Constant(1,physical_eta-reference);
            lm.Initialize(x); auto status=lm.Advance(x);
            EXPECT_GT(f.evaluations,1)<<f.failure;
            EXPECT_TRUE(x.allFinite());
            EXPECT_GT(x(0)+reference,physical_eta);
            EXPECT_NEAR(f.initial_radius,.1*f.sensitivity,1e-18);
            if(first_radius==0.0) {first_radius=f.initial_radius; first_eta=x(0)+reference;}
            EXPECT_NEAR(f.initial_radius,first_radius,1e-18);
            EXPECT_NEAR(x(0)+reference,first_eta,1e-8);
            EXPECT_LT(std::abs(x(0)+reference-f.target),std::abs(physical_eta-f.target));
            for(int i=0;i<30 && status==BoundedLmStatus::Running;++i)
                status=lm.Advance(x);
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
    EXPECT_EQ(search.bounded_lm_work.failed_subproblems,1u);
    EXPECT_EQ(search.stop_reason,"bounded-subproblem:injected-width-subproblem-failure");
    EXPECT_TRUE(search.eta.allFinite());
}

TEST(JointComponentSearchTest, SmallWeightedNormStillEvaluatesAMeaningfulPhysicalStep)
{
    ShiftedWidthResidual f; f.sensitivity=1e-16; f.background=1e-20;
    BoundedLmDriver lm(f,{.initial_radius_factor=.1});
    Vector x=Vector::Zero(1); lm.Initialize(x); lm.Advance(x);
    EXPECT_GT(f.evaluations,1)<<f.failure;
    EXPECT_TRUE(x.allFinite());
    EXPECT_GT(x(0),.09);
    EXPECT_LE(x(0),.1+1e-14);
    EXPECT_LT(f.sensitivity*std::abs(x(0)),64*std::numeric_limits<double>::epsilon());
}
