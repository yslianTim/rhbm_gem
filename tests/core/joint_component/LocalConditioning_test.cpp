#include <gtest/gtest.h>
#include "core/estimator/joint_component/GaussSeidelBcd.hpp"
#include "core/estimator/joint_component/LocalConditioning.hpp"
#include "core/estimator/joint_component/Problem.hpp"
#include "core/estimator/joint_component/TargetEvidence.hpp"
#include "core/estimator/joint_component/TiledDerivative.hpp"
#include "support/JointSyntheticWorkload.hpp"
#include "support/JointTestContext.hpp"
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>

namespace {
namespace c=rhbm_gem::core;
namespace n=c::joint_component;
namespace tc=joint_component_test::matched::joint_component;

n::LocalConditioningAssessment Assess(const n::Domain & domain,n::VectorRef y,const n::Vector & eta,
    const n::Vector & beta,const n::EvaluationContext & context)
{
    return n::AssessBcdLocalConditioning(domain,y,eta,beta,context,"test/block-1");
}

TEST(JointLocalConditioningTest, FullRankBlockUsesCompleteProfileDerivative)
{
    n::Domain domain(6,{{{0,0.0},{1,.25},{2,1.0},{3,2.25},{4,4.0},{5,6.25}}});
    n::Vector y(6); y<<.8,-.2,.5,1.1,-.6,.3;
    auto context=tc::MakeContext(y,1);
    n::Vector eta(1),beta(2); eta<<std::log(.8); beta<<1.2,.3;
    const auto result=Assess(domain,y,eta,beta,context);
    ASSERT_EQ(result.evidence.status,"assessed");
    EXPECT_EQ(result.evidence.free_ac_columns,2u);
    ASSERT_TRUE(result.evidence.free_ac_rank);
    EXPECT_EQ(*result.evidence.free_ac_rank,2u);
    ASSERT_TRUE(result.evidence.free_ac_condition);
    EXPECT_TRUE(std::isfinite(*result.evidence.free_ac_condition));
    EXPECT_EQ(result.evidence.free_ac_status,"full-rank");
    EXPECT_GT(*result.evidence.free_ac_rank_threshold,0.0);
    EXPECT_TRUE(result.evidence.free_ac_threshold_sensitive.has_value());
    EXPECT_TRUE(std::isfinite(*result.evidence.projected_width_sigma_min));
    EXPECT_TRUE(std::isfinite(*result.evidence.profile_sigma_max));
    EXPECT_EQ(result.evidence.conditioning_classification,"WellConditioned");
    ASSERT_TRUE(result.evidence.projected_width_rank);
    EXPECT_EQ(*result.evidence.projected_width_rank,1u);
    ASSERT_TRUE(result.evidence.projected_width_condition);
    ASSERT_TRUE(result.evidence.profile_condition);
    EXPECT_EQ(result.evidence.profile_rank,1u);

    const auto state=n::EvaluateState(domain,y,eta,beta,context);
    ASSERT_TRUE(state.valid);
    const auto width_kkt=n::ProfileGradientInfinityNorm(state,context);
    ASSERT_TRUE(result.evidence.local_ac_projected_kkt);
    ASSERT_TRUE(result.evidence.local_projected_width_kkt);
    ASSERT_TRUE(width_kkt);
    EXPECT_DOUBLE_EQ(*result.evidence.local_ac_projected_kkt,state.certificate.projected_kkt);
    EXPECT_DOUBLE_EQ(*result.evidence.local_projected_width_kkt,*width_kkt);
    const auto complete=n::PrepareDerivative(state,context.scale,&context);
    const auto projected=n::PrepareProjectedWidthDerivative(state,context.scale,&context);
    ASSERT_TRUE(complete.valid); ASSERT_TRUE(projected.valid);
    EXPECT_GT(complete.correction.norm(),1e-12);
    n::Matrix width_component,profile_jacobian;
    complete.Rows(0,domain.rows,width_component,profile_jacobian);
    EXPECT_GT((width_component-profile_jacobian).norm(),1e-12);
    EXPECT_GT(result.evidence.workspace_peak_memory_estimate_bytes.value_or(0u),0u);
}

TEST(JointLocalConditioningTest, ExactSingularityIsReportedWithoutInventingAConditionNumber)
{
    const std::vector<std::vector<n::Support>> supports{
        {{0,0.0},{1,.5},{2,1.0},{3,1.5},{4,2.0}},
        {{0,0.0},{1,.5},{2,1.0},{3,1.5},{4,2.0}}};
    n::Domain domain(5,supports);
    n::Vector y=n::Vector::Zero(5),eta=n::Vector::Constant(2,std::log(.7)),beta(4);
    beta<<1.0,.2,1.4,-.1;
    auto context=tc::MakeContext(y,2);
    const auto result=Assess(domain,y,eta,beta,context);
    ASSERT_EQ(result.evidence.status,"partial");
    ASSERT_TRUE(result.evidence.free_ac_rank);
    EXPECT_LT(*result.evidence.free_ac_rank,result.evidence.free_ac_columns);
    EXPECT_FALSE(result.evidence.free_ac_condition);
    EXPECT_EQ(result.evidence.free_ac_status,"rank-deficient");
    EXPECT_EQ(result.evidence.conditioning_classification,"Singular");
    ASSERT_TRUE(result.evidence.projected_width_rank);
    EXPECT_EQ(*result.evidence.projected_width_rank,2u);
    EXPECT_TRUE(result.evidence.projected_width_condition);
    EXPECT_NE(result.evidence.profile_status,"full-rank");
}

TEST(JointLocalConditioningTest, FullRankButIllConditionedBlockIsClassifiedSeparately)
{
    constexpr double offset=1e-8;
    n::Domain domain(5,{{{0,.1},{1,.5},{2,1.0},{3,1.5},{4,2.0}},
        {{0,.1+offset},{1,.5-offset},{2,1.0+offset},{3,1.5-offset},{4,2.0+offset}}});
    n::Vector y=n::Vector::Zero(5),eta=n::Vector::Constant(2,std::log(.7)),beta(4);
    beta<<1.0,.2,1.4,-.1;
    auto context=tc::MakeContext(y,2);
    const auto result=Assess(domain,y,eta,beta,context);
    ASSERT_TRUE(result.evidence.free_ac_rank);
    ASSERT_TRUE(result.evidence.free_ac_condition);
    EXPECT_EQ(*result.evidence.free_ac_rank,result.evidence.free_ac_columns);
    EXPECT_GT(*result.evidence.free_ac_condition,1e7);
    EXPECT_GE(*result.evidence.free_ac_condition*std::numeric_limits<double>::epsilon(),
        n::LocalConditioningRoundoffAmplificationWarning);
    EXPECT_EQ(result.evidence.free_ac_status,"full-rank");
    EXPECT_EQ(result.evidence.free_ac_conditioning,"IllConditioned");
    EXPECT_EQ(result.evidence.conditioning_classification,"IllConditioned");
}

TEST(JointLocalConditioningTest, IllConditioningWarnsOnPointQualityWithoutChangingConvergence)
{
    rhbm_gem::core::JointProblemInput input;
    input.atom_ids={"a","b"}; input.row_ids={"0","1","2","3","4"};
    input.observations.assign(5,0.0);
    input.support={{{0,.1},{1,.5},{2,1.0},{3,1.5},{4,2.0}},
        {{0,.1+1e-8},{1,.5-1e-8},{2,1.0+1e-8},{3,1.5-1e-8},{4,2.0+1e-8}}};
    input.selection_domain.emplace(); input.selection_domain->target_indices={0,1};
    const rhbm_gem::core::JointProblem problem(std::move(input));
    rhbm_gem::JointState state; state.ac={1.0,.2,1.4,-.1}; state.tau={.7,.7};
    state.log_tau={std::log(.7),std::log(.7)}; state.eta_gradient={0.,0.}; state.objective=0.;
    rhbm_gem::JointBlockConditioning conditioning; conditioning.block_id="component/block-1";
    conditioning.status="assessed"; conditioning.core_atoms={0,1};
    conditioning.conditioning_classification="IllConditioned";
    const auto evidence=n::AssessTargetPointState(problem,problem.ParameterLayout(),state,
        rhbm_gem::JointCheckStatus::Passed,{conditioning},rhbm_gem::JointEvidenceScope::AssembledGlobal,5);
    ASSERT_EQ(evidence.point_atoms.size(),2u);
    EXPECT_EQ(evidence.point_atoms[0].local_conditioning,"Warning");
    EXPECT_EQ(evidence.point_atoms[0].point_quality,"Warning");
    EXPECT_EQ(n::TargetConvergence(evidence,true),rhbm_gem::JointCheckStatus::Passed);
}

TEST(JointLocalConditioningTest, ActiveAmplitudeBoundaryIsKeptAsConditionalEvidence)
{
    n::Domain domain(5,{{{0,0.0},{1,.5},{2,1.0},{3,1.5},{4,2.0}}});
    n::Vector y=n::Vector::Zero(5),eta(1),beta(2); eta<<std::log(.7); beta<<0.,.2;
    auto context=tc::MakeContext(y,1);
    const auto result=Assess(domain,y,eta,beta,context);
    ASSERT_EQ(result.evidence.status,"assessed");
    EXPECT_EQ(result.evidence.active_amplitude_count,1u);
    EXPECT_EQ(result.evidence.free_ac_columns,1u);
    EXPECT_EQ(result.evidence.profile_status,"active-set-boundary");
    EXPECT_EQ(result.evidence.profile_reason,"profile-derivative-may-switch-active-face");
    EXPECT_EQ(result.evidence.projected_width_status,"active-face-full-rank");
    EXPECT_EQ(result.evidence.profile_conditioning,"ActiveFace");
    EXPECT_EQ(result.evidence.conditioning_classification,"ActiveFace");
}

TEST(JointLocalConditioningTest, ZeroAmplitudeWidthInformationDependsOnCharge)
{
    n::Domain domain(6,{{{0,0.0},{1,.25},{2,1.0},{3,2.25},{4,4.0},{5,6.25}}});
    const n::Vector eta=n::Vector::Constant(1,std::log(.7));
    for(double charge:{0.0,.2})
    {
        n::Vector beta(2); beta<<0.0,charge;
        n::Vector y(6);
        for(const auto & support:domain.atoms[0])
            y(support.row)=charge*n::EvaluateKernel(support.square,.7,2.5).charge;
        auto context=tc::MakeContext(y,1); context.eta_upper_bound=std::log(2.2);
        const auto state=n::EvaluateState(domain,y,eta,beta,context);
        ASSERT_TRUE(state.valid); EXPECT_TRUE(state.certificate.kkt_passed);
        EXPECT_LE(*n::ProfileGradientInfinityNorm(state,context),1e-9);
        const auto prepared=n::PrepareProjectedWidthDerivative(state,context.scale,&context);
        const auto projected=n::ReduceProjectedWidthDerivative(prepared);
        ASSERT_TRUE(prepared.valid); ASSERT_TRUE(projected.valid);
        const auto conditioning=Assess(domain,y,eta,beta,context).evidence;
        EXPECT_EQ(conditioning.active_amplitude_count,1u);
        EXPECT_EQ(conditioning.profile_status,"active-set-boundary");
        if(charge==0)
        {
            EXPECT_EQ(prepared.raw.norm(),0.0);
            EXPECT_EQ(projected.projected_norms(0),0.0);
            EXPECT_EQ(conditioning.conditioning_classification,"Singular");
        }
        else
        {
            EXPECT_GT(prepared.raw.norm(),.01);
            EXPECT_GT(projected.projected_norms(0),.001);
            EXPECT_EQ(conditioning.projected_width_rank,1u);
            EXPECT_EQ(conditioning.conditioning_classification,"ActiveFace");
        }
        auto varied=eta; varied(0)+=std::log(1.25);
        const auto profiled=n::EvaluateProfile(domain,y,varied,&context);
        ASSERT_TRUE(profiled.valid); EXPECT_TRUE(profiled.certificate.feasible);
        const n::Vector prediction=profiled.x*profiled.beta;
        if(charge==0) EXPECT_EQ(prediction.norm(),0.0);
        else EXPECT_GT((prediction-y).norm(),1e-4);
    }
}

TEST(JointLocalConditioningTest, NormalizedConditioningDoesNotMeasureAbsoluteWidthInformation)
{
    n::Domain domain(6,{{{0,0.0},{1,.25},{2,1.0},{3,2.25},{4,4.0},{5,6.25}}});
    const n::Vector eta=n::Vector::Constant(1,std::log(.8));
    double reference{};
    for(double scale:{1.0,1e-12})
    {
        n::Vector beta(2); beta<<1.2*scale,.3*scale;
        n::Vector y(6);
        for(const auto & s:domain.atoms[0])
        {
            const auto kernel=n::EvaluateKernel(s.square,.8,2.5);
            y(s.row)=beta(0)*kernel.gaussian+beta(1)*kernel.charge;
        }
        auto context=tc::MakeContext(y,1);
        const auto conditioning=Assess(domain,y,eta,beta,context).evidence;
        EXPECT_EQ(conditioning.conditioning_classification,"WellConditioned");
        ASSERT_EQ(conditioning.projected_width_column_norms.size(),1u);
        const auto norm=conditioning.projected_width_column_norms.front();
        if(scale==1) reference=norm;
        else EXPECT_LT(norm,reference*1e-10);
    }
}

TEST(JointLocalConditioningTest, ZeroSupportAndNonfiniteInputAreUnavailable)
{
    n::Vector empty_y(0); auto empty_context=tc::MakeContext(empty_y,1);
    n::Domain empty_domain(0,{{}});
    n::Vector eta(1),beta(2); eta<<0.; beta<<1.,.2;
    const auto empty=Assess(empty_domain,empty_y,eta,beta,empty_context);
    EXPECT_EQ(empty.evidence.status,"unavailable");
    EXPECT_EQ(empty.evidence.reason,"zero-support");

    n::Domain domain(5,{{{0,0.0},{1,.5},{2,1.0},{3,1.5},{4,2.0}}});
    n::Vector y=n::Vector::Zero(5); auto context=tc::MakeContext(y,1);
    eta(0)=std::numeric_limits<double>::quiet_NaN();
    const auto nonfinite=Assess(domain,y,eta,beta,context);
    EXPECT_EQ(nonfinite.evidence.status,"unavailable");
    EXPECT_EQ(nonfinite.evidence.reason,"nonfinite-local-state");

    eta(0)=0.; context.rank.rows=-1;
    const auto decomposition=Assess(domain,y,eta,beta,context);
    EXPECT_EQ(decomposition.evidence.status,"partial");
    EXPECT_EQ(decomposition.evidence.free_ac_status,"full-rank");
    EXPECT_EQ(decomposition.evidence.projected_width_status,"unavailable");
    EXPECT_EQ(decomposition.evidence.projected_width_reason,"nonfinite-derivative");
    EXPECT_EQ(decomposition.evidence.conditioning_classification,"Unavailable");
}

TEST(JointLocalConditioningTest, MultiBlockDiagnosticsLeaveTheAcceptedFitUnchanged)
{
    const c::JointProblem problem(joint_component_test::SyntheticJointWorkload("chain",2));
    const auto & data=c::JointProblemAccess::Get(problem);
    const auto run=[&](bool diagnostics) {
        n::BcdSearchPolicy policy; policy.core_atoms=1;
        policy.collect_diagnostics=diagnostics;
        policy.eta_upper_bound=std::log(2.2);
        const n::Vector eta=n::Vector::Constant(2,std::log(.55));
        return n::SearchBcdComponent(data.domain,data.layout,data.y,eta,data.context,policy);
    };
    const auto with_diagnostics=run(true),without_diagnostics=run(false);
    ASSERT_EQ(with_diagnostics.local_conditioning.size(),2u);
    EXPECT_EQ(without_diagnostics.local_conditioning.size(),2u);
    auto context=data.context;
    context.eta_upper_bound=std::log(2.2);
    const auto compensation=n::AssessAcCompensation(data.domain,data.y,with_diagnostics.state.eta,
        with_diagnostics.state.beta,context,0);
    EXPECT_TRUE(compensation.raw_correlation)<<compensation.reason;
    EXPECT_FALSE(with_diagnostics.blocks.empty());
    EXPECT_TRUE(without_diagnostics.blocks.empty());
    EXPECT_EQ(with_diagnostics.state.eta,without_diagnostics.state.eta);
    EXPECT_EQ(with_diagnostics.state.beta,without_diagnostics.state.beta);
    EXPECT_EQ(with_diagnostics.state.prediction,without_diagnostics.state.prediction);
    EXPECT_EQ(with_diagnostics.state.residual,without_diagnostics.state.residual);
    EXPECT_DOUBLE_EQ(with_diagnostics.state.objective,without_diagnostics.state.objective);
    EXPECT_EQ(with_diagnostics.search_converged,without_diagnostics.search_converged);
    EXPECT_EQ(with_diagnostics.reason,without_diagnostics.reason);
}

n::Domain CompensationDomain(double spacing,bool neighbor=false)
{
    std::vector<std::vector<n::Support>> supports(neighbor ? 2 : 1);
    for(Eigen::Index i=0;i<24;++i)
    {
        const double radius=static_cast<double>(i)*spacing;
        supports[0].push_back({i,radius*radius});
        if(neighbor) supports[1].push_back({i,std::pow(radius-.003,2)});
    }
    return n::Domain(24,std::move(supports));
}

void CheckCompensationPerturbation(const n::Domain & domain,const n::Vector & eta,n::Vector beta,
    const std::string & label)
{
    n::Vector y=n::Vector::Zero(domain.rows);
    auto context=tc::MakeContext(y,eta.size());
    const auto geometry=n::EvaluateState(domain,y,eta,beta,context);
    ASSERT_TRUE(geometry.valid);
    if(label=="collinear")
        beta(1)=-beta(0)*geometry.x.col(0).norm()/geometry.x.col(1).norm();
    y=geometry.x*beta;
    context=tc::MakeContext(y,eta.size());
    const auto endpoint=n::EvaluateState(domain,y,eta,beta,context);
    const auto d=n::AssessAcCompensation(domain,y,eta,beta,context,0);
    ASSERT_TRUE(d.reason.empty())<<d.reason;
    ASSERT_TRUE(d.raw_correlation && d.conditional_correlation && d.cancellation_ratio);
    ASSERT_TRUE(d.residualized_gaussian_norm && d.residualized_charge_norm);
    ASSERT_EQ(d.physical_weak_direction.size(),beta.size());
    const n::Vector delta=.01*d.physical_weak_direction;
    for(Eigen::Index j=0;j<beta.size();j+=2) ASSERT_GT(beta(j)+delta(j),0);
    const auto perturbed=n::EvaluateState(domain,y,eta,beta+delta,context);
    ASSERT_TRUE(perturbed.valid);
    const double sensitivity=(endpoint.x*delta).norm();
    const auto before=n::NormalizedProfileObjective(endpoint,context.scale);
    const auto after=n::NormalizedProfileObjective(perturbed,context.scale);
    ASSERT_TRUE(before && after);
    EXPECT_NEAR(sensitivity,.01*d.singular_values.tail(1)(0),3e-16);
    EXPECT_NEAR(*after-*before,.5*std::pow(sensitivity/context.scale,2),3e-16);
    if(label=="collinear")
    {
        EXPECT_GT(*d.raw_correlation,.99999);
        EXPECT_GT(*d.cancellation_ratio,1e3);
        EXPECT_GT(d.weak_cancellation_alignment.value_or(0),.999);
        EXPECT_GT(std::abs(delta(0)),.1);
        EXPECT_LT(sensitivity,3e-5);
    }
    else
    {
        EXPECT_LT(*d.raw_correlation,.9);
        EXPECT_LT(*d.cancellation_ratio,1.3);
        if(label=="isolated")
        {
            EXPECT_GT(sensitivity,.003);
            EXPECT_NEAR(*d.raw_correlation,*d.conditional_correlation,1e-15);
            EXPECT_NEAR(*d.residualized_gaussian_norm,*d.gaussian_norm,1e-14);
            EXPECT_NEAR(*d.residualized_charge_norm,*d.charge_norm,1e-14);
        }
        else
        {
            EXPECT_LT(sensitivity,5e-5);
            ASSERT_EQ(d.nuisance_rank,2);
            EXPECT_LT(*d.residualized_gaussian_norm,*d.gaussian_norm*.03);
            EXPECT_LT(*d.residualized_charge_norm,*d.charge_norm*.03);
            EXPECT_GT(d.weak_atom_shares.minCoeff(),.4);
            const n::Matrix x(endpoint.x);
            const Eigen::HouseholderQR<n::Matrix> qr(x.rightCols(2));
            const n::Matrix orthogonal=qr.householderQ().adjoint()*x.leftCols(2);
            const auto projected=orthogonal.bottomRows(domain.rows-2);
            EXPECT_NEAR(projected.col(0).norm(),*d.residualized_gaussian_norm,1e-13);
            EXPECT_NEAR(projected.col(1).norm(),*d.residualized_charge_norm,1e-13);
            EXPECT_NEAR(projected.col(0).normalized().dot(projected.col(1).normalized()),
                *d.conditional_correlation,1e-10);
        }
    }
    std::ostringstream metrics;
    metrics<<std::setprecision(12)<<"beta "<<beta.transpose()<<" rho "<<*d.raw_correlation
        <<" conditional "<<*d.conditional_correlation<<" residual norms "<<*d.residualized_gaussian_norm
        <<" "<<*d.residualized_charge_norm<<" cancellation "<<*d.cancellation_ratio
        <<" singular values "<<d.singular_values.transpose()<<" delta beta "<<delta.transpose()
        <<" delta prediction "<<sensitivity<<" delta objective "<<*after-*before;
    ::testing::Test::RecordProperty(label,metrics.str());
}

TEST(JointLocalConditioningTest, CompensationSeparatesStrongSignalFromCollinearityAndConfounding)
{
    const n::Vector isolated_eta=n::Vector::Constant(1,std::log(.4));
    n::Vector beta(2); beta<<120.,.35;
    CheckCompensationPerturbation(CompensationDomain(.1),isolated_eta,beta,"isolated");
    CheckCompensationPerturbation(CompensationDomain(.01),n::Vector::Constant(1,std::log(2.2)),beta,"collinear");
    n::Vector paired(4); paired<<120.,.35,50.,.2;
    CheckCompensationPerturbation(CompensationDomain(.1,true),n::Vector::Constant(2,std::log(.4)),paired,"neighbor");
}

TEST(JointLocalConditioningTest, CompensationRespectsActiveAmplitudeAndWidthFaces)
{
    const auto domain=CompensationDomain(.1);
    const n::Vector y=n::Vector::Zero(domain.rows),eta=n::Vector::Constant(1,std::log(2.2));
    auto context=tc::MakeContext(y,1); context.eta_upper_bound=eta(0);
    n::Vector beta(2); beta<<0.,-.3;
    const auto active=n::AssessAcCompensation(domain,y,eta,beta,context,0);
    EXPECT_TRUE(active.amplitude_active && active.width_upper_bound_active && active.conditional_on_active_face);
    EXPECT_TRUE(active.raw_correlation);
    EXPECT_FALSE(active.conditional_correlation || active.residualized_gaussian_norm);
    EXPECT_EQ(active.conditional_reason,"conditional-on-active-amplitude-face");
    EXPECT_EQ(active.residualized_charge_norm,active.charge_norm);
    EXPECT_EQ(active.cancellation_ratio,1.);
    ASSERT_EQ(active.physical_weak_direction.size(),2);
    EXPECT_EQ(active.physical_weak_direction(0),0.);
    beta<<120.,0.;
    const auto zero_charge=n::AssessAcCompensation(domain,y,eta,beta,context,0);
    EXPECT_FALSE(zero_charge.amplitude_active || zero_charge.conditional_on_active_face);
    EXPECT_TRUE(zero_charge.width_upper_bound_active);
    EXPECT_EQ(zero_charge.free_rank,2);
    EXPECT_EQ(zero_charge.cancellation_ratio,1.);
}

TEST(JointLocalConditioningTest, CompensationUsesTruncatedNuisanceSpanAndPreservesTheWeakSubspace)
{
    const std::vector<n::Support> support{{0,0.},{1,.5},{2,1.5}};
    const n::Domain domain(3,{support,support,support});
    const n::Vector y=n::Vector::Zero(3),eta=n::Vector::Constant(3,std::log(.7)),beta=n::Vector::Ones(6);
    const auto d=n::AssessAcCompensation(domain,y,eta,beta,tc::MakeContext(y,3),0);
    EXPECT_EQ(d.nuisance_rank,2);
    EXPECT_EQ(d.free_rank,2);
    EXPECT_EQ(d.weak_dimension,4);
    EXPECT_EQ(d.physical_weak_direction.size(),0);
    EXPECT_NEAR(d.weak_atom_shares.sum(),1.,1e-14);
    EXPECT_EQ(d.conditional_reason,"residualized-column-at-roundoff");
    EXPECT_FALSE(d.conditional_correlation);
    ASSERT_TRUE(d.residualized_gaussian_norm && d.residualized_charge_norm);
    EXPECT_LT(*d.residualized_gaussian_norm,1e-14);
    EXPECT_LT(*d.residualized_charge_norm,1e-14);
}

TEST(JointLocalConditioningTest, CompensationReportsZeroTinyAndInvalidEvidence)
{
    const n::Vector y=n::Vector::Zero(2),eta=n::Vector::Zero(1);
    const auto context=tc::MakeContext(y,1);
    n::Vector beta(2); beta<<1.,0.;
    const auto empty=n::AssessAcCompensation(n::Domain(2,{{}}),y,eta,beta,context,0);
    EXPECT_EQ(empty.reason,"empty-target-support");
    const auto zero=n::AssessAcCompensation(n::Domain(2,{{{0,9.},{1,10.}}}),y,eta,beta,context,0);
    EXPECT_EQ(zero.gaussian_norm,0.);
    EXPECT_FALSE(zero.raw_correlation || zero.cancellation_ratio);
    EXPECT_EQ(zero.free_rank,0);
    EXPECT_FALSE(zero.weak_reason.empty());
    const n::Domain domain(2,{{{0,0.},{1,0.}}});
    const auto tiny=n::AssessAcCompensation(domain,y,n::Vector::Constant(1,std::log(1e103)),beta,context,0);
    ASSERT_TRUE(tiny.gaussian_norm);
    EXPECT_GT(*tiny.gaussian_norm,0.);
    EXPECT_LT(*tiny.gaussian_norm,std::numeric_limits<double>::min());
    EXPECT_FALSE(tiny.raw_correlation || tiny.cancellation_ratio);
    EXPECT_TRUE(tiny.residualized_charge_norm);
    const auto kernel=n::EvaluateKernel(0.,1.,2.5);
    beta(1)=-kernel.gaussian/kernel.charge;
    for(double c:{beta(1),std::nextafter(beta(1),0.)})
    {
        beta(1)=c;
        const auto near_zero=n::AssessAcCompensation(domain,y,eta,beta,context,0);
        EXPECT_FALSE(near_zero.cancellation_ratio);
        EXPECT_EQ(near_zero.cancellation_reason,"net-prediction-at-roundoff");
        EXPECT_EQ(near_zero.free_rank,1);
        EXPECT_EQ(near_zero.weak_dimension,1);
    }
    const auto refused=n::AssessAcCompensation(domain,y,eta,beta,context,0,n::JointMemoryHardLimitBytes);
    EXPECT_EQ(refused.reason,"resource-limit:local-conditioning");
    EXPECT_FALSE(refused.raw_correlation);
    beta(0)=std::numeric_limits<double>::infinity();
    EXPECT_EQ(n::AssessAcCompensation(domain,y,eta,beta,context,0).reason,"invalid-local-state");
    beta(0)=1.;
    const n::Domain bad(2,{{{0,std::numeric_limits<double>::quiet_NaN()},{1,0.}}});
    EXPECT_EQ(n::AssessAcCompensation(bad,y,eta,beta,context,0).reason,"nonfinite-basis");
}
}
