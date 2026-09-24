#include <gtest/gtest.h>
#include "core/detail/joint_component/OperatorSearch.hpp"
#include "core/detail/joint_component/Problem.hpp"
#include "support/JointOperatorWorkload.hpp"
#include "support/JointDenseReference.hpp"
#include <map>

namespace {
namespace n=rhbm_gem::core::joint_component;
namespace c=rhbm_gem::core;
struct Sample
{
    c::JointProblem problem{second_stage_test::OperatorWorkload("chain",4)};
    const n::ProblemData & data=c::JointProblemAccess::Get(problem);
    n::Evaluation e=n::EvaluateProfile(data.domain,data.y,n::Vector::Constant(4,std::log(.55)),false,&data.context);
};
using Signature=std::map<std::string,std::pair<std::vector<std::string>,std::vector<std::string>>>;
Signature Describe(const n::PreconditionerPartition & p)
{
    Signature out;
    for(const auto & b:p.blocks)
    {
        auto & [core,overlap]=out[b.id];
        for(auto a:b.core_atoms) core.push_back(p.problem->atom_ids[static_cast<std::size_t>(a)]);
        for(auto a:b.overlap_atoms) overlap.push_back(p.problem->atom_ids[static_cast<std::size_t>(a)]);
    }
    return out;
}
}
TEST(JointOperatorSearchTest, PcgMatchesIndependentDenseSystemAndDetectsFailure)
{
    n::Matrix a(3,3); a<<4,1,0,1,3,1,0,1,2;
    n::Vector rhs(3); rhs<<1,-2,3; const n::Vector metric=n::Vector::LinSpaced(3,1,2);
    auto result=n::SolvePcg([&](n::VectorRef x)->n::Vector{return a*x;},[](n::VectorRef x)->n::Vector{return x;},rhs,metric);
    ASSERT_TRUE(result.valid)<<result.reason;
    EXPECT_LT((result.step-a.llt().solve(rhs)).norm(),1e-10); EXPECT_LE(result.relative_residual,1e-10);
    EXPECT_EQ(n::SolvePcg([&](n::VectorRef x)->n::Vector{return a*x;},[](n::VectorRef x)->n::Vector{return x;},rhs,metric,0).reason,"pcg-iteration-budget");
    EXPECT_EQ(n::SolvePcg([](n::VectorRef x)->n::Vector{return -x;},[](n::VectorRef x)->n::Vector{return x;},rhs,metric).reason,"pcg-nonpositive-curvature");
    EXPECT_EQ(n::SolvePcg([](n::VectorRef x)->n::Vector{return x;},[](n::VectorRef x)->n::Vector{return -x;},rhs,metric).reason,"pcg-nonpositive-preconditioner");
    EXPECT_EQ(n::SolvePcg([](n::VectorRef x)->n::Vector{return n::Vector::Constant(x.size(),n::unavailable);},[](n::VectorRef x)->n::Vector{return x;},rhs,metric).reason,"pcg-nonfinite");
}
TEST(JointOperatorSearchTest, FrozenTopologyDeterminismMappingsAndLimits)
{
    auto input=std::make_shared<c::JointProblemInput>(second_stage_test::OperatorWorkload("chain",8));
    n::PreconditionerLimits limits; limits.core_atoms=2;
    const auto original=n::BuildPreconditionerPartition(input,n::BuildParameterLayout(*input),limits);
    auto reordered=std::make_shared<c::JointProblemInput>(*input);
    std::reverse(reordered->atom_ids.begin(),reordered->atom_ids.end()); std::reverse(reordered->support.begin(),reordered->support.end());
    const auto permuted=n::BuildPreconditionerPartition(reordered,n::BuildParameterLayout(*reordered),limits);
    EXPECT_EQ(Describe(*original),Describe(*permuted)); ASSERT_EQ(original->blocks.size(),4);
    for(const auto & block:original->blocks) {EXPECT_EQ(block.core_atoms.size(),2); EXPECT_FALSE(block.overlap_atoms.empty());}
    auto w=n::WidthMapping(*original); n::Vector beta=n::Vector::Ones(16); beta(0)=0; beta(1)=-2;
    auto ac=n::FreeColumnMapping(*original,beta); EXPECT_EQ(ac.dimension,15); EXPECT_EQ(w.dimension,8);
    n::Vector summed=n::Vector::Zero(8),v=n::Vector::LinSpaced(8,1,2);
    for(const auto & b:w.blocks)
    {
        b.Scatter(b.Restrict(v),summed);
        for(std::size_t j=0;j<b.global.size();++j) EXPECT_EQ(b.LocalIndex(b.global[j]),j);
        EXPECT_EQ(b.LocalIndex(-1),-1);
    }
    EXPECT_LT((summed-v).norm(),1e-12);
    limits.block_atoms=2; EXPECT_THROW(n::BuildPreconditionerPartition(input,n::BuildParameterLayout(*input),limits),std::runtime_error);
    limits.block_atoms=512; limits.scratch_bytes=1; EXPECT_THROW(n::BuildPreconditionerPartition(input,n::BuildParameterLayout(*input),limits),std::runtime_error);
}
TEST(JointOperatorSearchTest, SchwarzSpdAndThreeStepSolversHaveSameGlobalSolution)
{
    Sample s; ASSERT_TRUE(s.e.valid); const n::ProfileJacobianOperator op(s.e,s.data.context); ASSERT_TRUE(op.Valid());
    const auto norms=n::WidthNorms(n::RawWidthDerivative(s.e),s.data.context.scale);
    const n::PreconditionerContext pc{op.Identity(),n::PreconditionerSpace::Width,n::WidthMetric(norms),.02};
    n::PreconditionerLimits limits; limits.core_atoms=2;
    const auto partition=n::SearchPartition(s.data.domain,s.data.context,limits);
    const n::SchwarzModel model(*partition,s.e,s.data.context.scale,pc,limits); const n::SchwarzPreconditioner schwarz(model,pc);
    const n::Vector u=n::Vector::LinSpaced(4,-1,2),v=n::Vector::LinSpaced(4,.3,-.8);
    EXPECT_GT(u.dot(schwarz.ApplyInverse(u,pc)),0);
    EXPECT_NEAR(u.dot(schwarz.ApplyInverse(v,pc)),v.dot(schwarz.ApplyInverse(u,pc)),1e-10*std::max(1.,schwarz.ApplyInverse(u,pc).norm()));
    n::Matrix j(op.Rows(),4); for(int k=0;k<4;++k) j.col(k)=op.Apply(n::Vector::Unit(4,k));
    n::Matrix h=j.transpose()*j; h.diagonal().array()+=pc.damping*pc.metric.array().square();
    const n::Vector gradient=op.ApplyAdjoint(s.e.residual/s.data.context.scale),expected=h.llt().solve(-gradient);
    for(int kind=0;kind<3;++kind)
    {
        const auto action=[&](n::VectorRef x)->n::Vector {
            if(kind==2) return schwarz.ApplyInverse(x,pc);
            if(kind==1) return x.array()/(norms.array().square()+pc.damping*pc.metric.array().square());
            return x;
        };
        const auto composed=n::SolvePcg([&](n::VectorRef x)->n::Vector {
            return op.ApplyAdjoint(op.Apply(x))+(pc.damping*pc.metric.array().square()*x.array()).matrix();
        },action,-gradient,pc.metric);
        ASSERT_TRUE(composed.valid)<<composed.reason;
        const auto result=n::WidthStepSolver(op,gradient,pc,action); ASSERT_TRUE(result.valid)<<result.reason;
        EXPECT_LT((result.step-composed.step).norm(),1e-10*(1+composed.step.norm()));
        EXPECT_LT((result.step-expected).norm(),1e-10*(1+expected.norm()));
        EXPECT_NEAR(result.predicted,-gradient.dot(expected)-.5*(j*expected).squaredNorm(),1e-12);
    }
    auto stale=pc; stale.damping*=2; EXPECT_THROW(schwarz.ApplyInverse(u,stale),std::logic_error);
    const n::SchwarzPreconditioner shifted(model,stale); EXPECT_GT(u.dot(shifted.ApplyInverse(u,stale)),0);
    stale.metric(0)*=2; EXPECT_THROW(n::SchwarzPreconditioner(model,stale),std::logic_error);
    stale=pc; stale.linearization=std::make_shared<const n::LinearizationIdentity>(); EXPECT_THROW(schwarz.ApplyInverse(u,stale),std::logic_error);
    EXPECT_FALSE(n::WidthStepSolver(op,gradient,stale,[](n::VectorRef x)->n::Vector{return x;}).valid);
    auto changed=s.e; changed.beta(0)=0;
    const n::SchwarzModel new_face(*partition,changed,s.data.context.scale,stale,limits);
    const n::SchwarzPreconditioner new_inverse(new_face,stale);
    EXPECT_EQ(n::FreeColumnMapping(*partition,changed.beta).dimension,7);
    EXPECT_GT(u.dot(new_inverse.ApplyInverse(u,stale)),0);
}
TEST(JointOperatorSearchTest, LocalRankFailureIsRegularizedWithoutGlobalRankWork)
{
    Sample s; s.e.x.col(1)=s.e.x.col(0); const auto before=n::SparseWorkForTesting();
    n::SearchWorkForTesting()={}; const n::PreconditionerContext pc{std::make_shared<const n::LinearizationIdentity>(),n::PreconditionerSpace::Width,n::Vector::Ones(4),1e-3};
    const auto partition=n::SearchPartition(s.data.domain,s.data.context);
    const bool audit=n::ResourceWorkForTesting().enabled; n::ResourceWorkForTesting().enabled=true;
    const n::SchwarzModel model(*partition,s.e,s.data.context.scale,pc); const n::SchwarzPreconditioner inverse(model,pc);
    n::ResourceWorkForTesting().enabled=audit;
    EXPECT_GT(n::Vector::Ones(4).dot(inverse.ApplyInverse(n::Vector::Ones(4),pc)),0);
    EXPECT_GT(n::SearchWorkForTesting().maximum_lambda,0); EXPECT_GT(n::SearchWorkForTesting().maximum_tau,0);
    const auto & records=n::SearchWorkForTesting().regularizations;
    ASSERT_EQ(records.size(),2); EXPECT_EQ(records[0].factor_build,0);
    EXPECT_DOUBLE_EQ(records[0].lambda,records[1].lambda); EXPECT_GT(records[1].tau,0);
    EXPECT_DOUBLE_EQ(records[1].damping,pc.damping); EXPECT_EQ(records[1].attempt,1);
    EXPECT_EQ(n::SparseWorkForTesting().numeric,before.numeric); EXPECT_EQ(n::SparseWorkForTesting().free_design_svds,before.free_design_svds);
    n::PreconditionerLimits limits; limits.storage_bytes=1;
    EXPECT_THROW(n::SchwarzModel(*partition,s.e,s.data.context.scale,pc,limits),std::runtime_error);
}
TEST(JointOperatorSearchTest, LocalSchurMatchesDenseEliminationAndMetricFloor)
{
    Sample s;
    const auto raw=n::RawWidthDerivative(s.e);
    const auto norms=n::WidthNorms(raw,s.data.context.scale);
    const n::PreconditionerContext pc{std::make_shared<const n::LinearizationIdentity>(),n::PreconditionerSpace::Width,n::WidthMetric(norms),.01};
    const auto partition=n::SearchPartition(s.data.domain,s.data.context);
    const n::SchwarzModel model(*partition,s.e,s.data.context.scale,pc);
    ASSERT_EQ(model.Matrices().size(),1);
    const auto & order=model.Mapping().blocks.front().global;
    n::Matrix z(s.e.x.rows(),s.e.x.cols()),d(raw.rows(),raw.cols());
    for(std::size_t j=0;j<order.size();++j)
    {
        const auto a=order[j],local=static_cast<Eigen::Index>(j); ASSERT_GT(s.e.beta(2*a),0);
        for(Eigen::Index k=0;k<2;++k) z.col(2*local+k)=n::Vector(s.e.x.col(2*a+k))/s.e.x.col(2*a+k).norm();
        d.col(local)=n::Vector(raw.col(a))/(s.data.context.scale*pc.metric(a));
    }
    n::Matrix a=z.transpose()*z;
    const double lambda=std::sqrt(std::numeric_limits<double>::epsilon())*std::max(1.,a.cwiseAbs().rowwise().sum().maxCoeff());
    a.diagonal().array()+=lambda;
    const n::Matrix cross=z.transpose()*d,expected=d.transpose()*d-cross.transpose()*a.llt().solve(cross);
    EXPECT_LT((model.Matrices().front()-expected).norm(),1e-12*(1+expected.norm()));
    EXPECT_EQ(n::WidthMetric(n::Vector::Zero(4)),n::Vector::Ones(4));
    n::Vector disparate(3); disparate<<0,1e-30,2;
    const auto floored=n::WidthMetric(disparate);
    EXPECT_DOUBLE_EQ(floored(0),2*std::sqrt(std::numeric_limits<double>::epsilon()));
    EXPECT_DOUBLE_EQ(floored(1),floored(0)); EXPECT_DOUBLE_EQ(floored(2),2);
}
TEST(JointOperatorSearchTest, ProfileRowsAndParentMappingsRemainStructural)
{
    auto input=second_stage_test::OperatorWorkload("chain",4);
    input.selection_domain.emplace(); input.selection_domain->target_indices={0,1,2};
    input.support[3]={{7,0}};
    const c::JointProblem problem(input); const auto & data=c::JointProblemAccess::Get(problem);
    const auto domain=n::ProfileDomain(data.domain,data.layout);
    const auto context=n::ProfileContext(data.context,data.layout,data.domain.rows);
    const auto partition=n::SearchPartition(domain,context);
    EXPECT_EQ(partition->problem.get(),data.input.get()); EXPECT_EQ(partition->layout.full_atoms,data.layout.full_atoms);
    for(const auto & block:partition->blocks)
    {
        EXPECT_EQ(std::count(block.informative_rows.begin(),block.informative_rows.end(),7),0);
        for(auto a:block.core_atoms) EXPECT_LT(a,3);
    }
    EXPECT_DOUBLE_EQ(context.scale,data.context.scale); EXPECT_EQ(context.rank.rows,data.domain.rows);
}
TEST(JointOperatorSearchTest, SearchUsesReplayWithoutGlobalReductionAndPreservesBudgets)
{
    Sample s; auto context=s.data.context; context.search.method=n::SearchMethod::OperatorPcg;
    for(auto kind:{n::PreconditionerKind::Identity,n::PreconditionerKind::Diagonal,n::PreconditionerKind::Schwarz})
    {
        context.search.preconditioner=kind; n::SparseWorkForTesting()={};
        const auto result=n::SearchProfile(s.data.domain,s.data.y,n::Vector::Constant(4,.55),context);
        EXPECT_TRUE(result.initial_accepted); EXPECT_GT(result.accepted,0)<<result.stop_reason; EXPECT_EQ(result.references,0);
        EXPECT_EQ(n::SparseWorkForTesting().derivative_preparations,0); EXPECT_EQ(n::SparseWorkForTesting().reference_solves,0);
        for(const auto & trial:result.trials) if(trial.accepted) {ASSERT_TRUE(trial.trust); EXPECT_TRUE(trial.trust->passed); EXPECT_FALSE(trial.trust->reference);}
        for(const auto & trial:result.trials) if(trial.lm)
            EXPECT_LE(trial.lm->diagonal.cwiseProduct(trial.lm->step).stableNorm(),trial.lm->radius);
    }
    context.update_budget=0; EXPECT_EQ(n::SearchProfile(s.data.domain,s.data.y,n::Vector::Constant(4,.55),context).stop_reason,"accepted-update-budget");
    context.update_budget=100; context.search.damping_trials=0;
    auto stopped=n::SearchProfile(s.data.domain,s.data.y,n::Vector::Constant(4,.55),context);
    EXPECT_TRUE(stopped.stopped); EXPECT_TRUE(stopped.initial_accepted); EXPECT_EQ(stopped.stop_reason,"damping-trial-budget");
    EXPECT_EQ(stopped.eta,stopped.initial.eta);
    context.search.damping_trials=20; context.search.pcg_iterations=0;
    EXPECT_EQ(n::SearchProfile(s.data.domain,s.data.y,n::Vector::Constant(4,.55),context).stop_reason,"pcg-iteration-budget");
    context.search.pcg_iterations=-1; context.profile_budget=1;
    const auto last=n::SearchProfile(s.data.domain,s.data.y,n::Vector::Constant(4,.55),context);
    EXPECT_EQ(last.stop_reason,"profile-budget"); EXPECT_TRUE(last.initial_accepted);
    EXPECT_EQ(last.evaluations,1); EXPECT_EQ(last.accepted,0); EXPECT_EQ(last.eta,last.initial.eta);
    const auto invalid=n::SearchProfile(s.data.domain,s.data.y,n::Vector::Zero(4),context);
    EXPECT_EQ(invalid.stop_reason,"inner-invalid-input"); EXPECT_FALSE(invalid.initial_accepted);
    context.profile_budget=0; EXPECT_EQ(n::SearchProfile(s.data.domain,s.data.y,n::Vector::Constant(4,.55),context).stop_reason,"profile-budget");
}

#include "core/detail/JointUncertainty.hpp"
#include "data/io/detail/JointResultJson.hpp"
TEST(JointOperatorSearchTest, FullRuntimeEvidenceUncertaintyAndPersistenceWithOperatorPolicy)
{
    for(bool halo:{false,true})
    {
        auto input=second_stage_test::OperatorWorkload("chain",4);
        if(halo)
        {
            input.selection_domain.emplace(); input.selection_domain->target_indices={0,1,2};
            input.support[3]={{7,0}};
        }
        const c::JointProblem problem(input); const std::vector<double> initial(4,.55);
        const auto legacy=c::FitJointComponents(problem,initial);
        n::SearchPolicy policy; policy.method=n::SearchMethod::OperatorPcg;
        const auto actual=n::FitWithSearchPolicy(problem,initial,policy);
        ASSERT_EQ(actual.assembled_state.has_value(),legacy.assembled_state.has_value());
        EXPECT_EQ(actual.RuntimeConvergence(),legacy.RuntimeConvergence());
        ASSERT_TRUE(actual.objective); ASSERT_TRUE(legacy.objective); EXPECT_NEAR(*actual.objective,*legacy.objective,1e-12);
        const auto capture=c::CaptureJointAnalysisResult(actual);
        const auto encoded=rhbm_gem::joint_result_io::Encode(capture);
        EXPECT_EQ(rhbm_gem::joint_result_io::Encode(rhbm_gem::joint_result_io::Decode(encoded)),encoded);
        const auto expected=c::detail::ComputeJointUncertainty(problem,c::CaptureJointAnalysisResult(legacy));
        const auto measured=c::detail::ComputeJointUncertainty(problem,capture);
        ASSERT_EQ(expected.size(),measured.size());
        for(const auto & [id,u]:expected)
        {
            const auto & v=measured.at(id); EXPECT_EQ(u.status,v.status); EXPECT_EQ(u.covariance.has_value(),v.covariance.has_value());
            if(u.covariance && v.covariance) EXPECT_LT((*u.covariance-*v.covariance).norm(),1e-10*(1+u.covariance->norm()));
        }
    }
}
