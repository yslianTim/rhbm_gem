#include <gtest/gtest.h>
#include "core/detail/joint_component/ProfileJacobianOperator.hpp"
#include "core/detail/joint_component/Preconditioner.hpp"
#include "core/detail/joint_component/Problem.hpp"
#include "support/JointDenseReference.hpp"
#include "support/JointOperatorWorkload.hpp"
#include <cmath>

namespace {
namespace n=rhbm_gem::core::joint_component;
namespace c=rhbm_gem::core;
namespace p=second_stage_test::matched::joint_abc;
struct Sample
{
    n::Domain domain{60,{}};
    n::Vector y,eta;
    n::EvaluationContext context;
    Sample()
    {
        std::vector<std::vector<n::Support>> support(3);
        y=n::Vector::Zero(60); eta=n::Vector::Constant(3,std::log(.65));
        for(Eigen::Index r=0;r<60;++r)
        {
            for(Eigen::Index a=0;a<3;++a)
            {
                const double x=-2.4+.08*static_cast<double>(r)-.6*static_cast<double>(a);
                const double square=x*x+.03*static_cast<double>(a);
                if(square<=6.25)
                {
                    support[static_cast<std::size_t>(a)].push_back({r,square});
                    const auto b=n::EvaluateKernel(square,.5+.06*static_cast<double>(a),2.5);
                    y(r)+=(2+.2*static_cast<double>(a))*b.gaussian+(a%2 ? -.2 : .3)*b.charge;
                }
            }
            y(r)+=1e-5*std::sin(static_cast<double>(r));
        }
        domain=n::Domain(60,std::move(support)); context=n::CreateContext(y,3);
    }
};
void Parity(const n::Evaluation & e,const n::EvaluationContext & context)
{
    const n::ProfileJacobianOperator op(e,context); ASSERT_TRUE(op.Valid())<<op.Reason();
    const auto dense=p::DenseDifferentiate(e,context.scale,&context); ASSERT_TRUE(dense.valid);
    n::Matrix j(op.Rows(),op.Columns());
    for(Eigen::Index k=0;k<op.Columns();++k) j.col(k)=op.Apply(n::Vector::Unit(op.Columns(),k));
    EXPECT_LE((j-dense.jacobian).norm(),1e-8*std::max(1e-12,dense.jacobian.norm()));
    const n::Vector v=n::Vector::LinSpaced(op.Columns(),-.3,.8),u=n::Vector::Ones(op.Columns());
    const n::Vector w=n::Vector::LinSpaced(op.Rows(),-.4,.7);
    EXPECT_LE((op.ApplyAdjoint(w)-dense.jacobian.transpose()*w).norm(),1e-8*std::max(1e-12,(dense.jacobian.transpose()*w).norm()));
    const auto jv=op.Apply(v); const auto jtw=op.ApplyAdjoint(w);
    EXPECT_NEAR(jv.dot(w),v.dot(jtw),1e-12*std::max({1.,jv.norm()*w.norm(),v.norm()*jtw.norm()}));
    EXPECT_LE((op.Apply(v+2*u)-jv-2*op.Apply(u)).norm(),1e-12*std::max(1.,jv.norm()));
    const auto q_before=n::SparseWorkForTesting().q_actions;
    const n::Vector normal=op.ApplyNormal(v);
    EXPECT_EQ(n::SparseWorkForTesting().q_actions-q_before,2);
    const auto composed_before=n::SparseWorkForTesting().q_actions;
    const n::Vector composed=op.ApplyAdjoint(op.Apply(v));
    EXPECT_EQ(n::SparseWorkForTesting().q_actions-composed_before,6);
    const n::Vector expected=dense.jacobian.transpose()*dense.jacobian*v;
    EXPECT_LE((normal-expected).norm(),1e-12+1e-8*expected.norm());
    EXPECT_LE((normal-composed).norm(),1e-12+1e-8*composed.norm());
    EXPECT_NEAR(u.dot(normal),v.dot(op.ApplyNormal(u)),1e-12*std::max({1.,u.norm()*normal.norm(),v.norm()*op.ApplyNormal(u).norm()}));
    EXPECT_NEAR(v.dot(normal),jv.squaredNorm(),1e-12*std::max(1.,jv.squaredNorm()));
    EXPECT_LE((op.ApplyNormal(v+2*u)-normal-2*op.ApplyNormal(u)).norm(),1e-12*std::max(1.,normal.norm()));
    const auto gradient=op.ApplyAdjoint(e.residual/context.scale);
    for(Eigen::Index k=0;k<gradient.size();++k) EXPECT_NEAR(gradient(k),e.gradient(k),1e-13+2e-9*std::abs(e.gradient(k)));
}
}
TEST(JointProfileOperatorTest, FullResidualCorrectionAndFixedStateParity)
{
    Sample s; const auto e=n::EvaluateProfile(s.domain,s.y,s.eta,false,&s.context); ASSERT_TRUE(e.valid);
    EXPECT_GT(e.residual.norm(),1e-4); Parity(e,s.context);
    const auto d=p::DenseDifferentiate(e,s.context.scale,&s.context);
    EXPECT_GT((d.jacobian-d.projected).norm(),1e-8);
    const n::Vector direction=n::Vector::LinSpaced(3,-.3,.4); const double h=1e-5;
    const auto plus=n::EvaluateProfile(s.domain,s.y,s.eta+h*direction,false,&s.context);
    const auto minus=n::EvaluateProfile(s.domain,s.y,s.eta-h*direction,false,&s.context);
    ASSERT_TRUE(plus.valid && minus.valid);
    ASSERT_EQ(plus.certificate.active_atoms,e.certificate.active_atoms);
    ASSERT_EQ(minus.certificate.active_atoms,e.certificate.active_atoms);
    EXPECT_LT(((plus.residual-minus.residual)/(2*h*s.context.scale)-n::ProfileJacobianOperator(e,s.context).Apply(direction)).norm(),1e-8);
}
TEST(JointProfileOperatorTest, ActiveFaceAndOwnership)
{
    Sample s; n::LinearWorkspace workspace;
    auto e=n::EvaluateProfile(s.domain,s.y,s.eta,false,&s.context,nullptr,&workspace); ASSERT_TRUE(e.valid);
    const n::ProfileJacobianOperator op(e,s.context); ASSERT_TRUE(op.Valid());
    const n::Vector v=n::Vector::Ones(3),before=op.Apply(v);
    n::EvaluateProfile(s.domain,s.y,s.eta.array()+.02,false,&s.context,nullptr,&workspace);
    const auto normal_before=op.ApplyNormal(v);
    e={}; EXPECT_EQ(op.Apply(v),before); EXPECT_EQ(op.ApplyNormal(v),normal_before);
    n::Vector beta(6); beta<<0,-.2,2,.3,1,-.1;
    const auto boundary=n::EvaluateState(s.domain,s.y,s.eta,beta,s.context);
    const n::ProfileJacobianOperator changed(boundary,s.context); ASSERT_TRUE(changed.Valid());
    EXPECT_EQ(changed.FreeColumns(),5); EXPECT_NE(changed.Identity(),op.Identity());
    const auto dense=p::DenseDifferentiate(boundary,s.context.scale,&s.context);
    EXPECT_LT((changed.Apply(v)-dense.jacobian*v).norm(),1e-10);
    const auto counts=n::SparseWorkForTesting(); const auto ranks=n::OperatorWorkForTesting().rank_checks;
    op.Apply(v); op.ApplyAdjoint(n::Vector::Ones(op.Rows())); op.ApplyNormal(v);
    EXPECT_EQ(n::SparseWorkForTesting().numeric,counts.numeric);
    EXPECT_EQ(n::OperatorWorkForTesting().rank_checks,ranks);
    EXPECT_THROW(op.Apply(n::Vector::Zero(2)),std::invalid_argument);
}
TEST(JointProfileOperatorTest, CancellationAndRankThreshold)
{
    Sample s; auto e=n::EvaluateProfile(s.domain,s.y,s.eta,false,&s.context); ASSERT_TRUE(e.valid);
    e.derivative=e.x;
    const n::ProfileJacobianOperator op(e,s.context); ASSERT_TRUE(op.Valid());
    const auto dense=p::DenseDifferentiate(e,s.context.scale,&s.context);
    EXPECT_LT((op.Apply(n::Vector::Ones(3))-dense.jacobian*n::Vector::Ones(3)).norm(),1e-12);
    EXPECT_LT((op.ApplyNormal(n::Vector::Ones(3))-dense.jacobian.transpose()*dense.jacobian*n::Vector::Ones(3)).norm(),1e-12);
    EXPECT_FALSE(n::ProfileJacobianOperator(e,s.context,100).Valid());
    e.x.col(1)=e.x.col(0);
    EXPECT_FALSE(n::ProfileJacobianOperator(e,s.context).Valid());
}
TEST(JointProfileOperatorTest, FactorAdjointsAndPermutation)
{
    n::Matrix a=n::Matrix::Zero(9,3); a(0,2)=2; a(1,0)=3; a(2,1)=4; a(5,0)=.3; a(7,2)=.2;
    n::LinearWorkspace workspace; auto factor=workspace.Factor(a.sparseView(),{4,1,3},0);
    const n::Vector v=n::Vector::LinSpaced(3,-.7,.6),w=n::Vector::LinSpaced(9,-.2,.5);
    const n::Matrix inverse=a.completeOrthogonalDecomposition().pseudoInverse();
    const auto fixed=n::FreeDesignFactor::Fixed(a.sparseView(),{4,1,3});
    EXPECT_LT((fixed->LeastSquares(w)-inverse*w).norm(),1e-12);
    EXPECT_LT((fixed->PseudoInverseTranspose(v)-inverse.transpose()*v).norm(),1e-12);
    EXPECT_LT((fixed->NormalSolve(v)-inverse*inverse.transpose()*v).norm(),1e-12);
    EXPECT_LT((fixed->Compact().transpose()*fixed->Compact()-a.transpose()*a).norm(),1e-12);
    EXPECT_LT((factor->LeastSquares(w)-inverse*w).norm(),1e-12);
    EXPECT_LT((factor->PseudoInverseTranspose(v)-inverse.transpose()*v).norm(),1e-12);
    EXPECT_LT((factor->ProjectComplement(w)-(w-a*inverse*w)).norm(),1e-12);
    workspace.Factor(a.sparseView(),{4,1,3},0);
    EXPECT_THROW(factor->LeastSquares(w),std::logic_error);
}
TEST(JointProfileOperatorTest, NearThresholdDecisionsUseOriginalRows)
{
    for(double delta:{1e-8,1e-12,1e-15,0.})
    {
        n::Evaluation e; e.valid=true; e.eta=n::Vector::Zero(2); e.beta=n::Vector::Ones(4);
        n::Matrix x=n::Matrix::Zero(9,4); x(0,0)=1; x(1,1)=1; x(2,2)=1; x(0,3)=1; x(3,3)=delta;
        e.x=x.sparseView(); e.derivative=e.x; e.residual=n::Vector::Zero(9);
        auto context=n::CreateContext(e.residual,2); context.rank.rows=1000;
        const auto dense=p::DenseDifferentiate(e,context.scale,&context);
        const n::ProfileJacobianOperator op(e,context);
        EXPECT_EQ(op.Valid(),dense.valid)<<delta;
        if(!dense.valid) EXPECT_EQ(op.Reason(),dense.reason);
    }
}
TEST(JointProfileOperatorTest, InstrumentationIsNeutralAndActionsDoNotReduceMatrices)
{
    Sample s; const auto e=n::EvaluateProfile(s.domain,s.y,s.eta,false,&s.context);
    const n::ProfileJacobianOperator control(e,s.context);
    auto & work=n::ResourceWorkForTesting(); work={}; work.enabled=true;
    const n::ProfileJacobianOperator op(e,s.context); ASSERT_TRUE(op.Valid());
    const n::Vector v=n::Vector::Ones(3);
    EXPECT_EQ(control.Apply(v),op.Apply(v)); op.ApplyAdjoint(e.residual); EXPECT_EQ(control.ApplyNormal(v),op.ApplyNormal(v));
    bool compact=false;
    for(const auto & shape:work.dense_shapes)
    {
        if(shape.phase=="operator-rank" && shape.role=="free-design-compact") compact=true;
        if(shape.phase=="operator-apply" || shape.phase=="operator-adjoint" || shape.phase=="operator-normal") EXPECT_EQ(shape.columns,1);
        EXPECT_NE(shape.role,"derivative-t"); EXPECT_NE(shape.role,"derivative-coefficients"); EXPECT_NE(shape.role,"derivative-correction");
    }
    EXPECT_TRUE(compact); work={};
}
TEST(JointProfileOperatorTest, ObservableRowsPreserveParentContext)
{
    c::JointProblemInput input; input.atom_ids={"target","halo"}; input.support.resize(2);
    input.selection_domain.emplace(); input.selection_domain->target_indices={0};
    for(std::size_t r=0;r<30;++r)
    {
        const double square=std::pow(.03+.07*static_cast<double>(r),2);
        input.row_ids.push_back(std::to_string(r)); input.support[0].push_back({r,square});
        const auto b=n::EvaluateKernel(square,.6,2.5);
        input.observations.push_back(2*b.gaussian-.1*b.charge+1e-5*std::sin(static_cast<double>(r)));
    }
    input.support[1]={{12,6.25}}; input.observations[12]+=.7;
    const c::JointProblem problem(std::move(input)); const auto & data=c::JointProblemAccess::Get(problem);
    const auto & layout=problem.ParameterLayout(); const auto domain=n::ProfileDomain(data.domain,layout);
    const auto context=n::ProfileContext(data.context,layout,30);
    const n::Indices rows(layout.informative_rows.begin(),layout.informative_rows.end());
    const auto y=n::SelectValues(data.y,rows);
    EXPECT_EQ(context.rank.rows,30); EXPECT_EQ(context.scale,data.context.scale); EXPECT_EQ(y.size(),29);
    const auto e=n::EvaluateProfile(domain,y,n::Vector::Constant(1,std::log(.7)),false,&context);
    ASSERT_TRUE(e.valid); Parity(e,context);
}
TEST(JointPreconditionerTest, SymmetricWeightsAndPositiveAction)
{
    const n::SolverBlockMapping map(n::PreconditionerSpace::Width,3,{{0,1},{1,2}});
    n::Matrix h(2,2); h<<2,.3,.3,1;
    const auto apply=[&](const n::Vector & v) {
        n::Vector out=n::Vector::Zero(3);
        for(const auto & block:map.blocks) block.Scatter(h.llt().solve(block.Restrict(v)),out);
        return out;
    };
    const n::Vector x=n::Vector::LinSpaced(3,-.4,.7),y=n::Vector::LinSpaced(3,.8,-.2);
    EXPECT_NEAR(x.dot(apply(y)),y.dot(apply(x)),1e-12); EXPECT_GT(x.dot(apply(x)),0);
    n::Vector unity=n::Vector::Zero(3);
    for(const auto & b:map.blocks) b.Scatter(b.Restrict(x),unity);
    EXPECT_LT((unity-x).norm(),1e-12);
    EXPECT_THROW((n::SolverBlockMapping(n::PreconditionerSpace::Width,3,{{0,1}})),std::invalid_argument);
    EXPECT_THROW((n::SolverBlockMapping(n::PreconditionerSpace::Width,3,{{0,1,1,2}})),std::invalid_argument);
}
TEST(JointPreconditionerTest, IdentityAndStaleContext)
{
    n::PreconditionerContext context{std::make_shared<const n::LinearizationIdentity>(),n::PreconditionerSpace::Width,n::Vector::Ones(3),.1};
    const n::IdentityPreconditioner identity(context); const n::Vector rhs=n::Vector::LinSpaced(3,-1,1);
    EXPECT_EQ(identity.ApplyInverse(rhs,context),rhs);
    auto stale=context; stale.damping=.2; EXPECT_THROW(identity.ApplyInverse(rhs,stale),std::logic_error);
    stale=context; stale.metric(0)=2; EXPECT_THROW(identity.ApplyInverse(rhs,stale),std::logic_error);
    stale=context; stale.space=n::PreconditionerSpace::FreeAC; EXPECT_THROW(identity.ApplyInverse(rhs,stale),std::logic_error);
    stale=context; stale.linearization=std::make_shared<const n::LinearizationIdentity>();
    EXPECT_THROW(identity.ApplyInverse(rhs,stale),std::logic_error);
}
TEST(JointPreconditionerTest, StructuralSnapshotAndSolverMappingsAreSeparate)
{
    auto input=std::make_shared<const c::JointProblemInput>(second_stage_test::OperatorWorkload("chain",3));
    const auto layout=n::BuildParameterLayout(*input);
    const n::PreconditionerPartition partition(input,layout,{{"left",{0,1},{2},{}},{"right",{2},{1},{}}});
    EXPECT_EQ(partition.problem,input); EXPECT_EQ(partition.atom_blocks[1],(std::vector<std::size_t>{0,1}));
    EXPECT_THROW((n::PreconditionerPartition(input,layout,{{"missing",{0},{1},{}}})),std::invalid_argument);
    const n::SolverBlockMapping widths(n::PreconditionerSpace::Width,3,{{0,1,2},{1,2}});
    const n::SolverBlockMapping free_ac(n::PreconditionerSpace::FreeAC,5,{{0,1,2,3,4},{2,3,4}});
    EXPECT_NE(widths.dimension,free_ac.dimension);
    EXPECT_THROW(widths.blocks[0].Restrict(n::Vector::Ones(5)),std::invalid_argument);
}
TEST(JointProfileOperatorTest, FrozenWorkloadHasCompleteDeterministicSupport)
{
    for(const auto * topology:{"chain","cube"})
    {
        const auto a=second_stage_test::OperatorWorkload(topology,9),b=second_stage_test::OperatorWorkload(topology,9);
        EXPECT_EQ(a.observations,b.observations); EXPECT_EQ(a.row_ids,b.row_ids);
        EXPECT_EQ(second_stage_test::OperatorWorkloadHash(a),second_stage_test::OperatorWorkloadHash(b));
        auto changed=b; changed.observations[0]+=.01;
        EXPECT_NE(second_stage_test::OperatorWorkloadHash(a),second_stage_test::OperatorWorkloadHash(changed));
        for(const auto & support:a.support)
        {
            ASSERT_EQ(support.size(),515);
            for(const auto & member:support) {EXPECT_LT(member.row,a.observations.size()); EXPECT_LE(member.squared_distance,6.25);}
        }
        const c::JointProblem problem(a); EXPECT_EQ(c::JointProblemAccess::Get(problem).partition.components.size(),1);
    }
}

TEST(JointProfileOperatorTest, NormalActionRetainsWeakDirectionsAndRejectsInvalidInputs)
{
    n::Evaluation e; e.valid=true; e.eta=n::Vector::Zero(2); e.beta=n::Vector::Ones(4);
    n::Matrix x=n::Matrix::Zero(9,4); x(0,0)=1; x(1,1)=1; x(2,2)=1; x(0,3)=1; x(3,3)=1e-8;
    e.x=x.sparseView(); x.row(8)<<.1,.2,.3,.4; e.derivative=x.sparseView();
    e.residual=n::Vector::Zero(9); e.residual(8)=.1;
    const auto context=n::CreateContext(e.residual,2);
    const n::ProfileJacobianOperator op(e,context); ASSERT_TRUE(op.Valid());
    const auto dense=p::DenseDifferentiate(e,context.scale,&context); ASSERT_TRUE(dense.valid);
    for(int k=0;k<2;++k)
    {
        const n::Vector v=n::Vector::Unit(2,k),expected=dense.jacobian.transpose()*dense.jacobian*v;
        EXPECT_LE((op.ApplyNormal(v)-expected).norm(),1e-12+1e-8*expected.norm());
    }
    EXPECT_THROW(op.ApplyNormal(n::Vector::Zero(3)),std::invalid_argument);
    EXPECT_THROW(op.ApplyNormal(n::Vector::Constant(2,n::unavailable)),std::invalid_argument);
    const auto work=n::SparseWorkForTesting(); const auto ranks=n::OperatorWorkForTesting().rank_checks;
    op.ApplyNormal(n::Vector::Ones(2));
    EXPECT_EQ(n::SparseWorkForTesting().numeric,work.numeric);
    EXPECT_EQ(n::SparseWorkForTesting().fixed_factorizations,work.fixed_factorizations);
    EXPECT_EQ(n::SparseWorkForTesting().compact_extractions,work.compact_extractions);
    EXPECT_EQ(n::OperatorWorkForTesting().rank_checks,ranks);
}
