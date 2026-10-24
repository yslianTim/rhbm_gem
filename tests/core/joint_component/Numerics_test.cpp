#include <gtest/gtest.h>
#include "core/detail/joint_component/TiledDerivative.hpp"
#include "core/detail/joint_component/ResourceWork.hpp"
#include "core/detail/joint_component/CompactSvd.hpp"
#include "core/detail/joint_component/Problem.hpp"
#include "support/JointSyntheticWorkload.hpp"
#include "support/JointTestNumerics.hpp"
#include <unsupported/Eigen/NonLinearOptimization>
#include <array>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <tuple>

namespace {
namespace m=second_stage_test::matched;
namespace p=m::joint_abc;
namespace n=p::runtime;
using Vector=Eigen::VectorXd;
using Matrix=Eigen::MatrixXd;
struct Sample
{
    m::unique_grid::Grid grid;
    std::vector<m::Atom> atoms{{{-.6,0,0},2,.42,-.3},{{.7,.2,0},1.5,.67,.4}};
    Vector y;
    Sample()
    {
        for (int x=-6;x<=6;++x) for (int z=-3;z<=3;++z) for (int y0=-4;y0<=4;++y0)
        {
            m::unique_grid::Voxel voxel; voxel.position={x*.4,y0*.4,z*.4}; voxel.index=grid.voxels.size();
            grid.voxels.push_back(voxel);
        }
        y.resize(static_cast<Eigen::Index>(grid.voxels.size()));
        for (Eigen::Index k=0;k<y.size();++k) y(k)=m::unique_grid::Direct(grid.voxels[static_cast<std::size_t>(k)].position,atoms,2.5);
    }
};
// Independent small direct joint solve, with numerical derivatives of A/C/log-B.
// This comparison uses an interior solution; boundary behavior is tested separately.
struct Direct
{
    const Sample & sample;
    int values() const {return static_cast<int>(sample.y.size());}
    int operator()(const Vector & x,Vector & residual) const
    {
        auto atoms=sample.atoms;
        for (std::size_t a=0;a<atoms.size();++a)
        {atoms[a].amplitude=x(static_cast<Eigen::Index>(3*a)); atoms[a].charge=x(static_cast<Eigen::Index>(3*a+1)); atoms[a].width=std::exp(x(static_cast<Eigen::Index>(3*a+2)));}
        residual.resize(sample.y.size());
        for (Eigen::Index k=0;k<residual.size();++k)
            residual(k)=m::unique_grid::Direct(sample.grid.voxels[static_cast<std::size_t>(k)].position,atoms,2.5)-sample.y(k);
        return 0;
    }
    int df(const Vector & x,Matrix & jacobian) const
    {
        jacobian.resize(sample.y.size(),x.size());
        for (Eigen::Index k=0;k<x.size();++k)
        {
            Vector plus=x,minus=x,rp,rm; plus(k)+=1e-5; minus(k)-=1e-5;
            (*this)(plus,rp); (*this)(minus,rm); jacobian.col(k)=(rp-rm)/2e-5;
        }
        return 0;
    }
};

;
;
}

TEST(JointTestNumericsTest, LogWidthBasisIncludesCenterNearZeroAndFixedCutoff)
{
    for (double r:{0.,.5e-5,1e-5,1.1e-5,.4,2.5,2.50000001})
    {
        const double b=.53,h=1e-5;
        const auto value=m::EvaluateBasis(r*r,b,2.5),plus=m::EvaluateBasis(r*r,b*std::exp(h),2.5),minus=m::EvaluateBasis(r*r,b*std::exp(-h),2.5);
        EXPECT_NEAR(value.gaussian_log_width,(plus.gaussian-minus.gaussian)/(2*h),1e-8);
        EXPECT_NEAR(value.charge_log_width,(plus.charge-minus.charge)/(2*h),1e-8);
        if (r>2.5) EXPECT_EQ(value.gaussian+value.charge+value.gaussian_log_width+value.charge_log_width,0);
    }
}

TEST(JointTestNumericsTest, FullJacobianAndEnvelopeGradientAtNonzeroResidual)
{
    Sample sample; const p::Domain domain(sample.grid,sample.atoms);
    for (Eigen::Index k=0;k<sample.y.size();++k) sample.y(k)+=.03*std::sin(static_cast<double>(k));
    const Vector eta=Eigen::Vector2d(.55,.51).array().log();
    const auto e=p::Evaluate(domain,sample.y,eta); ASSERT_TRUE(e.valid);
    const double scale=std::max(1.0,sample.y.norm()); const auto d=p::Differentiate(e,scale);
    ASSERT_TRUE(d.valid); EXPECT_GT((d.jacobian-d.projected).norm(),1e-3);
    for (Eigen::Index k=0;k<2;++k)
    {
        Vector plus=eta,minus=eta; plus(k)+=1e-5; minus(k)-=1e-5;
        const auto ep=p::Evaluate(domain,sample.y,plus),em=p::Evaluate(domain,sample.y,minus);
        ASSERT_TRUE(ep.valid && em.valid);
        const Vector finite=(ep.residual-em.residual)/(2e-5*scale);
        EXPECT_LT((finite-d.jacobian.col(k)).norm()/finite.norm(),1e-7);
        EXPECT_NEAR(e.gradient(k),(ep.residual.squaredNorm()-em.residual.squaredNorm())/(4e-5*scale*scale),1e-9);
    }
}

TEST(JointTestNumericsTest, BoundaryAndChangingActiveFaceRebuildDerivative)
{
    Sample sample; sample.atoms.resize(1); sample.atoms[0].amplitude=.05; sample.atoms[0].charge=-1;
    for (Eigen::Index k=0;k<sample.y.size();++k) sample.y(k)=m::unique_grid::Direct(sample.grid.voxels[static_cast<std::size_t>(k)].position,sample.atoms,2.5);
    const p::Domain domain(sample.grid,sample.atoms); bool active=false,interior=false;
    for (double width:{.25,.42,.8})
    {
        const Vector eta=Vector::Constant(1,std::log(width)); const auto e=p::Evaluate(domain,sample.y,eta);
        ASSERT_TRUE(e.valid); active |= e.beta(0)==0; interior |= e.beta(0)>0;
        const double scale=std::max(1.0,sample.y.norm()); const auto d=p::Differentiate(e,scale); ASSERT_TRUE(d.valid);
        const auto plus=p::Evaluate(domain,sample.y,eta.array()+1e-5),minus=p::Evaluate(domain,sample.y,eta.array()-1e-5);
        ASSERT_EQ(plus.certificate.at("active_atoms"),minus.certificate.at("active_atoms"));
        EXPECT_LT(((plus.residual-minus.residual)/(2e-5*scale)-d.jacobian.col(0)).norm(),1e-8);
    }
    EXPECT_TRUE(active && interior);
}

TEST(JointTestNumericsTest, InvalidBStillFailsAndUninformativeWidthCanConverge)
{
    Sample sample; auto duplicate=sample.atoms; duplicate[1].position=duplicate[0].position;
    const p::Domain bad(sample.grid,duplicate),good(sample.grid,sample.atoms);
    const Vector eta=Vector::Constant(2,std::log(.5));
    const auto invalid=p::Evaluate(bad,sample.y,eta); EXPECT_FALSE(invalid.valid);
    EXPECT_FALSE(p::Differentiate(invalid,1).valid);
    EXPECT_NE(p::Fit(bad,sample.y,Vector::Constant(2,.5)).at("runtime_convergence"),"passed");
    EXPECT_FALSE(p::Evaluate(good,sample.y,Vector::Constant(2,1000)).valid);
    EXPECT_FALSE(p::Evaluate(good,sample.y,Vector::Constant(2,-1000)).valid);
    const auto zero=p::Fit(good,Vector::Zero(sample.y.size()),Vector::Constant(2,.5));
    EXPECT_EQ(zero.at("runtime_convergence"),"passed");
    EXPECT_TRUE(zero.at("primary").at("kkt_passed").as_bool());
    EXPECT_EQ(zero.at("primary").at("rss"),0.0);
    EXPECT_TRUE(zero.at("numerical_convergence_passed").as_bool());
    // Geometry and observations are immutable across all evaluations.
    EXPECT_EQ(good.rows,sample.y.size()); EXPECT_EQ(good.atoms[0].size(),p::Domain(sample.grid,sample.atoms).atoms[0].size());
}

TEST(JointTestNumericsTest, TiledDerivativeMatchesDenseAcrossTileBoundaries)
{
    namespace n=p::runtime;
    Sample sample; const p::Domain base(sample.grid,sample.atoms);
    for(Eigen::Index rows:{8191,8192,8193})
    {
        std::vector<std::vector<p::Support>> support(base.atoms.size()); Vector y(rows);
        for(Eigen::Index first=0;first<rows;first+=base.rows)
            for(std::size_t a=0;a<support.size();++a) for(const auto & s:base.atoms[a])
                if(first+s.row<rows) support[a].push_back({first+s.row,s.square});
        for(Eigen::Index r=0;r<rows;++r) y(r)=sample.y(r%base.rows)+.03*std::sin(static_cast<double>(r));
        const p::Domain domain(rows,std::move(support)); const auto context=p::MakeContext(y,2);
        const auto e=p::Evaluate(domain,y,Eigen::Vector2d(.55,.51).array().log(),&context);
        ASSERT_TRUE(e.valid); const auto dense=p::DenseDifferentiate(e,context.scale,&context);
        ASSERT_TRUE(dense.valid); ASSERT_GT((dense.projected-dense.jacobian).norm(),1e-3);
        const auto correction=p::DenseLocalCorrection(e,dense,context);
        for(Eigen::Index tile:{17,8192})
        {
            n::DerivativeWorkForTesting()={};
            const auto d=n::PrepareDerivative(e,context.scale,&context,-1,tile);
            const auto reduced=n::ReduceDerivative(d,e.residual,tile);
            const auto width_derivative=n::PrepareProjectedWidthDerivative(e,context.scale,&context,-1,tile);
            const auto reduced_width=n::ReduceProjectedWidthDerivative(width_derivative,tile);
            ASSERT_TRUE(reduced.valid);
            ASSERT_TRUE(reduced_width.valid);
            EXPECT_LE(n::DerivativeWorkForTesting().maximum_generated_rows,tile);
            EXPECT_LE(n::DerivativeWorkForTesting().maximum_reduction_rows,tile+e.x.cols());
            const auto actual=p::MaterializeDerivative(e,context.scale,&context,-1,tile);
            EXPECT_LT((actual.projected-dense.projected).norm()/dense.projected.norm(),1e-8);
            EXPECT_LT((actual.jacobian-dense.jacobian).norm()/dense.jacobian.norm(),1e-8);
            for(const auto & matrices:{std::pair{dense.projected,reduced_width.projected},
                std::pair{dense.jacobian,reduced.jacobian}})
            {
                const auto left=n::ComputeSpectrum(matrices.first,context.rank,2,false);
                const auto right=n::ComputeSpectrum(matrices.second,context.rank,2,false);
                EXPECT_EQ(left.rank,right.rank);
                EXPECT_LE((left.singular_values-right.singular_values).lpNorm<Eigen::Infinity>()/left.singular_values(0),1e-10);
                EXPECT_NEAR(left.threshold,right.threshold,1e-10*left.threshold);
            }
            Eigen::JacobiSVD<Eigen::MatrixXd> svd(reduced.jacobian,Eigen::ComputeThinU|Eigen::ComputeThinV);
            svd.setThreshold(context.rank.Relative(2));
            const Vector compact=svd.solve(-reduced.response);
            EXPECT_LE(((correction-compact).array().abs()/(1+correction.array().abs().max(compact.array().abs()))).maxCoeff(),1e-10);
            EXPECT_LT((reduced_width.projected_norms-dense.projected.colwise().norm().transpose()).norm()/dense.projected.norm(),1e-10);
            EXPECT_LT((reduced.jacobian_norms-dense.jacobian.colwise().blueNorm().transpose()).norm()/dense.jacobian.norm(),1e-10);
        }
    }
}

TEST(JointComponentNumericsTest, SparseConstrainedSolvePreservesInteriorBoundaryAndReleaseFaces)
{
    namespace n=p::runtime;
    Matrix x=Matrix::Zero(7,4);
    x(0,0)=1; x(1,1)=1; x(2,2)=1; x(3,3)=1;
    const n::Sparse sparse=x.sparseView(); const Vector weights=Vector::Ones(x.rows());
    for(const Eigen::Vector4d truth:{Eigen::Vector4d(2,-.3,1.5,.4),Eigen::Vector4d(-.5,-.3,1.5,.4)})
    {
        const Vector y=x*truth; const auto solved=n::SolveLinear(sparse,y,weights);
        ASSERT_TRUE(solved.valid);
        const auto certificate=n::CertifyLinear(sparse,y,solved.beta,std::max(1.0,y.norm()));
        ASSERT_TRUE(certificate.available); EXPECT_TRUE(certificate.feasible); EXPECT_TRUE(certificate.kkt_passed);
        EXPECT_EQ(solved.releases,0);
        if(truth(0)>0) {EXPECT_TRUE(certificate.active_atoms.empty()); EXPECT_NEAR(solved.beta(0),truth(0),1e-12);}
        else {ASSERT_EQ(certificate.active_atoms,std::vector<Eigen::Index>{0}); EXPECT_EQ(solved.beta(0),0);}
    }

    Matrix release_design(7,4);
    release_design<<0,0,-1,2, 0,1,-1,-1, -2,0,-2,-2, 2,-2,-1,0,
        0,0,2,2, 0,2,-1,0, -1,1,2,-1;
    const Vector release_y=(Vector(7)<<0,1,2,1,1,2,0).finished();
    const auto release_sparse=n::Sparse(release_design.sparseView());
    const auto released=n::SolveLinear(release_sparse,release_y,Vector::Ones(release_y.size()));
    ASSERT_TRUE(released.valid); EXPECT_EQ(released.releases,1); EXPECT_GT(released.beta(0),0);
    EXPECT_EQ(released.beta(2),0);
    const auto released_certificate=n::CertifyLinear(release_sparse,release_y,released.beta,
        std::max(1.0,release_y.norm()));
    ASSERT_TRUE(released_certificate.available); EXPECT_TRUE(released_certificate.feasible);
    EXPECT_TRUE(released_certificate.kkt_passed);

    Matrix deficient=Matrix::Zero(7,4);
    deficient(0,0)=1; deficient(1,1)=1; deficient(0,2)=1; deficient(2,3)=1;
    const auto rejected=n::SolveLinear(n::Sparse(deficient.sparseView()),Vector::Ones(7),Vector::Ones(7));
    EXPECT_FALSE(rejected.valid); EXPECT_EQ(rejected.reason,"rank-deficient"); EXPECT_EQ(rejected.rank,3);
}

TEST(JointComponentNumericsTest, LocalLinearCertificateKeepsItsStrictKktLimit)
{
    n::Sparse design(2,2);
    design.insert(0,0)=1.0; design.insert(1,1)=1.0; design.makeCompressed();
    Vector beta=Vector::Ones(2),observations(2);
    observations<<1.0-2e-7,1.0;
    const auto certificate=n::CertifyLinear(design,observations,beta,observations.norm());
    ASSERT_TRUE(certificate.available);
    EXPECT_GT(certificate.projected_kkt,n::PracticalJointConvergenceCriteria.global_ac_kkt);
    EXPECT_FALSE(certificate.kkt_passed);
}

TEST(JointComponentNumericsTest, CompactSvdMatchesJacobiReferenceAtRankBoundaryAndAbsoluteOverride)
{
    namespace n=p::runtime;
    constexpr double relative=1e-9;
    for(double absolute:{-1.,1e-7}) for(double factor:{0.,.5,.999999,1.,1.000001,2.})
    {
        const double threshold=absolute<0 ? relative : absolute;
        Matrix x=Matrix::Identity(16,16); x(15,15)=threshold*factor;
        const Vector rhs=x*Vector::LinSpaced(16,.5,2.);
        Eigen::JacobiSVD<Matrix,Eigen::ComputeFullU|Eigen::ComputeFullV> reference(x);
        const double maximum=reference.singularValues()(0);
        reference.setThreshold(absolute>=0 && maximum>0 ? absolute/maximum : relative);
        const double expected_threshold=reference.threshold()*maximum;
        const Eigen::Index expected_rank=reference.rank();
        const Vector expected_solution=reference.solve(rhs);
        const auto actual=n::CompactSvd(x,relative,absolute,&rhs);
        ASSERT_TRUE(actual.valid); EXPECT_EQ(actual.rank,expected_rank);
        EXPECT_EQ(actual.threshold,expected_threshold);
        EXPECT_LT((actual.solution-expected_solution).norm(),1e-10);
        EXPECT_TRUE(actual.used_bdc);
        ASSERT_EQ(actual.solution.size(),rhs.size()); EXPECT_TRUE(actual.solution.allFinite());
        const bool near=std::abs(factor-1)*threshold<=64*std::numeric_limits<double>::epsilon()*16;
        EXPECT_EQ(actual.jacobi_retry,near);
    }
}

TEST(JointComponentNumericsTest, SharedRankContractPreservesBothBoundaryRules)
{
    namespace n=p::runtime;
    const double relative=n::RankPolicy{1000,0,0}.Relative(4);
    for(double absolute:{-1.,1e-7}) for(double multiple:{0.,.5,1.,2.})
    {
        const double threshold=absolute<0 ? relative : absolute;
        Matrix x=Matrix::Identity(4,4); x(3,3)=threshold*multiple;
        const auto old=n::CompactSvd(x,relative,absolute);
        const auto native=n::EvaluateRank(x,{{1000,8,4},4,absolute});
        const auto strict=n::EvaluateRank(x,{{1000,8,4},4,absolute,n::RankBoundary::StrictGreater});
        ASSERT_TRUE(native.valid && strict.valid);
        EXPECT_EQ(native.rank,old.rank); EXPECT_EQ(native.threshold,old.threshold);
        EXPECT_EQ(strict.rank,(old.singular_values.array()>old.threshold).count());
    }
    const auto global=n::EvaluateRank(Matrix::Identity(2,2),{{1000,4,2},4});
    ASSERT_TRUE(global.valid); EXPECT_EQ(global.threshold,n::RankPolicy({1000,4,2}).Relative(4));
    EXPECT_EQ(n::EvaluateRank(Matrix::Zero(4,4),{{1000,8,4},4}).rank,0);
    EXPECT_FALSE(n::EvaluateRank(Matrix::Constant(4,4,n::unavailable),{{1000,8,4},4}).valid);
}

TEST(JointComponentNumericsTest, MemoryPlannerRetainsNonRankStagesAndLocalBudgets)
{
    using rhbm_gem::core::JointProblemInput;
    std::size_t bytes{};
    EXPECT_TRUE(n::CheckedMultiplyBytes(8,16,bytes)); EXPECT_EQ(bytes,128u);
    EXPECT_FALSE(n::CheckedMultiplyBytes(std::numeric_limits<std::size_t>::max(),2,bytes));
    EXPECT_FALSE(n::CheckedAddBytes(std::numeric_limits<std::size_t>::max(),1,bytes));

    JointProblemInput input; input.observations={1,2}; input.atom_ids={"a","b"}; input.row_ids={"r0","r1"};
    input.support={{{0,0},{1,1}},{{1,2}}};
    ASSERT_EQ(n::CountSupportEntries(input),3u);
    ASSERT_TRUE(n::EstimateSnapshotBytes(input));
    ASSERT_TRUE(n::EstimateConstructionScratchBytes(input));

    n::MemoryShape large; large.component_atoms=37'406; large.component_observations=5'938'087;
    large.component_support_entries=2'000'000; large.component_blocks=3'118;
    large.affected_row_entries=3'000'000; large.component_free_columns=74'812;
    large.maximum_block_rows=8'000; large.maximum_block_support_entries=80'000;
    large.maximum_neighbor_atoms=12; large.component_width_columns=12;
    large.maximum_tile_rows=8192; large.resident_bytes=256'000'000;
    const auto construction=n::PlanMemory(large,n::MemoryStage::ProblemConstruction);
    const auto structural=n::PlanMemory(large,n::MemoryStage::StructuralConstruction);
    const auto preparation=n::PlanMemory(large,n::MemoryStage::FixedNeighborPreparation);
    const auto state=n::PlanMemory(large,n::MemoryStage::GlobalStateEvaluation);
    const auto endpoint=n::PlanMemory(large,n::MemoryStage::EndpointVerification);
    auto local=large; local.component_free_columns=24;
    const auto conditioning=n::PlanMemory(local,n::MemoryStage::LocalConditioning);
    local.component_stability_evaluations=2;
    const auto stability=n::PlanMemory(local,n::MemoryStage::ParameterStability);
    EXPECT_EQ(construction.status,n::MemoryStatus::Pass);
    EXPECT_EQ(structural.status,n::MemoryStatus::Pass);
    EXPECT_EQ(preparation.status,n::MemoryStatus::Pass);
    EXPECT_EQ(state.status,n::MemoryStatus::Pass);
    EXPECT_EQ(endpoint.status,n::MemoryStatus::Pass);
    EXPECT_EQ(conditioning.status,n::MemoryStatus::Pass);
    EXPECT_EQ(stability.status,n::MemoryStatus::Pass);
    EXPECT_GT(stability.estimated_peak_bytes,conditioning.estimated_peak_bytes);
    EXPECT_EQ(conditioning.factor_rows,24);
    EXPECT_EQ(conditioning.factor_columns,24);
    EXPECT_EQ(conditioning.budget_bytes,n::JointMemoryPlanningBudgetBytes);

    auto overflow=large; overflow.component_free_columns=std::numeric_limits<std::size_t>::max();
    EXPECT_EQ(n::PlanMemory(overflow,n::MemoryStage::LocalConditioning).status,n::MemoryStatus::Overflow);
    auto tile_overflow=large; tile_overflow.maximum_tile_rows=std::numeric_limits<std::size_t>::max();
    EXPECT_EQ(n::PlanMemory(tile_overflow,n::MemoryStage::LocalConditioning).status,n::MemoryStatus::Overflow);
    auto construction_overflow=large;
    construction_overflow.component_observations=std::numeric_limits<std::size_t>::max();
    EXPECT_EQ(n::PlanMemory(construction_overflow,n::MemoryStage::ProblemConstruction).status,n::MemoryStatus::Overflow);
}

TEST(JointComponentNumericsTest, StreamingStateMatchesSparseGlobalEvaluation)
{
    using rhbm_gem::core::JointProblem;
    using rhbm_gem::core::JointProblemAccess;
    JointProblem problem(second_stage_test::SyntheticJointWorkload("chain",8));
    const auto & data=JointProblemAccess::Get(problem);
    const Vector eta=Vector::Constant(8,std::log(.55));
    Vector beta(16);
    for(Eigen::Index atom=0;atom<8;++atom)
    {
        beta(2*atom)=atom%5==0 ? 0.0 : 1.7+.1*static_cast<double>(atom);
        beta(2*atom+1)=atom%2 ? -.4 : .3;
    }
    const auto sparse=n::EvaluateState(data.domain,data.y,eta,beta,data.context);
    const auto streaming=n::EvaluateStateStreaming(data.domain,data.y,eta,beta,data.context);
    ASSERT_TRUE(sparse.valid)<<sparse.reason;
    ASSERT_TRUE(streaming.valid)<<streaming.reason;
    ASSERT_TRUE(sparse.certificate.available && streaming.certificate.available);
    const Vector sparse_prediction=sparse.residual+data.y;
    const Vector streaming_prediction=streaming.residual+data.y;
    const double prediction_scale=std::max(1.0,sparse_prediction.cwiseAbs().maxCoeff());
    EXPECT_LE((streaming_prediction-sparse_prediction).lpNorm<Eigen::Infinity>(),2e-13*prediction_scale);
    EXPECT_LE((streaming.residual-sparse.residual).lpNorm<Eigen::Infinity>(),2e-13*prediction_scale);
    EXPECT_NEAR(streaming.certificate.objective,sparse.certificate.objective,2e-13*std::max(1.0,sparse.certificate.objective));
    EXPECT_NEAR(streaming.certificate.projected_kkt,sparse.certificate.projected_kkt,2e-12);
    EXPECT_EQ(streaming.certificate.active_atoms,sparse.certificate.active_atoms);
    EXPECT_EQ(streaming.certificate.feasible,sparse.certificate.feasible);
    EXPECT_LE((streaming.gradient-sparse.gradient).lpNorm<Eigen::Infinity>(),
        2e-12*std::max(1.0,sparse.gradient.lpNorm<Eigen::Infinity>()));
    EXPECT_EQ(streaming.x.size(),0);
    EXPECT_EQ(streaming.derivative.size(),0);
}

TEST(JointComponentNumericsTest, IndependentStreamingReplayMatchesLegacyMatrixReplay)
{
    for(const std::string topology:{"chain","cube"})
    {
        SCOPED_TRACE(topology);
        using rhbm_gem::core::JointProblem;
        using rhbm_gem::core::JointProblemAccess;
        JointProblem problem(second_stage_test::SyntheticJointWorkload(topology,8));
        const auto & data=JointProblemAccess::Get(problem);
        const Vector eta=Vector::Constant(8,std::log(.55));
        const auto primary=n::EvaluateProfile(data.domain,data.y,eta,&data.context);
        ASSERT_TRUE(primary.valid);
        const auto streaming=n::EvaluateStateStreaming(data.domain,data.y,primary.eta,primary.beta,data.context);
        ASSERT_TRUE(streaming.valid);
        const auto legacy=n::CheckReplay(data.domain,data.y,primary,data.context);
        const auto independent=n::CheckIndependentReplay(data.domain,data.y,streaming,data.context);
        ASSERT_TRUE(legacy.passed)<<legacy.reason;
        ASSERT_TRUE(independent.available && independent.valid)<<independent.reason;
        EXPECT_TRUE(independent.primary_kkt_passed);
        EXPECT_TRUE(independent.independent_kkt_passed);
        EXPECT_TRUE(independent.prediction_passed);
        EXPECT_TRUE(independent.residual_passed);
        EXPECT_TRUE(independent.objective_passed);
        EXPECT_TRUE(independent.kkt_agrees);
        EXPECT_TRUE(independent.gradient_agrees);
        EXPECT_TRUE(independent.passed)<<independent.reason;
        EXPECT_LE(independent.prediction_difference,2e-12+2e-13*(primary.x*primary.beta).cwiseAbs().maxCoeff());
        EXPECT_LE(independent.residual_difference,2e-12+2e-13*(primary.x*primary.beta).cwiseAbs().maxCoeff());
        EXPECT_LE(independent.objective_difference,1e-12+2e-12*primary.certificate.objective);
        EXPECT_LE(independent.kkt_difference,1e-13);
        EXPECT_LE(independent.gradient_difference,1e-13+2e-9*streaming.gradient.lpNorm<Eigen::Infinity>());
        EXPECT_TRUE(std::isfinite(independent.cancellation_ratio));
    }
}

TEST(JointComponentNumericsTest, PracticalStationarityCannotHidePerCoordinateReplayDisagreement)
{
    constexpr Eigen::Index rows=4;
    std::vector<std::vector<n::Support>> support(1);
    n::Vector y(rows),eta= n::Vector::Constant(1,std::log(.55)),beta(2);
    beta<<1.5,-.1;
    for(Eigen::Index row=0;row<rows;++row)
    {
        const double square=std::array<double,4>{0.0,.36,1.44,4.0}[static_cast<std::size_t>(row)];
        support[0].push_back({row,square});
        const auto basis=n::EvaluateKernel(square,.55,2.5);
        y(row)=beta(0)*basis.gaussian+beta(1)*basis.charge;
    }
    const n::Domain domain(rows,std::move(support));
    auto context=p::MakeContext(y,1);
    context.atom_ids=n::Identities({"target-42"});
    context.row_ids=n::Identities({"row-1","row-2","row-3","row-4"});
    auto primary=n::EvaluateStateStreaming(domain,y,eta,beta,context);
    ASSERT_TRUE(primary.valid)<<primary.reason;
    auto original=n::CheckIndependentReplay(domain,y,primary,context,nullptr,nullptr,
        n::PracticalJointConvergenceCriteria);
    ASSERT_TRUE(original.passed)<<original.reason;

    primary.gradient(0)+=5e-10;
    const auto replay=n::CheckIndependentReplay(domain,y,primary,context,nullptr,nullptr,
        n::PracticalJointConvergenceCriteria);
    ASSERT_TRUE(replay.valid)<<replay.reason;
    EXPECT_TRUE(replay.primary_kkt_passed);
    EXPECT_TRUE(replay.independent_kkt_passed);
    EXPECT_TRUE(replay.prediction_passed);
    EXPECT_TRUE(replay.residual_passed);
    EXPECT_TRUE(replay.objective_passed);
    EXPECT_TRUE(replay.kkt_agrees);
    EXPECT_FALSE(replay.gradient_agrees);
    EXPECT_FALSE(replay.passed);
    EXPECT_EQ(replay.reason,"per-coordinate-replay-enclosure");
    EXPECT_GT(replay.primary_width_gradient,n::StrictJointConvergenceCriteria.global_width_gradient);
    EXPECT_LE(replay.primary_width_gradient,n::PracticalJointConvergenceCriteria.global_width_gradient);
    EXPECT_EQ(replay.maximum_gradient_difference_atom_id,"target-42");
    EXPECT_EQ(replay.gradient_enclosure_failure_atom_id,"target-42");
    EXPECT_DOUBLE_EQ(replay.gradient_enclosure_failure_primary,primary.gradient(0));
    EXPECT_TRUE(std::isfinite(replay.gradient_enclosure_failure_independent));
    EXPECT_DOUBLE_EQ(replay.gradient_enclosure_failure_difference,
        std::abs(replay.gradient_enclosure_failure_independent-primary.gradient(0)));
    EXPECT_DOUBLE_EQ(replay.gradient_difference,replay.gradient_enclosure_failure_difference);
    EXPECT_DOUBLE_EQ(replay.maximum_gradient_difference_primary,primary.gradient(0));
    EXPECT_DOUBLE_EQ(replay.maximum_gradient_difference_independent,replay.gradient_enclosure_failure_independent);
    EXPECT_DOUBLE_EQ(replay.maximum_gradient_difference_bound,replay.gradient_enclosure_failure_bound);
    EXPECT_DOUBLE_EQ(replay.maximum_gradient_difference_ratio,replay.gradient_enclosure_failure_ratio);
    EXPECT_DOUBLE_EQ(replay.primary_width_gradient,std::abs(primary.gradient(0)));
    EXPECT_DOUBLE_EQ(replay.independent_width_gradient,std::abs(replay.gradient_enclosure_failure_independent));
    EXPECT_DOUBLE_EQ(replay.gradient_enclosure_failure_bound,
        1e-13+2e-9*std::abs(primary.gradient(0)));
    EXPECT_DOUBLE_EQ(replay.gradient_enclosure_failure_ratio,
        replay.gradient_enclosure_failure_difference/replay.gradient_enclosure_failure_bound);
    EXPECT_GT(replay.gradient_enclosure_failure_ratio,1.0);
    EXPECT_DOUBLE_EQ(replay.primary_objective,primary.certificate.objective);
    EXPECT_DOUBLE_EQ(replay.primary_residual_norm,primary.residual.norm());
    EXPECT_DOUBLE_EQ(replay.independent_residual_norm,std::sqrt(2.0*replay.independent_objective));

    const auto assessment=n::VerifyJointEndpoint(domain,y,primary,replay,context,
        n::PracticalJointConvergenceCriteria);
    EXPECT_FALSE(assessment.numerical_qualified);
    EXPECT_EQ(assessment.failure,"per-coordinate-replay-enclosure");
}

TEST(JointComponentNumericsTest, IndependentStreamingReplayPreservesActiveBoundaryAndCancellation)
{
    constexpr Eigen::Index rows=33;
    std::vector<std::vector<n::Support>> support(1);
    Vector y(rows),eta=Vector::Constant(1,std::log(.55)),beta(2);
    Vector gaussian(rows),charge(rows);
    for(Eigen::Index row=0;row<rows;++row)
    {
        const double square=6.25*static_cast<double>(row)/static_cast<double>(rows-1);
        support[0].push_back({row,square});
        const auto basis=n::EvaluateKernel(square,.55,2.5);
        gaussian(row)=basis.gaussian; charge(row)=basis.charge;
    }
    n::Domain domain(rows,std::move(support));
    y=-gaussian; beta<<0,charge.dot(y)/charge.squaredNorm();
    const auto context=p::MakeContext(y,1);
    const auto active=n::EvaluateStateStreaming(domain,y,eta,beta,context);
    ASSERT_TRUE(active.valid)<<active.reason;
    ASSERT_TRUE(active.certificate.kkt_passed);
    const auto active_replay=n::CheckIndependentReplay(domain,y,active,context);
    EXPECT_TRUE(active_replay.passed)<<active_replay.reason;
    EXPECT_TRUE(active_replay.feasible);
    EXPECT_DOUBLE_EQ(beta(0),0.0);

    std::vector<std::vector<n::Support>> duplicate_support(2);
    for(Eigen::Index row=0;row<rows;++row)
    {
        const double square=6.25*static_cast<double>(row)/static_cast<double>(rows-1);
        duplicate_support[0].push_back({row,square});
        duplicate_support[1].push_back({row,square});
    }
    n::Domain duplicate_domain(rows,std::move(duplicate_support));
    const Vector zero=Vector::Zero(rows),duplicate_eta=Vector::Constant(2,std::log(.55));
    Vector cancelling_beta(4); cancelling_beta<<0,1e12,0,-1e12;
    const auto duplicate_context=p::MakeContext(zero,2);
    const auto cancelling=n::EvaluateStateStreaming(duplicate_domain,zero,duplicate_eta,cancelling_beta,duplicate_context);
    ASSERT_TRUE(cancelling.valid)<<cancelling.reason;
    const auto cancellation_replay=n::CheckIndependentReplay(duplicate_domain,zero,cancelling,duplicate_context);
    EXPECT_TRUE(cancellation_replay.passed)<<cancellation_replay.reason;
    EXPECT_GT(cancellation_replay.cancellation_ratio,1e10);
    EXPECT_DOUBLE_EQ(cancellation_replay.objective_difference,0.0);
}

TEST(JointComponentNumericsTest, BoundedReplayUsesProjectedWidthKktAndChecksFeasibility)
{
    constexpr Eigen::Index rows=17;
    std::vector<std::vector<n::Support>> support(1);
    Vector eta=Vector::Constant(1,std::log(.55));
    for(Eigen::Index row=0;row<rows;++row)
        support[0].push_back({row,6.25*static_cast<double>(row)/static_cast<double>(rows-1)});
    const n::Domain domain(rows,std::move(support));
    auto evaluate=[&](double truth_width) {
        Vector y(rows);
        for(Eigen::Index row=0;row<rows;++row)
        {
            const double square=6.25*static_cast<double>(row)/static_cast<double>(rows-1);
            const auto basis=n::EvaluateKernel(square,truth_width,2.5);
            y(row)=2.0*basis.gaussian+.5*basis.charge;
        }
        auto context=p::MakeContext(y,1);
        context.eta_upper_bound=std::log(.55);
        const auto primary=n::EvaluateProfile(domain,y,eta,&context);
        EXPECT_TRUE(primary.valid)<<primary.reason;
        const auto replay=n::CheckIndependentReplay(domain,y,primary,context,nullptr,nullptr,
            n::StrictJointConvergenceCriteria);
        return std::tuple{y,context,primary,replay};
    };

    auto [outward_y,outward_context,outward,outward_replay]=evaluate(1.2);
    ASSERT_TRUE(outward.valid);
    ASSERT_TRUE(outward_replay.available && outward_replay.valid)<<outward_replay.reason;
    EXPECT_LT(outward.gradient(0),0.0);
    EXPECT_TRUE(outward_replay.width_bound_applied);
    EXPECT_TRUE(outward_replay.width_bound_feasible);
    EXPECT_TRUE(outward_replay.primary_width_kkt_passed);
    EXPECT_TRUE(outward_replay.independent_width_kkt_passed)<<outward_replay.reason;
    EXPECT_GT(outward_replay.primary_width_gradient,1e-12);
    EXPECT_DOUBLE_EQ(outward_replay.primary_width_projected_kkt,0.0);
    EXPECT_DOUBLE_EQ(outward_replay.independent_width_projected_kkt,0.0);
    EXPECT_DOUBLE_EQ(*n::ProfileGradientInfinityNorm(outward,outward_context),0.0);
    EXPECT_EQ(outward_replay.active_width_bounds,1u);
    EXPECT_GT(outward_replay.boundary_multiplier_max,0.0);
    EXPECT_LE(outward_replay.boundary_gradient_max,0.0);
    EXPECT_DOUBLE_EQ(outward_replay.complementarity_residual,0.0);
    const auto bounded_assessment=n::VerifyJointEndpoint(domain,outward_y,outward,outward_replay,outward_context);
    EXPECT_TRUE(bounded_assessment.numerical_qualified)<<bounded_assessment.failure;
    const auto bounded_evidence=n::AssessmentEvidence(bounded_assessment,
        rhbm_gem::JointEvidenceScope::AssembledGlobal);
    const auto width_stationarity=std::find_if(bounded_evidence.begin(),bounded_evidence.end(),[](const auto & check) {
        return check.name=="width-stationarity";
    });
    ASSERT_NE(width_stationarity,bounded_evidence.end());
    EXPECT_EQ(width_stationarity->status,rhbm_gem::JointCheckStatus::Passed);
    ASSERT_TRUE(width_stationarity->value.has_value());
    EXPECT_DOUBLE_EQ(*width_stationarity->value,0.0);
    const auto bound_feasibility=std::find_if(bounded_evidence.begin(),bounded_evidence.end(),[](const auto & check) {
        return check.name=="width-bound-feasibility";
    });
    ASSERT_NE(bound_feasibility,bounded_evidence.end());
    EXPECT_EQ(bound_feasibility->status,rhbm_gem::JointCheckStatus::Passed);

    auto [inward_y,inward_context,inward,inward_replay]=evaluate(.3);
    ASSERT_TRUE(inward.valid);
    ASSERT_TRUE(inward_replay.available && inward_replay.valid)<<inward_replay.reason;
    EXPECT_GT(inward.gradient(0),0.0);
    EXPECT_FALSE(inward_replay.primary_width_kkt_passed);
    EXPECT_FALSE(inward_replay.independent_width_kkt_passed);
    EXPECT_FALSE(inward_replay.passed);
    const auto rejected=n::VerifyJointEndpoint(domain,inward_y,inward,inward_replay,inward_context);
    EXPECT_FALSE(rejected.numerical_qualified);

    const Vector infeasible_eta=Vector::Constant(1,std::log(.55)+1e-8);
    const auto infeasible=n::EvaluateProfile(domain,outward_y,infeasible_eta,&outward_context);
    ASSERT_TRUE(infeasible.valid);
    const auto infeasible_replay=n::CheckIndependentReplay(domain,outward_y,infeasible,outward_context);
    ASSERT_TRUE(infeasible_replay.available && infeasible_replay.valid);
    EXPECT_FALSE(infeasible_replay.width_bound_feasible);
    EXPECT_FALSE(infeasible_replay.passed);
    EXPECT_NEAR(infeasible_replay.maximum_width_bound_violation,1e-8,1e-12);
    const auto infeasible_assessment=n::VerifyJointEndpoint(domain,outward_y,infeasible,
        infeasible_replay,outward_context);
    const auto infeasible_evidence=n::AssessmentEvidence(infeasible_assessment,
        rhbm_gem::JointEvidenceScope::AssembledGlobal);
    const auto infeasible_bound=std::find_if(infeasible_evidence.begin(),infeasible_evidence.end(),[](const auto & check) {
        return check.name=="width-bound-feasibility";
    });
    ASSERT_NE(infeasible_bound,infeasible_evidence.end());
    EXPECT_EQ(infeasible_bound->status,rhbm_gem::JointCheckStatus::Failed);
    const auto infeasible_trial_trust=n::CheckReplay(domain,outward_y,infeasible,outward_context);
    EXPECT_FALSE(infeasible_trial_trust.passed);
    EXPECT_EQ(infeasible_trial_trust.reason,"width-bound-infeasible");
}

TEST(JointComponentNumericsTest, WidthBoundKktUsesUpperFaceSignAndInteriorGradient)
{
    Vector eta(2),gradient(2);
    eta<<0.0,1.0;
    gradient<<0.0,-2.0;
    auto kkt=n::CalculateWidthBoundKkt(eta,gradient,1.0);
    EXPECT_TRUE(kkt.feasible);
    EXPECT_DOUBLE_EQ(kkt.projected,0.0);
    EXPECT_EQ(kkt.active,1u);
    EXPECT_DOUBLE_EQ(kkt.boundary_gradient,-2.0);
    EXPECT_DOUBLE_EQ(kkt.multiplier,2.0);
    EXPECT_DOUBLE_EQ(kkt.complementarity,0.0);

    gradient(1)=.25;
    kkt=n::CalculateWidthBoundKkt(eta,gradient,1.0);
    EXPECT_TRUE(kkt.feasible);
    EXPECT_DOUBLE_EQ(kkt.projected,.25);

    eta(0)=std::nextafter(1.0,0.0);
    eta(1)=1.0;
    gradient<<-.5,0.0;
    kkt=n::CalculateWidthBoundKkt(eta,gradient,1.0);
    EXPECT_TRUE(kkt.feasible);
    EXPECT_NEAR(kkt.projected,0.0,2e-16);

    eta(1)=1.0+1e-8;
    kkt=n::CalculateWidthBoundKkt(eta,gradient,1.0);
    EXPECT_FALSE(kkt.feasible);
    EXPECT_NEAR(kkt.violation,1e-8,1e-12);

    eta<<.2,.4;
    gradient<<-.3,.5;
    kkt=n::CalculateWidthBoundKkt(eta,gradient,1.0);
    EXPECT_TRUE(kkt.feasible);
    EXPECT_EQ(kkt.active,0u);
    EXPECT_DOUBLE_EQ(kkt.projected,gradient.lpNorm<Eigen::Infinity>());
}
