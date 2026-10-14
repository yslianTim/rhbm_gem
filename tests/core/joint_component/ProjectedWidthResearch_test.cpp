#include <gtest/gtest.h>
#include "support/ProjectedWidthResearch.hpp"
#include "support/MatrixFreeProjectedWidth.hpp"
#include "support/StructuredProjectionResearch.hpp"
#include "support/BoundedRankResearch.hpp"
#include "core/detail/joint_component/CompactSvd.hpp"
#include "core/detail/joint_component/ResourceWork.hpp"
#include "core/detail/joint_component/TiledDerivative.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <Eigen/QR>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>

namespace {
namespace pw=second_stage_test::projected_width;
namespace n=rhbm_gem::core::joint_component;

std::string Values(const Eigen::VectorXd & x)
{
    std::ostringstream out;
    out<<std::setprecision(17)<<"[";
    for(Eigen::Index k=0;k<x.size();++k) {if(k) out<<","; out<<x(k);}
    return out.str()+"]";
}

TEST(JointComponentProjectedWidthResearchTest, DenseOracleBaselinesAreReproducible)
{
    std::size_t full_rank{},deficient{},threshold_sensitive{},false_full_rank{};
    for(const auto & fixture:pw::Fixtures())
    {
        SCOPED_TRACE(fixture.name);
        const auto started=std::chrono::steady_clock::now();
        const auto oracle=pw::EvaluateDenseOracle(fixture);
        const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
        ASSERT_TRUE(oracle.valid)<<oracle.reason;
        EXPECT_EQ(oracle.free_columns,fixture.free_design.cols());
        EXPECT_EQ(oracle.width_columns,fixture.width_derivative.cols());
        EXPECT_EQ(oracle.design.full_rank,oracle.design.rank==oracle.free_columns);
        EXPECT_EQ(oracle.width.full_rank,oracle.width.rank==oracle.width_columns);
        EXPECT_GE(oracle.projection_residual,0);
        if(oracle.design.full_rank && oracle.width.full_rank) ++full_rank;
        if(!oracle.design.full_rank || !oracle.width.full_rank) ++deficient;
        if(oracle.width.threshold_sensitive) ++threshold_sensitive;
        if(oracle.width.full_rank && fixture.name=="exactly-rank-deficient") ++false_full_rank;
        std::cout<<"projected-width-baseline fixture="<<fixture.name
            <<" rows="<<oracle.rows<<" free_columns="<<oracle.free_columns
            <<" full_design_columns="<<fixture.full_design_columns
            <<" width_columns="<<oracle.width_columns
            <<" policy_rows="<<fixture.policy.rows
            <<" structural_nonzeros="<<oracle.structural_nonzeros
            <<" supported_rows="<<oracle.supported_rows
            <<" column_norms="<<Values(oracle.design_column_norms)
            <<" design_singular_values="<<Values(oracle.design.singular_values)
            <<" design_rank="<<oracle.design.rank<<" design_threshold="<<oracle.design.threshold
            <<" design_threshold_sensitive="<<oracle.design.threshold_sensitive
            <<" width_singular_values="<<Values(oracle.width.singular_values)
            <<" width_rank="<<oracle.width.rank<<" width_strict_rank="<<oracle.width.strict_rank
            <<" width_threshold="<<oracle.width.threshold
            <<" width_threshold_sensitive="<<oracle.width.threshold_sensitive
            <<" projection_residual="<<oracle.projection_residual
            <<" full_rank="<<(oracle.width.full_rank && oracle.design.full_rank)
            <<" target_columns="<<fixture.target_columns
            <<" active_amplitude_boundary="<<fixture.active_amplitude_boundary
            <<" seconds="<<seconds<<"\n";
    }
    EXPECT_GE(full_rank,3u);
    EXPECT_GE(deficient,1u);
    EXPECT_GE(threshold_sensitive,1u);
    EXPECT_EQ(false_full_rank,0u);
    const auto peak=n::ProcessPeakRssBytes();
    std::cout<<"projected-width-baseline-summary full_rank_controls="<<full_rank
        <<" deficient_controls="<<deficient<<" threshold_sensitive_controls="<<threshold_sensitive
        <<" false_full_rank="<<false_full_rank
        <<" process_peak_rss_bytes="<<(peak ? std::to_string(*peak) : "unavailable")<<"\n";
    RecordProperty("full_rank_controls",std::to_string(full_rank));
    RecordProperty("deficient_controls",std::to_string(deficient));
    RecordProperty("threshold_sensitive_controls",std::to_string(threshold_sensitive));
    RecordProperty("false_full_rank",std::to_string(false_full_rank));
}

TEST(JointComponentProjectedWidthResearchTest, ActiveFaceAndHaloDimensionsRemainExplicit)
{
    const auto fixtures=pw::Fixtures();
    const auto active=std::find_if(fixtures.begin(),fixtures.end(),[](const auto & f) {
        return f.name=="active-amplitude-boundary";
    });
    const auto mixed=std::find_if(fixtures.begin(),fixtures.end(),[](const auto & f) {
        return f.name=="mixed-target-halo";
    });
    ASSERT_NE(active,fixtures.end()); ASSERT_NE(mixed,fixtures.end());
    EXPECT_TRUE(active->active_amplitude_boundary);
    EXPECT_LT(active->free_design.cols(),active->full_design_columns);
    EXPECT_EQ(active->width_derivative.cols(),active->target_columns);
    EXPECT_EQ(mixed->target_columns,3);
    EXPECT_EQ(mixed->free_design.cols(),4);
}

TEST(JointComponentProjectedWidthResearchTest, ValuesOnlyAndRightVectorModesPreserveRankEvidence)
{
    for(const auto & fixture:pw::Fixtures())
    {
        SCOPED_TRACE(fixture.name);
        const auto oracle=pw::EvaluateDenseOracle(fixture);
        ASSERT_TRUE(oracle.valid)<<oracle.reason;
        const n::RankRequest request{fixture.policy,fixture.width_derivative.cols(),-1,n::RankBoundary::SvdNative};
        const auto started=std::chrono::steady_clock::now();
        const auto values=n::EvaluateRank(oracle.projected_width,request,nullptr,n::CompactSvdVectors::None);
        const double values_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
        const auto vectors_started=std::chrono::steady_clock::now();
        const auto with_vectors=n::EvaluateRank(oracle.projected_width,request,nullptr,n::CompactSvdVectors::Right);
        const double vectors_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-vectors_started).count();
        ASSERT_TRUE(values.valid && with_vectors.valid);
        EXPECT_EQ(values.rank,with_vectors.rank);
        if(!oracle.width.threshold_sensitive) EXPECT_EQ(values.rank,oracle.width.rank);
        else EXPECT_FALSE(n::CompactRankDecisionSafe(values));
        EXPECT_DOUBLE_EQ(values.threshold,with_vectors.threshold);
        EXPECT_EQ(n::CompactRankDecisionSafe(values),n::CompactRankDecisionSafe(with_vectors));
        EXPECT_LE((values.singular_values-with_vectors.singular_values).lpNorm<Eigen::Infinity>(),
            2e-12*std::max(1.0,with_vectors.singular_values(0)));
        EXPECT_EQ(values.right_vectors.size(),0);
        EXPECT_EQ(with_vectors.right_vectors.cols(),std::min(oracle.rows,oracle.width_columns));
        const auto vector_bytes=static_cast<std::size_t>(with_vectors.right_vectors.size())*sizeof(double);
        std::cout<<"projected-width-svd-mode fixture="<<fixture.name
            <<" rows="<<oracle.rows<<" width_columns="<<oracle.width_columns
            <<" rank="<<values.rank<<" threshold="<<values.threshold
            <<" threshold_sensitive="<<!n::CompactRankDecisionSafe(values)
            <<" right_vector_bytes="<<vector_bytes
            <<" values_only_seconds="<<values_seconds<<" right_vector_seconds="<<vectors_seconds
            <<" process_peak_rss_bytes="<<(n::ProcessPeakRssBytes() ?
                std::to_string(*n::ProcessPeakRssBytes()) : "unavailable")<<"\n";
    }
    constexpr Eigen::Index original_rows=1000;
    const n::RankPolicy policy{original_rows,3,3};
    const double threshold=policy.Relative(3);
    for(const auto & item:std::vector<std::pair<const char *,double>>{
        {"below",std::nextafter(threshold,0.0)},
        {"equal",threshold},
        {"above",std::nextafter(threshold,std::numeric_limits<double>::infinity())}})
    {
        pw::Matrix boundary=pw::Matrix::Zero(3,3);
        boundary(0,0)=1; boundary(1,1)=.5; boundary(2,2)=item.second;
        const n::RankRequest request{policy,3,-1,n::RankBoundary::StrictGreater};
        const auto values=n::EvaluateRank(boundary,request,nullptr,n::CompactSvdVectors::None);
        const auto with_vectors=n::EvaluateRank(boundary,request,nullptr,n::CompactSvdVectors::Right);
        ASSERT_TRUE(values.valid && with_vectors.valid)<<item.first;
        EXPECT_DOUBLE_EQ(values.threshold,threshold);
        EXPECT_EQ(values.rank,with_vectors.rank);
        EXPECT_EQ(values.rank,item.second>threshold ? 3 : 2);
        EXPECT_EQ(n::CompactRankDecisionSafe(values),n::CompactRankDecisionSafe(with_vectors));
        std::cout<<"projected-width-svd-boundary case="<<item.first
            <<" rows="<<original_rows<<" columns=3 threshold="<<values.threshold
            <<" minimum="<<values.singular_values(2)<<" rank="<<values.rank
            <<" threshold_sensitive="<<!n::CompactRankDecisionSafe(values)<<"\n";
    }
}

TEST(JointComponentProjectedWidthResearchTest, MatrixFreeActionsMatchDenseAndTiledQrOracles)
{
    for(const auto & fixture:pw::Fixtures())
    {
        SCOPED_TRACE(fixture.name);
        const auto started=std::chrono::steady_clock::now();
        const auto oracle=pw::EvaluateDenseOracle(fixture);
        ASSERT_TRUE(oracle.valid)<<oracle.reason;
        const n::Sparse design=fixture.free_design.sparseView();
        const n::Sparse derivative=fixture.width_derivative.sparseView();
        auto op=pw::BuildExactOracleOperator(design,derivative,fixture.policy);
        ASSERT_TRUE(op.available)<<op.reason;
        EXPECT_EQ(op.status,"Completed");

        const Eigen::VectorXd v=Eigen::VectorXd::LinSpaced(op.width_columns,1.,2.).normalized();
        const Eigen::VectorXd u=Eigen::VectorXd::LinSpaced(op.normalized_design.rows(),2.,4.).normalized();
        const Eigen::VectorXd forward=op.Apply(v);
        const Eigen::VectorXd direct=oracle.projected_width*v;
        const Eigen::VectorXd adjoint=op.ApplyAdjoint(u);
        const Eigen::VectorXd direct_adjoint=oracle.projected_width.transpose()*u;
        const double absolute_forward=(forward-direct).norm();
        const double absolute_adjoint=(adjoint-direct_adjoint).norm();
        const double scale=std::max({1.0,derivative.norm()*v.norm(),derivative.norm()*u.norm()});
        const double roundoff=64*std::numeric_limits<double>::epsilon()*
            static_cast<double>(std::max(op.normalized_design.rows(),op.normalized_design.cols()))*scale;
        const double relative_forward=absolute_forward/std::max(direct.norm(),std::numeric_limits<double>::min());
        const double relative_adjoint=absolute_adjoint/std::max(direct_adjoint.norm(),std::numeric_limits<double>::min());
        EXPECT_LE(absolute_forward,roundoff);
        EXPECT_LE(absolute_adjoint,roundoff);
        EXPECT_LE(op.OrthogonalityResidual(forward),roundoff);

        const Eigen::VectorXd projected=op.Project(u);
        const Eigen::VectorXd projected_twice=op.Project(projected);
        const double idempotence=(projected_twice-projected).norm();
        EXPECT_LE(idempotence,roundoff);
        const double lhs=forward.dot(u),rhs=v.dot(adjoint);
        EXPECT_LE(std::abs(lhs-rhs),roundoff);

        n::TiledDifferential differential;
        differential.free_design=op.normalized_design;
        differential.raw=derivative;
        differential.scale=1.; differential.valid=true; differential.projected_only=true;
        Eigen::JacobiSVD<pw::Matrix> design_svd(pw::Matrix(op.normalized_design),
            Eigen::ComputeThinU|Eigen::ComputeThinV);
        differential.coefficients=design_svd.solve(pw::Matrix(derivative));
        const auto reduced=n::ReduceProjectedWidthDerivative(differential,2);
        ASSERT_TRUE(reduced.valid)<<reduced.reason;
        const auto tiled_svd=n::EvaluateRank(reduced.projected,
            {fixture.policy,fixture.width_derivative.cols(),-1,n::RankBoundary::SvdNative});
        ASSERT_TRUE(tiled_svd.valid);
        const bool rank_parity=tiled_svd.rank==oracle.width.rank;
        if(!oracle.width.threshold_sensitive) EXPECT_TRUE(rank_parity);
        else EXPECT_FALSE(n::CompactRankDecisionSafe(tiled_svd));
        EXPECT_LE((tiled_svd.singular_values-oracle.width.singular_values).lpNorm<Eigen::Infinity>(),roundoff);

        std::cout<<"matrix-free-operator fixture="<<fixture.name
            <<" rows="<<op.normalized_design.rows()<<" free_columns="<<op.normalized_design.cols()
            <<" width_columns="<<op.width_columns<<" design_rank="<<op.design_rank
            <<" width_rank="<<tiled_svd.rank<<" rank_threshold="<<tiled_svd.threshold
            <<" oracle_rank_parity="<<rank_parity
            <<" threshold_sensitive="<<!n::CompactRankDecisionSafe(tiled_svd)
            <<" forward_absolute_error="<<absolute_forward<<" forward_relative_error="<<relative_forward
            <<" adjoint_absolute_error="<<absolute_adjoint<<" adjoint_relative_error="<<relative_adjoint
            <<" adjoint_inner_product_difference="<<std::abs(lhs-rhs)
            <<" orthogonality_residual="<<op.OrthogonalityResidual(forward)
            <<" projection_idempotence_error="<<idempotence
            <<" projection_residual="<<oracle.projection_residual
            <<" roundoff_over_rank_threshold="<<(roundoff/oracle.width.threshold)
            <<" oracle_width_singular_values="<<Values(oracle.width.singular_values)
            <<" tiled_width_singular_values="<<Values(tiled_svd.singular_values)
            <<" estimated_peak_bytes="<<op.estimated_peak_bytes
            <<" basis_bytes="<<op.basis_bytes
            <<" projection_calls="<<op.projection_calls
            <<" forward_calls="<<op.forward_calls<<" adjoint_calls="<<op.adjoint_calls
            <<" seconds="<<std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count()<<"\n";
    }
}

TEST(JointComponentProjectedWidthResearchTest, BoundedProjectionPathFailsClosedAtPreflight)
{
    const auto fixture=pw::Fixtures().front();
    const n::Sparse design=fixture.free_design.sparseView();
    const n::Sparse derivative=fixture.width_derivative.sparseView();
    const auto limited=pw::BuildExactOracleOperator(design,derivative,fixture.policy,0,1);
    EXPECT_FALSE(limited.available);
    EXPECT_EQ(limited.status,"ResourceLimit");
    EXPECT_EQ(limited.reason,"estimated-peak-exceeds-planning-budget");

    const auto bounded=pw::PreflightBoundedProjection(5'938'087,74'812,37'406,
        149'622,149'622,600'000'000);
    EXPECT_EQ(bounded.resource_status,"Pass");
    EXPECT_EQ(bounded.numerical_status,"Unavailable");
    EXPECT_EQ(bounded.reason,"no-bounded-global-projector-with-validated-rank-contract");
    std::cout<<"matrix-free-giant-action-preflight rows=5,938,087 free_columns=74,812 width_columns=37,406"
        <<" estimated_action_bytes="<<bounded.estimated_action_bytes
        <<" estimated_peak_bytes="<<bounded.estimated_peak_bytes
        <<" planning_budget_bytes="<<bounded.budget_bytes
        <<" resource_status="<<bounded.resource_status
        <<" numerical_status="<<bounded.numerical_status<<" reason="<<bounded.reason
        <<" projection_calls=0\n";
}

TEST(JointComponentStructuredProjectionResearchTest, SmallProjectionsMatchDenseOracleAndPreserveGlobalCoupling)
{
    for(const auto & fixture:pw::Fixtures())
    {
        SCOPED_TRACE(fixture.name);
        const auto oracle=pw::EvaluateDenseOracle(fixture);
        ASSERT_TRUE(oracle.valid)<<oracle.reason;
        const n::Sparse design=fixture.free_design.sparseView();
        const auto resident=n::CurrentRssBytes().value_or(0);
        const Eigen::VectorXd u=Eigen::VectorXd::LinSpaced(design.rows(),1.,2.).normalized();
        const auto result=pw::ProjectWithSparseNormalEquations(design,u,resident);
        ASSERT_TRUE(result.preflight_passed)<<result.reason;
        EXPECT_EQ(result.rank_status,"Unavailable");
        if(!oracle.design.full_rank)
        {
            std::cout<<"structured-projection-control fixture="<<fixture.name
                <<" rows="<<result.rows<<" columns="<<result.columns
                <<" exact_blocks="<<result.exact_blocks
                <<" separator_columns="<<result.separator_columns
                <<" maximum_interior_columns="<<result.maximum_interior_columns
                <<" estimated_peak_bytes="<<result.estimated_peak_bytes
                <<" rank_status="<<result.rank_status
                <<" status="<<result.status<<" reason="<<result.reason
                <<" seconds="<<result.seconds<<"\n";
            continue;
        }
        ASSERT_TRUE(result.valid)<<result.reason;
        {
            pw::Matrix normalized=fixture.free_design;
            for(Eigen::Index column=0;column<normalized.cols();++column)
                normalized.col(column)/=normalized.col(column).norm();
            Eigen::HouseholderQR<pw::Matrix> qr(normalized);
            const pw::Matrix q=qr.householderQ()*pw::Matrix::Identity(normalized.rows(),normalized.cols());
            const Eigen::VectorXd qr_projection=q*(q.transpose()*u);
            Eigen::JacobiSVD<pw::Matrix> dense_svd(normalized,Eigen::ComputeThinU);
            const auto design_rank=oracle.design.rank;
            const pw::Matrix dense_q=dense_svd.matrixU().leftCols(design_rank);
            const Eigen::VectorXd dense_projection=dense_q*(dense_q.transpose()*u);
            const double error=(result.projection-dense_projection).norm();
            const double relative=error/std::max(dense_projection.norm(),std::numeric_limits<double>::min());
            const double qr_error=(result.projection-qr_projection).norm();
            const double dense_svd_error=(result.projection-dense_projection).norm();
            const double objective_difference=std::abs((u-result.projection).squaredNorm()-
                (u-qr_projection).squaredNorm());
            const double scale=std::max(1.0,u.norm());
            const double roundoff=32*std::numeric_limits<double>::epsilon()*
                static_cast<double>(std::max(design.rows(),design.cols()))*
                std::max(1.0,result.diagonal_condition_estimate)*scale;
            EXPECT_LE(error,roundoff);
            EXPECT_LE(qr_error,roundoff);
            EXPECT_LE(dense_svd_error,roundoff);
            EXPECT_LE(objective_difference,roundoff*scale);
            EXPECT_LE(result.projection_residual,roundoff);
            EXPECT_LE(result.idempotence_error,roundoff);
            if(fixture.name=="weak-separator") EXPECT_GT(result.separator_columns,0);
            if(fixture.name=="multiple-components") EXPECT_EQ(result.exact_blocks,4);
            std::cout<<"structured-projection-control fixture="<<fixture.name
                <<" rows="<<result.rows<<" columns="<<result.columns
                <<" input_nonzeros="<<result.input_nonzeros<<" exact_blocks="<<result.exact_blocks
                <<" articulation_columns="<<result.articulation_columns
                <<" separator_columns="<<result.separator_columns
                <<" maximum_interior_columns="<<result.maximum_interior_columns
                <<" interface_dimension="<<result.interface_dimension
                <<" gram_nonzeros="<<result.gram_nonzeros
                <<" factor_nonzeros="<<result.factor_nonzeros
                <<" predicted_fill_ratio="<<(result.gram_nonzeros ?
                    static_cast<double>(result.factor_nonzeros)/static_cast<double>(result.gram_nonzeros) : 0.0)
                <<" maximum_factor_size="<<result.maximum_factor_size
                <<" diagonal_condition_estimate="<<result.diagonal_condition_estimate
                <<" projection_absolute_error="<<error<<" projection_relative_error="<<relative
                <<" dense_qr_projection_error="<<qr_error
                <<" dense_svd_projection_error="<<dense_svd_error
                <<" projection_objective_difference="<<objective_difference
                <<" projection_residual="<<result.projection_residual
                <<" idempotence_error="<<result.idempotence_error
                <<" roundoff_over_rank_threshold="<<(roundoff/oracle.width.threshold)
                <<" estimated_peak_bytes="<<result.estimated_peak_bytes
                <<" rss_before_bytes="<<(result.current_rss_before_bytes ?
                    std::to_string(*result.current_rss_before_bytes) : "unavailable")
                <<" process_peak_rss_bytes="<<(result.process_peak_rss_bytes ?
                    std::to_string(*result.process_peak_rss_bytes) : "unavailable")
                <<" seconds="<<result.seconds<<" rank_status="<<result.rank_status
                <<" status="<<result.status<<" reason="<<result.reason<<"\n";
        }
    }
}

TEST(JointComponentProjectedWidthResearchTest, RankCertificatesRequireStrictEnclosedBounds)
{
    constexpr Eigen::Index original_rows=1000,width_columns=3;
    const n::RankPolicy policy{original_rows,6,width_columns};
    const double threshold=policy.Relative(width_columns);
    const auto full=[&](double sigma_min,double sigma_max) {
        pw::ProjectedSpectrumEnvelope envelope;
        envelope.sigma_min_estimate=sigma_min; envelope.sigma_max_estimate=sigma_max;
        envelope.error_enclosure_valid=true;
        return pw::ClassifyProjectedRankEnvelope(policy,width_columns,envelope);
    };
    const auto full_rank=full(1.,1.);
    const auto ill_conditioned=full(4*threshold,1.);
    const auto exactly_deficient=full(0.,std::sqrt(3.));
    const auto equality=full(threshold,1.);
    const auto below=full(std::nextafter(threshold,0.),1.);
    const auto above=full(std::nextafter(threshold,std::numeric_limits<double>::infinity()),1.);
    EXPECT_EQ(full_rank.status,pw::RankCertificateStatus::FullRank);
    EXPECT_EQ(ill_conditioned.status,pw::RankCertificateStatus::FullRank);
    EXPECT_EQ(exactly_deficient.status,pw::RankCertificateStatus::RankDeficient);
    EXPECT_EQ(equality.status,pw::RankCertificateStatus::Unavailable);
    EXPECT_EQ(below.status,pw::RankCertificateStatus::Unavailable);
    EXPECT_EQ(above.status,pw::RankCertificateStatus::Unavailable);

    pw::ProjectedSpectrumEnvelope incomplete;
    incomplete.sigma_min_estimate=1.; incomplete.sigma_max_estimate=1.;
    const auto no_enclosure=pw::ClassifyProjectedRankEnvelope(policy,width_columns,incomplete);
    EXPECT_EQ(no_enclosure.status,pw::RankCertificateStatus::Unavailable);
    pw::ProjectedSpectrumEnvelope crossing;
    crossing.sigma_min_estimate=threshold*1.1; crossing.sigma_max_estimate=1.;
    crossing.projection_error_bound=threshold*.2;
    crossing.matrix_vector_error_bound=threshold*.2;
    crossing.spectral_error_bound=threshold*.1;
    crossing.floating_point_error_bound=threshold*.1;
    crossing.normalization_error_bound=threshold*.1;
    crossing.error_enclosure_valid=true;
    const auto threshold_overlap=pw::ClassifyProjectedRankEnvelope(policy,width_columns,crossing);
    EXPECT_EQ(threshold_overlap.status,pw::RankCertificateStatus::Unavailable);

    std::size_t oracle_full{},oracle_deficient{},unavailable_cases{},false_certifications{};
    for(const auto & fixture:pw::Fixtures())
    {
        const auto oracle=pw::EvaluateDenseOracle(fixture);
        ASSERT_TRUE(oracle.valid)<<fixture.name;
        if(oracle.width.full_rank) ++oracle_full; else ++oracle_deficient;

        // No experimental operator in this study supplies rigorous spectral and
        // floating-point enclosures for these fixture matrices.  Dense-oracle
        // agreement is a comparison, not a certificate.
        const pw::ProjectedSpectrumEnvelope envelope;
        const auto certificate=pw::ClassifyProjectedRankEnvelope(
            {fixture.policy.rows,fixture.policy.design_columns,fixture.width_derivative.cols()},
            fixture.width_derivative.cols(),envelope);
        if(certificate.status==pw::RankCertificateStatus::Unavailable) ++unavailable_cases;
        EXPECT_EQ(certificate.status,pw::RankCertificateStatus::Unavailable)<<fixture.name;
        std::cout<<"projected-width-certificate fixture="<<fixture.name
            <<" oracle_rank="<<oracle.width.rank<<" columns="<<oracle.width_columns
            <<" oracle_threshold="<<oracle.width.threshold
            <<" oracle_threshold_sensitive="<<oracle.width.threshold_sensitive
            <<" certificate="<<pw::RankCertificateStatusName(certificate.status)
            <<" reason="<<certificate.reason<<"\n";
    }
    EXPECT_GT(oracle_full,0u); EXPECT_GT(oracle_deficient,0u);
    EXPECT_EQ(false_certifications,0u);
    EXPECT_EQ(unavailable_cases,pw::Fixtures().size());
    std::cout<<"projected-width-certificate-summary oracle_full="<<oracle_full
        <<" oracle_deficient="<<oracle_deficient<<" proven_full=0 proven_deficient=0 unavailable_cases="
        <<unavailable_cases<<" unavailable_rate="
        <<static_cast<double>(unavailable_cases)/static_cast<double>(pw::Fixtures().size())
        <<" false_certifications="<<false_certifications<<"\n";
    RecordProperty("false_certifications",std::to_string(false_certifications));
    RecordProperty("unavailable_cases",std::to_string(unavailable_cases));
}

TEST(JointComponentStructuredProjectionResearchTest, GiantComponentPreflightsRejectBeforeAllocation)
{
    const auto overflow=pw::PreflightStructuredProjection(1,
        std::numeric_limits<std::size_t>::max(),1,0);
    EXPECT_EQ(overflow.status,"ResourceLimit");
    EXPECT_EQ(overflow.reason,"memory-estimate-overflow");

    constexpr std::size_t rows=5'938'087,atoms=37'406,design_columns=74'812;
    constexpr std::size_t support_entries=74'811,nonzeros=149'622,resident=600'000'000;
    n::MemoryShape shape;
    shape.component_atoms=atoms; shape.component_observations=rows;
    shape.component_support_entries=support_entries; shape.component_structural_nonzeros=nonzeros;
    shape.component_structural_supported_rows=atoms; shape.component_design_columns=design_columns;
    shape.component_width_columns=atoms; shape.component_target_atoms=atoms;
    shape.component_target_columns=3*atoms; shape.resident_bytes=resident;

    const auto projected=n::PlanMemory(shape,n::MemoryStage::ProjectedWidthAssessment);
    ASSERT_EQ(projected.status,n::MemoryStatus::Reject);
    const auto rank=n::PlanMemory(shape,n::MemoryStage::RankAssessment);
    EXPECT_EQ(rank.status,n::MemoryStatus::Reject);
    EXPECT_EQ(rank.estimated_peak_bytes,projected.estimated_peak_bytes);
    EXPECT_EQ(projected.free_design_factor_bytes,design_columns*design_columns*sizeof(double));
    EXPECT_EQ(projected.projection_coefficient_bytes,design_columns*atoms*sizeof(double));
    EXPECT_EQ(projected.projected_width_factor_bytes,atoms*atoms*sizeof(double));

    const auto matrix_free=pw::PreflightBoundedProjection(rows,design_columns,atoms,
        nonzeros,nonzeros,resident);
    EXPECT_EQ(matrix_free.resource_status,"Pass");
    EXPECT_EQ(matrix_free.numerical_status,"Unavailable");
    const auto structured=pw::PreflightStructuredProjection(rows,design_columns,nonzeros,
        resident);
    EXPECT_EQ(structured.status,"ResourceLimit");
    EXPECT_GT(structured.estimated_peak_bytes,n::JointMemoryPlanningBudgetBytes);

    std::cout<<"projected-width-giant-preflight rows="<<rows<<" atoms="<<atoms
        <<" components=1 free_columns="<<design_columns<<" width_columns="<<atoms
        <<" structural_nonzeros="<<nonzeros<<" structural_supported_rows="<<atoms
        <<" full_design_columns="<<design_columns<<" structural_rank_upper_bound="<<atoms
        <<" full_design_structural_status=Failed"
        <<" free_factor_bytes="<<projected.free_design_factor_bytes
        <<" free_svd_workspace_bytes="<<projected.free_design_svd_workspace_bytes
        <<" projection_coefficients_bytes="<<projected.projection_coefficient_bytes
        <<" width_factor_bytes="<<projected.projected_width_factor_bytes
        <<" width_svd_workspace_bytes="<<projected.projected_width_svd_workspace_bytes
        <<" sparse_workset_bytes="<<projected.sparse_design_bytes
        <<" qr_tile_workspace_bytes="<<projected.qr_tile_workspace_bytes
        <<" temporary_bytes="<<projected.temporary_bytes
        <<" production_estimated_peak_bytes="<<projected.estimated_peak_bytes
        <<" production_status="<<n::MemoryStatusName(projected.status)
        <<" action_estimated_peak_bytes="<<matrix_free.estimated_peak_bytes
        <<" action_resource_status="<<matrix_free.resource_status
        <<" action_numerical_status="<<matrix_free.numerical_status
        <<" structured_estimated_peak_bytes="<<structured.estimated_peak_bytes
        <<" structured_status="<<structured.status<<" structured_reason="<<structured.reason
        <<" planning_budget_bytes="<<n::JointMemoryPlanningBudgetBytes
        <<" observed_process_peak_bytes="<<(n::ProcessPeakRssBytes() ?
            std::to_string(*n::ProcessPeakRssBytes()) : "unavailable")
        <<" allocations_performed=0 rank_status=Unavailable\n";

    for(const std::size_t free_columns:{design_columns,std::size_t{56'109},atoms})
    {
        auto active_shape=shape; active_shape.component_free_columns=free_columns;
        const auto active_production=n::PlanMemory(active_shape,n::MemoryStage::ProjectedWidthAssessment);
        EXPECT_EQ(active_production.status,n::MemoryStatus::Reject);
        const auto active= pw::PreflightStructuredProjection(rows,free_columns,nonzeros,resident);
        EXPECT_EQ(active.status,"ResourceLimit");
        std::cout<<"production-active-face-preflight free_columns="<<free_columns
            <<" free_factor_bytes="<<active_production.free_design_factor_bytes
            <<" free_svd_workspace_bytes="<<active_production.free_design_svd_workspace_bytes
            <<" projection_coefficients_bytes="<<active_production.projection_coefficient_bytes
            <<" width_factor_bytes="<<active_production.projected_width_factor_bytes
            <<" width_svd_workspace_bytes="<<active_production.projected_width_svd_workspace_bytes
            <<" sparse_workset_bytes="<<active_production.sparse_design_bytes
            <<" qr_tile_workspace_bytes="<<active_production.qr_tile_workspace_bytes
            <<" temporary_bytes="<<active_production.temporary_bytes
            <<" estimated_peak_bytes="<<active_production.estimated_peak_bytes
            <<" status="<<n::MemoryStatusName(active_production.status)<<"\n";
        std::cout<<"structured-active-face-preflight free_columns="<<free_columns
            <<" estimated_peak_bytes="<<active.estimated_peak_bytes
            <<" status="<<active.status<<" reason="<<active.reason<<"\n";
    }
}

#define STRUCTURED_TOPOLOGY_CASE(test_name,topology_name,atom_count) \
TEST(JointComponentStructuredProjectionResearchTest,test_name) \
{ \
    const auto design=second_stage_test::bounded_rank::SyntheticDesign(topology_name,atom_count); \
    const auto resident=n::CurrentRssBytes().value_or(0); \
    const Eigen::VectorXd u=Eigen::VectorXd::LinSpaced(design.rows(),1.,2.).normalized(); \
    const auto result=pw::ProjectWithSparseNormalEquations(design,u,resident); \
    EXPECT_EQ(result.rank_status,"Unavailable"); \
    if(result.estimated_peak_bytes>n::JointMemoryPlanningBudgetBytes) \
    {EXPECT_FALSE(result.preflight_passed); EXPECT_EQ(result.status,"ResourceLimit");} \
    std::cout<<"structured-projection-topology topology="<<topology_name<<" atoms="<<atom_count \
        <<" rows="<<design.rows()<<" columns="<<design.cols()<<" input_nonzeros="<<design.nonZeros() \
        <<" exact_blocks="<<result.exact_blocks<<" separator_columns="<<result.separator_columns \
        <<" maximum_interior_columns="<<result.maximum_interior_columns \
        <<" interface_dimension="<<result.interface_dimension<<" gram_nonzeros="<<result.gram_nonzeros \
        <<" factor_nonzeros="<<result.factor_nonzeros<<" maximum_factor_size="<<result.maximum_factor_size \
        <<" estimated_peak_bytes="<<result.estimated_peak_bytes \
        <<" rss_before_bytes="<<(result.current_rss_before_bytes ? std::to_string(*result.current_rss_before_bytes) : "unavailable") \
        <<" rss_after_bytes="<<(result.current_rss_after_bytes ? std::to_string(*result.current_rss_after_bytes) : "unavailable") \
        <<" process_peak_rss_bytes="<<(result.process_peak_rss_bytes ? std::to_string(*result.process_peak_rss_bytes) : "unavailable") \
        <<" projection_residual="<<result.projection_residual<<" idempotence_error="<<result.idempotence_error \
        <<" diagonal_condition_estimate="<<result.diagonal_condition_estimate<<" seconds="<<result.seconds \
        <<" numerical_rank_status="<<result.rank_status<<" status="<<result.status<<" reason="<<result.reason<<"\n"; \
}
STRUCTURED_TOPOLOGY_CASE(Chain128,"chain",128)
STRUCTURED_TOPOLOGY_CASE(Cube128,"cube",128)
STRUCTURED_TOPOLOGY_CASE(Chain256,"chain",256)
STRUCTURED_TOPOLOGY_CASE(Cube256,"cube",256)
STRUCTURED_TOPOLOGY_CASE(Chain512,"chain",512)
STRUCTURED_TOPOLOGY_CASE(Cube512,"cube",512)
STRUCTURED_TOPOLOGY_CASE(Chain1024,"chain",1024)
STRUCTURED_TOPOLOGY_CASE(Cube1024,"cube",1024)
#undef STRUCTURED_TOPOLOGY_CASE
}
