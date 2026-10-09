#include <gtest/gtest.h>
#include "core/detail/joint_component/Numerics.hpp"
#include "core/detail/joint_component/TiledDerivative.hpp"
#include "core/detail/joint_component/CompactSvd.hpp"
#include "core/detail/joint_component/Problem.hpp"
#include "core/detail/joint_component/SparseFactor.hpp"
#include "support/JointOperatorWorkload.hpp"
#include "support/JointRuntimeJson.hpp"
#include <boost/json.hpp>
#include <cmath>

namespace {
namespace c=rhbm_gem::core;
namespace n=c::joint_component;
namespace j=boost::json;
using Matrix=Eigen::MatrixXd;

void CompareJson(const j::value & expected,const j::value & actual,const std::string & path={})
{
    SCOPED_TRACE(path);
    ASSERT_EQ(expected.kind(),actual.kind());
    if(expected.is_object())
    {
        const auto & left=expected.as_object(); const auto & right=actual.as_object();
        ASSERT_EQ(left.size(),right.size());
        for(const auto & [key,value]:left)
        {
            ASSERT_TRUE(right.contains(key));
            CompareJson(value,right.at(key),path+"/"+std::string(key));
        }
    }
    else if(expected.is_array())
    {
        const auto & left=expected.as_array(); const auto & right=actual.as_array();
        ASSERT_EQ(left.size(),right.size());
        for(std::size_t k=0;k<left.size();++k) CompareJson(left[k],right[k],path+"/"+std::to_string(k));
    }
    else if(expected.is_number() && actual.is_number())
    {
        if(expected.is_double() || actual.is_double())
        {
            const double a=j::value_to<double>(expected),b=j::value_to<double>(actual);
            EXPECT_LE(std::abs(a-b),1e-10*(1.+std::abs(a)));
        }
        else EXPECT_EQ(expected,actual);
    }
    else EXPECT_EQ(expected,actual);
}

j::array MatrixJson(const Matrix & value)
{
    j::array rows;
    for(Eigen::Index row=0;row<value.rows();++row)
    {
        j::array columns;
        for(Eigen::Index column=0;column<value.cols();++column) columns.push_back(value(row,column));
        rows.push_back(std::move(columns));
    }
    return rows;
}

void CompareAssessment(const n::Assessment & current,const n::Assessment & compact)
{
    auto a=j::value(second_stage_test::matched::runtime_json::Assessment(current)).as_object();
    auto b=j::value(second_stage_test::matched::runtime_json::Assessment(compact)).as_object();
    if(current.widths && compact.widths)
    {
        const Matrix current_projector=current.weak_directions*current.weak_directions.transpose();
        const Matrix compact_projector=compact.weak_directions*compact.weak_directions.transpose();
        a["width_spectrum"].as_object().erase("weak_directions");
        b["width_spectrum"].as_object().erase("weak_directions");
        a["width_spectrum"].as_object()["weak_direction_projector"]=MatrixJson(current_projector);
        b["width_spectrum"].as_object()["weak_direction_projector"]=MatrixJson(compact_projector);
    }
    CompareJson(a,b,"assessment");
}

void CompareEvidence(const n::Assessment & current,const n::Assessment & compact)
{
    const auto a=n::AssessmentEvidence(current,rhbm_gem::JointEvidenceScope::ComponentLocal);
    const auto b=n::AssessmentEvidence(compact,rhbm_gem::JointEvidenceScope::ComponentLocal);
    ASSERT_EQ(a.size(),b.size());
    for(std::size_t k=0;k<a.size();++k)
    {
        EXPECT_EQ(a[k].name,b[k].name); EXPECT_EQ(a[k].status,b[k].status);
        EXPECT_EQ(a[k].scope,b[k].scope); EXPECT_EQ(a[k].reason,b[k].reason);
        ASSERT_EQ(a[k].value.has_value(),b[k].value.has_value());
        if(a[k].value) EXPECT_LE(std::abs(*a[k].value-*b[k].value),1e-10*(1.+std::abs(*a[k].value)));
        ASSERT_EQ(a[k].threshold.has_value(),b[k].threshold.has_value());
        if(a[k].threshold) EXPECT_LE(std::abs(*a[k].threshold-*b[k].threshold),1e-10*(1.+std::abs(*a[k].threshold)));
    }
    EXPECT_EQ(n::ConvergenceStatus(a,rhbm_gem::JointEvidenceScope::ComponentLocal),
              n::ConvergenceStatus(b,rhbm_gem::JointEvidenceScope::ComponentLocal));
}

std::pair<n::Assessment,n::Assessment> CompareRoutes(const n::Domain & domain,n::VectorRef y,
    const n::Evaluation & endpoint,const n::Evaluation & reference,n::EvaluationContext context)
{
    n::JacobianReductionForTesting()=n::JacobianReductionKindForTesting::ObservationTsqr;
    n::AssessmentWorkForTesting()={};
    const auto current=n::AssessEvaluated(domain,y,endpoint,reference,context);
    n::JacobianReductionForTesting()=n::JacobianReductionKindForTesting::CompactStackQr;
    n::AssessmentWorkForTesting()={};
    const auto compact=n::AssessEvaluated(domain,y,endpoint,reference,context);
    CompareAssessment(current,compact);
    CompareEvidence(current,compact);
    return {current,compact};
}

std::pair<n::Assessment,n::Assessment> CompareProjectedRoutes(const n::Domain & domain,n::VectorRef y,
    const n::Evaluation & endpoint,const n::Evaluation & reference,n::EvaluationContext context,
    n::ProjectedReductionKindForTesting candidate_kind=n::ProjectedReductionKindForTesting::StructuredCompactQr)
{
    const auto saved=n::ProjectedReductionForTesting();
    n::JacobianReductionForTesting()=n::JacobianReductionKindForTesting::CompactStackQr;
    n::ProjectedReductionForTesting()=n::ProjectedReductionKindForTesting::ObservationTiledQr;
    n::AssessmentWorkForTesting()={}; n::DerivativeWorkForTesting()={};
    const auto current=n::AssessEvaluated(domain,y,endpoint,reference,context);
    n::ProjectedReductionForTesting()=candidate_kind;
    n::AssessmentWorkForTesting()={}; n::DerivativeWorkForTesting()={};
    const auto candidate=n::AssessEvaluated(domain,y,endpoint,reference,context);
    n::ProjectedReductionForTesting()=saved;
    CompareAssessment(current,candidate);
    CompareEvidence(current,candidate);
    return {current,candidate};
}

n::Evaluation BoundaryEndpoint(double minimum,bool active_face=false)
{
    n::Evaluation endpoint;
    endpoint.valid=true; endpoint.eta=n::Vector::Zero(2); endpoint.beta=n::Vector::Ones(4);
    endpoint.gradient=n::Vector::Zero(2); endpoint.residual=n::Vector::Zero(6);
    endpoint.x.resize(6,4);
    std::vector<Eigen::Triplet<double>> design;
    for(Eigen::Index k=0;k<4;++k) design.emplace_back(k,k,1.);
    endpoint.x.setFromTriplets(design.begin(),design.end());
    if(active_face) {endpoint.beta(0)=0.; endpoint.beta(2)=0.;}
    endpoint.derivative.resize(6,4);
    const Eigen::Index first_column=active_face ? 1 : 0;
    const Eigen::Index second_column=active_face ? 3 : 2;
    endpoint.derivative.insert(4,first_column)=1.;
    endpoint.derivative.insert(5,second_column)=minimum;
    return endpoint;
}

struct CompactMode
{
    n::JacobianReductionKindForTesting saved{n::JacobianReductionForTesting()};
    CompactMode() {n::JacobianReductionForTesting()=n::JacobianReductionKindForTesting::ObservationTsqr;}
    ~CompactMode() {n::JacobianReductionForTesting()=saved;}
};
struct EigenThreads
{
    int saved{Eigen::nbThreads()};
    EigenThreads() {Eigen::setNbThreads(1);}
    ~EigenThreads() {Eigen::setNbThreads(saved);}
};

void CompareComponentResult(const n::ComponentResult & current,const n::ComponentResult & compact)
{
    EXPECT_EQ(current.search.stop_reason,compact.search.stop_reason);
    EXPECT_EQ(current.search.evaluations,compact.search.evaluations);
    EXPECT_EQ(current.search.derivatives,compact.search.derivatives);
    EXPECT_EQ(current.search.accepted,compact.search.accepted);
    EXPECT_EQ(current.search.references,compact.search.references);
    EXPECT_EQ(current.search.stopped,compact.search.stopped);
    EXPECT_EQ(current.search.initial_accepted,compact.search.initial_accepted);
    EXPECT_EQ(current.search.eta.size(),compact.search.eta.size());
    EXPECT_LE((current.search.eta-compact.search.eta).norm(),1e-10*(1.+current.search.eta.norm()));
    EXPECT_EQ(current.search_success,compact.search_success);
    EXPECT_EQ(current.trusted_trial,compact.trusted_trial);
    ASSERT_EQ(current.trusted_state.has_value(),compact.trusted_state.has_value());
    if(current.trusted_state)
        CompareJson(second_stage_test::matched::runtime_json::Endpoint(*current.trusted_state),
                    second_stage_test::matched::runtime_json::Endpoint(*compact.trusted_state),"trusted_state");
    CompareAssessment(current.assessment,compact.assessment);
    CompareEvidence(current.assessment,compact.assessment);
    ASSERT_EQ(current.trusted_assessment.has_value(),compact.trusted_assessment.has_value());
    if(current.trusted_assessment)
    {
        CompareAssessment(*current.trusted_assessment,*compact.trusted_assessment);
        CompareEvidence(*current.trusted_assessment,*compact.trusted_assessment);
    }
    ASSERT_EQ(current.endpoint_trust.has_value(),compact.endpoint_trust.has_value());
    if(current.endpoint_trust)
        CompareJson(second_stage_test::matched::runtime_json::Trust(*current.endpoint_trust),
                    second_stage_test::matched::runtime_json::Trust(*compact.endpoint_trust),"endpoint_trust");
}

struct ProjectedMode
{
    n::ProjectedReductionKindForTesting saved{n::ProjectedReductionForTesting()};
    explicit ProjectedMode(n::ProjectedReductionKindForTesting mode) {n::ProjectedReductionForTesting()=mode;}
    ~ProjectedMode() {n::ProjectedReductionForTesting()=saved;}
};

n::TiledDifferential ProjectedFixture(const Matrix & z,const Matrix & raw,double scale)
{
    n::TiledDifferential out;
    out.free_design=z.sparseView(); out.raw=raw.sparseView(); out.coefficients=z.colPivHouseholderQr().solve(raw);
    out.correction=Matrix::Zero(z.cols(),raw.cols()); out.scale=scale; out.valid=true;
    Eigen::HouseholderQR<Matrix> qr(z);
    out.free_design_factor=qr.matrixQR().topRows(z.cols()).triangularView<Eigen::Upper>();
    out.free_design_response=n::Vector::Zero(z.cols());
    std::vector<Eigen::Index> columns(static_cast<std::size_t>(z.cols()));
    for(Eigen::Index column=0;column<z.cols();++column) columns[static_cast<std::size_t>(column)]=column;
    out.free_design_factor_for_testing=n::FreeDesignFactor::Fixed(out.free_design,columns);
    return out;
}

void CompareProjectedFactors(const std::string & topology,bool near_collinear,bool permute_rows,double scale)
{
    constexpr Eigen::Index rows=48,p=8,m=4;
    Matrix z=Matrix::Zero(rows,p),projected=Matrix::Zero(rows,m);
    for(Eigen::Index column=0;column<p;++column)
    {
        z(column,column)=1.+.03*static_cast<double>(column);
        if(topology=="chain") {z(column+1,column)=.3; z(8+column,column)=.1;}
        else {z((column+1)%p,column)=.08; z(8+column,column)=.27; z(16+column,column)=.11;}
    }
    for(Eigen::Index column=0;column<m;++column)
    {
        projected(24+column,column)=1.+.4*static_cast<double>(column);
        projected(28+column,column)=.2;
        projected(32+column,column)=.1;
    }
    if(near_collinear) projected.col(1)=projected.col(0)+1e-6*projected.col(1);
    Matrix coefficients(p,m);
    for(Eigen::Index row=0;row<p;++row) for(Eigen::Index column=0;column<m;++column)
        coefficients(row,column)=std::sin(.17*static_cast<double>(row+1)*static_cast<double>(column+2));
    const Matrix raw=z*coefficients+projected;
    n::Vector residual(rows);
    for(Eigen::Index row=0;row<rows;++row)
        residual(row)=std::cos(.13*static_cast<double>(row))+.01*static_cast<double>(row);
    Matrix ordered_z=z,ordered_raw=raw; n::Vector ordered_residual=residual;
    if(permute_rows)
    {
        for(Eigen::Index row=0;row<rows;++row)
        {
            const auto source=(row*17)%rows;
            ordered_z.row(row)=z.row(source); ordered_raw.row(row)=raw.row(source);
            ordered_residual(row)=residual(source);
        }
    }
    const auto differential=ProjectedFixture(ordered_z,ordered_raw,scale);
    n::ProjectedReductionForTesting()=n::ProjectedReductionKindForTesting::ObservationTiledQr;
    n::DerivativeWorkForTesting()={};
    const auto current=n::ReduceDerivativeForTesting(differential,ordered_residual,true,
        n::JacobianReductionKindForTesting::CompactStackQr,11);
    n::ProjectedReductionForTesting()=n::ProjectedReductionKindForTesting::StructuredCompactQr;
    n::DerivativeWorkForTesting()={};
    const auto candidate=n::ReduceDerivativeForTesting(differential,ordered_residual,true,
        n::JacobianReductionKindForTesting::CompactStackQr,11);
    ASSERT_TRUE(current.valid); ASSERT_TRUE(candidate.valid); ASSERT_TRUE(candidate.projected_candidate);
    ASSERT_EQ(current.projected.rows(),m); ASSERT_EQ(candidate.projected.rows(),m);
    const Matrix current_gram=current.projected.transpose()*current.projected;
    const Matrix candidate_gram=candidate.projected.transpose()*candidate.projected;
    EXPECT_LE((current_gram-candidate_gram).norm(),2e-11*(1.+current_gram.norm()));
    EXPECT_LE((current.projected_norms-candidate.projected_norms).norm(),2e-11*(1.+current.projected_norms.norm()));
    EXPECT_LE((current.projected.transpose()*current.projected_response-
        candidate.projected.transpose()*candidate.projected_response).norm(),
        2e-11*(1.+(current.projected.transpose()*current.projected_response).norm()));
    const n::RankRequest request{{rows,p,m},m,-1.,n::RankBoundary::SvdNative};
    const auto current_svd=n::EvaluateRank(current.projected,request,nullptr,n::CompactSvdVectors::Right);
    const auto candidate_svd=n::EvaluateRank(candidate.projected,request,nullptr,n::CompactSvdVectors::Right);
    ASSERT_TRUE(current_svd.valid); ASSERT_TRUE(candidate_svd.valid);
    EXPECT_EQ(current_svd.rank,candidate_svd.rank);
    EXPECT_LE((current_svd.singular_values-candidate_svd.singular_values).norm(),
        2e-11*(1.+current_svd.singular_values.norm()));
    Matrix current_normalized=current.projected,candidate_normalized=candidate.projected;
    for(Eigen::Index column=0;column<m;++column)
    {
        current_normalized.col(column)/=current.projected_norms(column);
        candidate_normalized.col(column)/=candidate.projected_norms(column);
    }
    const auto current_normalized_svd=n::EvaluateRank(current_normalized,request);
    const auto candidate_normalized_svd=n::EvaluateRank(candidate_normalized,request);
    ASSERT_TRUE(current_normalized_svd.valid); ASSERT_TRUE(candidate_normalized_svd.valid);
    EXPECT_EQ(current_normalized_svd.rank,candidate_normalized_svd.rank);
    EXPECT_LE((current_normalized_svd.singular_values-candidate_normalized_svd.singular_values).norm(),
        2e-11*(1.+current_normalized_svd.singular_values.norm()));
    const Matrix current_weak=current_svd.right_vectors.rightCols(2)*current_svd.right_vectors.rightCols(2).transpose();
    const Matrix candidate_weak=candidate_svd.right_vectors.rightCols(2)*candidate_svd.right_vectors.rightCols(2).transpose();
    EXPECT_LE((current_weak-candidate_weak).norm(),2e-9);

    n::ProjectedReductionForTesting()=n::ProjectedReductionKindForTesting::ProjectedTailQr;
    n::DerivativeWorkForTesting()={};
    const auto tail=n::ReduceDerivativeForTesting(differential,ordered_residual,true,
        n::JacobianReductionKindForTesting::CompactStackQr,11);
    ASSERT_TRUE(tail.valid);
    EXPECT_TRUE(tail.projected_candidate);
    const Matrix tail_gram=tail.projected.transpose()*tail.projected;
    EXPECT_LE((current_gram-tail_gram).norm(),2e-11*(1.+current_gram.norm()));
    EXPECT_LE((current.projected_norms-tail.projected_norms).norm(),
        2e-11*(1.+current.projected_norms.norm()));
    EXPECT_LE((current.projected.transpose()*current.projected_response-
        tail.projected.transpose()*tail.projected_response).norm(),
        2e-11*(1.+(current.projected.transpose()*current.projected_response).norm()));
    const auto tail_svd=n::EvaluateRank(tail.projected,request,nullptr,n::CompactSvdVectors::Right);
    ASSERT_TRUE(tail_svd.valid); EXPECT_EQ(current_svd.rank,tail_svd.rank);
    EXPECT_LE((current_svd.singular_values-tail_svd.singular_values).norm(),
        2e-11*(1.+current_svd.singular_values.norm()));
    Matrix tail_normalized=tail.projected;
    for(Eigen::Index column=0;column<m;++column) tail_normalized.col(column)/=tail.projected_norms(column);
    const auto tail_normalized_svd=n::EvaluateRank(tail_normalized,request);
    ASSERT_TRUE(tail_normalized_svd.valid); EXPECT_EQ(current_normalized_svd.rank,tail_normalized_svd.rank);
    EXPECT_LE((current_normalized_svd.singular_values-tail_normalized_svd.singular_values).norm(),
        2e-11*(1.+current_normalized_svd.singular_values.norm()));
    const Matrix tail_weak=tail_svd.right_vectors.rightCols(2)*tail_svd.right_vectors.rightCols(2).transpose();
    EXPECT_LE((current_weak-tail_weak).norm(),2e-9);
}

TEST(JointProjectedWidthTest, StructuredFactorMatchesObservationTiledQrFixtures)
{
    ProjectedMode mode(n::ProjectedReductionKindForTesting::ObservationTiledQr);
    for(const std::string topology:{"chain","cube"})
        for(const bool near_collinear:{false,true})
            for(const bool permute_rows:{false,true})
                for(const double scale:{.25,1.,4.})
                {
                    SCOPED_TRACE(topology+" near_collinear="+std::to_string(near_collinear)+
                        " permute_rows="+std::to_string(permute_rows)+" scale="+std::to_string(scale));
                    CompareProjectedFactors(topology,near_collinear,permute_rows,scale);
                }
}

TEST(JointProjectedWidthTest, FullAssessmentAndReturnedStateParityAcrossSmallLattices)
{
    CompactMode mode;
    ProjectedMode projected_mode(n::ProjectedReductionKindForTesting::ObservationTiledQr);
    EigenThreads threads;
    for(const std::string topology:{"chain","cube"}) for(const int atoms:{8,32})
    {
        SCOPED_TRACE(topology+"-"+std::to_string(atoms));
        const c::JointProblem problem(second_stage_test::OperatorWorkload(topology,atoms));
        const auto & data=c::JointProblemAccess::Get(problem);
        auto context=data.context;
        const auto search=n::SearchProfile(data.domain,data.y,n::Vector::Constant(atoms,.55),context);

        n::JacobianReductionForTesting()=n::JacobianReductionKindForTesting::CompactStackQr;
        n::ProjectedReductionForTesting()=n::ProjectedReductionKindForTesting::ObservationTiledQr;
        n::AssessmentWorkForTesting()={}; n::DerivativeWorkForTesting()={};
        const auto current=n::AssessComponentSearch(data.domain,data.y,context,search);
        n::ProjectedReductionForTesting()=n::ProjectedReductionKindForTesting::StructuredCompactQr;
        n::AssessmentWorkForTesting()={}; n::DerivativeWorkForTesting()={};
        const auto candidate=n::AssessComponentSearch(data.domain,data.y,context,search);
        CompareComponentResult(current,candidate);
        const auto & reduction=n::DerivativeWorkForTesting().projected_reduction;
        EXPECT_GT(reduction.attempts,0);
        EXPECT_GT(reduction.accepted,0);
        EXPECT_EQ(reduction.attempts,reduction.accepted+reduction.fallbacks);

        n::ProjectedReductionForTesting()=n::ProjectedReductionKindForTesting::ProjectedTailQr;
        n::AssessmentWorkForTesting()={}; n::DerivativeWorkForTesting()={};
        const auto tail=n::AssessComponentSearch(data.domain,data.y,context,search);
        CompareComponentResult(current,tail);
        const auto & tail_reduction=n::DerivativeWorkForTesting().projected_reduction;
        EXPECT_GT(tail_reduction.attempts,0);
        EXPECT_EQ(tail_reduction.attempts,tail_reduction.accepted+tail_reduction.fallbacks);
        EXPECT_GT(tail_reduction.accepted,0);
    }
}

TEST(JointProjectedWidthTest, BoundaryAndActiveFaceFallbackParity)
{
    CompactMode mode;
    EigenThreads threads;
    const n::Domain domain(6,{});
    const n::Vector y=n::Vector::Zero(6);
    auto context=n::CreateContext(y,2); context.rank={6,4,2};
    const double cutoff=context.rank.Relative(2);
    for(const double minimum:{std::nextafter(cutoff,0.),cutoff,
        std::nextafter(cutoff,std::numeric_limits<double>::infinity())})
    {
        SCOPED_TRACE("projected width boundary="+std::to_string(minimum));
        const auto endpoint=BoundaryEndpoint(minimum);
        const auto [current,candidate]=CompareProjectedRoutes(domain,y,endpoint,endpoint,context);
        EXPECT_EQ(current.jacobian->rank,candidate.jacobian->rank);
        const auto & reduction=n::DerivativeWorkForTesting().projected_reduction;
        EXPECT_EQ(reduction.attempts,1);
        EXPECT_EQ(reduction.accepted,0);
        EXPECT_EQ(reduction.fallbacks,1);
    }
    const auto endpoint=BoundaryEndpoint(1e-3,true);
    const auto [current,candidate]=CompareProjectedRoutes(domain,y,endpoint,endpoint,context);
    EXPECT_EQ(current.jacobian->rank,candidate.jacobian->rank);
    EXPECT_EQ(current.primary.beta(0),0.);
    EXPECT_EQ(current.primary.beta(2),0.);
    const auto & reduction=n::DerivativeWorkForTesting().projected_reduction;
    EXPECT_EQ(reduction.attempts,1);
    EXPECT_EQ(reduction.accepted,1);
    EXPECT_EQ(reduction.fallbacks,0);

    for(const double minimum:{std::nextafter(cutoff,0.),cutoff,
        std::nextafter(cutoff,std::numeric_limits<double>::infinity())})
    {
        const auto boundary_endpoint=BoundaryEndpoint(minimum);
        const auto [boundary_current,boundary_tail]=CompareProjectedRoutes(domain,y,boundary_endpoint,boundary_endpoint,context,
            n::ProjectedReductionKindForTesting::ProjectedTailQr);
        EXPECT_EQ(boundary_current.jacobian->rank,boundary_tail.jacobian->rank);
        const auto & tail_reduction=n::DerivativeWorkForTesting().projected_reduction;
        EXPECT_EQ(tail_reduction.attempts,1);
        EXPECT_EQ(tail_reduction.accepted,0);
        EXPECT_EQ(tail_reduction.fallbacks,1);
        EXPECT_EQ(tail_reduction.fallback_reason,"rank-decision-boundary");
    }
    const auto tail_active_endpoint=BoundaryEndpoint(1e-3,true);
    const auto [active_current,active_tail]=CompareProjectedRoutes(domain,y,tail_active_endpoint,
        tail_active_endpoint,context,n::ProjectedReductionKindForTesting::ProjectedTailQr);
    EXPECT_EQ(active_current.jacobian->rank,active_tail.jacobian->rank);
    const auto & active_tail_reduction=n::DerivativeWorkForTesting().projected_reduction;
    EXPECT_EQ(active_tail_reduction.attempts,1);
    EXPECT_EQ(active_tail_reduction.accepted,1);
    EXPECT_EQ(active_tail_reduction.fallbacks,0);
}

TEST(JointCompactAssessmentTest, FullAssessmentAndReturnedStateParityAcrossSmallLattices)
{
    CompactMode mode;
    ProjectedMode projected_mode(n::ProjectedReductionKindForTesting::ObservationTiledQr);
    EigenThreads threads;
    for(const std::string topology:{"chain","cube"}) for(const int atoms:{8,32})
    {
        SCOPED_TRACE(topology+"-"+std::to_string(atoms));
        const c::JointProblem problem(second_stage_test::OperatorWorkload(topology,atoms));
        const auto & data=c::JointProblemAccess::Get(problem);
        auto context=data.context;
        const auto search=n::SearchProfile(data.domain,data.y,n::Vector::Constant(atoms,.55),context);

        n::JacobianReductionForTesting()=n::JacobianReductionKindForTesting::ObservationTsqr;
        n::AssessmentWorkForTesting()={};
        const auto current=n::AssessComponentSearch(data.domain,data.y,context,search);
        n::JacobianReductionForTesting()=n::JacobianReductionKindForTesting::CompactStackQr;
        n::AssessmentWorkForTesting()={};
        const auto compact=n::AssessComponentSearch(data.domain,data.y,context,search);
        CompareComponentResult(current,compact);
        EXPECT_GT(n::AssessmentWorkForTesting().compact_attempts,0);
        EXPECT_EQ(n::AssessmentWorkForTesting().compact_attempts,
                  n::AssessmentWorkForTesting().compact_accepted+
                  n::AssessmentWorkForTesting().compact_boundary_fallbacks+
                  n::AssessmentWorkForTesting().compact_other_fallbacks);
    }
}

TEST(JointCompactAssessmentTest, FullAssessmentAndReturnedStateParityAtAcceptanceSizes)
{
    CompactMode mode;
    ProjectedMode projected_mode(n::ProjectedReductionKindForTesting::ObservationTiledQr);
    EigenThreads threads;
    for(const std::string topology:{"chain","cube"}) for(const int atoms:{128,256})
    {
        SCOPED_TRACE(topology+"-"+std::to_string(atoms));
        const c::JointProblem problem(second_stage_test::OperatorWorkload(topology,atoms));
        const auto & data=c::JointProblemAccess::Get(problem);
        auto context=data.context;
        const auto search=n::SearchProfile(data.domain,data.y,n::Vector::Constant(atoms,.55),context);

        n::JacobianReductionForTesting()=n::JacobianReductionKindForTesting::ObservationTsqr;
        n::AssessmentWorkForTesting()={};
        const auto current=n::AssessComponentSearch(data.domain,data.y,context,search);
        n::JacobianReductionForTesting()=n::JacobianReductionKindForTesting::CompactStackQr;
        n::AssessmentWorkForTesting()={};
        const auto compact=n::AssessComponentSearch(data.domain,data.y,context,search);
        CompareComponentResult(current,compact);
        EXPECT_GT(n::AssessmentWorkForTesting().compact_attempts,0);
        EXPECT_EQ(n::AssessmentWorkForTesting().compact_attempts,
                  n::AssessmentWorkForTesting().compact_accepted+
                  n::AssessmentWorkForTesting().compact_boundary_fallbacks+
                  n::AssessmentWorkForTesting().compact_other_fallbacks);
        n::JacobianReductionForTesting()=n::JacobianReductionKindForTesting::CompactStackQr;
        n::ProjectedReductionForTesting()=n::ProjectedReductionKindForTesting::ObservationTiledQr;
        n::AssessmentWorkForTesting()={}; n::DerivativeWorkForTesting()={};
        const auto projected_current=n::AssessComponentSearch(data.domain,data.y,context,search);
        n::ProjectedReductionForTesting()=n::ProjectedReductionKindForTesting::StructuredCompactQr;
        n::AssessmentWorkForTesting()={}; n::DerivativeWorkForTesting()={};
        const auto projected_candidate=n::AssessComponentSearch(data.domain,data.y,context,search);
        CompareComponentResult(projected_current,projected_candidate);
        const auto & projected=n::DerivativeWorkForTesting().projected_reduction;
        EXPECT_GT(projected.attempts,0);
        EXPECT_EQ(projected.attempts,projected.accepted+projected.fallbacks);
    }
}

TEST(JointCompactAssessmentTest, RankBoundaryAssessmentFallsBackForOneUlpBand)
{
    CompactMode mode;
    EigenThreads threads;
    const n::Domain domain(6,{});
    const n::Vector y=n::Vector::Zero(6);
    auto context=n::CreateContext(y,2);
    context.rank={6,4,2};
    const double cutoff=context.rank.Relative(2);
    const std::array<std::pair<const char *,double>,5> cases{{
        {"below",.5*cutoff},
        {"just-below",std::nextafter(cutoff,0.)},
        {"boundary",cutoff},
        {"just-above",std::nextafter(cutoff,std::numeric_limits<double>::infinity())},
        {"above",2.*cutoff}}};
    for(const auto & [name,minimum]:cases)
    {
        SCOPED_TRACE(name);
        const auto endpoint=BoundaryEndpoint(minimum);
        const auto [current,compact]=CompareRoutes(domain,y,endpoint,endpoint,context);
        ASSERT_TRUE(current.jacobian && compact.jacobian);
        EXPECT_EQ(current.jacobian->rank,compact.jacobian->rank);
        ASSERT_EQ(current.jacobian->singular_values.size(),compact.jacobian->singular_values.size());
        EXPECT_LE((current.jacobian->singular_values-compact.jacobian->singular_values).norm(),
            1e-14*std::max(1.,current.jacobian->singular_values.norm()));
        EXPECT_LE((current.correction-compact.correction).norm(),1e-12*(1.+current.correction.norm()));
        const bool near=std::string(name)=="just-below" || std::string(name)=="boundary" ||
            std::string(name)=="just-above";
        EXPECT_EQ(n::AssessmentWorkForTesting().compact_boundary_fallbacks,near ? 1 : 0);
        EXPECT_EQ(n::AssessmentWorkForTesting().compact_accepted,near ? 0 : 1);
        EXPECT_EQ(n::AssessmentWorkForTesting().compact_attempts,1);
    }
}

TEST(JointCompactAssessmentTest, CancellationAndActiveFaceAssessmentParity)
{
    CompactMode mode;
    EigenThreads threads;
    {
        const c::JointProblem problem(second_stage_test::OperatorWorkload("chain",8));
        const auto & data=c::JointProblemAccess::Get(problem);
        auto endpoint=n::EvaluateProfile(data.domain,data.y,n::Vector::Constant(8,std::log(.55)),false,&data.context);
        ASSERT_TRUE(endpoint.valid);
        endpoint.derivative=endpoint.x;
        const auto reference=n::EvaluateProfile(data.domain,data.y,endpoint.eta,true,&data.context);
        ASSERT_TRUE(reference.valid);
        const auto [current,compact]=CompareProjectedRoutes(data.domain,data.y,endpoint,reference,data.context);
        EXPECT_EQ(current.primary.beta,compact.primary.beta);
        EXPECT_EQ(n::DerivativeWorkForTesting().projected_reduction.attempts,0);
        const auto [tail_current,tail]=CompareProjectedRoutes(data.domain,data.y,endpoint,reference,data.context,
            n::ProjectedReductionKindForTesting::ProjectedTailQr);
        CompareAssessment(tail_current,tail);
        EXPECT_EQ(tail_current.primary.beta,tail.primary.beta);
        EXPECT_EQ(n::DerivativeWorkForTesting().projected_reduction.attempts,0);
    }
    {
        const n::Domain domain(6,{});
        const n::Vector y=n::Vector::Zero(6);
        auto context=n::CreateContext(y,2); context.rank={6,4,2};
        auto endpoint=BoundaryEndpoint(1e-3,true);
        const auto [current,compact]=CompareRoutes(domain,y,endpoint,endpoint,context);
        EXPECT_EQ(current.primary.beta(0),0.);
        EXPECT_EQ(current.primary.beta(2),0.);
        EXPECT_EQ(current.jacobian->rank,compact.jacobian->rank);
    }
}

} // namespace
