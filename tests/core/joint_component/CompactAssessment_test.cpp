#include <gtest/gtest.h>
#include "core/detail/joint_component/Numerics.hpp"
#include "core/detail/joint_component/TiledDerivative.hpp"
#include "core/detail/joint_component/Problem.hpp"
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

TEST(JointCompactAssessmentTest, FullAssessmentAndReturnedStateParityAcrossSmallLattices)
{
    CompactMode mode;
    EigenThreads threads;
    for(const std::string topology:{"chain","cube"}) for(const int atoms:{8,32})
    {
        SCOPED_TRACE(topology+"-"+std::to_string(atoms));
        const c::JointProblem problem(second_stage_test::OperatorWorkload(topology,atoms));
        const auto & data=c::JointProblemAccess::Get(problem);
        auto context=data.context; context.search.method=n::SearchMethod::OperatorPcg;
        context.search.preconditioner=n::PreconditionerKind::Schwarz;
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
        const auto [current,compact]=CompareRoutes(data.domain,data.y,endpoint,reference,data.context);
        EXPECT_EQ(current.primary.beta,compact.primary.beta);
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
