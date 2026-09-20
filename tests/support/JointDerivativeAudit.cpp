#include "support/JointDerivativeAudit.hpp"
#include "support/JointRuntimeJson.hpp"

namespace second_stage_test::matched::certification {
boost::json::object AssessWithDerivativeAudit(const joint_abc::Domain & domain,const Eigen::VectorXd & y,
    const Eigen::VectorXd & eta,const joint_abc::EvaluationContext & policy,const Eigen::VectorXd * beta)
{
    namespace j=boost::json;
    using namespace joint_abc::runtime;
    const auto * context=&policy; const double scale=policy.scale;
    const auto endpoint=beta ? EvaluateState(domain,y,eta,*beta,policy) : EvaluateProfile(domain,y,eta,false,context);
    const auto reference=EvaluateProfile(domain,y,eta,true,context);
    const auto assessment=AssessEvaluated(domain,y,endpoint,reference,policy,beta!=nullptr);
    auto out=runtime_json::Assessment(assessment);
    out.erase("runtime_convergence"); out.erase("runtime_checks"); out.erase("runtime_failure");
    out["joint_qualified"]=false; out["qualification_failure"]=assessment.failure;
    if(!assessment.jacobian) return out;
    const auto differential=DifferentiateProfile(endpoint,scale,context);
    j::array checks;
    std::vector<Vector> directions{Vector::Ones(eta.size()).normalized(),Vector(eta.size()),assessment.weak_directions.col(0)};
    for(Eigen::Index k=0;k<eta.size();++k) directions[1](k)=k%2 ? -1 : 1;
    directions[1].normalize();
    if(context->audit.directions.size())
    {directions.clear(); for(Eigen::Index k=0;k<context->audit.directions.cols();++k) directions.push_back(context->audit.directions.col(k));}
    bool verified=true;
    for(std::size_t k=0;k<directions.size();++k) for(double h:{1e-4,5e-5})
    {
        const auto plus=EvaluateProfile(domain,y,eta+h*directions[k],true,context),minus=EvaluateProfile(domain,y,eta-h*directions[k],true,context);
        bool passed=false,same_face=false; double error=std::numeric_limits<double>::infinity();
        if(plus.valid && minus.valid)
        {
            same_face=plus.certificate.active_atoms==endpoint.certificate.active_atoms && minus.certificate.active_atoms==endpoint.certificate.active_atoms;
            const Vector analytic=differential.jacobian*directions[k];
            const Vector finite=(plus.residual-minus.residual)/(2*h*scale);
            error=(finite-analytic).norm()/std::max({1e-12,finite.norm(),analytic.norm()});
            passed=same_face && error<=1e-6;
        }
        verified &= passed;
        checks.push_back(j::object{{"direction",k},{"h",h},{"relative_l2_difference",runtime_json::Number(error)},
            {"same_active_face",same_face},{"passed",passed},{"plus_valid",plus.valid},{"minus_valid",minus.valid}});
    }

    out["derivative_checks"]=checks; out["directional_evaluations"]=2*checks.size(); out["derivative_verified"]=verified;
    out["qualification_checks"]=j::object{{"inner",assessment.inner},{"b_gradient",assessment.gradient},
        {"local_correction",assessment.local},{"identified",assessment.identified},{"derivative",verified}};
    out["joint_qualified"]=assessment.inner && assessment.gradient && assessment.local && assessment.identified && verified;
    out["qualification_failure"]=!assessment.inner ? "inner-solve" : !verified ? "derivative-unverified" :
        !assessment.identified ? "width-unidentified" : !assessment.gradient || !assessment.local ? "b-not-stationary" : "none";
    return out;
}
}
