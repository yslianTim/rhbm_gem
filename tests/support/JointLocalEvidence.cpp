#include "support/JointLocalEvidence.hpp"
#include "support/JointDerivativeAudit.hpp"
#include <cmath>
#include <limits>

namespace second_stage_test::matched::certification {
namespace {
namespace j=boost::json;
using Vector=Eigen::VectorXd;
j::value Number(double v) {return std::isfinite(v) ? j::value(v) : j::value(nullptr);}
Vector Parse(const j::value & v) {Vector out(v.as_array().size()); for(Eigen::Index k=0;k<out.size();++k) out(k)=j::value_to<double>(v.at(static_cast<std::size_t>(k))); return out;}
}
LocalAuditState PrepareLocalAudit(const joint_abc::Domain & domain,const Vector & y,
    const j::object & search,const joint_abc::EvaluationContext & child)
{
    LocalAuditState out{child,search,{}};
    out.context.audit.directions.resize(0,0);
    out.scope={{"scope","component-local"},{"component_id",search.at("component_id")},
        {"parent_snapshot_sha256",child.snapshot_hash},{"weak_direction_available",false},
        {"direction_source","actual-component-endpoint"},{"state",nullptr}};
    if(!search.at("usable_state").as_bool())
    {
        out.fit["primary"]=j::object{{"valid",false},{"reason","missing-trusted-component-state"}};
        out.fit.erase("width_spectrum"); out.fit["joint_qualified"]=false;
        return out;
    }
    const auto & state=search.at("last_trusted_state");
    const Vector eta=Parse(state.at("eta")),beta=Parse(state.at("beta"));
    out.scope["state"]=j::object{{"eta",state.at("eta")},{"beta",state.at("beta")}};
    const auto assessment=AssessWithDerivativeAudit(domain,y,eta,out.context,&beta);
    for(const auto & field:assessment) out.fit[field.key()]=field.value();
    out.fit["assembled_state_preserved"]=true;
    const auto control=joint_abc::Evaluate(domain,y,eta,false,&out.context);
    const double difference=control.valid ? ((beta-control.beta).array().abs()/
        (1+beta.array().abs().max(control.beta.array().abs()))).maxCoeff() : std::numeric_limits<double>::infinity();
    out.fit["assembled_profile_agrees"]=difference<=1e-10;
    out.fit["assembled_profile_difference"]=Number(difference);
    out.fit["joint_qualified"]=out.fit.at("joint_qualified").as_bool() && difference<=1e-10 && search.at("search_success").as_bool();
    if(out.fit.contains("width_spectrum"))
    {
        auto & directions=out.context.audit.directions; directions.resize(eta.size(),3);
        directions.col(0)=Vector::Ones(eta.size()).normalized();
        for(Eigen::Index k=0;k<eta.size();++k) directions(k,1)=k%2 ? -1 : 1;
        directions.col(1).normalize();
        directions.col(2)=Parse(out.fit.at("width_spectrum").at("weak_directions").at(0));
        out.scope["weak_direction_available"]=true;
        const auto e=joint_abc::AtState(domain,y,eta,beta,out.context);
        const auto d=joint_abc::Differentiate(e,child.scale,&out.context);
        j::array diagnostics;
        for(Eigen::Index k=0;k<directions.cols();++k)
            diagnostics.push_back(j::object{{"direction",k},{"norm",directions.col(k).norm()},
                {"jacobian_direction_norm",d.valid ? Number((d.jacobian*directions.col(k)).norm()) : j::value(nullptr)}});
        out.scope["directions"]=diagnostics;
    }
    return out;
}

}
