#include "support/JointOfflineAudit.hpp"
#include "support/JointPrecisionAudit.hpp"
#include "core/command/detail/SimulationGeometry.hpp"
#include "core/command/detail/SimulationManifest.hpp"
#include <fstream>
#include <iomanip>
#include <chrono>
#include <iostream>
#include <sstream>
#include <map>
#include <sys/resource.h>

namespace second_stage_test::matched::certification {
namespace {
namespace j=boost::json;
namespace fs=std::filesystem;
namespace sim=rhbm_gem::core::simulation;
using Vector=Eigen::VectorXd;
using Matrix=Eigen::MatrixXd;
j::value Number(double v) {return std::isfinite(v) ? j::value(v) : j::value(nullptr);}
Vector Parse(const j::value & v) {Vector a(static_cast<Eigen::Index>(v.as_array().size())); for(Eigen::Index k=0;k<a.size();++k) a(k)=j::value_to<double>(v.at(static_cast<std::size_t>(k))); return a;}
void Write(const fs::path & path,const j::value & v) {std::ofstream out(path); out.exceptions(std::ios::failbit|std::ios::badbit); out<<j::serialize(v)<<'\n';}
double Relative(const Vector & a,const Vector & b) {return (a-b).norm()/std::max({1e-12,a.norm(),b.norm()});}
j::object EndpointEvidence(const joint_abc::Evaluation & e)
{
    return {{"valid",e.valid},{"reason",e.reason},{"certificate",e.certificate}};
}
}

void Audit(const joint_abc::Domain & domain,const Vector & y,const j::object & fit,const fs::path & output,const joint_abc::EvaluationContext * provided)
{
    const auto plan=joint_abc::RegisteredAudit(static_cast<Eigen::Index>(domain.atoms.size()),j::value_to<std::string>(fit.at("dataset")),j::value_to<std::string>(fit.at("case")));
    const auto fallback=provided ? joint_abc::EvaluationContext{} : joint_abc::MakeContext(y,static_cast<Eigen::Index>(domain.atoms.size()),"",&plan);
    const auto * context=provided ? provided : &fallback;
    const auto start=std::chrono::steady_clock::now();
    j::object audit{{"schema_version",2},{"case",fit.at("case")},{"dataset",fit.at("dataset")},
        {"legacy_joint_qualified",fit.at("joint_qualified")},{"derivative_status","unavailable"},
        {"derivative_verified",false},{"reason","invalid-endpoint"}};
    if(!fit.contains("primary")) {audit["reason"]="initialization-failed"; Write(output,audit); return;}
    const auto & primary=fit.at("primary");
    if (!primary.at("valid").as_bool()) {Write(output,audit); return;}
    const bool supplied=fit.if_contains("assembled_state_preserved") && fit.at("assembled_state_preserved").as_bool();
    const auto e=supplied ? joint_abc::AtState(domain,y,Parse(primary.at("eta")),Parse(primary.at("beta")),*context) :
        joint_abc::Evaluate(domain,y,Parse(primary.at("eta")),false,context);
    audit["endpoint_trust"]=joint_abc::Trust(domain,y,e,context);
    const double scale=context->scale;
    Vector norms(e.beta.size()); for(Eigen::Index k=0;k<norms.size();++k) norms(k)=e.x.col(k).norm();
    const Vector u=norms.array()*e.beta.array()/scale,dual=(e.x.transpose()*e.residual).array()/norms.array()/scale;
    j::array slack;
    for (Eigen::Index k=0;k<e.beta.size();k+=2) slack.push_back(j::object{{"atom",k/2},{"scaled_a",Number(u(k))},
        {"dual",Number(dual(k))},{"complementarity",Number(u(k)*dual(k))},{"exactly_active",e.beta(k)==0}});
    audit["constraint_evidence"]=slack;
    if (!fit.contains("width_spectrum")) {audit["reason"]="missing-width-spectrum"; Write(output,audit); return;}
    const auto d=joint_abc::Differentiate(e,scale,context);
    if(!d.valid) {audit["reason"]=d.reason; Write(output,audit); return;}
    Matrix directions(e.eta.size(),3); directions.col(0)=Vector::Ones(e.eta.size()).normalized();
    for(Eigen::Index k=0;k<e.eta.size();++k) directions(k,1)=k%2 ? -1 : 1;
    directions.col(1).normalize(); directions.col(2)=Parse(fit.at("width_spectrum").at("weak_directions").at(0));
    if(context->audit.directions.size()) directions=context->audit.directions;
    const bool precision=context->audit.precision;
    if(precision) audit["precision"]=PrecisionAudit(domain,y,e,directions,context,context->audit.block_precision);
    const bool old_verified=fit.at("derivative_verified").as_bool();
    bool verified=old_verified; bool crossed{},missing_interval{};
    j::array ladders;
    const bool expanded_audit=context->audit.expanded_if_unverified && (!old_verified || precision);
    if(expanded_audit)
    {
        verified=true;
        for(Eigen::Index direction=0;direction<directions.cols();++direction)
        {
            std::vector<Vector> finite; std::vector<double> noise; std::vector<bool> valid;
            j::array samples;
            for(int k=0;k<=16;++k)
            {
                const double h=std::ldexp(.01,-k);
                const auto plus=joint_abc::Evaluate(domain,y,e.eta+h*directions.col(direction),false,context);
                const auto minus=joint_abc::Evaluate(domain,y,e.eta-h*directions.col(direction),false,context);
                const auto rp=joint_abc::Evaluate(domain,y,e.eta+h*directions.col(direction),true,context);
                const auto rm=joint_abc::Evaluate(domain,y,e.eta-h*directions.col(direction),true,context);
                const bool solved=plus.valid && minus.valid && rp.valid && rm.valid;
                const bool face=solved && plus.certificate.at("active_atoms")==e.certificate.at("active_atoms") &&
                    minus.certificate.at("active_atoms")==e.certificate.at("active_atoms") &&
                    rp.certificate.at("active_atoms")==e.certificate.at("active_atoms") && rm.certificate.at("active_atoms")==e.certificate.at("active_atoms");
                crossed|=solved && !face; valid.push_back(face);
                finite.push_back(solved ? Vector((rp.residual-rm.residual)/(2*h*scale)) : Vector::Zero(y.size()));
                noise.push_back(solved ? ((plus.residual-rp.residual).norm()+(minus.residual-rm.residual).norm())/(2*h*scale) : std::numeric_limits<double>::infinity());
                const Vector analytic=d.jacobian*directions.col(direction);
                samples.push_back(j::object{{"h",h},{"same_face",face},{"plus",EndpointEvidence(rp)},{"minus",EndpointEvidence(rm)},
                    {"primary_plus",EndpointEvidence(plus)},{"primary_minus",EndpointEvidence(minus)},
                    {"plus_coefficient_difference",solved ? Number(((plus.beta-rp.beta).array().abs()/(1+plus.beta.array().abs().max(rp.beta.array().abs()))).maxCoeff()) : j::value(nullptr)},
                    {"minus_coefficient_difference",solved ? Number(((minus.beta-rm.beta).array().abs()/(1+minus.beta.array().abs().max(rm.beta.array().abs()))).maxCoeff()) : j::value(nullptr)},
                    {"cancellation_ratio",solved ? Number((rp.residual.norm()+rm.residual.norm())/std::max(std::numeric_limits<double>::min(),(rp.residual-rm.residual).norm())) : j::value(nullptr)},
                    {"absolute_error",solved ? Number((finite.back()-analytic).norm()) : j::value(nullptr)},
                    {"relative_error",solved ? Number(Relative(finite.back(),analytic)) : j::value(nullptr)},
                    {"inner_noise",Number(noise.back())},{"sensitivity",Number(analytic.norm())}});
            }
            // Selection uses only successive finite differences and inner disagreement.
            double best=std::numeric_limits<double>::infinity(); int selected=-1; Vector first,second;
            j::array candidates;
            for(int k=0;k<15;++k)
            {
                if(!valid[static_cast<std::size_t>(k)] || !valid[static_cast<std::size_t>(k+1)] || !valid[static_cast<std::size_t>(k+2)]) continue;
                const auto i=static_cast<std::size_t>(k);
                const Vector a=(4*finite[i+1]-finite[i])/3,b=(4*finite[i+2]-finite[i+1])/3;
                const double uncertainty=((a-b).norm()+(noise[i]+5*noise[i+1]+4*noise[i+2])/3)/std::max({1e-12,a.norm(),b.norm()});
                candidates.push_back(j::object{{"index",k},{"estimated_relative_error",Number(uncertainty)}});
                if(uncertainty<best) {best=uncertainty; selected=k; first=a; second=b;}
            }
            const Vector analytic=d.jacobian*directions.col(direction);
            const bool pass=selected>=0 && best<=1e-7 && Relative(first,analytic)<=1e-6 && Relative(second,analytic)<=1e-6;
            missing_interval|=selected<0;
            verified &= pass;
            ladders.push_back(j::object{{"direction",direction},{"samples",samples},{"candidates",candidates},
                {"selected",selected<0 ? j::value(nullptr) : j::value(selected)},{"estimated_relative_error",Number(best)},
                {"first_relative_error",selected>=0 ? Number(Relative(first,analytic)) : j::value(nullptr)},
                {"second_relative_error",selected>=0 ? Number(Relative(second,analytic)) : j::value(nullptr)},{"passed",pass}});
        }
        if(!precision || !audit.at("precision").at("agreement_passed").as_bool() ||
            !audit.at("precision").at("derivative_passed").as_bool()) verified=false;
    }
    audit["ladders"]=ladders; audit["derivative_verified"]=verified;
    audit["derivative_status"]=verified ? "verified" : crossed && missing_interval ? "transition-unverified" : "resolution-unverified";
    audit["reason"]=verified ? (expanded_audit ? "richardson-and-high-precision" : "legacy-two-step-evidence") : "insufficient-derivative-evidence";
    audit["seconds"]=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    struct rusage usage{}; getrusage(RUSAGE_SELF,&usage);
#ifdef __APPLE__
    audit["process_peak_rss_bytes"]=usage.ru_maxrss;
#else
    audit["process_peak_rss_bytes"]=usage.ru_maxrss*1024;
#endif
    Write(output,audit);
}

} // namespace second_stage_test::matched::certification
