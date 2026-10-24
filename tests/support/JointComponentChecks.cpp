#include "support/JointComponentChecks.hpp"
#include "support/JointRuntimeJson.hpp"
#include "core/detail/joint_component/FixedNeighborBlockCoordinate.hpp"
#include "core/detail/joint_component/Problem.hpp"
#include <sys/resource.h>
#include <algorithm>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <cmath>
#include <chrono>
#include <stdexcept>

namespace second_stage_test::matched::joint_abc {
namespace {
namespace j=boost::json;
std::size_t Index(Eigen::Index i) {return static_cast<std::size_t>(i);}
j::array Indices(const std::vector<Eigen::Index> & v)
{j::array out; for(auto i:v) out.push_back(i); return out;}
}
rhbm_gem::core::JointFitResult FitBoundedJointComponentsForTesting(
    const rhbm_gem::core::JointProblem & problem,const std::vector<double> & initial_widths,double maximum_width)
{
    if(!(std::isfinite(maximum_width) && maximum_width>0.0))
        throw std::invalid_argument("Bounded Joint test fit requires a finite positive Width upper bound.");
    runtime::FixedNeighborSearchPolicy policy;
    policy.eta_upper_bound=std::log(maximum_width);
    policy.configure_search_for_testing=[](runtime::FixedNeighborPolicy & search) {
        search.core_atoms=runtime::ProductionFixedNeighborCoreAtoms;
        search.maximum_sweeps=runtime::ProductionFixedNeighborMaximumSweeps;
        search.convergence_criteria=runtime::PracticalJointConvergenceCriteria;
    };
    return runtime::FitFixedNeighborComponents(problem,initial_widths,policy);
}
ComponentPartition BuildPartition(const Domain & d,const std::vector<std::string> & ids) {return runtime::Partition(d,ids);}
Eigen::VectorXd Select(const Eigen::VectorXd & v,const std::vector<Eigen::Index> & indices) {return runtime::SelectValues(v,indices);}
EvaluationContext ComponentContext(const EvaluationContext & p,const ComponentView & v,bool independent) {return runtime::ChildContext(p,v,independent);}
j::object Census(const Domain & domain,const ComponentPartition & partition,const EvaluationContext & context)
{
    j::array components; std::size_t memberships{},largest_atoms{},largest_rows{};
    for(const auto & c:partition.components)
    {
        std::size_t entries{}; for(const auto & a:c.domain.atoms) entries+=a.size(); memberships+=entries;
        largest_atoms=std::max(largest_atoms,c.atoms.size()); largest_rows=std::max(largest_rows,c.rows.size());
        std::vector<Eigen::Index> atom_map(partition.mappings->atom_component.size(),-1),row_map(partition.mappings->row_component.size(),-1);
        for(auto a:c.atoms) atom_map[Index(a)]=c.LocalAtom(a);
        for(auto r:c.rows) row_map[Index(r)]=c.LocalRow(r);
        components.push_back(j::object{{"id",c.id},{"atoms",Indices(c.atoms)},{"rows",Indices(c.rows)},
            {"parent_atom_to_local",Indices(atom_map)},{"parent_row_to_local",Indices(row_map)},
            {"atom_count",c.atoms.size()},{"row_count",c.rows.size()},{"memberships",entries}});
    }
    double constant{}; if(context.observations) for(auto r:partition.constant_rows) constant+=std::pow((*context.observations)(r),2)/2;
    return {{"schema_version",1},{"snapshot_sha256",context.snapshot_hash},{"atoms",domain.atoms.size()},
        {"rows",domain.rows},{"memberships",memberships},{"component_count",components.size()},{"components",components},
        {"largest_atom_count",largest_atoms},{"largest_row_count",largest_rows},{"atom_component",Indices(partition.mappings->atom_component)},
        {"row_component",Indices(partition.mappings->row_component)},{"constant_rows",Indices(partition.constant_rows)},
        {"constant_objective",constant},{"unobserved_atoms",Indices(partition.unobserved_atoms)}};
}

namespace {
using Vector=Eigen::VectorXd;
using Matrix=Eigen::MatrixXd;
bool IdentityPartition(const ComponentPartition & p)
{return p.components.size()==1 && p.constant_rows.empty() && p.unobserved_atoms.empty();}
std::vector<Eigen::Index> Coefficients(const ComponentView & view)
{std::vector<Eigen::Index> out; for(auto a:view.atoms) {out.push_back(2*a); out.push_back(2*a+1);} return out;}
Sparse Slice(const Sparse & x,const ComponentView & view,const std::vector<Eigen::Index> & columns)
{
    Sparse out(static_cast<Eigen::Index>(view.rows.size()),static_cast<Eigen::Index>(columns.size()));
    std::vector<Eigen::Triplet<double>> entries;
    for(std::size_t k=0;k<columns.size();++k) for(Sparse::InnerIterator e(x,columns[k]);e;++e)
    {
        const auto row=view.LocalRow(e.row());
        if(row<0) throw std::runtime_error("Cross-component numeric entry.");
        entries.emplace_back(row,static_cast<Eigen::Index>(k),e.value());
    }
    out.setFromTriplets(entries.begin(),entries.end()); return out;
}
Sparse Free(const Evaluation & e)
{
    Eigen::Index count{}; std::vector<Eigen::Triplet<double>> entries;
    for(Eigen::Index k=0;k<e.beta.size();++k) if(k%2 || e.beta(k)>0)
    {for(Sparse::InnerIterator v(e.x,k);v;++v) entries.emplace_back(v.row(),count,v.value()); ++count;}
    Sparse out(e.x.rows(),count); out.setFromTriplets(entries.begin(),entries.end()); return out;
}
double Relative(const Matrix & a,const Matrix & b)
{return (a-b).norm()/std::max({1e-12,a.norm(),b.norm()});}
double Scaled(const Vector & a,const Vector & b)
{return ((a-b).array().abs()/(1+a.array().abs().max(b.array().abs()))).maxCoeff();}
j::object State(const Evaluation & e)
{return {{"valid",e.valid},{"reason",e.reason},{"certificate",e.certificate}};}
j::array Values(const Vector & v)
{j::array out; for(double x:v) out.push_back(std::isfinite(x) ? j::value(x) : j::value(nullptr)); return out;}
Vector Parse(const j::value & v)
{Vector out(static_cast<Eigen::Index>(v.as_array().size())); for(Eigen::Index k=0;k<out.size();++k) out(k)=j::value_to<double>(v.at(Index(k))); return out;}
j::object FaceEvidence(const Domain & domain,const Vector & y,const Evaluation & e,const EvaluationContext & context)
{
    Matrix directions(e.eta.size(),2); directions.col(0)=Vector::Ones(e.eta.size()).normalized();
    for(Eigen::Index k=0;k<e.eta.size();++k) directions(k,1)=k%2 ? -1 : 1;
    directions.col(1).normalize();
    bool stable=true; j::array probes,saved;
    for(Eigen::Index k=0;k<directions.cols();++k)
    {
        j::array direction; for(double v:directions.col(k)) direction.push_back(v); saved.push_back(direction);
        for(double h:{1e-4,5e-5}) for(int sign:{-1,1})
        {
            const auto p=Evaluate(domain,y,e.eta+sign*h*directions.col(k),&context);
            const bool same=p.valid && p.certificate.at("active_atoms")==e.certificate.at("active_atoms");
            stable &= same;
            probes.push_back(j::object{{"direction",k},{"step",sign*h},{"same_face",same},
                {"primary",State(p)}});
        }
    }
    return {{"stable",stable},{"directions",saved},{"probes",probes},
        {"rule","Both declared step sizes must preserve the canonical active face in production solves."}};
}
}

Evaluation EvaluatePartitioned(const Domain & domain,const Vector & y,const Vector & eta,
    const ComponentPartition & partition,const EvaluationContext & context)
{
    if(IdentityPartition(partition)) return Evaluate(domain,y,eta,&context);
    Vector beta=Vector::Zero(2*eta.size());
    for(const auto & view:partition.components)
    {
        const auto local_y=Select(y,view.rows),local_eta=Select(eta,view.atoms);
        const auto child=ComponentContext(context,view,false);
        const auto local=runtime::EvaluateProfile(view.domain,local_y,local_eta,&child);
        if(!local.valid) return Evaluation(local);
        for(std::size_t atom=0;atom<view.atoms.size();++atom)
            beta.segment(2*view.atoms[atom],2)=local.beta.segment(2*static_cast<Eigen::Index>(atom),2);
    }
    return Evaluation(runtime::EvaluateState(domain,y,eta,beta,context));
}
Evaluation ComponentEvaluation(const Evaluation & e,const ComponentView & view)
{
    Evaluation out; out.valid=e.valid; out.reason=e.reason;
    const auto columns=Coefficients(view); out.x=Slice(e.x,view,columns); out.derivative=Slice(e.derivative,view,columns);
    out.eta=Select(e.eta,view.atoms); out.beta=Select(e.beta,columns);
    out.residual=Select(e.residual,view.rows); out.gradient=Select(e.gradient,view.atoms); return out;
}
Differential DifferentiatePartitioned(const Evaluation & e,const ComponentPartition & partition,const EvaluationContext & context)
{
    Differential out;
    if(!e.valid) {out.reason="invalid-inner"; return out;}
    if(IdentityPartition(partition)) return Differentiate(e,context.scale,&context);
    std::vector<Evaluation> local; double maximum{}; Eigen::Index columns{};
    for(const auto & view:partition.components)
    {
        local.push_back(ComponentEvaluation(e,view)); const auto free=Free(local.back()); columns+=free.cols();
        const auto spectrum=runtime::ComputeSpectrum(free,context.rank,free.cols(),true);
        if(spectrum.available && spectrum.singular_values.size())
            maximum=std::max(maximum,spectrum.singular_values(0));
    }
    const double threshold=context.rank.Relative(columns)*maximum;
    out.projected=Matrix::Zero(e.x.rows(),e.eta.size()); out.jacobian=out.projected;
    for(std::size_t c=0;c<local.size();++c)
    {
        const auto & view=partition.components[c]; const auto child=ComponentContext(context,view,false);
        const auto d=Differentiate(local[c],context.scale,&child,threshold);
        if(!d.valid) {out.reason=view.id+":"+d.reason; out.projected.resize(0,0); out.jacobian.resize(0,0); return out;}
        for(std::size_t r=0;r<view.rows.size();++r) for(std::size_t a=0;a<view.atoms.size();++a)
        {out.projected(view.rows[r],view.atoms[a])=d.projected(static_cast<Eigen::Index>(r),static_cast<Eigen::Index>(a));
         out.jacobian(view.rows[r],view.atoms[a])=d.jacobian(static_cast<Eigen::Index>(r),static_cast<Eigen::Index>(a));}
    }
    out.valid=true; out.reason="global-threshold-block-derivative"; return out;
}
j::object SameState(const Domain & domain,const Vector & y,const Vector & eta,const Vector & beta,
    const ComponentPartition & partition,const EvaluationContext & context)
{
    j::object out{{"passed",true},{"full_equivalence",false},{"context",ContextEvidence(context)}};
    const bool identity=IdentityPartition(partition); out["identity_backend"]=identity;
    bool passed=true;
    const auto raw=AtState(domain,y,eta,beta,context);
    if(!raw.valid) {out["raw"]=State(raw); out["reason"]="raw-state-unavailable"; return out;}
    Vector prediction=Vector::Zero(y.size()); double loss{}; double basis_error{},basis_derivative_error{};
    for(const auto & view:partition.components)
    {
        if(view.rows.empty()) continue;
        const auto child=ComponentContext(context,view,false);
        const auto e=AtState(view.domain,Select(y,view.rows),Select(eta,view.atoms),Select(beta,Coefficients(view)),child);
        if(!e.valid) {out["passed"]=false; out["reason"]="block-raw-unavailable"; return out;}
        const Vector p=e.x*e.beta; for(std::size_t r=0;r<view.rows.size();++r) prediction(view.rows[r])=p(static_cast<Eigen::Index>(r));
        loss+=e.residual.squaredNorm()/2;
        basis_error=std::max(basis_error,(e.x-Slice(raw.x,view,Coefficients(view))).norm());
        basis_derivative_error=std::max(basis_derivative_error,(e.derivative-Slice(raw.derivative,view,Coefficients(view))).norm());
    }
    double constant{}; for(auto r:partition.constant_rows) constant+=y(r)*y(r)/2; loss+=constant;
    const Vector original=raw.x*raw.beta; const double pred_error=(prediction-original).lpNorm<Eigen::Infinity>();
    const double objective_error=std::abs(loss-raw.residual.squaredNorm()/2)/context.scale/context.scale;
    const bool raw_pass=basis_error==0 && basis_derivative_error==0 &&
        ((prediction-original).array().abs()<=2e-12+2e-13*original.array().abs()).all() && objective_error<=1e-12;
    out["raw"]=j::object{{"available",true},{"passed",raw_pass},{"basis_difference",basis_error},
        {"basis_derivative_difference",basis_derivative_error},{"prediction_inf_difference",pred_error},
        {"rss",raw.residual.squaredNorm()},{"partitioned_rss",2*loss},{"constant_objective",constant},
        {"normalized_objective_difference",objective_error}}; passed &= raw_pass;
    j::array solves;
    {
        const auto a=Evaluate(domain,y,eta,&context),b=identity ? a : EvaluatePartitioned(domain,y,eta,partition,context);
        j::object pair{{"solver","production"},{"monolithic",State(a)},{"partitioned",State(b)},
            {"available",a.valid && b.valid},{"validity_agrees",a.valid==b.valid}};
        passed &= a.valid==b.valid;
        if(a.valid && b.valid)
        {
            const bool face=a.certificate.at("active_atoms")==b.certificate.at("active_atoms");
            const double coefficients=Scaled(a.beta,b.beta),objective=std::abs(a.residual.squaredNorm()-b.residual.squaredNorm())/2/context.scale/context.scale;
            const bool gradient=((a.gradient-b.gradient).array().abs()<=1e-13+2e-9*a.gradient.array().abs()).all();
            pair["same_face"]=face; pair["scaled_coefficient_difference"]=coefficients;
            pair["normalized_objective_difference"]=objective; pair["gradient_passed"]=gradient;
            pair["fixed_face_available"]=face;
            if(!face) pair["limitation"]="cross-active-face; fixed-face equivalence unavailable";
            pair["passed"]=coefficients<=1e-10 && objective<=1e-12 && gradient; passed &= pair.at("passed").as_bool();
            if(face)
            {
                // Isolate differentiation at exactly the same profiled A/C/B.
                // Independently solved coefficients are compared above, never
                // substituted silently into this fixed-state comparison.
                const auto da=Differentiate(a,context.scale,&context),db=identity ? da : DifferentiatePartitioned(a,partition,context);
                j::object derivative{{"available",da.valid && db.valid},{"monolithic_reason",da.reason},{"partitioned_reason",db.reason}};
                bool stable=true;
                if(domain.atoms.size()<=12)
                {
                    const auto evidence=FaceEvidence(domain,y,a,context); stable=evidence.at("stable").as_bool();
                    derivative["face_evidence"]=evidence;
                }
                if(!stable) {derivative["available"]=false; derivative["limitation"]="cross-active-face at registered derivative steps";}
                passed &= da.valid==db.valid;
                if(da.valid && db.valid && stable)
                {
                    const double proj=Relative(da.projected,db.projected),jac=Relative(da.jacobian,db.jacobian);
                    derivative["projected_relative_difference"]=proj; derivative["jacobian_relative_difference"]=jac;
                    derivative["passed"]=proj<=1e-8 && jac<=1e-8; passed &= derivative.at("passed").as_bool();
                    const auto ca=LocalCorrection(a,da,context),cb=identity ? ca : LocalCorrection(a,db,context);
                    derivative["local_correction_available"]=ca.size()>0 && cb.size()>0;
                    if(ca.size()>0 && cb.size()>0) derivative["local_correction_scaled_difference"]=Scaled(ca,cb);
                    out["full_equivalence"]=ca.size()>0 && cb.size()>0;
                    if(ca.size()==0 || cb.size()==0) derivative["limitation"]="local profile correction unavailable";
                }
                out["derivative"]=derivative;
            }
        }
        solves.push_back(pair);
    }
    for(const auto & pair:solves) if(pair.as_object().if_contains("same_face") && !pair.at("same_face").as_bool()) out["full_equivalence"]=false;
    out["profile_solves"]=solves; out["passed"]=passed;
    out["full_equivalence"]=passed && out.at("full_equivalence").as_bool(); return out;
}
runtime::FixedNeighborComponentResult SolveFixedNeighborComponentForTesting(const ComponentView & view,
    const Vector & initial_b,const EvaluationContext & parent)
{
    rhbm_gem::JointParameterLayout layout;
    layout.full_atoms.assign(view.atoms.begin(),view.atoms.end());
    layout.informative_rows.assign(view.rows.begin(),view.rows.end());
    const auto prepared=runtime::MakeComponentProblemView(layout,view);
    auto bounded_parent=parent;
    if(!bounded_parent.eta_upper_bound)
    {
        double maximum_initial_width{};
        for(const double width:initial_b)
            if(std::isfinite(width) && width>maximum_initial_width) maximum_initial_width=width;
        if(!(maximum_initial_width>0.0)) maximum_initial_width=1.0;
        bounded_parent.eta_upper_bound=std::log(4.0*maximum_initial_width);
    }
    runtime::FixedNeighborSearchPolicy policy;
    policy.eta_upper_bound=bounded_parent.eta_upper_bound;
    return runtime::SolveFixedNeighborComponent(prepared,initial_b,bounded_parent,policy);
}
FixedNeighborFitDiagnostic SolveFixedNeighborComponentForDiagnostics(const ComponentView & view,
    const Vector & initial_b,const EvaluationContext & parent)
{
    const auto fixed=SolveFixedNeighborComponentForTesting(view,initial_b,parent);
    FixedNeighborFitDiagnostic out;
    const auto & result=fixed.search;
    const bool endpoint_qualified=result.certification && result.certification->qualified();
    out.search.initial.eta.resize(static_cast<Eigen::Index>(view.atoms.size()));
    out.search.initial.beta=Vector::Zero(2*static_cast<Eigen::Index>(view.atoms.size()));
    for(std::size_t k=0;k<view.atoms.size();++k)
        out.search.initial.eta(static_cast<Eigen::Index>(k))=std::log(initial_b(view.atoms[k]));
    out.search.initial.valid=true;
    out.search.eta=result.state.eta;
    out.candidate_eta=result.state.eta;
    out.candidate_beta=result.state.beta;
    out.candidate_objective=result.state.objective;
    out.candidate_available=out.candidate_eta.size()==static_cast<Eigen::Index>(view.atoms.size()) &&
        out.candidate_beta.size()==2*static_cast<Eigen::Index>(view.atoms.size()) &&
        out.candidate_eta.allFinite() && out.candidate_beta.allFinite() && std::isfinite(out.candidate_objective);
    out.search.stopped=!result.search_converged || !endpoint_qualified;
    out.search.stop_reason=fixed.stop_reason;
    out.search.evaluations=static_cast<int>(result.total_profile_evaluations);
    out.search.seconds=result.search_seconds;
    out.search.accepted=static_cast<int>(result.total_accepted_local_updates);
    out.search.accepted_objective=runtime::LastFixedNeighborObjective(result);
    out.search.accepted_gradient_inf_norm=runtime::LastFixedNeighborGradient(result);
    if(result.certification) out.assessment=result.certification->assessment;
    if(endpoint_qualified)
    {
        out.endpoint_trust=result.certification->trust;
        out.trusted_state=result.certification->assessment.primary;
        out.trusted_assessment=result.certification->assessment;
    }
    out.search_success=fixed.succeeded;
    return out;
}
j::object FitComponent(const ComponentView & view,const Vector & initial_b,const EvaluationContext & parent)
{
    const auto context=ComponentContext(parent,view,true); const Vector start=Select(initial_b,view.atoms);
    const auto result=SolveFixedNeighborComponentForDiagnostics(view,initial_b,parent);
    auto fit=runtime_json::Search(result.search,context,view.domain.rows);
    auto assessment=runtime_json::Assessment(result.trusted_assessment ? *result.trusted_assessment : result.assessment);
    for(auto & field:assessment) fit[field.key()]=std::move(field.value());
    struct rusage usage{}; getrusage(RUSAGE_SELF,&usage);
#ifdef __APPLE__
    const auto bytes=usage.ru_maxrss;
#else
    const auto bytes=usage.ru_maxrss*1024;
#endif
    fit["resources"]=j::object{{"search",j::object{{"seconds",result.search.seconds},{"process_peak_rss_bytes",bytes}}},
        {"endpoint_audit",j::object{{"seconds",result.assessment_seconds},{"process_peak_rss_bytes",bytes}}}};
    fit["seconds"]=result.search.seconds+result.assessment_seconds;
    fit["component_id"]=view.id; fit["parent_atoms"]=Indices(view.atoms); fit["parent_rows"]=Indices(view.rows);
    fit["context"]=ContextEvidence(context); fit["initial_b"]=Values(start);
    fit["usable_state"]=result.trusted_state.has_value(); fit["last_trusted_state"]=nullptr;
    fit["candidate_available"]=result.candidate_available;
    if(result.candidate_available)
    {
        const auto & primary=result.assessment.primary;
        const double objective=primary.certificate.available ? primary.certificate.objective : result.candidate_objective;
        j::object candidate{{"beta",Values(result.candidate_beta)},
            {"eta",Values(result.candidate_eta)},
            {"b",Values(result.candidate_eta.array().exp().matrix())},
            {"objective",objective},{"valid",primary.valid}};
        if(primary.certificate.evaluated)
        {
            candidate["feasible"]=primary.certificate.feasible;
            candidate["kkt_passed"]=primary.certificate.kkt_passed;
            candidate["projected_kkt"]=runtime_json::Number(primary.certificate.projected_kkt);
            candidate["rss"]=runtime_json::Number(primary.certificate.rss);
            candidate["b_gradient"]=Values(primary.gradient);
        }
        fit["candidate_state"]=std::move(candidate);
    }
    else fit["candidate_state"]=nullptr;
    if(result.trusted_trial) fit["last_trusted_state"]=runtime_json::Trial(result.search.trials[*result.trusted_trial]);
    else if(result.trusted_state) fit["last_trusted_state"]=runtime_json::Endpoint(*result.trusted_state);
    if(result.endpoint_trust) fit["endpoint_trust"]=runtime_json::Trust(*result.endpoint_trust);
    fit["search_endpoint_eta"]=Values(result.search.eta);
    if(!result.trusted_state) {fit["runtime_convergence"]="unavailable"; fit["runtime_failure"]="missing-trusted-state";}
    fit["search_success"]=result.search_success; return fit;
}
j::object Assemble(const Domain & domain,const Vector & y,const ComponentPartition & partition,
    const EvaluationContext & context,const j::array & fits)
{
    const auto start=std::chrono::steady_clock::now();
    std::vector<rhbm_gem::core::JointComponentResult> children; int evaluations{},updates{}; j::array failed;
    for(const auto & view:partition.components)
    {
        const auto found=std::find_if(fits.begin(),fits.end(),[&](const auto & f){return f.at("component_id").as_string()==view.id;});
        if(found==fits.end()) throw std::invalid_argument("Missing component scientific record.");
        const auto & fit=*found; evaluations+=j::value_to<int>(fit.at("profile_evaluations"));
        updates+=j::value_to<int>(fit.at("accepted_updates"));
        rhbm_gem::core::JointComponentResult child; child.id=view.id;
        for(const auto atom:view.atoms) child.atoms.push_back(static_cast<std::size_t>(atom));
        for(const auto row:view.rows) child.rows.push_back(static_cast<std::size_t>(row));
        child.search_completed=fit.at("search_success").as_bool();
        child.stop_reason=j::value_to<std::string>(fit.at("stop_reason"));
        child.profile_evaluations=j::value_to<int>(fit.at("profile_evaluations"));
        child.accepted_updates=j::value_to<int>(fit.at("accepted_updates"));
        if(!child.search_completed) failed.push_back(j::object{{"id",view.id},{"reason",fit.at("stop_reason")}});
        if(fit.at("usable_state").as_bool())
        {
            const auto & endpoint=fit.at("last_trusted_state");
            const auto beta=Parse(endpoint.at("beta")),b=Parse(endpoint.at("b")),eta=Parse(endpoint.at("eta"));
            rhbm_gem::JointState state;
            state.ac.assign(beta.data(),beta.data()+beta.size());
            state.b.assign(b.data(),b.data()+b.size());
            state.log_b.assign(eta.data(),eta.data()+eta.size());
            const auto gradient=Parse(endpoint.at("b_gradient"));
            state.width_gradient.assign(gradient.data(),gradient.data()+gradient.size());
            child.state=std::move(state);
            const auto status=runtime_json::ParseStatus(fit.at("runtime_convergence"));
            for(const char * name:{"ac-feasibility","conditional-ac-optimality","independent-ac-optimality",
                "independent-kkt-consistency","independent-replay","width-gradient-consistency",
                "width-stationarity","endpoint-verification","numerical-qualification"})
                child.evidence.push_back({name,status,rhbm_gem::core::JointEvidenceScope::ComponentLocal,{},{},{}});
            child.evidence.push_back({"complete-sweep-stationarity",child.search_completed ?
                rhbm_gem::core::JointCheckStatus::Passed : rhbm_gem::core::JointCheckStatus::Unavailable,
                rhbm_gem::core::JointEvidenceScope::ComponentLocal,{},{},{}});
        }
        children.push_back(std::move(child));
    }
    const auto assembly=runtime::AssembleComponents(domain,y,partition,context,children);
    j::array mask; for(bool value:assembly.row_mask) mask.push_back(value);
    j::object out{{"schema_version",2},{"experiment","exact-component-equivalence"},{"variant","guarded-components"},
        {"execution_complete",true},{"search_stopped_without_convergence",!assembly.completed},{"prediction_available",assembly.available},
        {"objective_available",assembly.available},{"available_row_mask",mask},{"failed_components",failed},
        {"profile_evaluations",evaluations},{"accepted_updates",updates},
        {"runtime_convergence","unavailable"},{"context",ContextEvidence(context)},{"row_count",domain.rows},{"residual_scale",context.scale}};
    if(!assembly.available)
    {out["runtime_failure"]="missing-trusted-component-state"; out["prediction"]=nullptr; out["objective"]=nullptr; return out;}
    auto assessment=runtime_json::Assessment(assembly.assessment);
    for(auto & field:assessment) out[field.key()]=std::move(field.value());
    out["assembled_state"]=j::object{{"eta",Values(assembly.eta)},{"b",Values(assembly.eta.array().exp())},{"beta",Values(assembly.beta)},
        {"certificate",runtime_json::Certificate(assembly.raw.certificate)},{"valid",assembly.raw.valid}};
    out["assembled_state_preserved"]=true;
    out["assembled_rss"]=assembly.raw.valid ? j::value(assembly.raw.certificate.rss) : j::value(nullptr);
    out["assembly_mapping_valid"]=assembly.mapping_valid;
    out["assembled_consistency_passed"]=assembly.consistency_passed;
    out["independent_prediction_difference"]=runtime_json::Number(assembly.replay.prediction_difference);
    out["independent_objective_difference"]=runtime_json::Number(assembly.replay.objective_difference);
    const auto offline_profile=runtime::EvaluateProfile(domain,y,assembly.eta,&context);
    const bool profile_agrees=offline_profile.valid &&
        ((assembly.beta-offline_profile.beta).array().abs()<=1e-10*(1+assembly.beta.array().abs().max(offline_profile.beta.array().abs()))).all();
    out["offline_profile_control"]=j::object{{"valid",offline_profile.valid},{"reason",offline_profile.reason},
        {"certificate",runtime_json::Certificate(offline_profile.certificate)}};
    out["offline_profile_control_evaluations"]=1;
    out["assembled_profile_difference"]=offline_profile.valid ?
        j::value(((assembly.beta-offline_profile.beta).array().abs()/(1+assembly.beta.array().abs().max(offline_profile.beta.array().abs()))).maxCoeff()) : j::value(nullptr);
    out["assembled_profile_agrees"]=profile_agrees;
    auto evidence=runtime::AssessmentEvidence(assembly.assessment,rhbm_gem::core::JointEvidenceScope::AssembledGlobal);
    evidence.push_back({"assembled-consistency",assembly.consistency_passed ? rhbm_gem::core::JointCheckStatus::Passed : rhbm_gem::core::JointCheckStatus::Failed,
        rhbm_gem::core::JointEvidenceScope::AssembledGlobal,{},{},{}});
    auto convergence=runtime::ConvergenceStatus(evidence,rhbm_gem::core::JointEvidenceScope::AssembledGlobal,true);
    for(const auto & fit:fits) convergence=runtime::MergeConvergenceStatus(convergence,runtime_json::ParseStatus(fit.at("runtime_convergence")));
    out["runtime_convergence"]=runtime_json::Status(convergence);
    out["assembly_seconds"]=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count(); return out;
}
j::object FitComponents(std::shared_ptr<const rhbm_gem::core::JointProblemInput> input,
    const Vector & initial_b,const ComponentPartition & partition,const EvaluationContext & context)
{
    const Domain domain(input);
    const runtime::VectorMap y(input->observations.data(),static_cast<Eigen::Index>(input->observations.size()));
    j::array fits; for(const auto & view:partition.components) fits.push_back(FitComponent(view,initial_b,context));
    auto result=Assemble(domain,y,partition,context,fits); result["components"]=std::move(fits); return result;
}
} // namespace second_stage_test::matched::joint_abc
