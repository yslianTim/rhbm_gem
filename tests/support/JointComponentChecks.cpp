#include "support/JointComponentChecks.hpp"
#include "support/JointRuntimeJson.hpp"
#include <sys/resource.h>
#include <algorithm>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <cmath>
#include <chrono>

namespace second_stage_test::matched::joint_abc {
namespace {
namespace j=boost::json;
std::size_t Index(Eigen::Index i) {return static_cast<std::size_t>(i);}
j::array Indices(const std::vector<Eigen::Index> & v)
{j::array out; for(auto i:v) out.push_back(i); return out;}
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
double Maximum(const j::object & spectrum)
{return spectrum.at("available").as_bool() && !spectrum.at("singular_values").as_array().empty() ? j::value_to<double>(spectrum.at("singular_values").at(0)) : 0;}
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
    Matrix directions=context.audit.directions;
    if(!directions.size())
    {
        directions.resize(e.eta.size(),2); directions.col(0)=Vector::Ones(e.eta.size()).normalized();
        for(Eigen::Index k=0;k<e.eta.size();++k) directions(k,1)=k%2 ? -1 : 1;
        directions.col(1).normalize();
    }
    bool stable=true; j::array probes,saved;
    for(Eigen::Index k=0;k<directions.cols();++k)
    {
        j::array direction; for(double v:directions.col(k)) direction.push_back(v); saved.push_back(direction);
        for(double h:{1e-4,5e-5}) for(int sign:{-1,1})
        {
            const auto p=Evaluate(domain,y,e.eta+sign*h*directions.col(k),false,&context);
            const auto r=Evaluate(domain,y,e.eta+sign*h*directions.col(k),true,&context);
            const bool same=p.valid && r.valid && p.certificate.at("active_atoms")==e.certificate.at("active_atoms") &&
                r.certificate.at("active_atoms")==e.certificate.at("active_atoms");
            stable &= same;
            probes.push_back(j::object{{"direction",k},{"step",sign*h},{"same_face",same},
                {"primary",State(p)},{"reference",State(r)}});
        }
    }
    return {{"stable",stable},{"directions",saved},{"probes",probes},
        {"rule","Both declared legacy step sizes must preserve the canonical active face in primary and reference solves."}};
}
// Every spectrum is judged at the same current global threshold. Individual
// block ranks are deliberately ignored, including at repeated singular values.
j::object CompareSpectra(const j::object & global,const std::vector<j::object> & blocks)
{
    if(!global.at("available").as_bool()) return {{"available",false},{"reason","empty-global-matrix"}};
    std::vector<double> merged; j::array records;
    for(const auto & block:blocks)
    {
        if(!block.at("available").as_bool()) return {{"available",false},{"reason","empty-block-matrix"}};
        records.push_back(block);
        for(const auto & v:block.at("singular_values").as_array()) merged.push_back(j::value_to<double>(v));
    }
    const auto & expected=global.at("singular_values").as_array();
    std::sort(merged.begin(),merged.end(),std::greater<double>()); merged.resize(expected.size(),0);
    const double maximum=Maximum(global),threshold=j::value_to<double>(global.at("rank_threshold"));
    double difference{}; std::size_t rank{}; j::array values;
    for(std::size_t k=0;k<merged.size();++k)
    {difference=std::max(difference,std::abs(merged[k]-j::value_to<double>(expected[k]))); rank+=merged[k]>threshold; values.push_back(merged[k]);}
    const double error=maximum>0 ? difference/maximum : difference;
    const bool same_rank=rank==j::value_to<std::size_t>(global.at("rank"));
    return {{"available",true},{"passed",error<=1e-10 && same_rank},{"global",global},{"blocks",records},
        {"merged_singular_values",values},{"merged_rank",rank},{"shared_absolute_threshold",threshold},
        {"normalized_spectrum_difference",error},{"threshold_disagreement",!same_rank}};
}
}

Evaluation EvaluatePartitioned(const Domain & domain,const Vector & y,const Vector & eta,
    const ComponentPartition & partition,const EvaluationContext & context,bool reference)
{
    if(IdentityPartition(partition)) return Evaluate(domain,y,eta,reference,&context);
    std::vector<runtime::LinearBlock> blocks;
    for(const auto & view:partition.components) blocks.push_back({view.rows,Coefficients(view)});
    return Evaluate(domain,y,eta,reference,&context,&blocks);
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
        maximum=std::max(maximum,Maximum(MatrixSpectrum(free,context.rank,free.cols(),true)));
    }
    const double threshold=context.rank.Absolute(columns,maximum);
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
    std::vector<j::object> design_blocks; const auto design_spectrum=MatrixSpectrum(raw.x,context.rank,raw.x.cols(),true);
    for(const auto & view:partition.components)
    {
        if(view.rows.empty())
        {design_blocks.push_back(MatrixSpectrum(Slice(raw.x,view,Coefficients(view)),context.rank,raw.x.cols(),true)); continue;}
        const auto child=ComponentContext(context,view,false);
        const auto e=AtState(view.domain,Select(y,view.rows),Select(eta,view.atoms),Select(beta,Coefficients(view)),child);
        if(!e.valid) {out["passed"]=false; out["reason"]="block-raw-unavailable"; return out;}
        const Vector p=e.x*e.beta; for(std::size_t r=0;r<view.rows.size();++r) prediction(view.rows[r])=p(static_cast<Eigen::Index>(r));
        loss+=e.residual.squaredNorm()/2;
        basis_error=std::max(basis_error,(e.x-Slice(raw.x,view,Coefficients(view))).norm());
        basis_derivative_error=std::max(basis_derivative_error,(e.derivative-Slice(raw.derivative,view,Coefficients(view))).norm());
        design_blocks.push_back(identity ? design_spectrum : MatrixSpectrum(e.x,context.rank,raw.x.cols(),true));
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
    j::object spectra; spectra["column_normalized_design"]=CompareSpectra(design_spectrum,design_blocks);
    passed &= spectra.at("column_normalized_design").at("passed").as_bool();
    j::array solves;
    for(bool reference:{false,true})
    {
        const auto a=Evaluate(domain,y,eta,reference,&context),b=identity ? a : EvaluatePartitioned(domain,y,eta,partition,context,reference);
        j::object pair{{"reference",reference},{"monolithic",State(a)},{"partitioned",State(b)},
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
            if(!reference && face)
            {
                // Isolate differentiation at exactly the same profiled A/C/B.
                // Independently solved coefficients are compared above, never
                // substituted silently into this fixed-state comparison.
                const auto da=Differentiate(a,context.scale,&context),db=identity ? da : DifferentiatePartitioned(a,partition,context);
                j::object derivative{{"available",da.valid && db.valid},{"monolithic_reason",da.reason},{"partitioned_reason",db.reason}};
                bool stable=true;
                if(context.audit.expanded_if_unverified)
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
                    std::vector<j::object> frees,projected,normalized,jacobians;
                    const auto free=Free(a);
                    const auto free_spectrum=MatrixSpectrum(free,context.rank,free.cols(),true);
                    const auto projected_spectrum=MatrixSpectrum(da.projected,context.rank,eta.size(),false);
                    const auto normalized_spectrum=MatrixSpectrum(da.projected,context.rank,eta.size(),true);
                    const auto jacobian_spectrum=MatrixSpectrum(da.jacobian,context.rank,eta.size(),false);
                    for(const auto & view:partition.components)
                    {
                        if(identity) {frees.push_back(free_spectrum); projected.push_back(projected_spectrum);
                            normalized.push_back(normalized_spectrum); jacobians.push_back(jacobian_spectrum); continue;}
                        const auto e=ComponentEvaluation(a,view);
                        Matrix p(static_cast<Eigen::Index>(view.rows.size()),static_cast<Eigen::Index>(view.atoms.size())),jacobian=p;
                        for(std::size_t r=0;r<view.rows.size();++r) for(std::size_t k=0;k<view.atoms.size();++k)
                        {p(static_cast<Eigen::Index>(r),static_cast<Eigen::Index>(k))=db.projected(view.rows[r],view.atoms[k]);
                         jacobian(static_cast<Eigen::Index>(r),static_cast<Eigen::Index>(k))=db.jacobian(view.rows[r],view.atoms[k]);}
                        frees.push_back(MatrixSpectrum(Free(e),context.rank,free.cols(),true));
                        projected.push_back(MatrixSpectrum(p,context.rank,eta.size(),false));
                        normalized.push_back(MatrixSpectrum(p,context.rank,eta.size(),true));
                        jacobians.push_back(MatrixSpectrum(jacobian,context.rank,eta.size(),false));
                    }
                    spectra["free_design"]=CompareSpectra(free_spectrum,frees);
                    spectra["projected_width"]=CompareSpectra(projected_spectrum,projected);
                    spectra["normalized_projected_width"]=CompareSpectra(normalized_spectrum,normalized);
                    spectra["profile_jacobian"]=CompareSpectra(jacobian_spectrum,jacobians);
                    const auto ca=LocalCorrection(a,da,context),cb=identity ? ca : LocalCorrection(a,db,context,
                        j::value_to<double>(spectra.at("profile_jacobian").at("shared_absolute_threshold")));
                    derivative["local_correction_available"]=ca.size()>0 && cb.size()>0;
                    if(ca.size()>0 && cb.size()>0) derivative["local_correction_scaled_difference"]=Scaled(ca,cb);
                    out["full_equivalence"]=ca.size()>0 && cb.size()>0;
                    if(ca.size()==0 || cb.size()==0) derivative["limitation"]="rank-deficient profile Jacobian; local correction unavailable";
                }
                out["derivative"]=derivative;
            }
        }
        solves.push_back(pair);
    }
    for(const auto & entry:spectra) passed &= entry.value().at("available").as_bool() && entry.value().at("passed").as_bool();
    for(const auto & pair:solves) if(pair.as_object().if_contains("same_face") && !pair.at("same_face").as_bool()) out["full_equivalence"]=false;
    out["profile_solves"]=solves; out["spectra"]=spectra; out["passed"]=passed;
    out["full_equivalence"]=passed && out.at("full_equivalence").as_bool(); return out;
}
j::object FitComponent(const ComponentView & view,const Vector & y,const Vector & initial_b,const EvaluationContext & parent)
{
    const auto context=ComponentContext(parent,view,true); const Vector start=Select(initial_b,view.atoms);
    const auto result=runtime::SolveComponent(view,y,initial_b,parent);
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
    std::vector<runtime::ComponentResult> children; int evaluations{},references{},updates{}; j::array failed;
    for(const auto & view:partition.components)
    {
        const auto found=std::find_if(fits.begin(),fits.end(),[&](const auto & f){return f.at("component_id").as_string()==view.id;});
        if(found==fits.end()) throw std::invalid_argument("Missing component scientific record.");
        const auto & fit=*found; evaluations+=j::value_to<int>(fit.at("profile_evaluations"));
        references+=j::value_to<int>(fit.at("search_reference_evaluations")); updates+=j::value_to<int>(fit.at("accepted_updates"));
        runtime::ComponentResult child; child.search_success=fit.at("search_success").as_bool();
        if(!child.search_success) failed.push_back(j::object{{"id",view.id},{"reason",fit.at("stop_reason")}});
        if(fit.at("usable_state").as_bool())
        {
            runtime::Endpoint state; state.eta=Parse(fit.at("last_trusted_state").at("eta"));
            state.beta=Parse(fit.at("last_trusted_state").at("beta")); child.trusted_state=std::move(state);
        }
        children.push_back(std::move(child));
    }
    const auto assembly=runtime::AssembleComponents(domain,y,partition,context,children);
    j::array mask; for(bool value:assembly.row_mask) mask.push_back(value);
    j::object out{{"schema_version",2},{"experiment","exact-component-equivalence"},{"variant","guarded-components"},
        {"execution_complete",true},{"search_stopped_without_convergence",!assembly.completed},{"prediction_available",assembly.available},
        {"objective_available",assembly.available},{"available_row_mask",mask},{"failed_components",failed},
        {"profile_evaluations",evaluations},{"search_reference_evaluations",references},{"accepted_updates",updates},
        {"runtime_convergence","unavailable"},{"context",ContextEvidence(context)},{"row_count",domain.rows},{"residual_scale",context.scale}};
    if(!assembly.available)
    {out["runtime_failure"]="missing-trusted-component-state"; out["prediction"]=nullptr; out["objective"]=nullptr; return out;}
    auto assessment=runtime_json::Assessment(assembly.assessment);
    for(auto & field:assessment) out[field.key()]=std::move(field.value());
    out["assembled_state"]=j::object{{"eta",Values(assembly.eta)},{"b",Values(assembly.eta.array().exp())},{"beta",Values(assembly.beta)},
        {"certificate",runtime_json::Certificate(assembly.raw.certificate)},{"valid",assembly.raw.valid}};
    out["assembled_state_preserved"]=true;
    out["assembled_rss"]=assembly.raw.valid ? j::value(assembly.raw.certificate.rss) : j::value(nullptr);
    if(assembly.profile_evaluated)
    {
        out["profile_control"]=j::object{{"valid",assembly.profile_control.valid},{"reason",assembly.profile_control.reason},
            {"certificate",runtime_json::Certificate(assembly.profile_control.certificate)}};
        out["profile_control_evaluations"]=1;
        if(assembly.profile_control.valid) out["assembled_profile_difference"]=assembly.profile_difference;
    }
    out["assembled_profile_agrees"]=assembly.profile_agrees;
    auto evidence=runtime::AssessmentEvidence(assembly.assessment,rhbm_gem::core::JointEvidenceScope::AssembledGlobal);
    evidence.push_back({"assembled-profile",assembly.profile_agrees ? rhbm_gem::core::JointCheckStatus::Passed : rhbm_gem::core::JointCheckStatus::Failed,
        rhbm_gem::core::JointEvidenceScope::AssembledGlobal,{},{},{}});
    auto convergence=runtime::ConvergenceStatus(evidence,rhbm_gem::core::JointEvidenceScope::AssembledGlobal,true);
    for(const auto & fit:fits) convergence=runtime::MergeConvergenceStatus(convergence,runtime_json::ParseStatus(fit.at("runtime_convergence")));
    out["runtime_convergence"]=runtime_json::Status(convergence);
    out["assembly_seconds"]=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count(); return out;
}
j::object FitComponents(const Domain & domain,const Vector & y,const Vector & initial_b,
    const ComponentPartition & partition,const EvaluationContext & context)
{
    j::array fits; for(const auto & view:partition.components) fits.push_back(FitComponent(view,y,initial_b,context));
    auto result=Assemble(domain,y,partition,context,fits); result["components"]=std::move(fits); return result;
}
} // namespace second_stage_test::matched::joint_abc
