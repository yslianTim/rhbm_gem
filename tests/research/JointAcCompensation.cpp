#include "core/estimator/InitialStage.hpp"
#include "core/estimator/joint_component/LocalConditioning.hpp"
#include "core/estimator/joint_component/Problem.hpp"
#include "core/estimator/joint_component/WidthCalibration.hpp"
#include "data/io/detail/JointResultJson.hpp"
#include <rhbm_gem/data/io/ModelMapFileIO.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <boost/json.hpp>
#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>
#include <stdexcept>

namespace {
using namespace rhbm_gem;
using namespace rhbm_gem::core;
namespace n=joint_component;
namespace json=boost::json;

json::array Values(n::VectorRef values)
{
    json::array out;
    for(double value:values) out.push_back(value);
    return out;
}

json::object Diagnose(const n::ProblemData & data,const JointFitResult & fit,const std::string & id,
    const n::EvaluationContext & context)
{
    const auto found=std::find(data.input->atom_ids.begin(),data.input->atom_ids.end(),id);
    if(found==data.input->atom_ids.end()) return {{"reason","unknown-fitting-target"}};
    const auto atom=static_cast<std::size_t>(found-data.input->atom_ids.begin());
    const auto & endpoint=*fit.assembled_state;
    for(const auto & component:fit.components) for(const auto & block:component.local_conditioning)
    {
        const auto member=std::find(block.core_atoms.begin(),block.core_atoms.end(),atom);
        if(member==block.core_atoms.end()) continue;
        const auto target=static_cast<Eigen::Index>(member-block.core_atoms.begin());
        n::Indices core(block.core_atoms.begin(),block.core_atoms.end()),rows;
        n::Vector eta(static_cast<Eigen::Index>(core.size())),beta(2*eta.size());
        json::array core_ids;
        for(Eigen::Index j=0;j<eta.size();++j)
        {
            const auto index=static_cast<std::size_t>(core[static_cast<std::size_t>(j)]);
            core_ids.push_back(json::value(data.input->atom_ids[index]));
            const auto p=std::find(data.layout.full_atoms.begin(),data.layout.full_atoms.end(),index);
            if(p==data.layout.full_atoms.end()) return {{"reason","missing-parameter-mapping"}};
            const auto parameter=static_cast<std::size_t>(p-data.layout.full_atoms.begin());
            eta(j)=endpoint.log_tau.at(parameter);
            beta(2*j)=endpoint.ac.at(2*parameter); beta(2*j+1)=endpoint.ac.at(2*parameter+1);
            for(const auto & support:data.domain.atoms[index])
                if(std::binary_search(data.layout.informative_rows.begin(),data.layout.informative_rows.end(),
                    static_cast<std::size_t>(support.row))) rows.push_back(support.row);
        }
        std::sort(rows.begin(),rows.end());
        rows.erase(std::unique(rows.begin(),rows.end()),rows.end());
        if(rows.size()!=block.affected_row_count) return {{"reason","block-row-count-mismatch"}};
        const auto domain=data.domain.Select(core,static_cast<Eigen::Index>(rows.size()),
            n::CompactRowMapping::FromParentRows(rows));
        n::Vector residual(domain.rows);
        for(Eigen::Index j=0;j<domain.rows;++j)
        {
            const auto row=rows[static_cast<std::size_t>(j)];
            residual(j)=fit.prediction->at(static_cast<std::size_t>(row))-data.y(row);
        }
        const auto core_prediction=n::EvaluateStateStreaming(domain,n::Vector::Zero(domain.rows),eta,beta,context);
        const n::Vector effective_response=core_prediction.residual-residual;
        const auto started=std::chrono::steady_clock::now();
        const auto d=n::AssessAcCompensation(domain,effective_response,eta,beta,context,target,
            data.construction_memory.resident_bytes);
        const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
        json::object out{
            {"component",component.id},{"block",block.block_id},{"core_ids",std::move(core_ids)},{"rows",domain.rows},
            {"target_support",domain.atoms[static_cast<std::size_t>(target)].size()},{"A",beta(2*target)},
            {"C",beta(2*target+1)},{"width",std::exp(eta(target))},{"width_upper_bound",d.width_upper_bound_active},
            {"active_face",d.conditional_on_active_face},{"amplitude_active",d.amplitude_active},{"seconds",seconds},
            {"temporary_bytes",d.memory_estimate.temporary_bytes},{"reason",d.reason},
            {"raw",json::value_from(d.raw_correlation)},{"raw_reason",d.raw_reason},
            {"conditional",json::value_from(d.conditional_correlation)},{"conditional_reason",d.conditional_reason},
            {"g_norm",json::value_from(d.gaussian_norm)},{"k_norm",json::value_from(d.charge_norm)},
            {"g_residual_norm",json::value_from(d.residualized_gaussian_norm)},
            {"k_residual_norm",json::value_from(d.residualized_charge_norm)},
            {"cancellation",json::value_from(d.cancellation_ratio)},{"cancellation_reason",d.cancellation_reason},
            {"free_rank",json::value_from(d.free_rank)},{"nuisance_rank",json::value_from(d.nuisance_rank)},
            {"free_rank_threshold_sensitive",d.free_rank_threshold_sensitive},
            {"nuisance_rank_threshold_sensitive",d.nuisance_rank_threshold_sensitive},
            {"free_ac_condition",json::value_from(block.free_ac_condition)},
            {"projected_width_rank",json::value_from(block.projected_width_rank)},
            {"projected_width_condition",json::value_from(block.projected_width_condition)},
            {"profile_rank",json::value_from(block.profile_rank)},{"profile_condition",json::value_from(block.profile_condition)},
            {"profile_status",block.profile_status},{"singular_values",Values(d.singular_values)},
            {"weak_dimension",d.weak_dimension},{"weak_atom_shares",Values(d.weak_atom_shares)},
            {"weak_alignment",json::value_from(d.weak_cancellation_alignment)},{"weak_reason",d.weak_reason}
        };
        const auto & points=fit.target_evidence->point_atoms;
        const auto point=std::find_if(points.begin(),points.end(),[&](const auto & p){return p.atom==atom;});
        if(point!=points.end()) out["point_quality"]=point->point_quality;
        if(d.physical_weak_direction.size())
        {
            double step=0.01; // ||D delta_beta|| before feasibility clipping.
            for(Eigen::Index j=0;j<beta.size();j+=2)
                if(d.physical_weak_direction(j)!=0)
                    step=std::min(step,0.25*beta(j)/std::abs(d.physical_weak_direction(j)));
            const n::Vector delta=step*d.physical_weak_direction;
            const auto variation=n::EvaluateStateStreaming(domain,n::Vector::Zero(domain.rows),eta,delta,context);
            const auto & prediction=variation.residual;
            const double delta_objective=(residual.dot(prediction)+.5*prediction.squaredNorm())/context.scale/context.scale;
            out["perturbation"]={{"normalized_step",step},{"delta_A",delta(2*target)},{"delta_C",delta(2*target+1)},
                {"parameter_norm",delta.stableNorm()},{"prediction_norm",prediction.stableNorm()},
                {"delta_objective",delta_objective}};
        }
        return out;
    }
    return {{"reason","missing-conditioning-block"}};
}
}

int main(int argc,char ** argv)
try
{
    if(argc<8)
    {
        std::cerr<<"Usage: joint_ac_compensation MODEL MAP CHAIN FIRST_RESIDUE LAST_RESIDUE OUTPUT ATOM_ID...\n";
        return 2;
    }
    auto model=ReadModel(argv[1]);
    auto map=ReadMap(argv[2]);
    const double divisor=map->MapValueArrayNormalization().value_or(1);
    const int first=std::stoi(argv[4]),last=std::stoi(argv[5]);
    model->SelectAtoms([&](const AtomObject & atom) {
        return atom.GetElement()!=Element::HYDROGEN && atom.GetChainID()==argv[3] &&
            atom.GetSequenceID()>=first && atom.GetSequenceID()<=last;
    });
    const auto problem=BuildPureHaloSupportExclusionProblem(*map,*model);
    model->EditAnalysis().InitializeFromSelection();
    const auto workset=detail::BuildInitialStageWorkset(*model,problem);
    FitOptions options; options.quiet_mode=true; options.joint_initialization_thread_size=1;
    const auto init=detail::RunInitialStage(*map,*model,workset,options);
    const auto fit=FitJointComponentsInitialized(problem,init);
    auto context=JointProblemAccess::Get(problem).context;
    const auto calibration=n::CalibrateGlobalWidth(init);
    context.width_bounds.upper=calibration.eta_max4;
    json::object checks,targets;
    for(const auto & check:fit.evidence)
        checks[check.name]={{"status",joint_result_io::StatusText(check.status)},
            {"value",json::value_from(check.value)},{"reason",check.reason}};
    const bool qualified=fit.TargetRuntimeConvergence()==JointCheckStatus::Passed && fit.assembled_state &&
        fit.prediction && fit.target_evidence;
    if(qualified)
        for(int i=7;i<argc;++i) targets[argv[i]]=Diagnose(JointProblemAccess::Get(problem),fit,argv[i],context);
    std::ofstream output(argv[6]);
    output<<json::serialize(json::object{
        {"target_count",init.atoms.size()},{"qualified",qualified},{"map_divisor",divisor},
        {"tau_max",json::value_from(calibration.tau_max4)},{"scale",fit.observation_scale},
        {"objective",json::value_from(fit.objective)},{"checks",std::move(checks)},{"targets",std::move(targets)},
        {"row_ids",json::value_from(problem.Input().row_ids)},{"prediction",json::value_from(fit.prediction)}
    })<<'\n';
    if(!output) throw std::runtime_error("Could not write research output.");
    return qualified ? 0 : 1;
}
catch(const std::exception & error)
{
    std::cerr<<error.what()<<'\n';
    return 2;
}
