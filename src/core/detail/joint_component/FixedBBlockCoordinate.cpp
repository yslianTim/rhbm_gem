#include "FixedBBlockCoordinate.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numeric>

namespace rhbm_gem::core::joint_component {
namespace {
using Clock=std::chrono::steady_clock;
double Seconds(Clock::time_point start) {return std::chrono::duration<double>(Clock::now()-start).count();}
struct ReplayState {Vector prediction,residual; double objective{};};
ReplayState Replay(const JointProblemInput & input,const JointParameterLayout & layout,VectorRef observations,
    VectorRef eta,VectorRef beta,double scale)
{
    ReplayState out; out.prediction=Vector::Zero(observations.size()); out.residual=Vector::Zero(observations.size());
    std::vector<bool> informative(static_cast<std::size_t>(observations.size()),false);
    for(auto row:layout.informative_rows) informative.at(row)=true;
    for(auto atom:layout.full_atoms)
    {
        const auto a=static_cast<Eigen::Index>(atom); const double width=std::exp(eta(a));
        for(const auto & support:input.support.at(atom)) if(informative.at(support.row))
        {
            const auto basis=EvaluateKernel(support.squared_distance,width,2.5);
            out.prediction(static_cast<Eigen::Index>(support.row))+=beta(2*a)*basis.gaussian+beta(2*a+1)*basis.charge;
        }
    }
    double squared{};
    for(auto row:layout.informative_rows)
    {
        const auto r=static_cast<Eigen::Index>(row); out.residual(r)=out.prediction(r)-observations(r);
        squared+=out.residual(r)*out.residual(r);
    }
    out.objective=.5*squared/(scale*scale); return out;
}
struct BlockDesign {Sparse design; Indices rows;};
BlockDesign BuildBlockDesign(const JointProblemInput & input,const StructuralCore & core,
    VectorRef eta,const JointParameterLayout & layout)
{
    BlockDesign out; out.rows=core.affected_rows;
    std::vector<Eigen::Index> row_to_local(input.observations.size(),-1);
    for(std::size_t k=0;k<out.rows.size();++k) row_to_local.at(static_cast<std::size_t>(out.rows[k]))=static_cast<Eigen::Index>(k);
    out.design.resize(static_cast<Eigen::Index>(out.rows.size()),2*static_cast<Eigen::Index>(core.atoms.size()));
    std::vector<std::uint8_t> informative(input.observations.size());
    for(auto row:layout.informative_rows) informative.at(row)=1;
    std::vector<Eigen::Triplet<double>> entries;
    for(std::size_t k=0;k<core.atoms.size();++k)
    {
        const auto atom=static_cast<std::size_t>(core.atoms[k]); const auto local=static_cast<Eigen::Index>(k);
        const double width=std::exp(eta(static_cast<Eigen::Index>(atom)));
        for(const auto & support:input.support.at(atom))
        {
            if(!informative[support.row]) continue;
            const auto row=row_to_local.at(support.row);
            if(row<0) throw std::logic_error("block-invalid-partition");
            const auto basis=EvaluateKernel(support.squared_distance,width,2.5);
            if(basis.gaussian!=0) entries.emplace_back(row,2*local,basis.gaussian);
            if(basis.charge!=0) entries.emplace_back(row,2*local+1,basis.charge);
        }
    }
    out.design.setFromTriplets(entries.begin(),entries.end());
    return out;
}
struct GlobalDiagnostics {double feasibility{},kkt{};};
GlobalDiagnostics DiagnoseGlobalAC(const JointProblemInput & input,const JointParameterLayout & layout,
    const BlockCoordinateState & state,double scale)
{
    const auto columns=2*static_cast<Eigen::Index>(input.atom_ids.size());
    Vector norms=Vector::Zero(columns),gradient=Vector::Zero(columns);
    std::vector<std::uint8_t> informative(input.observations.size());
    for(auto row:layout.informative_rows) informative.at(row)=1;
    for(auto atom:layout.full_atoms)
    {
        const auto a=static_cast<Eigen::Index>(atom); const double width=std::exp(state.eta(a));
        for(const auto & support:input.support.at(atom)) if(informative[support.row])
        {
            const auto row=static_cast<Eigen::Index>(support.row); const auto basis=EvaluateKernel(support.squared_distance,width,2.5);
            for(Eigen::Index kind=0;kind<2;++kind)
            {
                const auto column=2*a+kind; const double value=kind==0 ? basis.gaussian : basis.charge;
                norms(column)+=value*value; gradient(column)+=value*state.residual(row);
            }
        }
    }
    norms=norms.cwiseSqrt(); GlobalDiagnostics out; double feasibility{};
    for(auto atom:layout.full_atoms)
    {
        const auto a=static_cast<Eigen::Index>(atom); feasibility=std::max(feasibility,std::max(0.0,-state.beta(2*a)));
        for(Eigen::Index kind=0;kind<2;++kind)
        {
            const auto column=2*a+kind; if(!(norms(column)>0)) {out.kkt=std::numeric_limits<double>::infinity(); continue;}
            const double scaled=norms(column)*state.beta(column)/scale;
            double projected=scaled-gradient(column)/norms(column)/scale;
            if(kind==0) projected=std::max(0.0,projected);
            out.kkt=std::max(out.kkt,std::abs(scaled-projected));
        }
    }
    out.feasibility=feasibility; return out;
}
bool WithinPredictionReplay(double error,VectorRef reference)
{return error<=2e-12+2e-13*std::max(1.0,reference.cwiseAbs().maxCoeff());}
bool WithinObjectiveReplay(double error,double reference)
{return error<=1e-12+2e-12*std::abs(reference);}
}
FixedBBlockResult SearchFixedBBlocks(const JointProblemInput & input,const JointParameterLayout & layout,VectorRef observations,
    VectorRef eta,const EvaluationContext & context,const FixedBBlockPolicy & policy)
{
    FixedBBlockResult out; out.reason="block-search-failed";
    if(observations.size()!=static_cast<Eigen::Index>(input.observations.size()) || eta.size()!=static_cast<Eigen::Index>(input.atom_ids.size()) ||
        !observations.allFinite() || !eta.allFinite() || !(context.scale>0) || !std::isfinite(context.scale) ||
        context.scale!=std::max(1.0,observations.norm()) || policy.core_atoms==0)
    {out.reason="block-invalid-partition"; return out;}
    StructuralBlockPartition partition;
    try {partition=BuildStructuralBlockPartition(input,layout,policy.core_atoms);}
    catch(const std::exception &) {out.reason="block-invalid-partition"; return out;}
    if(partition.cores.empty()) {out.reason="block-invalid-partition"; return out;}
    out.state.eta=eta; out.state.beta=Vector::Zero(2*static_cast<Eigen::Index>(input.atom_ids.size()));
    const auto initial=Replay(input,layout,observations,out.state.eta,out.state.beta,context.scale);
    out.state.prediction=initial.prediction; out.state.residual=initial.residual; out.state.objective=initial.objective;
    std::vector<std::size_t> order(partition.cores.size());
    std::iota(order.begin(),order.end(),0);
    if(policy.order==FixedBBlockOrder::Reverse) std::reverse(order.begin(),order.end());
    for(std::size_t sweep_index=0;sweep_index<policy.maximum_sweeps;++sweep_index)
    {
        const auto sweep_started=Clock::now(); FixedBBlockSweep sweep;
        sweep.objective_before=out.state.objective;
        for(auto block_index:order)
        {
            const auto & core=partition.cores[block_index]; FixedBBlockRecord record;
            record.sweep=sweep_index+1; record.block=block_index+1; record.core_atoms=core.atoms;
            record.affected_rows=core.affected_rows; record.neighbor_atoms=core.neighbor_atoms;
            record.diagnostics_captured=policy.capture_diagnostics;
            if(policy.capture_diagnostics)
            {
                record.global_kkt_before=DiagnoseGlobalAC(input,layout,out.state,context.scale).kkt;
                record.global_kkt_after=record.global_kkt_before;
            }
            record.objective_before=out.state.objective;
            const auto local=BuildBlockDesign(input,core,out.state.eta,layout);
            Vector local_y(static_cast<Eigen::Index>(local.rows.size())),old_beta(2*static_cast<Eigen::Index>(core.atoms.size()));
            for(std::size_t a=0;a<core.atoms.size();++a)
            {
                const auto atom=static_cast<Eigen::Index>(core.atoms[a]);
                old_beta.segment<2>(2*static_cast<Eigen::Index>(a))=out.state.beta.segment<2>(2*atom);
            }
            const auto old_core=(local.design*old_beta).eval();
            for(std::size_t k=0;k<local.rows.size();++k)
            {
                const auto row=local.rows[k];
                local_y(static_cast<Eigen::Index>(k))=old_core(static_cast<Eigen::Index>(k))-out.state.residual(row);
            }
            if(!local_y.allFinite())
            {record.status="failed"; record.reason="block-invalid-effective-response"; ++sweep.failed_blocks; out.blocks.push_back(std::move(record)); out.reason="block-invalid-effective-response"; return out;}
            const auto local_before=.5*(local.design*old_beta-local_y).squaredNorm()/(context.scale*context.scale);
            if(policy.capture_diagnostics) record.local_objective_before=local_before;
            const auto factor_started=Clock::now();
            const auto solved=SolveLinear(local.design,local_y,Vector::Ones(local_y.size()),false,true,nullptr,&context.linear);
            record.factor_seconds=Seconds(factor_started); record.local_iterations=solved.solves;
            record.factorizations=solved.block_factorizations;
            if(policy.capture_diagnostics)
            {
                record.linear_solves=solved.solves;
                record.rank=solved.rank;
                if(solved.beta.size()==old_beta.size())
                {
                    record.solver_beta_equal_old=(solved.beta.array()==old_beta.array()).all();
                    record.solver_beta_norm_difference=(solved.beta-old_beta).norm();
                    record.raw_ac_change=record.solver_beta_norm_difference;
                }
            }
            record.active_A=0; record.free_columns=0;
            if(!solved.valid)
            {
                record.status="failed"; record.reason=solved.reason=="rank-deficient" ? "block-inner-rank-failed" : "block-inner-invalid";
                ++sweep.failed_blocks; out.reason=record.reason; out.blocks.push_back(std::move(record)); return out;
            }
            const auto certificate=CertifyLinear(local.design,local_y,solved.beta,context.scale);
            if(policy.capture_diagnostics)
            {
                record.local_kkt=certificate.projected_kkt;
                record.local_feasibility=certificate.feasible ? 0.0 : 1.0;
                record.active_atoms=certificate.active_atoms.size();
                record.local_certificate_available=certificate.available;
                record.local_certificate_feasible=certificate.feasible;
                record.local_certificate_kkt_passed=certificate.kkt_passed;
                record.local_certificate_projected_kkt=certificate.projected_kkt;
                record.local_certificate_rss=certificate.rss;
                record.local_certificate_objective=certificate.objective;
            }
            if(!certificate.available || !certificate.feasible || !certificate.kkt_passed)
            {record.status="failed"; record.reason="block-inner-invalid"; ++sweep.failed_blocks; out.blocks.push_back(std::move(record)); out.reason="block-inner-invalid"; return out;}
            for(Eigen::Index k=0;k<solved.beta.size();++k)
            {
                if(k%2==0) record.active_A+=solved.beta(k)==0;
                if(k%2!=0 || solved.beta(k)>0) ++record.free_columns;
            }
            const auto local_after=.5*(local.design*solved.beta-local_y).squaredNorm()/(context.scale*context.scale);
            if(policy.capture_diagnostics)
            {
                record.local_objective_after=local_after;
                record.local_objective_reduction=local_before-local_after;
            }
            Vector candidate_beta=out.state.beta;
            for(std::size_t k=0;k<core.atoms.size();++k)
                candidate_beta.segment<2>(2*static_cast<Eigen::Index>(core.atoms[k]))=solved.beta.segment<2>(2*static_cast<Eigen::Index>(k));
            const auto replay=Replay(input,layout,observations,out.state.eta,candidate_beta,context.scale);
            record.objective_after=replay.objective; record.objective_reduction=record.objective_before-record.objective_after;
            if(policy.capture_diagnostics)
            {
                record.global_replay_delta=replay.objective-out.state.objective;
                record.objective_replay_enclosure=1e-12+2e-12*std::max(std::abs(replay.objective),std::abs(out.state.objective));
            }
            double local_global_error=std::abs((replay.objective-out.state.objective)-(local_after-local_before));
            if(policy.capture_diagnostics) record.local_global_delta_error=local_global_error;
            const double replay_reference=std::max(std::abs(replay.objective),std::abs(out.state.objective));
            const bool local_global_delta_enclosed=WithinObjectiveReplay(local_global_error,replay_reference);
            const double replay_enclosure=1e-12+2e-12*replay_reference;
            if(!local_global_delta_enclosed)
            {record.status="failed"; record.reason="block-objective-replay-failed"; ++sweep.failed_blocks; out.blocks.push_back(std::move(record)); out.reason="block-objective-replay-failed"; return out;}
            const double global_replay_delta=replay.objective-out.state.objective;
            if(global_replay_delta>0 && !(local_after<=local_before && local_global_delta_enclosed &&
                global_replay_delta<=replay_enclosure))
            {
                record.status="unchanged"; record.reason="block-objective-increase"; ++sweep.unchanged_blocks;
                sweep.objective_after=out.state.objective; out.blocks.push_back(std::move(record)); continue;
            }
            const auto new_core=(local.design*solved.beta).eval();
            double old_squared{},new_squared{};
            for(std::size_t k=0;k<local.rows.size();++k)
            {
                const auto row=local.rows[k]; const double before=out.state.residual(row);
                const double change=new_core(static_cast<Eigen::Index>(k))-old_core(static_cast<Eigen::Index>(k));
                out.state.prediction(row)+=change; out.state.residual(row)+=change;
                old_squared+=before*before; new_squared+=out.state.residual(row)*out.state.residual(row);
            }
            out.state.objective+=(new_squared-old_squared)/(2*context.scale*context.scale);
            out.state.beta=std::move(candidate_beta);
            for(Eigen::Index k=0;k<solved.beta.size();++k)
                record.max_scaled_ac_change=std::max(record.max_scaled_ac_change,
                    std::abs(solved.beta(k)-old_beta(k))*local.design.col(k).norm()/context.scale);
            sweep.max_scaled_ac_change=std::max(sweep.max_scaled_ac_change,record.max_scaled_ac_change);
            const bool unchanged=(solved.beta.array()==old_beta.array()).all();
            record.status=unchanged ? "unchanged" : "accepted"; record.reason=unchanged ? "block-unchanged" : "";
            if(policy.capture_diagnostics)
            {
                record.accepted=!unchanged;
                if(!unchanged) record.global_kkt_after=DiagnoseGlobalAC(input,layout,out.state,context.scale).kkt;
            }
            if(unchanged) ++sweep.unchanged_blocks; else ++sweep.accepted_blocks;
            sweep.objective_after=out.state.objective; out.blocks.push_back(std::move(record));
        }
        sweep.objective_after=out.state.objective;
        sweep.relative_objective_change=std::abs(sweep.objective_after-sweep.objective_before)/std::max(1.0,std::abs(sweep.objective_before));
        const auto replay=Replay(input,layout,observations,out.state.eta,out.state.beta,context.scale);
        sweep.cache_replay_error=(out.state.prediction-replay.prediction).lpNorm<Eigen::Infinity>();
        sweep.objective_replay_error=std::abs(out.state.objective-replay.objective);
        if(!WithinPredictionReplay(sweep.cache_replay_error,replay.prediction) ||
            !WithinObjectiveReplay(sweep.objective_replay_error,replay.objective))
        {sweep.failed_blocks=1; out.reason="block-cache-replay-failed";}
        const auto diagnostics=DiagnoseGlobalAC(input,layout,out.state,context.scale);
        sweep.global_a_feasibility=diagnostics.feasibility; sweep.global_ac_kkt=diagnostics.kkt;
        sweep.wall_seconds=Seconds(sweep_started); out.sweeps.push_back(sweep);
        if(sweep.failed_blocks) return out;
        if(sweep.global_ac_kkt<=1e-10)
        {out.success=true; out.reason="block-converged"; out.sweeps_to_global_kkt=sweep_index+1; return out;}
    }
    out.reason="block-sweep-budget"; return out;
}
}
