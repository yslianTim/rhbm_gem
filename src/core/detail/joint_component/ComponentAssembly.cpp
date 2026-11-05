#include "Numerics.hpp"

#include <algorithm>
#include <chrono>
#include <utility>

namespace rhbm_gem::core::joint_component {
namespace {
using Clock=std::chrono::steady_clock;
double Seconds(Clock::time_point start) {return std::chrono::duration<double>(Clock::now()-start).count();}
}
AssemblyResult AssembleComponents(const Domain & domain,VectorRef y,const ComponentPartition & partition,
    const EvaluationContext & context,const std::vector<JointComponentResult> & components,
    const JointConvergenceCriteria & convergence_criteria)
{
    AssemblyResult out; out.assessment.convergence_criteria=convergence_criteria;
    if(components.size()!=partition.components.size())
    {out.mapping_reason="missing-component-result"; return out;}
    out.available=true;
    out.row_mask.assign(static_cast<std::size_t>(domain.rows),false);
    out.eta=Vector::Constant(static_cast<Eigen::Index>(domain.atoms.size()),unavailable);
    out.beta=Vector::Constant(2*out.eta.size(),unavailable);
    if(domain.rows<0 || y.size()!=domain.rows || !partition.mappings ||
        partition.mappings->atom_component.size()!=domain.atoms.size() ||
        partition.mappings->atom_to_local.size()!=domain.atoms.size() ||
        partition.mappings->row_component.size()!=static_cast<std::size_t>(domain.rows) ||
        partition.mappings->row_to_local.size()!=static_cast<std::size_t>(domain.rows))
    {out.available=false; out.mapping_reason="invalid-partition-mapping"; return out;}
    std::vector<bool> mapped_atoms(domain.atoms.size(),false),mapped_rows(static_cast<std::size_t>(domain.rows),false);
    for(auto r:partition.constant_rows)
    {
        if(r<0 || r>=domain.rows || mapped_rows[static_cast<std::size_t>(r)] ||
            partition.mappings->row_component[static_cast<std::size_t>(r)]!=-1)
        {out.available=false; out.mapping_reason="invalid-constant-row-mapping"; return out;}
        mapped_rows[static_cast<std::size_t>(r)]=true;
        out.row_mask[static_cast<std::size_t>(r)]=true;
    }
    for(std::size_t k=0;k<components.size();++k)
    {
        const auto & component=components[k]; const auto & view=partition.components[k];
        out.completed &= component.search_completed;
        const bool same_atoms=component.atoms.size()==view.atoms.size() &&
            std::equal(component.atoms.begin(),component.atoms.end(),view.atoms.begin(),[](std::size_t a,Eigen::Index b) {
                return b>=0 && a==static_cast<std::size_t>(b);
            });
        const bool same_rows=component.rows.size()==view.rows.size() &&
            std::equal(component.rows.begin(),component.rows.end(),view.rows.begin(),[](std::size_t a,Eigen::Index b) {
                return b>=0 && a==static_cast<std::size_t>(b);
            });
        if(component.id!=view.id || view.component_index!=static_cast<Eigen::Index>(k) || !same_atoms || !same_rows)
        {
            out.available=false; out.mapping_reason="component-identity-mapping-mismatch"; return out;
        }
        for(std::size_t atom_index=0;atom_index<view.atoms.size();++atom_index)
        {
            const auto atom=view.atoms[atom_index];
            if(atom<0 || atom>=static_cast<Eigen::Index>(domain.atoms.size()) ||
                partition.mappings->atom_component[static_cast<std::size_t>(atom)]!=static_cast<Eigen::Index>(k) ||
                partition.mappings->atom_to_local[static_cast<std::size_t>(atom)]!=static_cast<Eigen::Index>(atom_index) ||
                mapped_atoms[static_cast<std::size_t>(atom)])
            {out.available=false; out.mapping_reason="component-atom-mapping-mismatch"; return out;}
            mapped_atoms[static_cast<std::size_t>(atom)]=true;
            for(const auto & support:domain.atoms[static_cast<std::size_t>(atom)])
                if(support.row<0 || support.row>=domain.rows ||
                    partition.mappings->row_component[static_cast<std::size_t>(support.row)]!=static_cast<Eigen::Index>(k))
                {out.available=false; out.mapping_reason="component-support-row-mapping-mismatch"; return out;}
        }
        for(std::size_t row_index=0;row_index<view.rows.size();++row_index)
        {
            const auto r=view.rows[row_index];
            if(r<0 || r>=domain.rows || mapped_rows[static_cast<std::size_t>(r)] ||
                partition.mappings->row_component[static_cast<std::size_t>(r)]!=static_cast<Eigen::Index>(k) ||
                partition.mappings->row_to_local[static_cast<std::size_t>(r)]!=static_cast<Eigen::Index>(row_index))
            {out.available=false; out.mapping_reason="component-row-mapping-mismatch"; return out;}
            mapped_rows[static_cast<std::size_t>(r)]=true;
        }
    }
    for(std::size_t atom=0;atom<domain.atoms.size();++atom)
        if(!mapped_atoms[atom] || partition.mappings->atom_component[atom]<0)
        {out.available=false; out.mapping_reason="unmapped-component-atom"; return out;}
    if(!std::all_of(mapped_rows.begin(),mapped_rows.end(),[](bool mapped){return mapped;}))
    {out.available=false; out.mapping_reason="unmapped-observation-row"; return out;}
    for(const auto & support:domain.atoms)
        for(const auto & entry:support)
            if(entry.row<0 || entry.row>=domain.rows || !mapped_rows[static_cast<std::size_t>(entry.row)] ||
                std::binary_search(partition.constant_rows.begin(),partition.constant_rows.end(),entry.row))
            {out.available=false; out.mapping_reason="constant-row-mapping-mismatch"; return out;}
    out.mapping_valid=true;
    for(std::size_t k=0;k<components.size();++k)
    {
        const auto & component=components[k]; const auto & view=partition.components[k];
        if(!component.state)
        {
            out.available=false;
            if(out.mapping_reason.empty()) out.mapping_reason="missing-component-state";
            continue;
        }
        if(component.RuntimeConvergence()!=JointCheckStatus::Passed)
        {
            out.available=false;
            if(out.mapping_reason.empty()) out.mapping_reason="component-not-converged";
            continue;
        }
        const auto & state=*component.state;
        if(state.ac.size()!=2*view.atoms.size() || state.tau.size()!=view.atoms.size() ||
            state.log_tau.size()!=view.atoms.size() || state.eta_gradient.size()!=view.atoms.size())
        {
            out.available=false;
            if(out.mapping_reason.empty()) out.mapping_reason="invalid-component-state-dimensions";
            continue;
        }
        for(std::size_t a=0;a<view.atoms.size();++a)
        {
            out.eta(view.atoms[a])=state.log_tau[a];
            out.beta(2*view.atoms[a])=state.ac[2*a];
            out.beta(2*view.atoms[a]+1)=state.ac[2*a+1];
        }
        for(const auto row:view.rows) out.row_mask[static_cast<std::size_t>(row)]=true;
    }
    if(!out.available) return out;
    out.row_mask=std::move(mapped_rows);
    if(!out.eta.allFinite() || !out.beta.allFinite())
    {out.available=false; out.mapping_reason="invalid-assembled-parameter-state"; return out;}
    out.verification_memory_estimate=EndpointVerificationMemoryEstimate(domain,y,context);
    if(out.verification_memory_estimate.status!=MemoryStatus::Pass)
    {
        LogMemoryEstimate(out.verification_memory_estimate);
        out.assessment.failure="resource-limit:endpoint-verification:"+out.verification_memory_estimate.reason;
        out.assessment.verification_memory_estimate=out.verification_memory_estimate;
        out.assessment.primary.eta=out.eta; out.assessment.primary.reason=out.assessment.failure;
        out.replay.reason=out.assessment.failure;
        return out;
    }
    if(out.verification_memory_estimate.estimated_peak_bytes>out.verification_memory_estimate.budget_bytes/2)
        LogMemoryEstimate(out.verification_memory_estimate);
    const auto primary_started=Clock::now();
    const auto raw=EvaluateStateStreaming(domain,y,out.eta,out.beta,context);
    const double primary_seconds=Seconds(primary_started);
    out.raw=raw;
    const auto replay_started=Clock::now();
    out.replay=CheckIndependentReplay(domain,y,raw,context,nullptr,nullptr,convergence_criteria);
    const double replay_seconds=Seconds(replay_started);
    out.consistency_passed=out.replay.passed;
    out.assessment=VerifyJointEndpoint(domain,y,raw,out.replay,context,convergence_criteria);
    out.assessment.verification_memory_estimate=out.verification_memory_estimate;
    out.assessment.primary_endpoint_seconds=primary_seconds;
    out.assessment.independent_replay_seconds=replay_seconds;
    if(raw.valid) {out.prediction=raw.residual+y; out.objective=raw.certificate.objective;}
    return out;
}
}
