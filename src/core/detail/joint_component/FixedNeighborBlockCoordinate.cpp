#include "FixedNeighborBlockCoordinate.hpp"
#include "Problem.hpp"
#include "ResourceWork.hpp"
#include "StructuralPartition.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <numeric>

namespace rhbm_gem::core::joint_component {
double BlockObjectiveReplayEnclosure(double reference)
{return 1e-12+2e-12*std::abs(reference);}
bool WithinBlockObjectiveReplay(double error,double reference)
{return error<=BlockObjectiveReplayEnclosure(reference);}

namespace {
using Clock=std::chrono::steady_clock;
constexpr double EtaChangeConfirmationThreshold=1e-10;
constexpr std::size_t EndpointCertificationRetryInterval=10;
double Seconds(Clock::time_point start) {return std::chrono::duration<double>(Clock::now()-start).count();}
struct OptionalPhaseTimer
{
    double * target{};
    Clock::time_point started{Clock::now()};
    explicit OptionalPhaseTimer(double * target):target(target) {}
    ~OptionalPhaseTimer() {if(target) *target+=Seconds(started);}
};
OptionalPhaseTimer Measure(bool enabled,double & seconds)
{return OptionalPhaseTimer{enabled ? &seconds : nullptr};}
Vector Select(VectorRef values,const Indices & indices)
{
    Vector out(static_cast<Eigen::Index>(indices.size()));
    for(std::size_t k=0;k<indices.size();++k) out(static_cast<Eigen::Index>(k))=values(indices[k]);
    return out;
}
Indices IndicesOf(const std::vector<std::size_t> & values)
{return {values.begin(),values.end()};}
bool AddMemoryBytes(std::size_t & total,std::size_t count,std::size_t width)
{
    std::size_t bytes{}; return CheckedMultiplyBytes(count,width,bytes) && CheckedAddBytes(total,bytes,total);
}
bool PartitionStorageBytes(const StructuralBlockPartition & partition,bool include_metadata,std::size_t & bytes)
{
    bytes=0;
    if(!AddMemoryBytes(bytes,partition.cores.capacity(),sizeof(StructuralCore))) return false;
    for(const auto & core:partition.cores)
    {
        std::size_t id_bytes{};
        if(!CheckedAddBytes(core.id.capacity(),1,id_bytes) || !CheckedAddBytes(id_bytes,32,id_bytes) ||
            !CheckedAddBytes(bytes,id_bytes,bytes) || !AddMemoryBytes(bytes,core.atoms.capacity(),sizeof(Eigen::Index)) ||
            !AddMemoryBytes(bytes,1,4*32)) return false;
        if(include_metadata && (!AddMemoryBytes(bytes,core.affected_rows.capacity(),sizeof(Eigen::Index)) ||
            !AddMemoryBytes(bytes,core.neighbor_atoms.capacity(),sizeof(Eigen::Index)) ||
            !AddMemoryBytes(bytes,core.context_atoms.capacity(),sizeof(Eigen::Index)))) return false;
    }
    return true;
}
Vector SelectBeta(VectorRef beta,const std::vector<std::size_t> & atoms)
{
    Vector out(2*static_cast<Eigen::Index>(atoms.size()));
    for(std::size_t k=0;k<atoms.size();++k)
        out.segment<2>(2*static_cast<Eigen::Index>(k))=beta.segment<2>(2*static_cast<Eigen::Index>(atoms[k]));
    return out;
}
struct ReplayState {Vector prediction,residual; double objective{};};
ReplayState Replay(const Domain & domain,const JointParameterLayout & layout,VectorRef observations,
    VectorRef eta,VectorRef beta,double scale)
{
    ReplayState out; out.prediction=Vector::Zero(observations.size()); out.residual=Vector::Zero(observations.size());
    std::vector<std::uint8_t> informative(static_cast<std::size_t>(domain.rows));
    for(auto row:layout.informative_rows) informative.at(row)=1;
    for(auto atom:layout.full_atoms)
    {
        const auto a=static_cast<Eigen::Index>(atom); const double width=std::exp(eta(a));
        for(const auto & support:domain.atoms[atom]) if(informative.at(static_cast<std::size_t>(support.row)))
        {
            const auto basis=EvaluateKernel(support.square,width,2.5);
            out.prediction(support.row)+=beta(2*a)*basis.gaussian+beta(2*a+1)*basis.charge;
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
Sparse Design(const Domain & domain,VectorRef eta)
{
    Sparse out(domain.rows,2*eta.size()); std::vector<Eigen::Triplet<double>> entries;
    for(Eigen::Index a=0;a<eta.size();++a)
    {
        const double width=std::exp(eta(a));
        for(const auto & support:domain.atoms[static_cast<std::size_t>(a)])
        {
            const auto basis=EvaluateKernel(support.square,width,2.5);
            if(basis.gaussian!=0) entries.emplace_back(support.row,2*a,basis.gaussian);
            if(basis.charge!=0) entries.emplace_back(support.row,2*a+1,basis.charge);
        }
    }
    out.setFromTriplets(entries.begin(),entries.end()); return out;
}
double PredictionDifference(VectorRef a,VectorRef b)
{return (a-b).lpNorm<Eigen::Infinity>();}
double ScaledCoefficientDifference(VectorRef lhs,VectorRef rhs)
{
    if(lhs.size()!=rhs.size() || lhs.size()==0) return std::numeric_limits<double>::infinity();
    return ((lhs-rhs).array().abs()/(1.0+lhs.array().abs().max(rhs.array().abs()))).maxCoeff();
}
}
bool IsFixedNeighborEtaChangeConfirmed(double eta_change_inf,bool has_previous_complete_sweep)
{
    return has_previous_complete_sweep && std::isfinite(eta_change_inf) && eta_change_inf>=0.0 &&
        eta_change_inf<=EtaChangeConfirmationThreshold;
}
FixedNeighborStationarityState ClassifyFixedNeighborStationarity(
    double ac_kkt,double width_gradient_inf,double eta_change_inf,bool has_previous_complete_sweep)
{
    if(!(ac_kkt<=1e-10) || !(width_gradient_inf<=1e-12))
        return FixedNeighborStationarityState::NotStationary;
    return IsFixedNeighborEtaChangeConfirmed(eta_change_inf,has_previous_complete_sweep) ?
        FixedNeighborStationarityState::ConfirmedStationary : FixedNeighborStationarityState::CandidateStationary;
}
bool ShouldAttemptFixedNeighborCertification(FixedNeighborStationarityState stationarity,std::size_t sweep,
    std::size_t & next_retry_sweep)
{
    if(stationarity!=FixedNeighborStationarityState::ConfirmedStationary)
    {next_retry_sweep=0; return false;}
    if(next_retry_sweep!=0 && sweep<next_retry_sweep) return false;
    next_retry_sweep=sweep+EndpointCertificationRetryInterval;
    return true;
}
FixedNeighborSearchResult SearchFixedNeighborComponent(
    const Domain & domain,const JointParameterLayout & layout,VectorRef y,VectorRef initial_eta,const EvaluationContext & context,
    const FixedNeighborPolicy & policy,Vector initial_beta)
{
    FixedNeighborSearchResult out;
    const bool has_initial_beta=initial_beta.size()!=0;
    if(initial_eta.size()!=static_cast<Eigen::Index>(domain.atoms.size()) || !initial_eta.allFinite() ||
        y.size()!=domain.rows || !y.allFinite() || !(context.scale>0) || !std::isfinite(context.scale) ||
        (has_initial_beta && (initial_beta.size()!=2*static_cast<Eigen::Index>(domain.atoms.size()) ||
            !initial_beta.allFinite())) ||
        policy.core_atoms==0 || layout.groups.size()!=0 || layout.full_atoms.size()!=domain.atoms.size() ||
        layout.informative_rows.size()!=static_cast<std::size_t>(domain.rows) ||
        context.atom_ids.size()!=domain.atoms.size() || context.row_ids.size()!=static_cast<std::size_t>(domain.rows))
    {out.reason="block-invalid-partition"; return out;}
    for(Eigen::Index a=0;a<initial_eta.size();++a)
        if(!(std::exp(initial_eta(a))>0) || !std::isfinite(std::exp(initial_eta(a))))
        {out.reason="block-invalid-partition"; return out;}

    auto structural_timer=Measure(policy.collect_telemetry,out.work.structural_preparation_seconds);
    MemoryShape memory_shape; memory_shape.component_atoms=domain.atoms.size();
    memory_shape.component_observations=static_cast<std::size_t>(domain.rows);
    memory_shape.component_target_columns=context.target_columns;
    if(!CheckedMultiplyBytes(memory_shape.component_atoms,2,memory_shape.component_free_columns))
    {out.memory_estimate.stage=MemoryStage::StructuralConstruction; out.memory_estimate.status=MemoryStatus::Overflow;
        out.memory_estimate.reason="size-overflow"; out.reason="resource-limit:structural-construction:size-overflow"; return out;}
    for(std::size_t a=0;a<domain.atoms.size();++a)
        if(!CheckedAddBytes(memory_shape.component_support_entries,domain.atoms[a].size(),memory_shape.component_support_entries))
        {out.memory_estimate.stage=MemoryStage::StructuralConstruction; out.memory_estimate.status=MemoryStatus::Overflow;
            out.memory_estimate.reason="size-overflow"; out.reason="resource-limit:structural-construction:size-overflow"; return out;}
    const auto * input=domain.atoms.Storage();
    const auto snapshot_bytes=input ? EstimateSnapshotBytes(*input) : std::optional<std::size_t>{0};
    std::size_t base_resident=snapshot_bytes.value_or(0),root_mappings{},root_entries{},selector_entries{},twice_rows{};
    const auto local_support_bytes=domain.atoms.LocalSupportIndexBytes();
    if(!snapshot_bytes || local_support_bytes==std::numeric_limits<std::size_t>::max() ||
        (input && (!CheckedAddBytes(input->atom_ids.size(),input->observations.size(),root_entries) ||
            !AddMemoryBytes(root_mappings,root_entries,4*sizeof(Eigen::Index)) ||
            !CheckedAddBytes(base_resident,root_mappings,base_resident))) ||
        !AddMemoryBytes(base_resident,local_support_bytes,1) ||
        !CheckedMultiplyBytes(memory_shape.component_observations,2,twice_rows) ||
        !CheckedAddBytes(memory_shape.component_atoms,twice_rows,selector_entries) ||
        !AddMemoryBytes(base_resident,selector_entries,sizeof(Eigen::Index)))
    {out.memory_estimate.shape=memory_shape; out.memory_estimate.stage=MemoryStage::StructuralConstruction;
        out.memory_estimate.status=MemoryStatus::Overflow; out.memory_estimate.reason="size-overflow";
        out.reason="resource-limit:structural-construction:size-overflow"; return out;}
    memory_shape.resident_bytes=base_resident;
    memory_shape.component_blocks=memory_shape.component_atoms;
    out.memory_estimate=PlanMemory(memory_shape,MemoryStage::StructuralConstruction);
    if(out.memory_estimate.status!=MemoryStatus::Pass)
    {LogMemoryEstimate(out.memory_estimate); out.reason="resource-limit:structural-construction:"+out.memory_estimate.reason; return out;}

    StructuralBlockPartition partition;
    StructuralBlockCensus census;
    std::size_t partition_bytes{};
    try
    {
        StructuralTopology topology(domain,context.atom_ids,layout);
        partition=topology.BuildCores(policy.core_atoms,false);
        if(partition.cores.empty()) {out.reason="block-invalid-partition"; return out;}
        census=topology.Census(partition);
        memory_shape.component_blocks=partition.cores.size();
        memory_shape.affected_row_entries=census.affected_row_entries;
        memory_shape.component_neighbor_entries=census.neighbor_atom_entries;
        memory_shape.maximum_block_rows=census.maximum_block_rows;
        memory_shape.maximum_block_support_entries=census.maximum_block_support_entries;
        memory_shape.total_block_support_entries=census.total_block_support_entries;
        memory_shape.maximum_neighbor_atoms=census.maximum_neighbor_atoms;
        std::size_t topology_bytes{},core_bytes{},planned_metadata{},planned_entries{};
        if(!CheckedAddBytes(census.affected_row_entries,census.neighbor_atom_entries,planned_entries) ||
            !CheckedMultiplyBytes(planned_entries,2*sizeof(Eigen::Index),planned_metadata) ||
            !AddMemoryBytes(planned_metadata,partition.cores.size(),4*32) ||
            !PartitionStorageBytes(partition,false,core_bytes) ||
            (topology_bytes=topology.StorageBytes())==std::numeric_limits<std::size_t>::max() ||
            !CheckedAddBytes(base_resident,topology_bytes,memory_shape.resident_bytes) ||
            !CheckedAddBytes(memory_shape.resident_bytes,core_bytes,memory_shape.resident_bytes) ||
            !CheckedAddBytes(memory_shape.resident_bytes,planned_metadata,memory_shape.resident_bytes))
        {out.memory_estimate.shape=memory_shape; out.memory_estimate.stage=MemoryStage::FixedNeighborPreparation;
            out.memory_estimate.status=MemoryStatus::Overflow; out.memory_estimate.reason="size-overflow";
            out.reason="resource-limit:fixed-neighbor-preparation:size-overflow"; return out;}
        const auto structural_estimate=PlanMemory(memory_shape,MemoryStage::FixedNeighborPreparation);
        if(structural_estimate.status!=MemoryStatus::Pass)
        {out.memory_estimate=structural_estimate; LogMemoryEstimate(out.memory_estimate);
            out.reason="resource-limit:fixed-neighbor-preparation:"+out.memory_estimate.reason; return out;}
        topology.PopulateMetadata(partition);
        if(!PartitionStorageBytes(partition,true,partition_bytes))
        {out.memory_estimate.shape=memory_shape; out.memory_estimate.stage=MemoryStage::FixedNeighborPreparation;
            out.memory_estimate.status=MemoryStatus::Overflow; out.memory_estimate.reason="size-overflow";
            out.reason="resource-limit:fixed-neighbor-preparation:size-overflow"; return out;}
    }
    catch(const std::overflow_error &)
    {out.memory_estimate.shape=memory_shape; out.memory_estimate.stage=MemoryStage::FixedNeighborPreparation;
        out.memory_estimate.status=MemoryStatus::Overflow; out.memory_estimate.reason="size-overflow";
        out.reason="resource-limit:fixed-neighbor-preparation:size-overflow"; return out;}
    catch(const std::exception &) {out.reason="block-invalid-partition"; return out;}
    std::size_t resident{};
    if(!CheckedAddBytes(base_resident,partition_bytes,resident))
    {out.memory_estimate.shape=memory_shape; out.memory_estimate.stage=MemoryStage::FixedNeighborPreparation;
        out.memory_estimate.status=MemoryStatus::Overflow; out.memory_estimate.reason="size-overflow";
        out.reason="resource-limit:fixed-neighbor-preparation:size-overflow"; return out;}
    memory_shape.resident_bytes=resident;
    out.memory_estimate=PlanMemory(memory_shape,MemoryStage::FixedNeighborPreparation);
    if(out.memory_estimate.status!=MemoryStatus::Pass ||
        out.memory_estimate.estimated_peak_bytes>out.memory_estimate.budget_bytes/2)
        LogMemoryEstimate(out.memory_estimate);
    if(out.memory_estimate.status!=MemoryStatus::Pass)
    {out.reason="resource-limit:fixed-neighbor-preparation:"+out.memory_estimate.reason; return out;}
    if(structural_timer.target)
    {
        out.work.structural_preparation_seconds+=Seconds(structural_timer.started);
        structural_timer.target=nullptr;
    }
    auto global_shape=memory_shape;
    auto global_resident=resident;
    if(!AddMemoryBytes(global_resident,memory_shape.component_support_entries,sizeof(LocalSupportIndex)) ||
        !AddMemoryBytes(global_resident,memory_shape.affected_row_entries,5*sizeof(Eigen::Index)) ||
        !AddMemoryBytes(global_resident,memory_shape.component_atoms,4*sizeof(Eigen::Index)+3*sizeof(double)+
            sizeof(std::vector<LocalSupportIndex>)) ||
        !AddMemoryBytes(global_resident,memory_shape.component_blocks,
            sizeof(PreparedFixedNeighborBlock)+sizeof(CompactRowMapping)+sizeof(std::size_t)) ||
        !AddMemoryBytes(global_resident,memory_shape.component_observations,3*sizeof(double)+sizeof(Eigen::Index)))
    {
        out.global_state_memory_estimate.shape=memory_shape;
        out.global_state_memory_estimate.stage=MemoryStage::GlobalStateEvaluation;
        out.global_state_memory_estimate.status=MemoryStatus::Overflow;
        out.global_state_memory_estimate.reason="size-overflow";
        out.memory_estimate=out.global_state_memory_estimate;
        out.reason="resource-limit:global-state-evaluation:size-overflow"; return out;
    }
    global_shape.resident_bytes=global_resident;
    out.global_state_memory_estimate=PlanMemory(global_shape,MemoryStage::GlobalStateEvaluation);
    if(out.global_state_memory_estimate.status!=MemoryStatus::Pass ||
        out.global_state_memory_estimate.estimated_peak_bytes>out.global_state_memory_estimate.budget_bytes/2)
        LogMemoryEstimate(out.global_state_memory_estimate);
    if(out.global_state_memory_estimate.status!=MemoryStatus::Pass)
    {
        out.memory_estimate=out.global_state_memory_estimate;
        out.reason="resource-limit:global-state-evaluation:"+out.memory_estimate.reason; return out;
    }
    auto * const profile_work=policy.collect_telemetry ? &out.work.local_profile_work : nullptr;
    const auto certify_endpoint=[&] {
        const auto assessment_timer=Measure(policy.collect_telemetry,out.work.reference_assessment_seconds);
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        if(policy.certification_for_testing)
            return policy.certification_for_testing(domain,y,out.state,context,profile_work);
#endif
        auto certification_resident=out.global_state_memory_estimate.shape.resident_bytes;
        if(!AddMemoryBytes(certification_resident,memory_shape.component_observations,3*sizeof(double)) ||
            !AddMemoryBytes(certification_resident,memory_shape.component_atoms,8*sizeof(double)))
            certification_resident=std::numeric_limits<std::size_t>::max();
        return CertifyFixedNeighborEndpoint(domain,y,out.state,context,profile_work,
            certification_resident);
    };

    out.state.eta=initial_eta;
    out.state.beta=has_initial_beta ? std::move(initial_beta) :
        Vector::Zero(2*static_cast<Eigen::Index>(domain.atoms.size()));
    auto initial=Replay(domain,layout,y,out.state.eta,out.state.beta,context.scale);
    out.state.prediction=std::move(initial.prediction); out.state.residual=std::move(initial.residual);
    out.state.objective=initial.objective;
    auto mapping_timer=Measure(policy.collect_telemetry,out.work.structural_preparation_seconds);
    Indices atom_position(domain.atoms.size(),-1),row_position(static_cast<std::size_t>(domain.rows),-1);
    for(std::size_t k=0;k<layout.full_atoms.size();++k) atom_position.at(layout.full_atoms[k])=static_cast<Eigen::Index>(k);
    for(std::size_t k=0;k<layout.informative_rows.size();++k) row_position.at(layout.informative_rows[k])=static_cast<Eigen::Index>(k);
    std::vector<PreparedFixedNeighborBlock> prepared_blocks;
    prepared_blocks.reserve(partition.cores.size());
    for(std::size_t block_index=0;block_index<partition.cores.size();++block_index)
    {
        const auto & core=partition.cores[block_index]; PreparedFixedNeighborBlock prepared;
        prepared.block=block_index; prepared.core_atoms=core.atoms; prepared.affected_rows=core.affected_rows;
        prepared.profile_atoms.reserve(core.atoms.size()); prepared.profile_rows.reserve(core.affected_rows.size());
        for(auto atom:core.atoms)
        {
            const auto position=atom_position.at(static_cast<std::size_t>(atom));
            if(position<0) {out.reason="block-invalid-partition"; return out;}
            prepared.profile_atoms.push_back(position);
        }
        for(auto row:core.affected_rows)
        {
            const auto position=row_position.at(static_cast<std::size_t>(row));
            if(position<0) {out.reason="block-invalid-partition"; return out;}
            prepared.profile_rows.push_back(position);
        }
        prepared.row_mapping=CompactRowMapping::FromParentRows(prepared.profile_rows);
        prepared.domain=domain.Select(prepared.profile_atoms,
            static_cast<Eigen::Index>(prepared.profile_rows.size()),prepared.row_mapping);
        std::size_t mapping_entries{},mapping_bytes{};
        const auto support_entries=prepared.domain.atoms.LocalSupportIndexEntries();
        const auto support_bytes=prepared.domain.atoms.LocalSupportIndexBytes();
        if(!CheckedAddBytes(prepared.row_mapping->parent_rows.size(),prepared.row_mapping->local_rows.size(),mapping_entries) ||
            !CheckedAddBytes(mapping_entries,support_entries,mapping_entries) ||
            !CheckedAddBytes(prepared.row_mapping->StorageBytes(),support_bytes,mapping_bytes) ||
            !CheckedAddBytes(out.prepared_mapping_entries,mapping_entries,out.prepared_mapping_entries) ||
            !CheckedAddBytes(out.prepared_mapping_bytes,mapping_bytes,out.prepared_mapping_bytes))
        {out.memory_estimate.status=MemoryStatus::Overflow; out.memory_estimate.reason="size-overflow";
            out.reason="resource-limit:fixed-neighbor-mapping:size-overflow"; return out;}
        prepared.context=context;
        prepared.context.atom_ids=context.atom_ids.Select(prepared.profile_atoms);
        prepared.context.row_ids=context.row_ids.Select(prepared.profile_rows);
        prepared.context.independent_search=true;
        prepared.local_atoms=static_cast<Eigen::Index>(prepared.profile_atoms.size());
        prepared.context.rank={prepared.domain.rows,2*prepared.local_atoms,prepared.local_atoms};
        prepared.context.linear.rank_relative=prepared.context.rank.Relative(2*prepared.local_atoms);
        // FixedNeighbor performs at most one trusted accepted local update per block visit.
        prepared.context.update_budget=1;
        prepared_blocks.push_back(std::move(prepared));
    }
    out.block_preparations=prepared_blocks.size();
    out.domain_preparations=prepared_blocks.size();
    out.mapping_preparations=prepared_blocks.size();
    std::vector<std::size_t> order(partition.cores.size()); std::iota(order.begin(),order.end(),0);
    if(policy.order==FixedNeighborBlockOrder::Reverse) std::reverse(order.begin(),order.end());
    if(mapping_timer.target)
    {
        out.work.structural_preparation_seconds+=Seconds(mapping_timer.started);
        mapping_timer.target=nullptr;
    }
    if(policy.preparation_observer) policy.preparation_observer();
    const bool retain_diagnostics=policy.collect_diagnostics || policy.capture_local_trajectory ||
        static_cast<bool>(policy.state_observer);
    bool stop=false;
    std::size_t next_certification_retry_sweep{};
    for(std::size_t sweep_index=0;sweep_index<policy.maximum_sweeps && !stop;++sweep_index)
    {
        const auto sweep_started=Clock::now(); FixedNeighborBlockSweep sweep; sweep.sweep=sweep_index+1;
        const Vector eta_before=out.state.eta,beta_before=out.state.beta;
        sweep.objective_before=out.state.objective;
        auto finish_record=[&](FixedNeighborBlockRecord record,Clock::time_point block_started) {
            record.block_seconds=Seconds(block_started);
            if(policy.collect_telemetry)
            {
                ++out.work.block_visits;
                out.work.maximum_block_seconds=std::max(out.work.maximum_block_seconds,record.block_seconds);
            }
            if(record.accepted)
            {
                ++sweep.accepted_blocks;
                sweep.accepted_local_updates+=static_cast<std::size_t>(record.accepted_updates);
                out.total_accepted_local_updates+=static_cast<std::size_t>(record.accepted_updates);
            }
            else if(record.status=="unchanged") ++sweep.unchanged_blocks;
            if(policy.diagnostic_sink && policy.diagnostic_sink->block) policy.diagnostic_sink->block(record);
            if(retain_diagnostics) out.blocks.push_back(std::move(record));
        };
        for(auto block_index:order)
        {
            const auto & prepared=prepared_blocks.at(block_index);
            const auto & core_atoms=prepared.core_atoms; const auto & affected_rows=prepared.affected_rows;
            const auto & local_domain=prepared.domain; const auto & local_context=prepared.context;
            const auto local_atoms=prepared.local_atoms; const auto block_started=Clock::now(); FixedNeighborBlockRecord record;
            record.sweep=sweep_index+1; record.block=block_index+1; record.core_atoms=core_atoms;
            record.affected_rows=affected_rows.size(); record.objective_before=out.state.objective;

            Vector old_eta(local_atoms),old_beta(2*local_atoms),old_core;
            {
                const auto old_core_timer=Measure(policy.collect_telemetry,out.work.old_core_seconds);
                if(policy.collect_telemetry) ++out.work.old_core_basis_builds;
                for(std::size_t k=0;k<core_atoms.size();++k)
                {
                    const auto atom=static_cast<Eigen::Index>(core_atoms[k]);
                    old_eta(static_cast<Eigen::Index>(k))=out.state.eta(atom);
                    old_beta.segment<2>(2*static_cast<Eigen::Index>(k))=out.state.beta.segment<2>(2*atom);
                }
                const Sparse old_design=Design(local_domain,old_eta);
                old_core=(old_design*old_beta).eval();
                (void)old_core_timer;
            }

            Vector local_y(static_cast<Eigen::Index>(affected_rows.size()));
            {
                const auto response_timer=Measure(policy.collect_telemetry,out.work.effective_response_seconds);
                double local_before_squared{};
                for(std::size_t k=0;k<affected_rows.size();++k)
                {
                    const auto row=static_cast<Eigen::Index>(affected_rows[k]);
                    local_y(static_cast<Eigen::Index>(k))=old_core(static_cast<Eigen::Index>(k))-out.state.residual(row);
                    local_before_squared+=out.state.residual(row)*out.state.residual(row);
                }
                record.local_objective_before=.5*local_before_squared/(context.scale*context.scale);
                (void)response_timer;
            }

            const auto search_started=Clock::now();
            if(policy.collect_telemetry)
            {
                const auto atoms=static_cast<std::size_t>(local_atoms);
                ++out.work.local_problem_count;
                out.work.total_local_atoms+=static_cast<double>(atoms);
                if(out.work.local_problem_count==1)
                    out.work.minimum_local_atoms=out.work.maximum_local_atoms=atoms;
                else
                {
                    out.work.minimum_local_atoms=std::min(out.work.minimum_local_atoms,atoms);
                    out.work.maximum_local_atoms=std::max(out.work.maximum_local_atoms,atoms);
                }
            }
            const auto old_widths=old_eta.array().exp().eval();
            const auto profile_search=SearchProfile(local_domain,local_y,old_widths,local_context,{},nullptr,profile_work);
            if(policy.collect_telemetry) out.work.local_profile_work.Merge(profile_search.profile_work);
            record.search_seconds=Seconds(search_started); record.profile_evaluations=profile_search.evaluations;
            if(policy.collect_telemetry) out.work.local_search_seconds+=record.search_seconds;
            record.accepted_updates=profile_search.accepted; record.local_search_stop_reason=profile_search.stop_reason;
            if(policy.capture_local_trajectory)
            {
                record.profile_trials.reserve(profile_search.trials.size());
                double accepted_objective=record.local_objective_before;
                Vector accepted_eta=old_eta; double cumulative_factor_seconds{};
                for(std::size_t trial_index=0;trial_index<profile_search.trials.size();++trial_index)
                {
                    const auto & trial=profile_search.trials[trial_index];
                    const double objective=trial.endpoint.certificate.evaluated && trial.endpoint.certificate.available ?
                        trial.endpoint.certificate.objective/(context.scale*context.scale) : unavailable;
                    FixedNeighborProfileTrial telemetry;
                    telemetry.trial_index=trial_index; telemetry.profile_evaluation=trial.evaluation;
                    telemetry.accepted=trial.accepted; telemetry.accepted_update=trial.accepted_update;
                    telemetry.local_objective_before=accepted_objective;
                    telemetry.local_objective_after=objective;
                    telemetry.objective_reduction=std::isfinite(accepted_objective) && std::isfinite(objective) ?
                        accepted_objective-objective : unavailable;
                    telemetry.eta_change_inf=trial.endpoint.eta.size()==accepted_eta.size() ?
                        (trial.endpoint.eta-accepted_eta).lpNorm<Eigen::Infinity>() : unavailable;
                    telemetry.gradient_inf_norm=trial.endpoint.gradient.size() ?
                        trial.endpoint.gradient.lpNorm<Eigen::Infinity>() : unavailable;
                    telemetry.profile_seconds=trial.seconds; telemetry.factor_seconds=trial.factor_seconds;
                    cumulative_factor_seconds+=trial.factor_seconds;
                    telemetry.cumulative_factor_seconds=cumulative_factor_seconds;
                    record.profile_trials.push_back(std::move(telemetry));
                    if(trial.accepted && std::isfinite(objective))
                    {
                        accepted_objective=objective;
                        if(trial.endpoint.eta.size()==accepted_eta.size()) accepted_eta=trial.endpoint.eta;
                    }
                }
                record.profile_factor_seconds=cumulative_factor_seconds;
            }
            ++sweep.block_solves; ++out.total_block_solves;
            if(policy.collect_telemetry) ++out.work.block_solves;
            sweep.profile_evaluations+=static_cast<std::size_t>(profile_search.evaluations);
            out.total_profile_evaluations+=static_cast<std::size_t>(profile_search.evaluations);
            sweep.maximum_block_rows=std::max(sweep.maximum_block_rows,static_cast<std::size_t>(local_domain.rows));
            sweep.maximum_block_columns=std::max(sweep.maximum_block_columns,static_cast<std::size_t>(2*local_atoms));
            const Trial * accepted=nullptr;
            for(auto it=profile_search.trials.rbegin();it!=profile_search.trials.rend();++it)
                if(it->accepted && it->endpoint.valid && it->trust && it->trust->passed) {accepted=&*it; break;}
            if(!accepted)
            {
                record.status="failed"; record.reason="block-search-failed"; out.reason=record.reason;
                finish_record(std::move(record),block_started); stop=true; break;
            }
            const auto local_state_started=Clock::now();
            const auto local_state=EvaluateState(local_domain,local_y,accepted->endpoint.eta,accepted->endpoint.beta,local_context,
                ProfileEvaluationRole::AcceptedEndpoint,profile_work);
            if(policy.collect_telemetry) out.work.local_state_seconds+=Seconds(local_state_started);
            if(!local_state.valid || !local_state.certificate.available || !local_state.certificate.feasible ||
                !local_state.certificate.kkt_passed)
            {
                record.status="failed"; record.reason="block-inner-invalid"; out.reason=record.reason;
                finish_record(std::move(record),block_started); stop=true; break;
            }
            record.local_final_gradient_inf_norm=local_state.gradient.lpNorm<Eigen::Infinity>();
            record.local_final_ac_kkt=local_state.certificate.projected_kkt;
            record.local_objective_after=local_state.certificate.objective/(context.scale*context.scale);
            if(record.local_objective_after>record.local_objective_before)
            {
                record.status="unchanged"; record.reason="block-objective-increase";
                record.objective_after=record.objective_before; finish_record(std::move(record),block_started); continue;
            }
            Vector candidate_eta,candidate_beta;
            {
                const auto candidate_copy_timer=Measure(policy.collect_telemetry,out.work.candidate_copy_seconds);
                if(policy.collect_telemetry) ++out.work.candidate_state_full_copies;
                candidate_eta=out.state.eta; candidate_beta=out.state.beta;
                for(std::size_t k=0;k<core_atoms.size();++k)
                {
                    const auto atom=static_cast<Eigen::Index>(core_atoms[k]);
                    candidate_eta(atom)=accepted->endpoint.eta(static_cast<Eigen::Index>(k));
                    candidate_beta.segment<2>(2*atom)=accepted->endpoint.beta.segment<2>(2*static_cast<Eigen::Index>(k));
                }
                (void)candidate_copy_timer;
            }
            ReplayState replay;
            {
                const auto candidate_replay_timer=Measure(policy.collect_telemetry,out.work.candidate_replay_seconds);
                if(policy.collect_telemetry) ++out.work.full_candidate_replays;
                replay=Replay(domain,layout,y,candidate_eta,candidate_beta,context.scale);
                (void)candidate_replay_timer;
            }
            record.objective_after=replay.objective; record.global_replay_delta=replay.objective-out.state.objective;
            const auto replay_reference=std::max(std::abs(replay.objective),std::abs(out.state.objective));
            record.objective_replay_enclosure=BlockObjectiveReplayEnclosure(replay_reference);
            record.local_global_delta_error=std::abs(record.global_replay_delta+
                record.local_objective_before-record.local_objective_after);
            const bool delta_enclosed=WithinBlockObjectiveReplay(record.local_global_delta_error,replay_reference);
            if(!delta_enclosed)
            {
                record.status="failed"; record.reason="block-objective-replay-failed"; out.reason=record.reason;
                finish_record(std::move(record),block_started); stop=true; break;
            }
            if(record.global_replay_delta>0 && !(record.local_objective_after<=record.local_objective_before &&
                record.global_replay_delta<=record.objective_replay_enclosure))
            {
                record.status="unchanged"; record.reason="block-objective-increase";
                finish_record(std::move(record),block_started); continue;
            }
            const bool unchanged=(out.state.eta.array()==candidate_eta.array()).all() &&
                (out.state.beta.array()==candidate_beta.array()).all();
            const auto cache_started=Clock::now();
            const Vector new_core=(local_state.x*local_state.beta).eval(); double old_squared{},new_squared{};
            for(std::size_t k=0;k<affected_rows.size();++k)
            {
                const auto row=static_cast<Eigen::Index>(affected_rows[k]); const double before=out.state.residual(row);
                const double change=new_core(static_cast<Eigen::Index>(k))-old_core(static_cast<Eigen::Index>(k));
                out.state.prediction(row)+=change; out.state.residual(row)+=change;
                old_squared+=before*before; new_squared+=out.state.residual(row)*out.state.residual(row);
                if(policy.collect_telemetry) ++out.work.affected_row_updates;
            }
            out.state.objective+=(new_squared-old_squared)/(2*context.scale*context.scale);
            out.state.eta=std::move(candidate_eta); out.state.beta=std::move(candidate_beta);
            if(policy.collect_telemetry) out.work.cache_update_seconds+=Seconds(cache_started);
            record.accepted=!unchanged; record.status=unchanged ? "unchanged" : "accepted";
            record.reason=unchanged ? "block-unchanged" : "";
            finish_record(std::move(record),block_started);
        }
        if(stop) break;
        sweep.objective_after=out.state.objective;
        ReplayState replay;
        {
            const auto timer=Measure(policy.collect_telemetry,out.work.sweep_replay_seconds);
            replay=Replay(domain,layout,y,out.state.eta,out.state.beta,context.scale);
            (void)timer;
        }
        sweep.cache_replay_error=PredictionDifference(out.state.prediction,replay.prediction);
        sweep.objective_replay_error=std::abs(out.state.objective-replay.objective);
        const Vector eta=Select(out.state.eta,IndicesOf(layout.full_atoms));
        const Vector beta=SelectBeta(out.state.beta,layout.full_atoms);
        const auto global_started=Clock::now();
        const auto global=EvaluateStateStreaming(domain,y,eta,beta,context,ProfileEvaluationRole::AcceptedEndpoint,profile_work);
        if(policy.collect_telemetry) out.work.sweep_global_state_seconds+=Seconds(global_started);
        if(!global.valid || !global.certificate.available)
        {out.reason="block-inner-invalid"; stop=true; break;}
        sweep.global_a_feasibility=global.certificate.feasible ? 0.0 : 1.0;
        sweep.global_ac_kkt=global.certificate.projected_kkt;
        sweep.global_width_gradient_inf_norm=global.gradient.lpNorm<Eigen::Infinity>();
        sweep.eta_change_inf=(out.state.eta-eta_before).lpNorm<Eigen::Infinity>();
        sweep.beta_scaled_change=ScaledCoefficientDifference(out.state.beta,beta_before);
        sweep.coordinate_confirmation_available=sweep_index>0;
        sweep.wall_seconds=Seconds(sweep_started); out.search_seconds+=sweep.wall_seconds;
        ++out.sweep_count; out.final_sweep=sweep;
        if(policy.diagnostic_sink && policy.diagnostic_sink->sweep) policy.diagnostic_sink->sweep(sweep);
        if(policy.production_progress) policy.production_progress(sweep);
        if(retain_diagnostics) out.sweeps.push_back(sweep);
        if(policy.sweep_observer)
        {
            const auto & observed=retain_diagnostics ? out.sweeps.back() : sweep;
            policy.sweep_observer(observed);
        }
        if(policy.state_observer) policy.state_observer(sweep_index+1,out.state,out.sweeps.back(),out.blocks);
        if(sweep.cache_replay_error>2e-12+2e-13*std::max(1.0,replay.prediction.cwiseAbs().maxCoeff()) ||
            !WithinBlockObjectiveReplay(sweep.objective_replay_error,replay.objective))
        {out.reason="block-cache-replay-failed"; stop=true; break;}
        const auto stationarity=ClassifyFixedNeighborStationarity(sweep.global_ac_kkt,
            sweep.global_width_gradient_inf_norm,sweep.eta_change_inf,sweep.coordinate_confirmation_available);
        if(stationarity!=FixedNeighborStationarityState::NotStationary && out.first_order_stationarity_sweep==0)
            out.first_order_stationarity_sweep=sweep_index+1;
        const bool certification_retry_due=ShouldAttemptFixedNeighborCertification(
            stationarity,sweep_index+1,next_certification_retry_sweep);
        if(stationarity==FixedNeighborStationarityState::ConfirmedStationary)
        {
            if(out.confirmed_stationarity_sweep==0) out.confirmed_stationarity_sweep=sweep_index+1;
            if(!policy.stop_after_stationarity || !policy.certify_endpoint)
            {
                out.search_converged=true; out.reason="block-stationary";
                if(policy.stop_after_stationarity) break;
                continue;
            }
            const bool retry_due=!out.certification || out.certification->passed() || certification_retry_due;
            if(!retry_due) continue;
            out.certification=certify_endpoint();
            out.certification_sweep=sweep_index+1;
            if(out.certification->passed())
            {
                out.search_converged=true; out.reason="block-stationary";
                break;
            }
            if(!out.certification->trust.primary_valid ||
                out.certification->trust.reason.starts_with("resource-limit:reference-solve"))
            {
                out.reason=out.certification->trust.reason.empty() ?
                    "endpoint-certification-invalid-primary" : out.certification->trust.reason;
                stop=true;
            }
        }
    }
    if(!out.search_converged && !stop)
    {
        out.reason="block-sweep-budget";
        if(policy.stop_after_stationarity && policy.certify_endpoint && out.sweep_count &&
            out.certification_sweep!=out.sweep_count)
        {
            out.certification=certify_endpoint();
            out.certification_sweep=out.sweep_count;
            if(!out.certification->trust.primary_valid ||
                out.certification->trust.reason.starts_with("resource-limit:reference-solve"))
                out.reason=out.certification->trust.reason.empty() ? "endpoint-certification-invalid-primary" :
                    out.certification->trust.reason;
        }
    }
    return out;
}

FixedNeighborEndpointCertification CertifyFixedNeighborEndpoint(
    const Domain & domain,VectorRef y,const BlockCoordinateState & state,const EvaluationContext & context,
    ProfileSearchWork * profile_work,std::size_t resident_override)
{
    FixedNeighborEndpointCertification out;
    const auto atom_count=static_cast<Eigen::Index>(domain.atoms.size());
    if(state.eta.size()!=atom_count || state.beta.size()!=2*atom_count || y.size()!=domain.rows)
        return out;
    out.reference_memory_estimate=ReferenceSolveMemoryEstimate(domain,y,context,resident_override);
    if(out.reference_memory_estimate.status!=MemoryStatus::Pass)
    {
        LogMemoryEstimate(out.reference_memory_estimate);
        const auto reason=std::string("resource-limit:reference-solve");
        out.assessment.primary=EvaluateStateStreaming(domain,y,state.eta,state.beta,context,
            ProfileEvaluationRole::AcceptedEndpoint,profile_work);
        out.assessment.reference.eta=state.eta; out.assessment.reference.reason=reason;
        out.assessment.failure=reason;
        out.trust.primary_valid=out.assessment.primary.valid; out.trust.reason=reason;
        return out;
    }
    if(out.reference_memory_estimate.estimated_peak_bytes>out.reference_memory_estimate.budget_bytes/2)
        LogMemoryEstimate(out.reference_memory_estimate);
    const auto endpoint=EvaluateState(domain,y,state.eta,state.beta,context,
        ProfileEvaluationRole::AcceptedEndpoint,profile_work);
    const auto reference=EvaluateProfile(domain,y,state.eta,true,&context,nullptr,
        ProfileEvaluationRole::Reference,profile_work);
    out.assessment=AssessEvaluated(domain,y,endpoint,reference,context,true);
    out.trust=CheckTrust(domain,y,endpoint,context,reference);
    return out;
}

std::optional<double> LastFixedNeighborObjective(const FixedNeighborSearchResult & result)
{
    if(result.final_sweep) return result.final_sweep->objective_after;
    return std::nullopt;
}

std::optional<double> LastFixedNeighborGradient(const FixedNeighborSearchResult & result)
{
    if(result.final_sweep) return result.final_sweep->global_width_gradient_inf_norm;
    return std::nullopt;
}

FixedNeighborComponentResult SolveFixedNeighborComponent(
    const ComponentProblemView & prepared,VectorRef initial_b,const EvaluationContext & parent_context,
    const FixedNeighborSearchPolicy & production_policy,const JointProgressObserver & observer,
    const JointProgressComponent * progress_component)
{
    FixedNeighborComponentResult out;
    const auto & local_layout=prepared.local_layout;
    if(local_layout.full_atoms.empty() || local_layout.informative_rows.empty())
    {
        out.stop_reason=local_layout.full_atoms.empty() ? "analytic-nuisance-only" : "unobserved-full-parameters";
        return out;
    }
    for(const auto atom:prepared.local_to_parent_atoms)
        if(atom<0 || atom>=initial_b.size() || !std::isfinite(initial_b(atom)) || !(initial_b(atom)>0))
        {
            out.stop_reason="invalid-initial-widths"; return out;
        }

    Vector local_y(static_cast<Eigen::Index>(prepared.local_to_parent_rows.size()));
    for(std::size_t k=0;k<prepared.local_to_parent_rows.size();++k)
        local_y(static_cast<Eigen::Index>(k))=(*parent_context.observations)(prepared.local_to_parent_rows[k]);
    auto context=parent_context;
    context.atom_ids=prepared.atom_ids; context.row_ids=prepared.row_ids;
    context.independent_search=true;
    context.target_columns=0;
    if(const auto * input=prepared.domain.atoms.Storage(); input && input->selection_domain)
    {
        std::size_t targets{};
        for(const auto atom:prepared.local_to_parent_atoms)
            targets+=std::binary_search(input->selection_domain->target_indices.begin(),
                input->selection_domain->target_indices.end(),static_cast<std::size_t>(atom));
        if(!CheckedMultiplyBytes(targets,3,context.target_columns))
            context.target_columns=std::numeric_limits<std::size_t>::max();
    }
    const auto atom_count=static_cast<Eigen::Index>(local_layout.full_atoms.size());
    context.rank={prepared.parent_row_count,2*atom_count,atom_count};
    context.linear.rank_relative=context.rank.Relative(2*atom_count);

    Vector initial_eta(atom_count);
    for(std::size_t k=0;k<prepared.local_to_parent_atoms.size();++k)
        initial_eta(static_cast<Eigen::Index>(k))=std::log(initial_b(prepared.local_to_parent_atoms[k]));

    FixedNeighborPolicy policy;
    policy.core_atoms=production_policy.core_atoms;
    policy.maximum_sweeps=production_policy.maximum_sweeps;
    policy.order=production_policy.order;
    if(progress_component && observer)
    {
        policy.production_progress=[&observer,progress_component,&production_policy,context](
            const FixedNeighborBlockSweep & sweep) {
            auto event=MakeJointProgressEvent(JointProgressPhase::SearchProgress,*progress_component);
            event.profile_evaluations=static_cast<int>(sweep.profile_evaluations);
            event.profile_budget=context.profile_budget;
            event.accepted_updates=static_cast<int>(sweep.accepted_blocks);
            event.update_budget=static_cast<int>(sweep.block_solves);
            event.elapsed_seconds=sweep.wall_seconds;
            event.accepted_objective=sweep.objective_after;
            event.accepted_gradient_inf_norm=sweep.global_width_gradient_inf_norm;
            event.fixed_neighbor=JointFixedNeighborProgress{
                sweep.sweep,production_policy.maximum_sweeps,sweep.block_solves,sweep.accepted_blocks,
                sweep.accepted_local_updates,
                sweep.objective_after,sweep.global_ac_kkt,sweep.global_width_gradient_inf_norm,sweep.eta_change_inf};
            observer(event);
        };
    }
    out.search=SearchFixedNeighborComponent(prepared.domain,local_layout,local_y,initial_eta,context,policy);
    const bool certified=out.search.certification && out.search.certification->passed();
    out.succeeded=out.search.search_converged && certified;
    out.stop_reason=out.search.reason;
    if(progress_component)
        NotifyJointProgress(observer,JointProgressPhase::CertificationStarted,*progress_component,
            static_cast<int>(out.search.total_profile_evaluations),parent_context.profile_budget,
            static_cast<int>(out.search.total_accepted_local_updates),parent_context.update_budget,
            out.search.search_seconds,out.stop_reason,certified,
            LastFixedNeighborObjective(out.search),LastFixedNeighborGradient(out.search));
    return out;
}

}
