#include "ResourceWork.hpp"
#include "Numerics.hpp"
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <algorithm>
#include <cstring>
#include <iomanip>
#include <limits>
#include <sstream>

namespace rhbm_gem::core::joint_component {
bool CheckedAddBytes(std::size_t a,std::size_t b,std::size_t & out)
{
    if(b>std::numeric_limits<std::size_t>::max()-a) return false;
    out=a+b; return true;
}
bool CheckedMultiplyBytes(std::size_t a,std::size_t b,std::size_t & out)
{
    if(a && b>std::numeric_limits<std::size_t>::max()/a) return false;
    out=a*b; return true;
}
namespace {
bool AddProduct(std::size_t & total,std::size_t count,std::size_t width)
{
    std::size_t bytes{}; return CheckedMultiplyBytes(count,width,bytes) && CheckedAddBytes(total,bytes,total);
}
bool AddProduct2D(std::size_t & total,std::size_t rows,std::size_t columns,std::size_t width)
{
    std::size_t count{}; return CheckedMultiplyBytes(rows,columns,count) && AddProduct(total,count,width);
}
bool AddSquare(std::size_t & total,std::size_t n,std::size_t multiplier=1)
{
    std::size_t square{},bytes{};
    return CheckedMultiplyBytes(n,n,square) && CheckedMultiplyBytes(square,sizeof(double),bytes) &&
        CheckedMultiplyBytes(bytes,multiplier,bytes) && CheckedAddBytes(total,bytes,total);
}
std::string Bytes(std::size_t value)
{
    std::ostringstream out; out<<value<<" bytes ("<<std::fixed<<std::setprecision(3)
        <<static_cast<double>(value)/1'000'000'000.0<<" GB; "
        <<static_cast<double>(value)/(1024.0*1024.0*1024.0)<<" GiB)"; return out.str();
}
void Overflow(MemoryEstimate & estimate)
{estimate.status=MemoryStatus::Overflow; estimate.reason="size-overflow";}
}
std::optional<std::size_t> EstimateSnapshotBytes(const JointProblemInput & input)
{
    std::size_t total{};
    if(!AddProduct(total,input.observations.capacity(),sizeof(double)) ||
        !AddProduct(total,input.support.capacity(),sizeof(std::vector<JointSupport>)) ||
        !AddProduct(total,input.atom_ids.capacity(),sizeof(std::string)) ||
        !AddProduct(total,input.row_ids.capacity(),sizeof(std::string))) return {};
    for(const auto & atom:input.support)
        if(!AddProduct(total,atom.capacity(),sizeof(JointSupport))) return {};
    for(const auto & id:input.atom_ids)
    {
        std::size_t bytes{}; if(!CheckedAddBytes(id.capacity(),1,bytes) || !AddProduct(total,bytes,sizeof(char))) return {};
    }
    for(const auto & id:input.row_ids)
    {
        std::size_t bytes{}; if(!CheckedAddBytes(id.capacity(),1,bytes) || !AddProduct(total,bytes,sizeof(char))) return {};
    }
    if(input.selection_domain && !AddProduct(total,input.selection_domain->target_indices.capacity(),sizeof(std::size_t))) return {};
    return total;
}
std::optional<std::size_t> CountSupportEntries(const JointProblemInput & input)
{
    std::size_t total{};
    for(const auto & atom:input.support) if(!CheckedAddBytes(total,atom.size(),total)) return {};
    return total;
}
const char * MemoryStageName(MemoryStage stage)
{
    switch(stage)
    {
    case MemoryStage::ProblemConstruction: return "problem-construction";
    case MemoryStage::FixedNeighborPreparation: return "fixed-neighbor-preparation";
    case MemoryStage::GlobalStateEvaluation: return "global-state-evaluation";
    case MemoryStage::ReferenceSolve: return "reference-solve";
    case MemoryStage::RankAssessment: return "rank-assessment";
    case MemoryStage::TargetIdentifiability: return "target-identifiability";
    }
    return "unknown";
}
const char * MemoryStatusName(MemoryStatus status)
{
    switch(status)
    {
    case MemoryStatus::Pass: return "PASS";
    case MemoryStatus::Reject: return "REJECT";
    case MemoryStatus::Overflow: return "OVERFLOW";
    }
    return "UNKNOWN";
}
MemoryEstimate PlanMemory(const MemoryShape & shape,MemoryStage stage,std::size_t budget)
{
    MemoryEstimate out; out.shape=shape; out.stage=stage; out.resident_bytes=shape.resident_bytes; out.budget_bytes=budget;
    std::size_t temporary{};
    const auto rows=shape.component_observations,atoms=shape.component_atoms,entries=shape.component_support_entries;
    const auto blocks=shape.component_blocks,affected=shape.affected_row_entries,free=shape.component_free_columns;
    bool ok=true; std::size_t rows_plus_atoms{},rows_plus_one{};
    switch(stage)
    {
    case MemoryStage::ProblemConstruction:
        // Current partition construction holds row incidence, four identity/local maps,
        // union-find state, and structural row offsets concurrently.
        ok=CheckedAddBytes(rows,atoms,rows_plus_atoms) && CheckedAddBytes(rows,1,rows_plus_one) &&
            AddProduct(temporary,rows,sizeof(std::vector<Eigen::Index>)) &&
            AddProduct(temporary,entries,4*sizeof(Eigen::Index)) &&
            AddProduct(temporary,rows_plus_atoms,4*sizeof(Eigen::Index)) &&
            AddProduct(temporary,rows_plus_one,sizeof(std::size_t)) && AddProduct(temporary,rows,sizeof(std::uint8_t));
        break;
    case MemoryStage::FixedNeighborPreparation:
        // Each mapping stores sorted parent/local row pairs. The view can also retain
        // one support index per original membership; neither term scales as rows*blocks.
        ok=AddProduct(temporary,affected,6*sizeof(Eigen::Index)) &&
            AddProduct(temporary,entries,2*sizeof(Eigen::Index)) &&
            AddProduct(temporary,atoms,2*sizeof(Eigen::Index)) &&
            AddProduct(temporary,blocks,6*sizeof(std::vector<Eigen::Index>)) &&
            AddProduct(temporary,rows,6*sizeof(double)) &&
            AddProduct(temporary,atoms,3*sizeof(double)) &&
            AddProduct(temporary,shape.maximum_block_rows,4*sizeof(double)) &&
            AddProduct(temporary,shape.maximum_block_support_entries,4*(sizeof(double)+sizeof(Eigen::Index)));
        break;
    case MemoryStage::GlobalStateEvaluation:
        ok=AddProduct(temporary,rows,2*sizeof(double)+sizeof(std::uint8_t)) &&
            AddProduct(temporary,atoms,24*sizeof(double)+sizeof(Eigen::Index));
        break;
    case MemoryStage::ReferenceSolve:
    {
        // ReferenceQR retains a dense R and forms a dense tiled Householder workspace.
        // Include the sparse design/value copies while the factorization runs.
        std::size_t sparse_bytes{},factor_bytes{},workspace{};
        ok=AddProduct(sparse_bytes,entries,9*(sizeof(double)+sizeof(Eigen::Index))) &&
            CheckedMultiplyBytes(free,free,factor_bytes) && CheckedMultiplyBytes(factor_bytes,sizeof(double),factor_bytes) &&
            CheckedMultiplyBytes(factor_bytes,3,workspace) &&
            AddProduct(workspace,free,2*8192*sizeof(double)) &&
            CheckedAddBytes(sparse_bytes,factor_bytes,temporary) && CheckedAddBytes(temporary,workspace,temporary);
        break;
    }
    case MemoryStage::RankAssessment:
    {
        // TiledDerivative currently holds three free-column by width dense arrays,
        // then compact projected/Jacobian factors. This is an upper-bound workspace.
        const auto widths=atoms;
        ok=AddProduct2D(temporary,free,widths,3*sizeof(double)) &&
            AddProduct2D(temporary,widths,widths,2*sizeof(double)) &&
            AddProduct2D(temporary,2*8192,widths,sizeof(double));
        break;
    }
    case MemoryStage::TargetIdentifiability:
    {
        const auto columns=shape.component_target_columns;
        ok=AddProduct(temporary,entries,6*(sizeof(double)+sizeof(Eigen::Index))) &&
            AddSquare(temporary,columns,3) && AddProduct(temporary,columns,2*8192*sizeof(double));
        break;
    }
    }
    if(!ok) {Overflow(out); return out;}
    out.temporary_bytes=temporary;
    if(!CheckedAddBytes(out.resident_bytes,out.temporary_bytes,out.estimated_peak_bytes)) {Overflow(out); return out;}
    if(out.estimated_peak_bytes>budget) {out.status=MemoryStatus::Reject; out.reason="estimated-peak-exceeds-budget";}
    return out;
}
ReferenceFeasibility AssessReferenceFeasibility(const MemoryShape & shape,std::size_t budget)
{
    ReferenceFeasibility out;
    std::size_t temporary{},factor{},workspace{},matrix{};
    bool ok=CheckedMultiplyBytes(shape.component_free_columns,shape.component_free_columns,factor) &&
        CheckedMultiplyBytes(factor,sizeof(double),factor) &&
        AddProduct(matrix,shape.component_support_entries,9*(sizeof(double)+sizeof(Eigen::Index))) &&
        CheckedMultiplyBytes(factor,3,workspace) &&
        AddProduct(workspace,shape.component_free_columns,2*8192*sizeof(double)) &&
        CheckedAddBytes(matrix,factor,temporary) && CheckedAddBytes(temporary,workspace,temporary) &&
        CheckedAddBytes(shape.resident_bytes,temporary,out.estimated_peak_bytes);
    if(!ok) {out.budget_status=MemoryStatus::Overflow; out.reason="size-overflow"; return out;}
    out.estimated_matrix_bytes=matrix; out.estimated_factor_bytes=factor; out.estimated_workspace_bytes=workspace;
    if(out.estimated_peak_bytes>budget)
    {out.budget_status=MemoryStatus::Reject; out.numerical_status="unavailable"; out.reason="resource-limit-reference-solve";}
    else out.numerical_status="not-validated";
    return out;
}
std::string FormatMemoryEstimate(const MemoryEstimate & estimate)
{
    const auto & s=estimate.shape; std::ostringstream out;
    out<<"[Joint][Memory] stage="<<MemoryStageName(estimate.stage)<<" scope="<<(s.problem_scope ? "problem" : "component")
       <<" component_atoms="<<s.component_atoms
       <<" component_rows="<<s.component_observations<<" support_entries="<<s.component_support_entries
       <<" prepared_blocks="<<s.component_blocks<<" affected_row_entries="<<s.affected_row_entries
       <<" neighbor_entries="<<s.component_neighbor_entries
       <<" max_block_rows="<<s.maximum_block_rows<<" max_block_support_entries="<<s.maximum_block_support_entries
       <<" free_columns="<<s.component_free_columns<<" target_columns="<<s.component_target_columns
       <<" resident="<<Bytes(estimate.resident_bytes)<<" temporary="<<Bytes(estimate.temporary_bytes)
       <<" estimated_peak="<<Bytes(estimate.estimated_peak_bytes)<<" budget="<<Bytes(estimate.budget_bytes)
       <<" preflight="<<MemoryStatusName(estimate.status);
    if(!estimate.reason.empty()) out<<" reason="<<estimate.reason;
    return out.str();
}
void LogMemoryEstimate(const MemoryEstimate & estimate)
{Logger::Log(estimate.status==MemoryStatus::Pass ? LogLevel::Info : LogLevel::Warning,FormatMemoryEstimate(estimate));}

NumericsWork & NumericsWorkForTesting() {static thread_local NumericsWork work; return work;}
ResourceWork & ResourceWorkForTesting() {static thread_local ResourceWork work; return work;}
ProfileEvaluationRole & ProfileEvaluationRoleForTesting()
{static thread_local auto role=ProfileEvaluationRole::Unspecified; return role;}
ProfileEvaluationRoleScopeForTesting::ProfileEvaluationRoleScopeForTesting(ProfileEvaluationRole role)
    :previous_(ProfileEvaluationRoleForTesting())
{ProfileEvaluationRoleForTesting()=role;}
ProfileEvaluationRoleScopeForTesting::~ProfileEvaluationRoleScopeForTesting()
{ProfileEvaluationRoleForTesting()=previous_;}
namespace {
ResourceStageObserverForTesting & ResourceStageObserver()
{static thread_local ResourceStageObserverForTesting observer{}; return observer;}
void * & ResourceStageObserverContext()
{static thread_local void * context{}; return context;}
SearchStageWork & SearchStage(ResourceWork & work,const char * name)
{
    const auto found=std::find_if(work.search_stages.begin(),work.search_stages.end(),
        [&](const auto & stage){return stage.stage==name;});
    if(found!=work.search_stages.end()) return *found;
    work.search_stages.push_back({name}); return work.search_stages.back();
}
void NotifyResourceStageObserver()
{if(const auto observer=ResourceStageObserver()) observer(ResourceWorkForTesting(),ResourceStageObserverContext());}
}
void SetResourceStageObserverForTesting(ResourceStageObserverForTesting observer,void * context)
{ResourceStageObserver()=observer; ResourceStageObserverContext()=context;}
void RecordDenseShape(const char * role,Eigen::Index rows,Eigen::Index columns)
{
    auto & w=ResourceWorkForTesting(); if(!w.enabled) return;
    const auto bytes=static_cast<std::size_t>(rows)*static_cast<std::size_t>(columns)*sizeof(double);
    for(auto & record:w.dense_shapes) if(record.phase==w.phase && record.role==role)
    {
        ++record.observations;
        if(bytes>record.maximum_bytes) {record.rows=rows; record.columns=columns; record.maximum_bytes=bytes;}
        return;
    }
    w.dense_shapes.push_back({w.phase,role,rows,columns,1,bytes});
}
void RecordSparseShape(const char * role,Eigen::Index rows,Eigen::Index columns,std::size_t nonzeros)
{
    auto & w=ResourceWorkForTesting(); if(!w.enabled) return;
    for(auto & record:w.search_stages) if(record.stage==w.active_search_stage)
    {
        record.rows=std::max(record.rows,rows); record.columns=std::max(record.columns,columns);
        record.nonzeros=std::max(record.nonzeros,nonzeros); break;
    }
    // Preserve the phase/role record separately from stage aggregation.
    const std::string phase=w.phase, name=role;
    auto & shapes=w.sparse_shapes;
    for(auto & record:shapes) if(record.phase==phase && record.role==name)
    {
        ++record.observations;
        if(nonzeros>record.maximum_nonzeros) {record.rows=rows; record.columns=columns; record.maximum_nonzeros=nonzeros;}
        return;
    }
    shapes.push_back({phase,name,rows,columns,1,nonzeros});
}
ResourcePhase::ResourcePhase(const char * name,bool search_stage,Eigen::Index rows,Eigen::Index columns,std::size_t nonzeros)
    :enabled_(ResourceWorkForTesting().enabled),name_(name)
{
    if(!enabled_) return;
    auto & w=ResourceWorkForTesting(); previous_=w.phase; w.phase=name;
    search_scope_=std::strcmp(name,"search")==0;
    if(search_scope_) ++w.search_depth;
    search_stage_=search_scope_ || (search_stage && w.search_depth>0);
    if(search_stage_)
    {
        previous_search_stage_=w.active_search_stage;
        w.active_search_stage=name;
        auto & stage=SearchStage(w,name); ++stage.calls;
        if(rows>=0) stage.rows=std::max(stage.rows,rows);
        if(columns>=0) stage.columns=std::max(stage.columns,columns);
        stage.nonzeros=std::max(stage.nonzeros,nonzeros);
        NotifyResourceStageObserver();
    }
    started_=std::chrono::steady_clock::now();
}
ResourcePhase::~ResourcePhase()
{
    if(!enabled_) return;
    auto & w=ResourceWorkForTesting();
    const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-started_).count();
    bool found=false;
    for(auto & record:w.phases) if(record.phase==name_) {++record.calls; record.inclusive_seconds+=seconds; found=true; break;}
    if(!found) w.phases.push_back({name_,1,seconds});
    if(search_stage_)
    {
        auto & stage=SearchStage(w,name_); ++stage.completed_calls; stage.seconds+=seconds;
        w.last_completed_search_stage=name_;
        if(std::find(w.completed_search_stages.begin(),w.completed_search_stages.end(),name_)==w.completed_search_stages.end())
            w.completed_search_stages.push_back(name_);
        w.active_search_stage=previous_search_stage_;
    }
    if(search_scope_ && w.search_depth>0) --w.search_depth;
    w.phase=previous_;
    if(search_stage_) NotifyResourceStageObserver();
}
}
