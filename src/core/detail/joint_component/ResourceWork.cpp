#include "ResourceWork.hpp"
#include "Numerics.hpp"
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#if defined(__APPLE__)
#include <mach/mach.h>
#include <mach/task_info.h>
#elif defined(__linux__)
#include <unistd.h>
#endif

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
    if(input.selection_domain)
    {
        if(!AddProduct(total,input.selection_domain->target_indices.capacity(),sizeof(std::size_t))) return {};
    }
    return total;
}
std::optional<std::size_t> EstimateConstructionScratchBytes(const JointProblemInput & input)
{
    std::size_t total{},characters{},maximum_atom_support{};
    const auto node_bytes=sizeof(std::string)+4*sizeof(void *)+32;
    for(const auto & id:input.row_ids)
    {
        std::size_t bytes{};
        if(!CheckedAddBytes(id.capacity(),1,bytes) || !CheckedAddBytes(characters,bytes,characters)) return {};
    }
    if(!AddProduct(total,input.row_ids.size(),node_bytes) ||
        !AddProduct(total,characters,sizeof(char))) return {};
    characters=0;
    for(const auto & id:input.atom_ids)
    {
        std::size_t bytes{};
        if(!CheckedAddBytes(id.capacity(),1,bytes) || !CheckedAddBytes(characters,bytes,characters)) return {};
    }
    if(!AddProduct(total,input.atom_ids.size(),2*sizeof(std::string)+node_bytes) ||
        !AddProduct(total,characters,2*sizeof(char))) return {};
    for(const auto & atom:input.support) maximum_atom_support=std::max(maximum_atom_support,atom.size());
    if(!AddProduct(total,maximum_atom_support,sizeof(std::size_t)+4*sizeof(void *)+32)) return {};
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
    case MemoryStage::StructuralConstruction: return "structural-construction";
    case MemoryStage::FixedNeighborPreparation: return "fixed-neighbor-preparation";
    case MemoryStage::GlobalStateEvaluation: return "global-state-evaluation";
    case MemoryStage::EndpointVerification: return "endpoint-verification";
    case MemoryStage::LocalConditioning: return "local-conditioning";
    case MemoryStage::ParameterStability: return "parameter-stability";
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
    const auto blocks=shape.component_blocks,affected=shape.affected_row_entries;
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
            AddProduct(temporary,rows_plus_one,sizeof(std::size_t)) && AddProduct(temporary,rows,sizeof(std::uint8_t)) &&
            CheckedAddBytes(temporary,shape.identity_validation_bytes,temporary);
        break;
    case MemoryStage::StructuralConstruction:
    {
        // Coarse upper bound before StructuralTopology and core traversal allocate
        // row incidence, census marks and traversal state.
        ok=CheckedAddBytes(rows,1,rows_plus_one) &&
            AddProduct(temporary,rows_plus_one,2*sizeof(std::size_t)) &&
            AddProduct(temporary,rows,sizeof(std::uint8_t)+sizeof(std::size_t)) &&
            AddProduct(temporary,entries,2*sizeof(Eigen::Index)) &&
            AddProduct(temporary,atoms,sizeof(Eigen::Index)+5*sizeof(std::size_t)) &&
            AddProduct(temporary,blocks,3*sizeof(std::vector<Eigen::Index>)+64) &&
            AddProduct(temporary,atoms,sizeof(Eigen::Index));
        break;
    }
    case MemoryStage::FixedNeighborPreparation:
        // Each mapping stores sorted parent/local row pairs. The view can also retain
        // one support index per original membership; neither term scales as rows*blocks.
        ok=AddProduct(temporary,affected,6*sizeof(Eigen::Index)) &&
            AddProduct(temporary,entries,2*sizeof(Eigen::Index)) &&
            AddProduct(temporary,atoms,2*sizeof(Eigen::Index)) &&
            AddProduct(temporary,atoms,sizeof(std::vector<LocalSupportIndex>)) &&
            AddProduct(temporary,blocks,6*sizeof(std::vector<Eigen::Index>)) &&
            AddProduct(temporary,rows,6*sizeof(double)) &&
            AddProduct(temporary,atoms,3*sizeof(double)) &&
            AddProduct(temporary,shape.maximum_block_rows,4*sizeof(double)) &&
            AddProduct(temporary,shape.maximum_block_support_entries,4*(sizeof(double)+sizeof(Eigen::Index)));
        break;
    case MemoryStage::GlobalStateEvaluation:
        ok=AddProduct(temporary,rows,3*sizeof(double)+sizeof(std::uint8_t)) &&
            AddProduct(temporary,atoms,24*sizeof(double)+sizeof(Eigen::Index));
        break;
    case MemoryStage::EndpointVerification:
        // The primary residual remains live while the independent replay holds
        // compensated prediction accumulators and O(P) gradient/norm vectors.
        // Neither path forms an N-by-2P or N-by-3P design.
        ok=AddProduct(temporary,rows,4*sizeof(double)) &&
            AddProduct(temporary,atoms,24*sizeof(double)+3*sizeof(Eigen::Index));
        break;
    case MemoryStage::LocalConditioning:
    {
        const auto local_rows=shape.maximum_block_rows ? shape.maximum_block_rows : rows;
        const auto local_entries=shape.maximum_block_support_entries ? shape.maximum_block_support_entries : entries;
        const auto local_atoms=shape.maximum_neighbor_atoms ? shape.maximum_neighbor_atoms : atoms;
        std::size_t free_columns=shape.component_free_columns,two_widths{},tile_columns{},free_square{},projection{},
            width_square{},free_square_bytes{},projection_bytes{},width_square_bytes{},tile_bytes{};
        const auto width_columns=shape.component_width_columns ? shape.component_width_columns : local_atoms;
        if(!free_columns) ok=CheckedMultiplyBytes(local_atoms,2,free_columns);
        ok=ok && CheckedMultiplyBytes(width_columns,2,two_widths) &&
            CheckedAddBytes(free_columns,two_widths,tile_columns) &&
            AddProduct(temporary,local_rows,5*sizeof(double)) &&
            AddProduct(temporary,local_entries,12*(sizeof(double)+sizeof(Eigen::Index))) &&
            AddProduct(temporary,local_atoms,48*sizeof(double)+8*sizeof(Eigen::Index)) &&
            AddProduct2D(temporary,free_columns,free_columns,8*sizeof(double)) &&
            AddProduct2D(temporary,free_columns,width_columns,4*sizeof(double)) &&
            AddProduct2D(temporary,width_columns,width_columns,9*sizeof(double)) &&
            AddProduct2D(temporary,shape.maximum_tile_rows,tile_columns,6*sizeof(double)) &&
            CheckedMultiplyBytes(free_columns,free_columns,free_square) &&
            CheckedMultiplyBytes(free_columns,width_columns,projection) &&
            CheckedMultiplyBytes(width_columns,width_columns,width_square) &&
            CheckedMultiplyBytes(free_square,sizeof(double),free_square_bytes) &&
            CheckedMultiplyBytes(projection,sizeof(double),projection_bytes) &&
            CheckedMultiplyBytes(width_square,sizeof(double),width_square_bytes);
        if(ok)
        {
            out.factor_rows=free_columns; out.factor_columns=free_columns;
            ok=CheckedAddBytes(free_square_bytes,projection_bytes,out.dense_factor_bytes) &&
                CheckedAddBytes(out.dense_factor_bytes,width_square_bytes,out.dense_factor_bytes) &&
                CheckedMultiplyBytes(shape.maximum_tile_rows,tile_columns,tile_bytes) &&
                CheckedMultiplyBytes(tile_bytes,6*sizeof(double),out.factor_workspace_bytes) &&
                CheckedMultiplyBytes(local_entries,sizeof(double)+sizeof(Eigen::Index),out.sparse_design_bytes);
        }
        break;
    }
    case MemoryStage::ParameterStability:
    {
        const auto local=PlanMemory(shape,MemoryStage::LocalConditioning,budget);
        if(local.status==MemoryStatus::Overflow) {ok=false; break;}
        temporary=local.temporary_bytes;
        out.sparse_design_bytes=local.sparse_design_bytes;
        out.dense_factor_bytes=local.dense_factor_bytes;
        out.factor_workspace_bytes=local.factor_workspace_bytes;
        out.factor_rows=local.factor_rows; out.factor_columns=local.factor_columns;
        const auto local_atoms=shape.maximum_neighbor_atoms ? shape.maximum_neighbor_atoms : atoms;
        ok=AddProduct2D(temporary,shape.component_stability_evaluations,local_atoms,
            8*sizeof(double)+2*sizeof(Eigen::Index));
        break;
    }
    }
    if(!ok) {Overflow(out); return out;}
    out.temporary_bytes=temporary;
    if(!CheckedAddBytes(out.resident_bytes,out.temporary_bytes,out.estimated_peak_bytes)) {Overflow(out); return out;}
    if(out.estimated_peak_bytes>budget) {out.status=MemoryStatus::Reject; out.reason="estimated-peak-exceeds-budget";}
    return out;
}
std::optional<std::size_t> CurrentRssBytes()
{
#if defined(__APPLE__)
    mach_task_basic_info_data_t info{}; mach_msg_type_number_t count=MACH_TASK_BASIC_INFO_COUNT;
    if(task_info(mach_task_self(),MACH_TASK_BASIC_INFO,reinterpret_cast<task_info_t>(&info),&count)!=KERN_SUCCESS) return {};
    return static_cast<std::size_t>(info.resident_size);
#elif defined(__linux__)
    std::ifstream statm("/proc/self/statm"); std::size_t total_pages{},resident_pages{};
    const auto page_bytes=sysconf(_SC_PAGESIZE);
    if(!(statm>>total_pages>>resident_pages) || page_bytes<=0 ||
        resident_pages>std::numeric_limits<std::size_t>::max()/static_cast<std::size_t>(page_bytes)) return {};
    return resident_pages*static_cast<std::size_t>(page_bytes);
#else
    return {};
#endif
}
std::string FormatMemoryEstimate(const MemoryEstimate & estimate)
{
    const auto & s=estimate.shape; std::ostringstream out;
    out<<"[Joint][Memory]\n"
       <<"  stage="<<MemoryStageName(estimate.stage)<<'\n'
       <<"  scope="<<(s.problem_scope ? "problem" : "component")<<'\n'
       <<"  component_atoms="<<s.component_atoms<<'\n'
       <<"  component_rows="<<s.component_observations<<'\n'
       <<"  support_entries="<<s.component_support_entries<<'\n'
       <<"  prepared_blocks="<<s.component_blocks<<'\n'
       <<"  affected_row_entries="<<s.affected_row_entries<<'\n'
       <<"  neighbor_entries="<<s.component_neighbor_entries<<'\n'
       <<"  max_block_rows="<<s.maximum_block_rows<<'\n'
       <<"  max_block_support_entries="<<s.maximum_block_support_entries<<'\n'
       <<"  total_block_support_entries="<<s.total_block_support_entries<<'\n'
       <<"  maximum_neighbor_atoms="<<s.maximum_neighbor_atoms<<'\n'
       <<"  free_columns="<<s.component_free_columns<<'\n'
       <<"  width_columns="<<s.component_width_columns<<'\n'
       <<"  identity_validation_scratch="<<Bytes(s.identity_validation_bytes)<<'\n'
       <<"  resident="<<Bytes(estimate.resident_bytes)<<'\n'
       <<"  temporary="<<Bytes(estimate.temporary_bytes)<<'\n'
       <<"  estimated_peak="<<Bytes(estimate.estimated_peak_bytes)<<'\n'
       <<"  sparse_working_set="<<Bytes(estimate.sparse_design_bytes)<<'\n'
       <<"  dense_factor_workspace="<<Bytes(estimate.dense_factor_bytes)<<'\n'
       <<"  factor_workspace="<<Bytes(estimate.factor_workspace_bytes)<<'\n'
       <<"  budget="<<Bytes(estimate.budget_bytes)<<'\n'
       <<"  preflight="<<MemoryStatusName(estimate.status);
    if(!estimate.reason.empty()) out<<'\n'<<"  reason="<<estimate.reason;
    return out.str();
}
void LogMemoryEstimate(const MemoryEstimate & estimate)
{
    Logger::FinishProgressLine();
    if (estimate.status!=MemoryStatus::Pass)
    {
        Logger::Log(LogLevel::Warning,FormatMemoryEstimate(estimate));
        return;
    }
    if (estimate.estimated_peak_bytes>=estimate.budget_bytes-estimate.budget_bytes/5)
    {
        Logger::Log(LogLevel::Info,"[Joint][Memory] " + std::string(MemoryStageName(estimate.stage))+
            " | estimated_peak=" + Bytes(estimate.estimated_peak_bytes)+
            " | budget=" + Bytes(estimate.budget_bytes));
    }
    Logger::Log(LogLevel::Debug,FormatMemoryEstimate(estimate));
}


}
