#include "ResourceWork.hpp"
#include "Numerics.hpp"
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <algorithm>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#if defined(__APPLE__) || defined(__linux__)
#include <sys/resource.h>
#endif
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
bool AddTileWorkspace(std::size_t & total,std::size_t tile_rows,std::size_t columns)
{
    std::size_t rows{}; return CheckedMultiplyBytes(tile_rows,2,rows) &&
        AddProduct2D(total,rows,columns,sizeof(double));
}
std::string Bytes(std::size_t value)
{
    std::ostringstream out; out<<value<<" bytes ("<<std::fixed<<std::setprecision(3)
        <<static_cast<double>(value)/1'000'000'000.0<<" GB; "
        <<static_cast<double>(value)/(1024.0*1024.0*1024.0)<<" GiB)"; return out.str();
}
void Overflow(MemoryEstimate & estimate)
{estimate.status=MemoryStatus::Overflow; estimate.reason="size-overflow";}
MemoryStageEstimate Attribute(const MemoryEstimate & estimate)
{
    MemoryStageEstimate out;
    out.stage=estimate.stage; out.resident_bytes=estimate.resident_bytes; out.temporary_bytes=estimate.temporary_bytes;
    out.estimated_peak_bytes=estimate.estimated_peak_bytes; out.budget_bytes=estimate.budget_bytes;
    out.sparse_design_bytes=estimate.sparse_design_bytes; out.dense_factor_bytes=estimate.dense_factor_bytes;
    out.factor_workspace_bytes=estimate.factor_workspace_bytes;
    out.right_singular_vector_bytes=estimate.right_singular_vector_bytes;
    out.maximum_tile_rows=estimate.shape.maximum_tile_rows;
    out.factor_rows=estimate.factor_rows; out.factor_columns=estimate.factor_columns;
    out.free_design_factor_bytes=estimate.free_design_factor_bytes;
    out.free_design_svd_workspace_bytes=estimate.free_design_svd_workspace_bytes;
    out.projection_coefficient_bytes=estimate.projection_coefficient_bytes;
    out.projected_width_factor_bytes=estimate.projected_width_factor_bytes;
    out.projected_width_svd_workspace_bytes=estimate.projected_width_svd_workspace_bytes;
    out.qr_tile_workspace_bytes=estimate.qr_tile_workspace_bytes;
    out.status=estimate.status; out.reason=estimate.reason;
    return out;
}
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
    case MemoryStage::ReferenceSolve: return "reference-solve";
    case MemoryStage::RankAssessment: return "rank-assessment";
    case MemoryStage::StructuralRank: return "structural-rank";
    case MemoryStage::StreamingSpectralBounds: return "streaming-spectral-bounds";
    case MemoryStage::SparseQRResearch: return "sparseqr-research";
    case MemoryStage::DesignSpectrum: return "design-spectrum";
    case MemoryStage::ProjectedWidthAssessment: return "projected-width-assessment";
    case MemoryStage::TargetIdentifiability: return "target-identifiability";
    case MemoryStage::TargetCovariance: return "target-covariance";
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
        const auto design=PlanMemory(shape,MemoryStage::DesignSpectrum,budget);
        const auto projected=PlanMemory(shape,MemoryStage::ProjectedWidthAssessment,budget);
        if(design.status==MemoryStatus::Overflow || projected.status==MemoryStatus::Overflow) {ok=false; break;}
        temporary=std::max(design.temporary_bytes,projected.temporary_bytes);
        out.attributed_stages={Attribute(design),Attribute(projected)};
        out.sparse_design_bytes=projected.sparse_design_bytes;
        out.free_design_factor_bytes=projected.free_design_factor_bytes;
        out.free_design_svd_workspace_bytes=projected.free_design_svd_workspace_bytes;
        out.projection_coefficient_bytes=projected.projection_coefficient_bytes;
        out.projected_width_factor_bytes=projected.projected_width_factor_bytes;
        out.projected_width_svd_workspace_bytes=projected.projected_width_svd_workspace_bytes;
        out.qr_tile_workspace_bytes=projected.qr_tile_workspace_bytes;
        break;
    }
    case MemoryStage::StructuralRank:
    {
        const auto supported=shape.component_structural_supported_rows ?
            shape.component_structural_supported_rows : std::min(rows,shape.component_structural_nonzeros);
        std::size_t offsets{};
        ok=AddProduct(temporary,shape.component_structural_nonzeros,2*sizeof(Eigen::Index)) &&
            AddProduct(temporary,supported,sizeof(Eigen::Index)) &&
            AddProduct(temporary,shape.component_design_columns,sizeof(std::uint8_t));
        if(ok && shape.component_structural_matching)
        {
            ok=CheckedAddBytes(shape.component_design_columns,1,offsets) &&
                AddProduct(temporary,offsets,sizeof(std::size_t)) &&
                AddProduct(temporary,shape.component_design_columns,7*sizeof(Eigen::Index)+sizeof(std::size_t)) &&
                AddProduct(temporary,supported,sizeof(Eigen::Index));
        }
        break;
    }
    case MemoryStage::StreamingSpectralBounds:
    {
        std::size_t columns=shape.component_design_columns;
        if(!columns) columns=shape.component_free_columns;
        ok=AddProduct(temporary,rows,2*sizeof(double)) &&
            AddProduct(temporary,columns,4*sizeof(double)+2*sizeof(Eigen::Index)+sizeof(std::uint8_t));
        break;
    }
    case MemoryStage::SparseQRResearch:
    {
        const auto columns=shape.component_design_columns;
        // Conservative experimental allowance for column-sized pivot data, dense
        // worst-case R/work arrays, dense row-column workspace, sparse copies and
        // a 5% allocation margin for Eigen/allocator overhead.
        ok=AddProduct2D(temporary,columns,columns,3*sizeof(double)+6*sizeof(Eigen::Index)) &&
            AddProduct2D(temporary,rows,columns,sizeof(double)) &&
            AddProduct(temporary,shape.component_structural_nonzeros,6*(sizeof(double)+sizeof(Eigen::Index))) &&
            AddProduct(temporary,columns,16*(sizeof(double)+sizeof(Eigen::Index)));
        if(ok) ok=CheckedAddBytes(temporary,temporary/20,temporary);
        break;
    }
    case MemoryStage::DesignSpectrum:
    {
        // DesignSpectrum calls the dense tiled ReferenceQR rank primitive. These are
        // sparse copies and dense factor work, not the removed production reference solve.
        std::size_t columns{},factor{},workspace{},sparse{};
        columns=shape.component_design_columns;
        if(!columns) ok=CheckedMultiplyBytes(atoms,2,columns);
        ok=ok && CheckedMultiplyBytes(columns,columns,factor) &&
            CheckedMultiplyBytes(factor,sizeof(double),factor) &&
            CheckedMultiplyBytes(factor,3,workspace) &&
            AddTileWorkspace(workspace,shape.maximum_tile_rows,columns) &&
            AddProduct(sparse,entries,4*(sizeof(double)+sizeof(Eigen::Index))) &&
            CheckedAddBytes(factor,workspace,temporary) && CheckedAddBytes(temporary,sparse,temporary);
        if(ok)
        {
            out.sparse_design_bytes=sparse; out.dense_factor_bytes=factor;
            out.factor_workspace_bytes=workspace; out.factor_rows=columns; out.factor_columns=columns;
        }
        break;
    }
    case MemoryStage::ProjectedWidthAssessment:
    {
        const auto widths=shape.component_width_columns ? shape.component_width_columns : atoms;
        std::size_t design_columns=shape.component_design_columns,free_columns=free;
        if(!design_columns) ok=CheckedMultiplyBytes(atoms,2,design_columns);
        if(ok && !free_columns) free_columns=design_columns;
        std::size_t free_factor{},free_svd_workspace{},projection_coefficients{},width_factor{},width_svd_workspace{};
        std::size_t free_entries{},raw_entries{},triplet_entries{},matrix_entries{},outer_rows{},sparse_matrices{},
            triplets{},sparse_peak{},free_index_bytes{},free_common{},tile_rows{},free_assembled_rows{},
            width_assembled_rows{},free_qr_workspace{},free_qr_peak{},free_svd_peak{},
            free_phase_peak{},width_qr_workspace{},width_qr_peak{},width_svd_peak{},width_phase_peak{};
        const std::size_t index_bytes=sizeof(Eigen::Index),triplet_bytes=sizeof(double)+2*sizeof(Eigen::Index);
        const auto add_product_bytes=[](std::size_t rows,std::size_t columns,std::size_t width,std::size_t & bytes) {
            std::size_t count{}; return CheckedMultiplyBytes(rows,columns,count) && CheckedMultiplyBytes(count,width,bytes);
        };
        ok=ok && design_columns>0 && free_columns<=design_columns && widths>0 &&
            add_product_bytes(free_columns,free_columns,sizeof(double),free_factor) &&
            CheckedMultiplyBytes(free_factor,3,free_svd_workspace) &&
            add_product_bytes(free_columns,widths,sizeof(double),projection_coefficients) &&
            add_product_bytes(widths,widths,sizeof(double),width_factor) &&
            CheckedMultiplyBytes(width_factor,3,width_svd_workspace);
        if(ok && entries)
        {
            std::size_t scaled_entries{},rounded_entries{},twice_entries{};
            ok=CheckedMultiplyBytes(entries,free_columns,scaled_entries) &&
                CheckedAddBytes(scaled_entries,design_columns-1,rounded_entries) &&
                CheckedMultiplyBytes(entries,2,twice_entries);
            if(ok)
            {
                free_entries=rounded_entries/design_columns; raw_entries=entries;
                ok=CheckedAddBytes(free_entries,raw_entries,matrix_entries) &&
                    CheckedMultiplyBytes(matrix_entries,sizeof(double)+index_bytes,sparse_matrices) &&
                    CheckedAddBytes(rows,1,outer_rows) &&
                    AddProduct(sparse_matrices,outer_rows,2*index_bytes) &&
                    CheckedAddBytes(free_entries,twice_entries,triplet_entries) &&
                    CheckedMultiplyBytes(triplet_entries,triplet_bytes,triplets) &&
                    CheckedMultiplyBytes(triplets,2,triplets);
            }
        }
        if(ok) ok=CheckedAddBytes(sparse_matrices,triplets,sparse_peak) &&
            CheckedMultiplyBytes(free_columns,index_bytes,free_index_bytes) &&
            CheckedAddBytes(sparse_peak,free_index_bytes,free_common) &&
            CheckedMultiplyBytes(shape.maximum_tile_rows,1,tile_rows) &&
            CheckedAddBytes(std::min(free_columns,rows),tile_rows,free_assembled_rows) &&
            CheckedAddBytes(std::min(widths,rows),tile_rows,width_assembled_rows) &&
            AddProduct(free_qr_workspace,2,free_factor) &&
            AddProduct(free_qr_workspace,2,projection_coefficients) &&
            AddProduct2D(free_qr_workspace,tile_rows,free_columns,sizeof(double)) &&
            AddProduct2D(free_qr_workspace,tile_rows,widths,sizeof(double)) &&
            AddProduct2D(free_qr_workspace,free_assembled_rows,free_columns,sizeof(double)) &&
            AddProduct2D(free_qr_workspace,free_assembled_rows,widths,sizeof(double)) &&
            AddProduct2D(free_qr_workspace,free_assembled_rows,free_columns,sizeof(double)) &&
            AddProduct(free_qr_workspace,free_assembled_rows,sizeof(double));
        if(ok) ok=CheckedAddBytes(free_common,free_qr_workspace,free_qr_peak) &&
            AddProduct(free_svd_peak,4,free_factor) &&
            CheckedAddBytes(free_svd_peak,projection_coefficients,free_svd_peak) &&
            AddProduct(free_svd_peak,free_columns,sizeof(double)) &&
            CheckedAddBytes(free_common,free_svd_peak,free_svd_peak) &&
            CheckedAddBytes(free_common,free_qr_peak,free_qr_peak) &&
            AddProduct(width_qr_workspace,2,width_factor) &&
            AddProduct2D(width_qr_workspace,width_assembled_rows,widths,sizeof(double)) &&
            AddProduct2D(width_qr_workspace,tile_rows,widths,3*sizeof(double)) &&
            AddProduct(width_qr_workspace,free_entries,sizeof(double)+index_bytes) &&
            AddProduct(width_qr_workspace,width_assembled_rows,sizeof(double));
        if(ok) ok=CheckedAddBytes(sparse_matrices,projection_coefficients,width_qr_peak) &&
            CheckedAddBytes(width_qr_peak,width_qr_workspace,width_qr_peak) &&
            AddProduct(width_svd_peak,5,width_factor) &&
            CheckedAddBytes(width_svd_peak,projection_coefficients,width_svd_peak) &&
            AddProduct(width_svd_peak,widths,sizeof(double)) &&
            CheckedAddBytes(sparse_matrices,width_svd_peak,width_svd_peak);
        if(ok)
        {
            free_phase_peak=std::max(free_qr_peak,free_svd_peak);
            width_phase_peak=std::max(width_qr_peak,width_svd_peak);
            temporary=std::max(free_phase_peak,width_phase_peak);
            out.sparse_design_bytes=sparse_peak;
            out.free_design_factor_bytes=free_factor;
            out.free_design_svd_workspace_bytes=free_svd_workspace;
            out.projection_coefficient_bytes=projection_coefficients;
            out.projected_width_factor_bytes=width_factor;
            out.projected_width_svd_workspace_bytes=width_svd_workspace;
            out.qr_tile_workspace_bytes=std::max(free_qr_workspace,width_qr_workspace);
            out.dense_factor_bytes=free_factor;
            out.factor_workspace_bytes=std::max({free_svd_workspace,width_svd_workspace,out.qr_tile_workspace_bytes});
            out.right_singular_vector_bytes=width_factor;
            out.factor_rows=free_columns; out.factor_columns=free_columns;
        }
        break;
    }
    case MemoryStage::TargetIdentifiability:
    {
        std::size_t evidence_atoms{},rank_columns{},checks_bytes{};
        ok=CheckedAddBytes(shape.component_target_atoms,shape.component_halo_atoms,evidence_atoms) &&
            CheckedAddBytes(evidence_atoms,shape.component_nuisance_atoms,evidence_atoms) &&
            CheckedAddBytes(shape.component_design_columns,shape.component_width_columns,rank_columns) &&
            AddProduct(temporary,evidence_atoms,sizeof(JointAtomIdentifiability)+32) &&
            AddProduct(temporary,rank_columns,sizeof(double)) &&
            AddProduct(temporary,2,sizeof(JointRankEvidence)+64) &&
            CheckedMultiplyBytes(8,sizeof(JointCheck)+32,checks_bytes) &&
            CheckedAddBytes(temporary,checks_bytes,temporary);
        break;
    }
    case MemoryStage::TargetCovariance:
    {
        const auto columns=shape.component_covariance_columns;
        std::size_t factor{},workspace{},tile{},tile_rows{},tile_count{},factor_count{},factor_workspace{};
        ok=columns>0 && columns<=static_cast<std::size_t>(std::numeric_limits<Eigen::Index>::max()) &&
            shape.maximum_tile_rows>0 &&
            CheckedMultiplyBytes(columns,columns,factor_count) &&
            CheckedMultiplyBytes(factor_count,sizeof(double),factor) &&
            CheckedMultiplyBytes(factor,2,workspace) &&
            CheckedMultiplyBytes(shape.maximum_tile_rows,2,tile_rows) &&
            AddProduct2D(tile,tile_rows,columns,sizeof(double)) &&
            CheckedAddBytes(temporary,factor,temporary) && CheckedAddBytes(temporary,factor,temporary) &&
            AddProduct(temporary,shape.component_support_entries,3*sizeof(std::size_t)) &&
            AddProduct(temporary,rows,4*sizeof(double)+2*sizeof(Eigen::Index)) &&
            AddProduct(temporary,columns,4*sizeof(double)) &&
            CheckedAddBytes(rows/shape.maximum_tile_rows,rows%shape.maximum_tile_rows ? 1 : 0,tile_count) &&
            AddProduct(temporary,tile_count,sizeof(std::vector<std::size_t>)) &&
            CheckedAddBytes(workspace,tile,factor_workspace) &&
            CheckedAddBytes(temporary,workspace,temporary) && CheckedAddBytes(temporary,tile,temporary);
        if(ok)
        {
            out.dense_factor_bytes=factor; out.right_singular_vector_bytes=factor;
            out.factor_workspace_bytes=factor_workspace;
            out.factor_rows=columns;
            out.factor_columns=columns;
        }
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
std::optional<std::size_t> ProcessPeakRssBytes()
{
#if defined(__APPLE__) || defined(__linux__)
    rusage usage{}; if(getrusage(RUSAGE_SELF,&usage)!=0 || usage.ru_maxrss<0) return {};
#if defined(__APPLE__)
    return static_cast<std::size_t>(usage.ru_maxrss);
#else
    std::size_t bytes{};
    if(!CheckedMultiplyBytes(static_cast<std::size_t>(usage.ru_maxrss),1024,bytes)) return {};
    return bytes;
#endif
#else
    return {};
#endif
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
       <<"  target_columns="<<s.component_target_columns<<'\n'
       <<"  design_columns="<<s.component_design_columns<<'\n'
       <<"  width_columns="<<s.component_width_columns<<'\n'
       <<"  target_atoms="<<s.component_target_atoms<<'\n'
       <<"  halo_atoms="<<s.component_halo_atoms<<'\n'
       <<"  nuisance_atoms="<<s.component_nuisance_atoms<<'\n'
       <<"  covariance_columns="<<s.component_covariance_columns<<'\n'
       <<"  identity_validation_scratch="<<Bytes(s.identity_validation_bytes)<<'\n'
       <<"  resident="<<Bytes(estimate.resident_bytes)<<'\n'
       <<"  temporary="<<Bytes(estimate.temporary_bytes)<<'\n'
       <<"  estimated_peak="<<Bytes(estimate.estimated_peak_bytes)<<'\n'
       <<"  sparse_working_set="<<Bytes(estimate.sparse_design_bytes)<<'\n'
       <<"  free_design_factor="<<Bytes(estimate.free_design_factor_bytes)<<'\n'
       <<"  free_design_svd_workspace="<<Bytes(estimate.free_design_svd_workspace_bytes)<<'\n'
       <<"  projection_coefficients="<<Bytes(estimate.projection_coefficient_bytes)<<'\n'
       <<"  projected_width_factor="<<Bytes(estimate.projected_width_factor_bytes)<<'\n'
       <<"  projected_width_svd_workspace="<<Bytes(estimate.projected_width_svd_workspace_bytes)<<'\n'
       <<"  qr_tile_workspace="<<Bytes(estimate.qr_tile_workspace_bytes)<<'\n'
       <<"  budget="<<Bytes(estimate.budget_bytes)<<'\n'
       <<"  preflight="<<MemoryStatusName(estimate.status);
    if(!estimate.reason.empty()) out<<'\n'<<"  reason="<<estimate.reason;
    for(const auto & stage:estimate.attributed_stages)
    {
        out<<'\n'<<"  attributed_stage="<<MemoryStageName(stage.stage)
           <<" status="<<MemoryStatusName(stage.status)
           <<" resident="<<Bytes(stage.resident_bytes)
           <<" temporary="<<Bytes(stage.temporary_bytes)
           <<" estimated_peak="<<Bytes(stage.estimated_peak_bytes)
           <<" factor="<<stage.factor_rows<<'x'<<stage.factor_columns
           <<" factor_bytes="<<Bytes(stage.dense_factor_bytes)
           <<" workspace_bytes="<<Bytes(stage.factor_workspace_bytes)
           <<" right_singular_vector_bytes="<<Bytes(stage.right_singular_vector_bytes)
           <<" free_design_factor_bytes="<<Bytes(stage.free_design_factor_bytes)
           <<" free_design_svd_workspace_bytes="<<Bytes(stage.free_design_svd_workspace_bytes)
           <<" projection_coefficient_bytes="<<Bytes(stage.projection_coefficient_bytes)
           <<" projected_width_factor_bytes="<<Bytes(stage.projected_width_factor_bytes)
           <<" projected_width_svd_workspace_bytes="<<Bytes(stage.projected_width_svd_workspace_bytes)
           <<" qr_tile_workspace_bytes="<<Bytes(stage.qr_tile_workspace_bytes)
           <<" sparse_design_bytes="<<Bytes(stage.sparse_design_bytes)
           <<" maximum_tile_rows="<<stage.maximum_tile_rows;
        if(!stage.reason.empty()) out<<" reason="<<stage.reason;
    }
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
    std::size_t bytes{};
    if(rows<0 || columns<0 || !CheckedMultiplyBytes(static_cast<std::size_t>(rows),static_cast<std::size_t>(columns),bytes) ||
        !CheckedMultiplyBytes(bytes,sizeof(double),bytes)) bytes=std::numeric_limits<std::size_t>::max();
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
