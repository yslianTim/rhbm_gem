#pragma once
#include <Eigen/Core>
#include <rhbm_gem/core/JointComponentEstimator.hpp>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace rhbm_gem::core::joint_component {
enum class ProfileEvaluationRole;

struct DenseShapeWork
{
    std::string phase,role;
    Eigen::Index rows{},columns{};
    std::size_t observations{},maximum_bytes{};
};
struct SparseShapeWork
{
    std::string phase,role;
    Eigen::Index rows{},columns{};
    std::size_t observations{},maximum_nonzeros{};
};
struct PhaseWork {std::string phase; std::size_t calls{}; double inclusive_seconds{};};
struct SearchStageWork
{
    std::string stage;
    std::size_t calls{},completed_calls{},nonzeros{};
    Eigen::Index rows{},columns{};
    double seconds{};
};

inline constexpr std::size_t JointMemoryHardLimitBytes=4'000'000'000;
inline constexpr std::size_t JointMemoryPlanningBudgetBytes=3'200'000'000;
enum class MemoryStage {ProblemConstruction,StructuralConstruction,FixedNeighborPreparation,GlobalStateEvaluation,
    EndpointVerification,ReferenceSolve,RankAssessment,StructuralRank,StreamingSpectralBounds,SparseQRResearch,DesignSpectrum,ProjectedWidthAssessment,
    ProjectedWidthDiagnostics,TargetIdentifiability};
enum class MemoryStatus {Pass,Reject,Overflow};
struct MemoryShape
{
    std::size_t component_atoms{},component_observations{},component_support_entries{},component_blocks{},
        affected_row_entries{},component_neighbor_entries{},maximum_block_rows{},maximum_block_support_entries{},
        total_block_support_entries{},maximum_neighbor_atoms{},component_free_columns{},component_target_columns{},
        identity_validation_bytes{},resident_bytes{},component_design_columns{},component_width_columns{},
        component_target_atoms{},component_halo_atoms{},
        maximum_tile_rows{8192},
        component_structural_nonzeros{},component_structural_supported_rows{};
    bool problem_scope{},component_structural_matching{};
};
struct MemoryStageEstimate
{
    MemoryStage stage{MemoryStage::ProblemConstruction};
    std::size_t resident_bytes{},temporary_bytes{},estimated_peak_bytes{},budget_bytes{JointMemoryPlanningBudgetBytes},
        sparse_design_bytes{},dense_factor_bytes{},factor_workspace_bytes{},right_singular_vector_bytes{},
        maximum_tile_rows{};
    std::size_t factor_rows{},factor_columns{};
    std::size_t free_design_factor_bytes{},free_design_svd_workspace_bytes{},projection_coefficient_bytes{},
        projected_width_factor_bytes{},projected_width_svd_workspace_bytes{},qr_tile_workspace_bytes{};
    MemoryStatus status{MemoryStatus::Pass};
    std::string reason;
};
struct MemoryEstimate
{
    MemoryShape shape;
    MemoryStage stage{MemoryStage::ProblemConstruction};
    std::size_t resident_bytes{},temporary_bytes{},estimated_peak_bytes{},budget_bytes{JointMemoryPlanningBudgetBytes};
    MemoryStatus status{MemoryStatus::Pass};
    std::string reason;
    std::size_t sparse_design_bytes{},dense_factor_bytes{},factor_workspace_bytes{},right_singular_vector_bytes{},
        factor_rows{},factor_columns{};
    std::size_t free_design_factor_bytes{},free_design_svd_workspace_bytes{},projection_coefficient_bytes{},
        projected_width_factor_bytes{},projected_width_svd_workspace_bytes{},qr_tile_workspace_bytes{};
    std::vector<MemoryStageEstimate> attributed_stages;
};
struct ReferenceFeasibility
{
    std::size_t estimated_matrix_bytes{},estimated_factor_bytes{},estimated_workspace_bytes{},estimated_peak_bytes{};
    MemoryStatus budget_status{MemoryStatus::Pass};
    std::string numerical_status{"not-run"};
    std::string reason;
};
bool CheckedAddBytes(std::size_t,std::size_t,std::size_t &);
bool CheckedMultiplyBytes(std::size_t,std::size_t,std::size_t &);
std::optional<std::size_t> EstimateSnapshotBytes(const JointProblemInput &);
std::optional<std::size_t> EstimateConstructionScratchBytes(const JointProblemInput &);
std::optional<std::size_t> CountSupportEntries(const JointProblemInput &);
MemoryEstimate PlanMemory(const MemoryShape &,MemoryStage,std::size_t=JointMemoryPlanningBudgetBytes);
ReferenceFeasibility AssessReferenceFeasibility(const MemoryShape &,std::size_t=JointMemoryPlanningBudgetBytes);
const char * MemoryStageName(MemoryStage);
const char * MemoryStatusName(MemoryStatus);
std::optional<std::size_t> ProcessPeakRssBytes();
std::optional<std::size_t> CurrentRssBytes();
std::string FormatMemoryEstimate(const MemoryEstimate &);
void LogMemoryEstimate(const MemoryEstimate &);

struct ResourceWork
{
    bool enabled{};
    const char * phase{"unscoped"};
    std::vector<DenseShapeWork> dense_shapes;
    std::vector<SparseShapeWork> sparse_shapes;
    std::vector<PhaseWork> phases;
    std::string active_search_stage,last_completed_search_stage;
    std::vector<std::string> completed_search_stages;
    std::vector<SearchStageWork> search_stages;
    std::size_t search_depth{};
};

struct NumericsWork
{
    std::size_t reference{},cancellation_reductions{},derivative_preparations{};
    std::size_t free_design_svds{},reference_svds{},reference_solves{},bdc_svds{},jacobi_retries{},design_spectrum_calls{};
    double derivative_seconds{},jacobi_retry_seconds{};
};
NumericsWork & NumericsWorkForTesting();
struct WorkTimer
{
    double & seconds;
    std::chrono::steady_clock::time_point started{std::chrono::steady_clock::now()};
    explicit WorkTimer(double & value):seconds(value) {}
    ~WorkTimer() {seconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();}
};
ProfileEvaluationRole & ProfileEvaluationRoleForTesting();
class ProfileEvaluationRoleScopeForTesting
{
    ProfileEvaluationRole previous_;
public:
    explicit ProfileEvaluationRoleScopeForTesting(ProfileEvaluationRole);
    ~ProfileEvaluationRoleScopeForTesting();
    ProfileEvaluationRoleScopeForTesting(const ProfileEvaluationRoleScopeForTesting &)=delete;
    ProfileEvaluationRoleScopeForTesting & operator=(const ProfileEvaluationRoleScopeForTesting &)=delete;
};

ResourceWork & ResourceWorkForTesting();
using ResourceStageObserverForTesting=void (*)(const ResourceWork &,void *);
void SetResourceStageObserverForTesting(ResourceStageObserverForTesting,void *);
// Shape probes describe known matrices, not every allocation inside Eigen/BLAS.
// Aggregate by role and phase so diagnostic storage cannot grow with iterations.
void RecordDenseShape(const char *,Eigen::Index,Eigen::Index);
void RecordSparseShape(const char *,Eigen::Index,Eigen::Index,std::size_t);
class ResourcePhase
{
    bool enabled_,search_stage_{},search_scope_{};
    const char * previous_{};
    const char * name_;
    std::string previous_search_stage_;
    std::chrono::steady_clock::time_point started_;
public:
    explicit ResourcePhase(const char *,bool=false,Eigen::Index=-1,Eigen::Index=-1,std::size_t=0);
    ~ResourcePhase();
};
}
