#pragma once
#include <Eigen/Core>
#include <chrono>
#include <cstddef>
#include <string>
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
    std::size_t reference{},factor_nonzeros{},cancellation_reductions{};
    double matrix_preparation_seconds{},reference_seconds{},reference_svd_seconds{},derivative_seconds{};
    std::size_t derivative_preparations{},derivative_compacts{},reference_compacts{},free_design_svds{},reference_svds{},reference_solves{},bdc_svds{},jacobi_retries{};
    double derivative_compact_seconds{},reference_compact_seconds{},free_design_svd_seconds{},reference_solve_seconds{},cancellation_seconds{},jacobi_retry_seconds{};
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
