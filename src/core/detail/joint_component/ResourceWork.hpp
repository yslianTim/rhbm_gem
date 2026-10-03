#pragma once
#include <Eigen/Core>
#include <chrono>
#include <cstddef>
#include <string>
#include <vector>

namespace rhbm_gem::core::joint_component {
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
