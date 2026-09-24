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
struct PhaseWork {std::string phase; std::size_t calls{}; double inclusive_seconds{};};
struct ResourceWork
{
    bool enabled{};
    const char * phase{"unscoped"};
    std::vector<DenseShapeWork> dense_shapes;
    std::vector<PhaseWork> phases;
};
ResourceWork & ResourceWorkForTesting();
// Shape probes describe known matrices, not every allocation inside Eigen/BLAS.
// Aggregate by role and phase so diagnostic storage cannot grow with iterations.
void RecordDenseShape(const char *,Eigen::Index,Eigen::Index);
class ResourcePhase
{
    bool enabled_;
    const char * previous_{};
    const char * name_;
    std::chrono::steady_clock::time_point started_;
public:
    explicit ResourcePhase(const char *);
    ~ResourcePhase();
};
}
