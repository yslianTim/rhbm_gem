#pragma once

#include "MapInterpolation.hpp"
#include "QScoreSampling.hpp"
#include <string_view>

namespace rhbm_gem::core::detail {
enum class JointSeedSamplingMethod { FibonacciShells, QScoreSpiralShells };
struct JointSeedSamplingPolicy
{
    JointSeedSamplingMethod sampling{JointSeedSamplingMethod::FibonacciShells};
    MapInterpolationMethod interpolation{MapInterpolationMethod::Tricubic};
};
std::string_view JointSeedSamplingName(JointSeedSamplingMethod);
std::string_view MapInterpolationName(MapInterpolationMethod);

struct JointSeedAtomSamplingDiagnostics
{
    int atom_id{};
    std::vector<QScoreShellDiagnostics> shells;
    std::size_t raw_count{}, selected_count{}, positive_count{}, fitting_count{}, positive_fitting_count{};
};
// Per-run research evidence; no public options, global flags or persistence schema changes.
struct JointSeedSamplingReport
{
    JointSeedSamplingPolicy policy;
    std::vector<JointSeedAtomSamplingDiagnostics> atoms;
};
void RunJointSeedSamplingWorkflow(
    MapObject &, ModelObject &, const std::vector<AtomObject *> &, int thread_count,
    std::string_view progress_label, bool quiet_mode, const JointSeedSamplingPolicy &,
    JointSeedSamplingReport * = nullptr);
} // namespace rhbm_gem::core::detail
