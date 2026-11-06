#pragma once

#include <rhbm_gem/utils/domain/SamplingTypes.hpp>
#include <array>
#include <cstddef>
#include <vector>

namespace rhbm_gem { class AtomObject; class ModelObject; }

namespace rhbm_gem::core::detail {
using RadialPoint = std::array<double, 3>;
std::vector<RadialPoint> GenerateSpiralSpherePoints(
    const RadialPoint & center, double radius, std::size_t point_count);

struct QScoreShellDiagnostics
{
    double radius{};
    std::size_t requested_count{}, attempt_count{}, candidate_count{}, accepted_count{};
    bool successful{};
};
struct QScoreShellSamples
{
    SamplingPointList points;
    QScoreShellDiagnostics diagnostics;
};
QScoreShellSamples SampleQScoreShell(const AtomObject &, const ModelObject &, double radius, int num_points);

struct QScoreSeedSamples
{
    SamplingPointList points;
    std::vector<QScoreShellDiagnostics> shells;
};
// QScore-inspired FirstStage observations: one center, no repeated-center weight.
QScoreSeedSamples BuildQScoreSeedSamples(const AtomObject &, const ModelObject &);
} // namespace rhbm_gem::core::detail
