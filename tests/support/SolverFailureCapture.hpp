#pragma once
#include <string>
#include <rhbm_gem/utils/hrl/RHBMTypes.hpp>
#include <rhbm_gem/utils/algorithm/WeightedRidgeSolver.hpp>

namespace rhbm_gem::core::detail { struct JointOffsetSolveResult; }
namespace second_stage_test {
// Available only in BUILD_TESTING binaries. The directory is supplied by the test runner.
void CaptureShapeFailure(const rhbm_gem::RHBMMemberDataset &, double,
    const rhbm_gem::RHBMExecutionOptions &, const rhbm_gem::RHBMBetaEstimateResult &) noexcept;
void CaptureJointFailure(const rhbm_gem::algorithm::WeightedRidgeSystem &,
    const rhbm_gem::core::detail::JointOffsetSolveResult &) noexcept;
bool ReplaySolverFailure(const std::string & path);
}
