#pragma once
#include <string>
#include <memory>
#include <boost/json.hpp>
#include <rhbm_gem/utils/hrl/RHBMTypes.hpp>
#include <rhbm_gem/utils/algorithm/WeightedRidgeSolver.hpp>
#include "core/detail/second_stage/SecondStageState.hpp"

namespace rhbm_gem::core::detail { struct JointOffsetSolveResult; }
namespace second_stage_test {
// Available only in BUILD_TESTING binaries. The directory is supplied by the test runner.
void CaptureShapeFailure(const rhbm_gem::RHBMMemberDataset &, double,
    const rhbm_gem::RHBMExecutionOptions &, const rhbm_gem::RHBMBetaEstimateResult &) noexcept;
void CaptureJointFailure(const rhbm_gem::algorithm::WeightedRidgeSystem &,
    const rhbm_gem::core::detail::JointOffsetSolveResult &) noexcept;
bool ReplaySolverFailure(const std::string & path);
struct SolverCapturePhase
{
    std::string name;
    std::size_t attempt{};
};
class ScopedSolverCapturePhase
{
    SolverCapturePhase previous;
public:
    explicit ScopedSolverCapturePhase(std::string name, std::size_t attempt = 0);
    ~ScopedSolverCapturePhase();
};
using SolverCaptureContext = std::shared_ptr<const boost::json::object>;
SolverCaptureContext CaptureOperatorContext(const rhbm_gem::core::detail::SecondStageContext &,
    const rhbm_gem::core::detail::FitState &,
    const std::vector<rhbm_gem::core::detail::ClusterKey> &) noexcept;
class ScopedSolverCaptureMember
{
    boost::json::object previous;
public:
    ScopedSolverCaptureMember(const SolverCaptureContext &, const std::vector<std::size_t> &,
        const char * role) noexcept;
    ~ScopedSolverCaptureMember();
};
void CaptureShapeResponse(const std::vector<double> &, const rhbm_gem::GaussianModel3D &) noexcept;
void CaptureShapeSolve(const rhbm_gem::RHBMBetaEstimateResult &) noexcept;
boost::json::object CurrentSolverCaptureMember();
}
