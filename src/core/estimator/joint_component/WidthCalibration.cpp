#include "WidthCalibration.hpp"

#include <rhbm_gem/utils/math/ArrayHelper.hpp>
#include <algorithm>
#include <cmath>
#include <limits>

namespace rhbm_gem::core::joint_component {
namespace {
constexpr std::size_t MinimumQualifiedWidthDonors=20;
constexpr double MinimumQualifiedDonorCoverage=0.95;

std::optional<double> BoundFrom(double reference_width, double multiplier)
{
    if (!(reference_width>0.0) || !std::isfinite(reference_width) ||
        reference_width>std::numeric_limits<double>::max()/multiplier)
        return std::nullopt;
    const double bound=multiplier*reference_width;
    if (!(bound>0.0) || !std::isfinite(bound) || !std::isfinite(std::log(bound)))
        return std::nullopt;
    return bound;
}

} // namespace

WidthBoundCalibration CalibrateGlobalWidth(const JointInitialization & initialization)
{
    WidthBoundCalibration result;
    if(initialization.atoms.size()!=initialization.tau.size())
    {
        result.unavailable_reason="initialization-record-size-mismatch";
        return result;
    }

    for(std::size_t i=0;i<initialization.atoms.size();++i)
        if(initialization.atoms[i].target &&
            (!std::isfinite(initialization.tau[i]) || !(initialization.tau[i]>0.0)))
        {
            result.unavailable_reason="invalid-target-initial-width";
            return result;
        }

    std::vector<double> donors;
    donors.reserve(initialization.atoms.size());
    for(const auto & atom:initialization.atoms)
    {
        if(!atom.target) continue;
        ++result.target_count;
        const double tau=atom.mdpde[1];
        if(atom.seed_source!="median-fallback" && !atom.reason.starts_with("initialization-exception") &&
            std::isfinite(tau) && tau>0.0 &&
            IsQualifiedInitialStageEstimate(atom.native_status,atom.initial_stage_fit_qualified,
                atom.initial_stage_qualification,atom.initial_stage_recovery))
            donors.push_back(tau);
    }

    result.qualified_donor_count=donors.size();
    result.qualified_donor_coverage=result.target_count ?
        static_cast<double>(result.qualified_donor_count)/static_cast<double>(result.target_count) : 0.0;
    if(donors.empty())
    {
        result.unavailable_reason="no-qualified-initial-stage-target-widths";
        return result;
    }
    if(donors.size()<MinimumQualifiedWidthDonors)
    {
        result.unavailable_reason="too-few-qualified-initial-stage-target-widths";
        return result;
    }
    if(result.qualified_donor_coverage<MinimumQualifiedDonorCoverage)
    {
        result.unavailable_reason="qualified-initial-stage-target-coverage-below-0.95";
        return result;
    }

    result.median_width=array_helper::ComputeMedian(donors);
    result.tau_max4=BoundFrom(*result.median_width,4.0);
    if(!result.tau_max4)
    {
        result.unavailable_reason="global-width-bound-overflow-or-invalid-log-bound";
        result.median_width.reset();
        return result;
    }
    result.eta_max4=std::log(*result.tau_max4);
    result.status="Available";
    return result;
}

BoundedInitialization MakeBoundedInitialization(
    const JointInitialization & initialization, const WidthBoundCalibration & calibration)
{
    BoundedInitialization result;
    if(initialization.atoms.size()!=initialization.tau.size())
    {
        result.reason="initialization-record-size-mismatch";
        return result;
    }
    if(calibration.status!="Available" || !calibration.tau_max4)
    {
        result.reason=calibration.unavailable_reason.empty() ? "calibration-unavailable" : calibration.unavailable_reason;
        return result;
    }
    result.tau_max=*calibration.tau_max4;
    result.tau_values=initialization.tau;
    for(std::size_t i=0;i<initialization.atoms.size();++i)
    {
        const auto & atom=initialization.atoms[i];
        if(!atom.target) continue;
        const double initial=initialization.tau[i];
        double bounded=initial;
        if(std::isfinite(initial) && initial>result.tau_max) bounded=result.tau_max;
        if(bounded!=initial)
        {
            ++result.projection_count;
            result.maximum_projection_magnitude=std::max(result.maximum_projection_magnitude,initial-bounded);
            result.tau_values[i]=bounded;
        }
    }
    result.available=true;
    return result;
}

} // namespace rhbm_gem::core::joint_component
