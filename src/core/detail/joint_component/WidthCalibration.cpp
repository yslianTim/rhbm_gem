#include "WidthCalibration.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>

namespace rhbm_gem::core::joint_component {
namespace {
constexpr std::size_t MinimumQualifiedWidthDonors=20;
constexpr double MinimumQualifiedDonorCoverage=0.95;

double QuantileType7(const std::vector<double> & sorted, double probability)
{
    const double position=static_cast<double>(sorted.size()-1)*probability;
    const auto lower=static_cast<std::size_t>(std::floor(position));
    const auto upper=static_cast<std::size_t>(std::ceil(position));
    const double fraction=position-static_cast<double>(lower);
    return sorted[lower]+(sorted[upper]-sorted[lower])*fraction;
}

bool HasInitializationFailure(const JointInitializationAtom & atom)
{
    return atom.reason.rfind("initialization-exception",0)==0;
}

std::optional<double> BoundFrom(double p95, double multiplier)
{
    if (!(p95>0.0) || !std::isfinite(p95) || p95>std::numeric_limits<double>::max()/multiplier)
        return std::nullopt;
    const double bound=multiplier*p95;
    if (!(bound>0.0) || !std::isfinite(bound) || !std::isfinite(std::log(bound)))
        return std::nullopt;
    return bound;
}
}

WidthBoundCalibration CalibrateGlobalWidth(
    const JointInitialization & initialization, std::string dataset_identity)
{
    WidthBoundCalibration result;
    result.dataset_identity=std::move(dataset_identity);
    if(initialization.atoms.size()!=initialization.b.size())
    {
        result.unavailable_reason="initialization-record-size-mismatch";
        return result;
    }

    for(std::size_t i=0;i<initialization.atoms.size();++i)
        if(initialization.atoms[i].target &&
            (!std::isfinite(initialization.b[i]) || !(initialization.b[i]>0.0)))
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
        const auto qualification=FirstStageQualification(atom.native_status,
            atom.first_stage_fit_qualified,atom.first_stage_qualification,atom.first_stage_recovery);
        const bool qualified=IsQualifiedFirstStageEstimate(atom.native_status,
            atom.first_stage_fit_qualified,atom.first_stage_qualification,atom.first_stage_recovery);
        if(qualified && qualification==RHBMSolveQualification::NativeSuccess) ++result.native_success_count;
        else if(qualified && qualification==RHBMSolveQualification::RefinedSuccess) ++result.refined_success_count;
        else ++result.unqualified_count;

        const double width=atom.mdpde[1];
        const bool valid_width=std::isfinite(width) && width>0.0;
        if(!valid_width) ++result.invalid_width_count;
        const bool fallback=atom.seed_source=="median-fallback";
        if(fallback) ++result.fallback_count;
        const bool initialization_failure=HasInitializationFailure(atom);
        if(initialization_failure) ++result.initialization_failure_count;
        std::string exclusion;
        if(fallback) exclusion="median-fallback";
        else if(initialization_failure) exclusion="initialization-failure";
        else if(!qualified) exclusion="unqualified-first-stage";
        else if(!std::isfinite(width)) exclusion="non-finite-width";
        else if(!(width>0.0)) exclusion="non-positive-width";
        if(exclusion.empty()) donors.push_back(width);
        else ++result.excluded_donor_reasons[exclusion];
    }

    result.qualified_donor_count=donors.size();
    result.qualified_donor_coverage=result.target_count ?
        static_cast<double>(result.qualified_donor_count)/static_cast<double>(result.target_count) : 0.0;
    if(donors.empty())
    {
        result.unavailable_reason="no-qualified-first-stage-target-widths";
        return result;
    }
    if(donors.size()<MinimumQualifiedWidthDonors)
    {
        result.unavailable_reason="too-few-qualified-first-stage-target-widths";
        return result;
    }
    if(result.qualified_donor_coverage<MinimumQualifiedDonorCoverage)
    {
        result.unavailable_reason="qualified-first-stage-target-coverage-below-0.95";
        return result;
    }

    std::sort(donors.begin(),donors.end());
    result.minimum_width=donors.front();
    result.maximum_width=donors.back();
    result.median_width=QuantileType7(donors,0.5);
    result.p90=QuantileType7(donors,0.90);
    result.p95=QuantileType7(donors,0.95);
    result.p99=QuantileType7(donors,0.99);
    result.bmax2=BoundFrom(*result.p95,2.0);
    result.bmax4=BoundFrom(*result.p95,4.0);
    if(!result.bmax2 || !result.bmax4)
    {
        result.unavailable_reason="global-width-bound-overflow-or-invalid-log-bound";
        result.p95.reset();
        result.p99.reset();
        result.median_width.reset();
        result.p90.reset();
        result.minimum_width.reset();
        result.maximum_width.reset();
        return result;
    }
    result.eta_max2=std::log(*result.bmax2);
    result.eta_max4=std::log(*result.bmax4);
    for(std::size_t i=0;i<initialization.atoms.size();++i)
    {
        const auto & atom=initialization.atoms[i];
        if(!atom.target || !std::isfinite(initialization.b[i]) || !(initialization.b[i]>0.0)) continue;
        if(initialization.b[i]>*result.bmax2) ++result.initial_violations2;
        if(initialization.b[i]>*result.bmax4) ++result.initial_violations4;
    }
    result.status="Available";
    return result;
}

BoundedInitialization MakeBoundedInitialization(
    const JointInitialization & initialization, const WidthBoundCalibration & calibration,
    GlobalWidthBoundMultiplier multiplier)
{
    BoundedInitialization result;
    if(initialization.atoms.size()!=initialization.b.size())
    {
        result.reason="initialization-record-size-mismatch";
        return result;
    }
    const auto selected_bound=multiplier==GlobalWidthBoundMultiplier::TwoP95 ? calibration.bmax2 : calibration.bmax4;
    if(calibration.status!="Available" || !selected_bound)
    {
        result.reason=calibration.unavailable_reason.empty() ? "calibration-unavailable" : calibration.unavailable_reason;
        return result;
    }
    result.bmax=*selected_bound;
    result.widths=initialization.b;
    result.projections.reserve(calibration.target_count);
    for(std::size_t i=0;i<initialization.atoms.size();++i)
    {
        const auto & atom=initialization.atoms[i];
        if(!atom.target) continue;
        const double initial=initialization.b[i];
        double bounded=initial;
        if(std::isfinite(initial) && initial>result.bmax) bounded=result.bmax;
        const bool projected=bounded!=initial;
        if(projected)
        {
            ++result.projection_count;
            result.maximum_projection_magnitude=std::max(result.maximum_projection_magnitude,initial-bounded);
            result.widths[i]=bounded;
        }
        result.projections.push_back({i,atom.id,atom.original_b,initial,bounded,projected});
    }
    result.available=true;
    return result;
}

} // namespace rhbm_gem::core::joint_component
