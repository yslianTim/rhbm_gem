#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>

namespace rhbm_gem::core::joint_component {

enum class WidthBoundaryStatus { Interior, LowerBound, UpperBound, Unavailable };

struct WidthBoundaryAssessment
{
    double upper_slack{}, lower_slack{};
    bool feasible{};
    WidthBoundaryStatus active{WidthBoundaryStatus::Unavailable};
};

struct WidthBounds
{
    std::optional<double> upper, lower;

    static double RoundoffTolerance(double eta)
    {
        return 64.0*std::numeric_limits<double>::epsilon()*std::max(1.0,std::abs(eta));
    }

    bool IsValid() const
    {
        return upper && std::isfinite(*upper) &&
            (!lower || (std::isfinite(*lower) && *lower < *upper));
    }

    bool Contains(double eta) const
    {
        return upper && eta <= *upper && (!lower || eta >= *lower);
    }

    bool ContainsWithRoundoff(double eta) const
    {
        return upper && *upper-eta>=-RoundoffTolerance(*upper) &&
            (!lower || eta-*lower>=-RoundoffTolerance(*lower));
    }

    WidthBoundaryAssessment Assess(double eta) const
    {
        if (!IsValid() || !std::isfinite(eta)) return {};
        const double upper_slack=*upper-eta;
        const double lower_slack=lower ? eta-*lower : std::numeric_limits<double>::infinity();
        const double upper_tolerance=RoundoffTolerance(*upper);
        const double lower_tolerance=lower ? RoundoffTolerance(*lower) : 0.0;
        if (lower_slack<=lower_tolerance && (upper_slack>upper_tolerance || lower_slack<upper_slack))
            return {upper_slack,lower_slack,ContainsWithRoundoff(eta),WidthBoundaryStatus::LowerBound};
        return {upper_slack,lower_slack,ContainsWithRoundoff(eta),upper_slack<=upper_tolerance ?
            WidthBoundaryStatus::UpperBound : WidthBoundaryStatus::Interior};
    }

    WidthBoundaryStatus Classify(double eta) const
    {
        const auto assessment=Assess(eta);
        return assessment.feasible ? assessment.active : WidthBoundaryStatus::Unavailable;
    }
};

} // namespace rhbm_gem::core::joint_component
