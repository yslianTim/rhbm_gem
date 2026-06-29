#pragma once

namespace rhbm_gem::algorithm {

struct RobustSlopeOptions
{
    int maximum_iterations{ 50 };
    double tolerance{ 1.0e-8 };
    double scale_multiplier{ 1.4826 };
    double scale_min{ 1.0e-12 };
    double cutoff_multiplier{ 1.345 };
    double regularization_prior_scale{ 0.0 };
};

} // namespace rhbm_gem::algorithm
