#pragma once

#include <optional>
#include <vector>

#include <rhbm_gem/utils/hrl/RHBMTypes.hpp>
#include <rhbm_gem/utils/math/GaussianModel3D.hpp>

namespace rhbm_gem {

enum class LocalGaussianFitModel
{
    LogQuadratic,
    OffsetTauGrid
};

struct LocalPotentialAnnotation
{
    GaussianModel3DWithUncertainty gaussian{
        GaussianModel3D{ 0.0, 0.0 },
        GaussianModel3DUncertainty{}
    };
    bool is_outlier{ false };
    double statistical_distance{ 0.0 };
};

struct LocalGaussianResult
{
    double alpha_r{ 0.0 };
    GaussianModel3DWithUncertainty ols{
        GaussianModel3D{ 0.0, 0.0 },
        GaussianModel3DUncertainty{}
    };
    GaussianModel3DWithUncertainty mdpde{
        GaussianModel3D{ 0.0, 0.0 },
        GaussianModel3DUncertainty{}
    };
    bool is_outlier{ false };
    double statistical_distance{ 0.0 };
    LocalGaussianFitModel fit_model{ LocalGaussianFitModel::LogQuadratic };
    std::optional<RHBMBetaEstimateResult> fit_result{};
};

struct GroupGaussianResult
{
    double alpha_g{ 0.0 };
    GaussianModel3D mean{ 0.0, 0.0 };
    GaussianModel3D mdpde{ 0.0, 0.0 };
    GaussianModel3DWithUncertainty prior{
        GaussianModel3D{ 0.0, 0.0 },
        GaussianModel3DUncertainty{}
    };
    std::vector<LocalGaussianResult> member_results{};
};

} // namespace rhbm_gem
