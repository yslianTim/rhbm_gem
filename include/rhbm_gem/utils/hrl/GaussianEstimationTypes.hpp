#pragma once

#include <optional>
#include <vector>

#include <rhbm_gem/utils/domain/SamplingTypes.hpp>
#include <rhbm_gem/utils/hrl/RHBMTypes.hpp>
#include <rhbm_gem/utils/math/GaussianModel3D.hpp>

namespace rhbm_gem {

enum class FittingStage
{
    First,
    Second,
    Third
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
    std::optional<RHBMBetaEstimateResult> fit_result{};
};

struct GroupGaussianMemberResult
{
    GaussianModel3DWithUncertainty posterior{};
    bool is_outlier{ false };
    double statistical_distance{ 0.0 };
};

struct GroupGaussianMemberInput
{
    LocalPotentialSampleList sample_entries{};
    double alpha_r{ 0.0 };
    GaussianModel3D local_model{};
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
    std::vector<GroupGaussianMemberResult> member_results{};
};

} // namespace rhbm_gem
