#pragma once

#include <cmath>
#include <optional>

#include <rhbm_gem/utils/math/GaussianModel3D.hpp>

#include "core/detail/LocalFittingTransformedChange.hpp"

namespace rhbm_gem::core::detail {

enum class SecondStageSeedSource
{
    GroupPosterior,
    GroupPrior,
    GroupMedian,
    GlobalMedian
};

struct SecondStageSeedCandidates
{
    std::optional<GaussianModel3DWithUncertainty> group_posterior{};
    std::optional<GaussianModel3DWithUncertainty> group_prior{};
    std::optional<GaussianModel3DWithUncertainty> group_median{};
    std::optional<GaussianModel3DWithUncertainty> global_median{};
};

struct SecondStageSeedSelection
{
    SecondStageSeedSource source{ SecondStageSeedSource::GlobalMedian };
    GaussianModel3DWithUncertainty model{};
};

inline bool IsValidSecondStageGaussianModel(const GaussianModel3D & model)
{
    return std::isfinite(model.GetAmplitude()) && model.GetAmplitude() > 0.0 &&
        std::isfinite(model.GetWidth()) && model.GetWidth() > 0.0 &&
        std::isfinite(model.GetOffset()) &&
        EncodeLocalFittingTransformedCoordinates(model).has_value();
}

inline std::optional<SecondStageSeedSelection> SelectSecondStageSeed(
    const SecondStageSeedCandidates & candidates)
{
    const auto select = [](
        SecondStageSeedSource source,
        const std::optional<GaussianModel3DWithUncertainty> & candidate)
        -> std::optional<SecondStageSeedSelection>
    {
        if (!candidate.has_value() ||
            !IsValidSecondStageGaussianModel(candidate->GetModel()))
        {
            return std::nullopt;
        }
        return SecondStageSeedSelection{ source, *candidate };
    };

    if (const auto selected{
            select(SecondStageSeedSource::GroupPosterior, candidates.group_posterior) })
    {
        return selected;
    }
    if (const auto selected{
            select(SecondStageSeedSource::GroupPrior, candidates.group_prior) })
    {
        return selected;
    }
    if (const auto selected{
            select(SecondStageSeedSource::GroupMedian, candidates.group_median) })
    {
        return selected;
    }
    return select(SecondStageSeedSource::GlobalMedian, candidates.global_median);
}

} // namespace rhbm_gem::core::detail
