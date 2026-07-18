#pragma once

#include <cmath>
#include <optional>

#include <rhbm_gem/utils/math/GaussianModel3D.hpp>

#include "core/detail/LocalFittingTransformedChange.hpp"

namespace rhbm_gem::core::detail {

enum class SecondStageSeedRepairSource
{
    GroupPosterior,
    GroupPrior,
    LocalOls,
    GroupMedian,
    GlobalMedian
};

struct SecondStageSeedRepairCandidates
{
    std::optional<GaussianModel3DWithUncertainty> group_posterior{};
    std::optional<GaussianModel3DWithUncertainty> group_prior{};
    std::optional<GaussianModel3DWithUncertainty> local_ols{};
    std::optional<GaussianModel3DWithUncertainty> group_median{};
    std::optional<GaussianModel3DWithUncertainty> global_median{};
};

struct SecondStageSeedRepairSelection
{
    SecondStageSeedRepairSource source{ SecondStageSeedRepairSource::GlobalMedian };
    GaussianModel3DWithUncertainty model{};
};

inline bool IsValidSecondStageGaussianModel(const GaussianModel3D & model)
{
    return std::isfinite(model.GetAmplitude()) && model.GetAmplitude() > 0.0 &&
        std::isfinite(model.GetWidth()) && model.GetWidth() > 0.0 &&
        std::isfinite(model.GetOffset()) &&
        EncodeLocalFittingTransformedCoordinates(model).has_value();
}

inline std::optional<SecondStageSeedRepairSelection> SelectSecondStageSeedRepair(
    const SecondStageSeedRepairCandidates & candidates)
{
    const auto select = [](
        SecondStageSeedRepairSource source,
        const std::optional<GaussianModel3DWithUncertainty> & candidate)
        -> std::optional<SecondStageSeedRepairSelection>
    {
        if (!candidate.has_value() ||
            !IsValidSecondStageGaussianModel(candidate->GetModel()))
        {
            return std::nullopt;
        }
        return SecondStageSeedRepairSelection{ source, *candidate };
    };

    if (const auto selected{
            select(SecondStageSeedRepairSource::GroupPosterior, candidates.group_posterior) })
    {
        return selected;
    }
    if (const auto selected{
            select(SecondStageSeedRepairSource::GroupPrior, candidates.group_prior) })
    {
        return selected;
    }
    if (const auto selected{
            select(SecondStageSeedRepairSource::LocalOls, candidates.local_ols) })
    {
        return selected;
    }
    if (const auto selected{
            select(SecondStageSeedRepairSource::GroupMedian, candidates.group_median) })
    {
        return selected;
    }
    return select(SecondStageSeedRepairSource::GlobalMedian, candidates.global_median);
}

inline GaussianModel3DWithUncertainty BuildRepairedSecondStageSeed(
    const GaussianModel3D & original_model,
    const SecondStageSeedRepairSelection & selection)
{
    const auto & fallback_model{ selection.model.GetModel() };
    return GaussianModel3DWithUncertainty{
        fallback_model.WithOffset(
            std::isfinite(original_model.GetOffset()) ?
                original_model.GetOffset() : fallback_model.GetOffset()),
        selection.model.GetStandardDeviationModel()
    };
}

} // namespace rhbm_gem::core::detail
