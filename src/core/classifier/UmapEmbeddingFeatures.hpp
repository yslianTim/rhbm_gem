#pragma once

#include "core/classifier/LocalFittingFeatures.hpp"
#include "core/classifier/UmapEmbeddingTypes.hpp"

#include <rhbm_gem/data/object/JointAnalysisResult.hpp>
#include <rhbm_gem/utils/domain/GlobalEnumClass.hpp>

#include <optional>
#include <string>
#include <string_view>

namespace rhbm_gem::core::classifier {

inline constexpr double kMadNormalConsistencyFactor{ 1.482602218505602 };

struct PreparedUmapInput
{
    std::vector<LocalFittingFeatureRow> feature_rows;
    std::vector<std::optional<std::size_t>> configured_spot_indices;
    std::vector<Element> elements;
    UmapFeatureMatrix features;
    std::vector<std::string_view> constant_features;
    int effective_neighbors{ 0 };
    std::map<std::string, std::size_t> exclusions;
    std::optional<JointMapNormalization> map_normalization;
    std::size_t excluded_rows{};
};

std::optional<PreparedUmapInput> BuildAndStandardizeInput(
    const ModelObject & model_object,
    std::vector<LocalFittingFeatureRow> source_rows,
    std::string_view model_key_tag,
    int requested_neighbors,
    std::string & error_message);

void PopulateUmapBinFeatureMedians(
    const std::vector<LocalFittingFeatureRow> & rows,
    UmapBinning & binning);

UmapFeatureMatrix BuildSecondUmapFeatures(
    const std::vector<LocalFittingFeatureRow> & rows,
    const UmapBinning & first_binning,
    std::string_view model_key_tag);

std::vector<std::string_view> GetUmapFeatureNames();
std::vector<std::string_view> GetSecondUmapFeatureNames();

} // namespace rhbm_gem::core::classifier
