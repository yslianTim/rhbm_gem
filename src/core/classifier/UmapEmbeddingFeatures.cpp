#include "core/classifier/UmapEmbeddingFeatures.hpp"
#include "core/classifier/UmapSpotCategories.hpp"

#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelAnalysisView.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <span>
#include <stdexcept>
#include <utility>

namespace rhbm_gem::core::classifier {

namespace {

struct UmapFeatureDefinition
{
    std::string_view name;
    bool include_in_umap;
};

// Keep these entries in local fitting feature order. Toggle include_in_umap to
// choose the features that are standardized and passed to UMAP, then rebuild.
constexpr std::array<UmapFeatureDefinition, kLocalFittingFeatureCount>
    kFeatureDefinitions{
        UmapFeatureDefinition{ kLocalFittingFeatureNames[0], true },
        UmapFeatureDefinition{ kLocalFittingFeatureNames[1], true },
        UmapFeatureDefinition{ kLocalFittingFeatureNames[2], false },
        UmapFeatureDefinition{ kLocalFittingFeatureNames[3], true },
        UmapFeatureDefinition{ kLocalFittingFeatureNames[4], false },
        UmapFeatureDefinition{ kLocalFittingFeatureNames[5], false },
        UmapFeatureDefinition{ kLocalFittingFeatureNames[6], false },
    };
constexpr std::size_t kInputFeatureCount{ kFeatureDefinitions.size() };
constexpr std::size_t kSelectedFeatureCount{ static_cast<std::size_t>(std::count_if(
    kFeatureDefinitions.begin(),
    kFeatureDefinitions.end(),
    [](const UmapFeatureDefinition & feature) { return feature.include_in_umap; }))
};
static_assert(
    kSelectedFeatureCount > 0,
    "At least one feature must be selected in kFeatureDefinitions for UMAP.");

constexpr std::array<std::string_view, kOutputDimensionCount> kSecondSpatialFeatureNames{
    "first bin center x (diagonal-normalized)",
    "first bin center y (diagonal-normalized)",
};
constexpr std::size_t kSecondSpatialFeatureCount{
    kIncludeSecondUmapSpatialFeatures ? kSecondSpatialFeatureNames.size() : 0
};

// Keep these entries in local fitting feature order.
constexpr std::array<UmapFeatureDefinition, kLocalFittingFeatureCount> kSecondUmapFeatureDefinitions{
    UmapFeatureDefinition{ kLocalFittingFeatureNames[0], true },
    UmapFeatureDefinition{ kLocalFittingFeatureNames[1], false },
    UmapFeatureDefinition{ kLocalFittingFeatureNames[2], true },
    UmapFeatureDefinition{ kLocalFittingFeatureNames[3], false },
    UmapFeatureDefinition{ kLocalFittingFeatureNames[4], false },
    UmapFeatureDefinition{ kLocalFittingFeatureNames[5], true },
    UmapFeatureDefinition{ kLocalFittingFeatureNames[6], true },
};
constexpr std::size_t kSecondSelectedFeatureCount{ kSecondSpatialFeatureCount + static_cast<std::size_t>(std::count_if(
    kSecondUmapFeatureDefinitions.begin(),
    kSecondUmapFeatureDefinitions.end(),
    [](const UmapFeatureDefinition & feature) { return feature.include_in_umap; }))
};
static_assert(
    kSecondSelectedFeatureCount > 0,
    "At least one feature must be selected for the second UMAP.");

std::string FormatFeatureLocation(
    std::string_view model_key_tag,
    int serial_id,
    std::string_view column_name)
{
    return "model key '" + std::string(model_key_tag)
        + "', atom serial " + std::to_string(serial_id)
        + ", column '" + std::string(column_name) + "'";
}

std::vector<std::string_view> CollectSelectedFeatureNames(
    std::span<const UmapFeatureDefinition> definitions)
{
    std::vector<std::string_view> names;
    for (const auto & feature : definitions)
    {
        if (feature.include_in_umap) names.push_back(feature.name);
    }
    return names;
}

} // namespace

std::vector<std::string_view> GetUmapFeatureNames()
{
    return CollectSelectedFeatureNames(kFeatureDefinitions);
}

std::vector<std::string_view> GetSecondUmapFeatureNames()
{
    auto names{ CollectSelectedFeatureNames(kSecondUmapFeatureDefinitions) };
    if constexpr (kIncludeSecondUmapSpatialFeatures)
        names.insert(names.begin(), kSecondSpatialFeatureNames.begin(), kSecondSpatialFeatureNames.end());
    return names;
}

std::optional<PreparedUmapInput> BuildAndStandardizeInput(
    const ModelObject & model_object,
    std::vector<LocalFittingFeatureRow> source_rows,
    std::string_view model_key_tag,
    int requested_neighbors,
    std::string & error_message)
{
    std::vector<LocalFittingFeatureRow> feature_rows;
    std::vector<std::optional<std::size_t>> configured_spot_indices;
    std::vector<Element> elements;
    feature_rows.reserve(source_rows.size());
    configured_spot_indices.reserve(source_rows.size());
    elements.reserve(source_rows.size());
    std::size_t excluded_rows{ 0 };
    std::map<std::string, std::size_t> exclusions;
    std::array<long double, kInputFeatureCount> means{};
    std::array<long double, kInputFeatureCount> sum_squared_differences{};

    for (auto & row : source_rows)
    {
        const auto configured_spot_index{
            FindConfiguredUmapSpotIndex(row.residue, row.spot)
        };
        bool valid{ true };
        for (std::size_t feature = 0; feature < kInputFeatureCount; ++feature)
        {
            if (kFeatureDefinitions[feature].include_in_umap && !std::isfinite(row.features[feature]))
            {
                std::string reason{ "unavailable: " + std::string(kFeatureDefinitions[feature].name) };
                if (feature == 3 && !row.tail_peeling_ratio_assessment.reason.empty())
                    reason += " (" + row.tail_peeling_ratio_assessment.reason + ")";
                ++exclusions[reason];
                valid = false;
                Logger::Log(LogLevel::Warning,
                    FormatFeatureLocation(model_key_tag, row.serial_id, kFeatureDefinitions[feature].name)
                        + ": excluded; required feature unavailable.");
            }
        }
        if (!valid)
        {
            ++excluded_rows;
            continue;
        }

        elements.push_back(model_object.FindAtomPtr(row.serial_id)->GetElement());
        feature_rows.emplace_back(std::move(row));
        configured_spot_indices.push_back(configured_spot_index);
        const auto & stored_row{ feature_rows.back() };
        const long double count{ static_cast<long double>(feature_rows.size()) };
        for (std::size_t feature = 0; feature < kInputFeatureCount; ++feature)
        {
            if (!kFeatureDefinitions[feature].include_in_umap) continue;
            const long double value{ static_cast<long double>(stored_row.features[feature]) };
            const long double delta{ value - means[feature] };
            means[feature] += delta / count;
            const long double updated_delta{ value - means[feature] };
            sum_squared_differences[feature] += delta * updated_delta;
        }
    }
    if (feature_rows.size() < 3)
    {
        error_message = "UMAP embedding requires at least 3 data rows; found "
            + std::to_string(feature_rows.size()) + ".";
        return std::nullopt;
    }
    if (feature_rows.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()))
    {
        error_message = "UMAP embedding input exceeds the supported number of rows.";
        return std::nullopt;
    }

    std::array<long double, kInputFeatureCount> standard_deviations{};
    std::vector<std::string_view> constant_features;
    for (std::size_t feature = 0; feature < kInputFeatureCount; ++feature)
    {
        if (!kFeatureDefinitions[feature].include_in_umap) continue;
        const long double variance{
            sum_squared_differences[feature]
                / static_cast<long double>(feature_rows.size())
        };
        if (!std::isfinite(means[feature]) || !std::isfinite(variance) || variance < 0)
        {
            error_message = "column '" + std::string(kFeatureDefinitions[feature].name)
                + "' could not be standardized to finite values.";
            return std::nullopt;
        }
        standard_deviations[feature] = std::sqrt(variance);
        if (standard_deviations[feature] == 0)
        {
            constant_features.push_back(kFeatureDefinitions[feature].name);
        }
    }
    if (constant_features.size() == kSelectedFeatureCount)
    {
        error_message = "All " + std::to_string(kSelectedFeatureCount)
            + " selected UMAP feature columns are constant; no embedding can be created.";
        return std::nullopt;
    }

    std::vector<double> standardized_features(
        feature_rows.size() * kSelectedFeatureCount);
    for (std::size_t observation = 0; observation < feature_rows.size(); ++observation)
    {
        std::size_t selected_feature{ 0 };
        for (std::size_t feature = 0; feature < kInputFeatureCount; ++feature)
        {
            if (!kFeatureDefinitions[feature].include_in_umap) continue;
            double standardized_value{ 0 };
            if (standard_deviations[feature] != 0)
            {
                const long double value{
                    (static_cast<long double>(feature_rows[observation].features[feature])
                        - means[feature])
                        / standard_deviations[feature]
                };
                standardized_value = static_cast<double>(value);
            }
            if (!std::isfinite(standardized_value))
            {
                error_message = FormatFeatureLocation(
                    model_key_tag,
                    feature_rows[observation].serial_id,
                    kFeatureDefinitions[feature].name)
                    + ": Z-score is not finite.";
                return std::nullopt;
            }
            standardized_features[
                observation * kSelectedFeatureCount + selected_feature] = standardized_value;
            ++selected_feature;
        }
    }

    PreparedUmapInput prepared;
    prepared.feature_rows = std::move(feature_rows);
    prepared.configured_spot_indices = std::move(configured_spot_indices);
    prepared.elements = std::move(elements);
    prepared.features = { std::move(standardized_features), kSelectedFeatureCount };
    prepared.constant_features = std::move(constant_features);
    prepared.effective_neighbors = std::min(
        requested_neighbors,
        static_cast<int>(prepared.feature_rows.size() - 1));
    prepared.exclusions = std::move(exclusions);
    prepared.excluded_rows = excluded_rows;
    const auto & joint{ model_object.GetAnalysisView().GetJointResult() };
    if (joint) prepared.map_normalization = joint->metadata.map_normalization;
    return prepared;
}

double MedianOfValues(std::vector<double> & values)
{
    std::sort(values.begin(), values.end());
    const auto middle{ values.size() / 2 };
    return std::midpoint(values[(values.size() - 1) / 2], values[middle]);
}

UmapFeatureMatrix BuildSecondUmapFeatures(
    const std::vector<LocalFittingFeatureRow> & rows,
    const UmapBinning & first_binning)
{
    constexpr auto feature_count{ kSecondSelectedFeatureCount };
    if (rows.empty()) return { {}, feature_count };
    std::vector<double> features(rows.size() * feature_count);
    bool has_variable_feature{};
    if constexpr (kIncludeSecondUmapSpatialFeatures)
    {
        if (first_binning.scale.diagonal != 0.0)
        {
            // Evaluate the normalized center directly in grid units: lambda * (center - origin) / R.
            constexpr double center_scale{ kSecondUmapSpatialWeight * kUmapBinWidthRatio };
            for (const auto & [index, bin] : first_binning.bins)
            {
                const double x{ center_scale * (static_cast<double>(index[0]) + 0.5) };
                const double y{ center_scale * (static_cast<double>(index[1]) + 0.5) };
                for (const auto point : bin.points)
                {
                    features[point * feature_count] = x;
                    features[point * feature_count + 1] = y;
                }
            }
            // Each occupied grid cell contributes a distinct center.
            has_variable_feature = first_binning.bins.size() > 1;
        }
    }
    std::vector<double> finite_values;
    finite_values.reserve(rows.size());
    std::size_t selected_feature{ kSecondSpatialFeatureCount };
    for (std::size_t feature = 0; feature < kSecondUmapFeatureDefinitions.size(); ++feature)
    {
        const auto & definition{ kSecondUmapFeatureDefinitions[feature] };
        if (!definition.include_in_umap) continue;
        const auto column{ selected_feature++ };
        const std::string name{ definition.name };
        finite_values.clear();
        for (std::size_t point = 0; point < rows.size(); ++point)
        {
            const double value{ rows[point].features[feature] };
            features[point * feature_count + column] = value;
            if (std::isfinite(value)) finite_values.push_back(value);
        }
        double fill{};
        if (finite_values.empty())
        {
            Logger::Log(LogLevel::Warning,
                "Second UMAP feature '" + name + "' has no finite values. Using a zero column.");
        }
        else if (finite_values.size() < rows.size())
        {
            fill = MedianOfValues(finite_values);
            Logger::Log(LogLevel::Warning, "Second UMAP feature '" + name + "' filled "
                + std::to_string(rows.size() - finite_values.size()) + " missing values with its median.");
        }

        long double mean{};
        long double sum_squared_differences{};
        for (std::size_t point = 0; point < rows.size(); ++point)
        {
            auto & value{ features[point * feature_count + column] };
            if (!std::isfinite(value)) value = fill;
            const long double delta{ value - mean };
            mean += delta / static_cast<long double>(point + 1);
            sum_squared_differences += delta * (value - mean);
        }
        const auto variance{ sum_squared_differences / static_cast<long double>(rows.size()) };
        if (!std::isfinite(variance) || variance < 0)
        {
            throw std::runtime_error("Column '" + name + "' could not be standardized to finite values.");
        }
        const auto deviation{ std::sqrt(variance) };
        has_variable_feature |= deviation > 0;
        for (std::size_t point = 0; point < rows.size(); ++point)
        {
            auto & value{ features[point * feature_count + column] };
            value = deviation > 0 ? static_cast<double>((value - mean) / deviation) : 0.0;
            if (!std::isfinite(value))
            {
                throw std::runtime_error("Column '" + name + "' has a non-finite Z-score.");
            }
        }
    }
    if (!has_variable_feature) features.clear();
    return { std::move(features), feature_count };
}

} // namespace rhbm_gem::core::classifier
