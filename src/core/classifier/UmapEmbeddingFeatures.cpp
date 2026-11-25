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

enum class SecondUmapFeatureSource { BinMedian, Atom };

struct SecondUmapFeatureDefinition
{
    std::string_view name;
    SecondUmapFeatureSource source;
    std::size_t source_index;
    bool include_in_umap;
};

constexpr std::array kSecondUmapFeatureDefinitions{
    SecondUmapFeatureDefinition{ "neighbor distance sum in 2A median of bin", SecondUmapFeatureSource::BinMedian, 0, false },
    SecondUmapFeatureDefinition{ "distance to closest neighbor median of bin", SecondUmapFeatureSource::BinMedian, 1, false },
    SecondUmapFeatureDefinition{ "signal peeling ratio median of bin", SecondUmapFeatureSource::BinMedian, 2, false },
    SecondUmapFeatureDefinition{ "tail peeling ratio median of bin", SecondUmapFeatureSource::BinMedian, 3, false },
    SecondUmapFeatureDefinition{ "amplitude joint median of bin", SecondUmapFeatureSource::BinMedian, 4, true },
    SecondUmapFeatureDefinition{ "width joint median of bin", SecondUmapFeatureSource::BinMedian, 5, true },
    SecondUmapFeatureDefinition{ "offset joint median of bin", SecondUmapFeatureSource::BinMedian, 6, true },
    SecondUmapFeatureDefinition{ kLocalFittingFeatureNames[4], SecondUmapFeatureSource::Atom, 4, true },
    SecondUmapFeatureDefinition{ kLocalFittingFeatureNames[5], SecondUmapFeatureSource::Atom, 5, true },
    SecondUmapFeatureDefinition{ kLocalFittingFeatureNames[6], SecondUmapFeatureSource::Atom, 6, true },
};
constexpr std::size_t kSecondSelectedFeatureCount{ static_cast<std::size_t>(std::count_if(
    kSecondUmapFeatureDefinitions.begin(),
    kSecondUmapFeatureDefinitions.end(),
    [](const SecondUmapFeatureDefinition & feature) { return feature.include_in_umap; }))
};
static_assert(kSecondSelectedFeatureCount > 0,
    "At least one feature must be selected in kSecondUmapFeatureDefinitions for UMAP.");

std::string FormatFeatureLocation(
    std::string_view model_key_tag,
    int serial_id,
    std::string_view column_name)
{
    return "model key '" + std::string(model_key_tag)
        + "', atom serial " + std::to_string(serial_id)
        + ", column '" + std::string(column_name) + "'";
}

std::vector<std::string_view> CollectSelectedFeatureNames(const auto & definitions)
{
    std::vector<std::string_view> names;
    for (const auto & feature : definitions)
    {
        if (feature.include_in_umap) names.push_back(feature.name);
    }
    return names;
}

struct PopulationZScore
{
    int exponent;
    long double reference;
    long double mean;
    long double deviation;

    double Standardize(double value) const
    {
        if (deviation == 0) return 0.0;
        const auto scaled{ std::scalbn(static_cast<long double>(value), -exponent) };
        return static_cast<double>(((scaled - reference) - mean) / deviation);
    }
};

std::optional<PopulationZScore> CalculatePopulationZScore(
    std::size_t row_count, const auto & value_at)
{
    double minimum{ value_at(0) };
    double maximum{ minimum };
    for (std::size_t row = 1; row < row_count; ++row)
    {
        const double value{ value_at(row) };
        minimum = std::min(minimum, value);
        maximum = std::max(maximum, value);
    }
    const double magnitude{ std::max(std::abs(minimum), std::abs(maximum)) };
    const int exponent{ magnitude == 0 ? 0 : std::ilogb(magnitude) };
    const auto scale{ [exponent](double value)
    {
        return std::scalbn(static_cast<long double>(value), -exponent);
    } };
    // Scale before subtracting or squaring. Center at the point in range nearest zero
    // to retain differences among nearby values without subtracting a large offset from small ones.
    const long double reference{ scale(std::clamp(0.0, minimum, maximum)) };
    long double mean{};
    long double sum_squared_differences{};
    for (std::size_t row = 0; row < row_count; ++row)
    {
        const long double centered{ scale(value_at(row)) - reference };
        const long double delta{ centered - mean };
        mean += delta / static_cast<long double>(row + 1);
        sum_squared_differences += delta * (centered - mean);
    }
    const auto variance{ sum_squared_differences / static_cast<long double>(row_count) };
    if (!std::isfinite(variance) || variance < 0) return std::nullopt;
    return PopulationZScore{ exponent, reference, mean, std::sqrt(variance) };
}

} // namespace

std::vector<std::string_view> GetUmapFeatureNames()
{
    return CollectSelectedFeatureNames(kFeatureDefinitions);
}

std::vector<std::string_view> GetSecondUmapFeatureNames()
{
    return CollectSelectedFeatureNames(kSecondUmapFeatureDefinitions);
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

    std::array<PopulationZScore, kInputFeatureCount> standardizers{};
    std::vector<std::string_view> constant_features;
    for (std::size_t feature = 0; feature < kInputFeatureCount; ++feature)
    {
        if (!kFeatureDefinitions[feature].include_in_umap) continue;
        const auto standardizer{ CalculatePopulationZScore(feature_rows.size(),
            [&](std::size_t row) { return feature_rows[row].features[feature]; }) };
        if (!standardizer)
        {
            error_message = "column '" + std::string(kFeatureDefinitions[feature].name)
                + "' could not be standardized to finite values.";
            return std::nullopt;
        }
        standardizers[feature] = *standardizer;
        if (standardizer->deviation == 0)
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
            const double standardized_value{ standardizers[feature].Standardize(
                feature_rows[observation].features[feature]) };
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

void PopulateUmapBinFeatureMedians(
    const std::vector<LocalFittingFeatureRow> & rows,
    UmapBinning & binning)
{
    std::vector<double> values;
    for (auto & [index, bin] : binning.bins)
    {
        values.reserve(bin.points.size());
        for (std::size_t feature = 0; feature < bin.feature_medians.size(); ++feature)
        {
            values.clear();
            for (const auto point : bin.points)
            {
                const double value{ rows[point].features[feature] };
                if (std::isfinite(value)) values.push_back(value);
            }
            bin.feature_medians[feature] = values.empty() ? std::nullopt : std::optional{ MedianOfValues(values) };
        }
    }
}

UmapFeatureMatrix BuildSecondUmapFeatures(
    const std::vector<LocalFittingFeatureRow> & rows,
    const UmapBinning & first_binning)
{
    constexpr auto feature_count{ kSecondSelectedFeatureCount };
    constexpr auto missing{ std::numeric_limits<double>::quiet_NaN() };
    if (rows.empty()) return { {}, feature_count };
    std::vector<double> features(rows.size() * feature_count, missing);
    bool has_variable_feature{};
    std::vector<double> finite_values;
    finite_values.reserve(rows.size());
    std::size_t column{};
    for (const auto & definition : kSecondUmapFeatureDefinitions)
    {
        if (!definition.include_in_umap) continue;
        const std::string name{ definition.name };
        if (definition.source == SecondUmapFeatureSource::BinMedian)
        {
            for (const auto & [index, bin] : first_binning.bins)
            {
                const double value{ bin.feature_medians[definition.source_index].value_or(missing) };
                for (const auto point : bin.points)
                    features[point * feature_count + column] = value;
            }
        }
        else
        {
            for (std::size_t point = 0; point < rows.size(); ++point)
                features[point * feature_count + column] = rows[point].features[definition.source_index];
        }

        const auto value_at{ [&](std::size_t point) -> double &
        {
            return features[point * feature_count + column];
        } };
        finite_values.clear();
        for (std::size_t point = 0; point < rows.size(); ++point)
        {
            const double value{ value_at(point) };
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

        if (finite_values.size() < rows.size())
        {
            for (std::size_t point = 0; point < rows.size(); ++point)
            {
                auto & value{ value_at(point) };
                if (!std::isfinite(value)) value = fill;
            }
        }
        const auto standardizer{ CalculatePopulationZScore(rows.size(), value_at) };
        if (!standardizer)
        {
            throw std::runtime_error("Column '" + name + "' could not be standardized to finite values.");
        }
        has_variable_feature |= standardizer->deviation > 0;
        for (std::size_t point = 0; point < rows.size(); ++point)
        {
            auto & value{ value_at(point) };
            value = standardizer->Standardize(value);
            if (!std::isfinite(value))
            {
                throw std::runtime_error("Column '" + name + "' has a non-finite Z-score.");
            }
        }
        ++column;
    }
    if (!has_variable_feature) features.clear();
    return { std::move(features), feature_count };
}

} // namespace rhbm_gem::core::classifier
