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
#include <tuple>
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
        UmapFeatureDefinition{ kLocalFittingFeatureNames[2], true },
        UmapFeatureDefinition{ kLocalFittingFeatureNames[3], true },
        UmapFeatureDefinition{ kLocalFittingFeatureNames[4], false },
        UmapFeatureDefinition{ kLocalFittingFeatureNames[5], false },
        UmapFeatureDefinition{ kLocalFittingFeatureNames[6], false },
        UmapFeatureDefinition{ kLocalFittingFeatureNames[7], false },
        UmapFeatureDefinition{ kLocalFittingFeatureNames[8], false },
        UmapFeatureDefinition{ kLocalFittingFeatureNames[9], false },
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
    SecondUmapFeatureDefinition{ "neighbor distance sum in 2A median of bin", SecondUmapFeatureSource::BinMedian, 0, true },
    SecondUmapFeatureDefinition{ "distance to closest neighbor median of bin", SecondUmapFeatureSource::BinMedian, 1, true },
    SecondUmapFeatureDefinition{ "signal peeling ratio median of bin", SecondUmapFeatureSource::BinMedian, 2, true },
    SecondUmapFeatureDefinition{ "tail peeling ratio median of bin", SecondUmapFeatureSource::BinMedian, 3, true },
    SecondUmapFeatureDefinition{ "amplitude joint median of bin", SecondUmapFeatureSource::BinMedian, 4, false },
    SecondUmapFeatureDefinition{ "width joint median of bin", SecondUmapFeatureSource::BinMedian, 5, false },
    SecondUmapFeatureDefinition{ "offset joint median of bin", SecondUmapFeatureSource::BinMedian, 6, false },
    SecondUmapFeatureDefinition{ "joint model center height median of bin", SecondUmapFeatureSource::BinMedian, 7, false },
    SecondUmapFeatureDefinition{ "joint signal gaussian integral median of bin", SecondUmapFeatureSource::BinMedian, 8, false },
    SecondUmapFeatureDefinition{ "joint signal charge integral median of bin", SecondUmapFeatureSource::BinMedian, 9, false },
    SecondUmapFeatureDefinition{ kLocalFittingFeatureNames[4], SecondUmapFeatureSource::Atom, 4, true },
    SecondUmapFeatureDefinition{ kLocalFittingFeatureNames[5], SecondUmapFeatureSource::Atom, 5, true },
    SecondUmapFeatureDefinition{ kLocalFittingFeatureNames[6], SecondUmapFeatureSource::Atom, 6, true },
    SecondUmapFeatureDefinition{ kLocalFittingFeatureNames[7], SecondUmapFeatureSource::Atom, 7, true },
    SecondUmapFeatureDefinition{ kLocalFittingFeatureNames[8], SecondUmapFeatureSource::Atom, 8, true },
    SecondUmapFeatureDefinition{ kLocalFittingFeatureNames[9], SecondUmapFeatureSource::Atom, 9, true },
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

double MedianOfValues(std::vector<double> & values)
{
    std::sort(values.begin(), values.end());
    const auto middle{ values.size() / 2 };
    return std::midpoint(values[(values.size() - 1) / 2], values[middle]);
}

auto FindMiddleValues(std::span<long double> values)
{
    const auto middle{ values.begin() + static_cast<std::ptrdiff_t>(values.size() / 2) };
    std::nth_element(values.begin(), middle, values.end());
    const auto lower{ values.size() % 2 == 0
        ? *std::max_element(values.begin(), middle) : *middle };
    return std::pair{ lower, *middle };
}

class MadRobustScaler
{
public:
    MadRobustScaler(
        std::span<const LocalFittingFeatureRow> rows,
        std::string_view stage,
        std::string_view model_key_tag)
        : rows_{ rows }, stage_{ stage }, model_key_tag_{ model_key_tag }, scratch_(rows.size())
    {
    }

    // Returns whether the original (or imputed) column varies.
    bool ScaleColumn(
        std::string_view name,
        std::size_t column,
        UmapFeatureMatrix & matrix,
        const auto & input_at)
    {
        const auto output_at{ [&](std::size_t row) -> double &
        {
            return matrix.values[row * matrix.column_count + column];
        } };
        const auto column_location{ [&]
        {
            return "model key '" + std::string(model_key_tag_) + "', column '" + std::string(name) + "'";
        } };
        const auto error_at{ [&](std::optional<std::size_t> row, std::string_view reason)
        {
            return std::runtime_error(std::string(stage_) + " MAD scaling for "
                + (row ? FormatFeatureLocation(model_key_tag_, rows_[*row].serial_id, name) : column_location())
                + ": " + std::string(reason) + ".");
        } };

        for (std::size_t row = 0; row < rows_.size(); ++row) scratch_[row] = input_at(row);
        const auto [minimum, maximum]{ std::minmax_element(scratch_.begin(), scratch_.end()) };
        if (*minimum == *maximum)
        {
            for (std::size_t row = 0; row < rows_.size(); ++row) output_at(row) = 0.0;
            return false;
        }

        // Normalize before subtracting. Only a subnormal (or zero) scaled value
        // can lose bits. Check those values even when long double has binary64 range.
        const int exponent{ std::ilogb(std::max(std::abs(*minimum), std::abs(*maximum))) };
        const long double multiplier{ -exponent < std::numeric_limits<long double>::max_exponent
            ? std::scalbn(1.0L, -exponent) : 0.0L };
        const auto scale{ [=](long double value)
        {
            return multiplier == 0 ? std::scalbn(value, -exponent) : value * multiplier;
        } };
        for (std::size_t row = 0; row < rows_.size(); ++row)
        {
            const auto raw{ scratch_[row] };
            const auto scaled{ scale(raw) };
            if (raw != 0 && std::abs(scaled) < std::numeric_limits<long double>::min()
                && std::scalbn(scaled, exponent) != raw)
            {
                throw error_at(row, "column exponent scaling loses a nonzero value or its precision");
            }
            scratch_[row] = scaled;
        }

        // Keep the midpoint as two endpoints plus half their gap. Rounding
        // the center to a single value can erase half-ULP differences at large offsets.
        const auto [lower, upper]{ FindMiddleValues(scratch_) };
        const auto gap{ upper - lower };
        const auto half_gap{ gap * 0.5L };
        if (half_gap * 2 != gap)
            throw error_at(std::nullopt, "median midpoint underflows in scaled coordinates");
        const auto centered{ [=](long double value)
        {
            return value <= lower ? (value - lower) - half_gap : (value - upper) + half_gap;
        } };
        for (auto & value : scratch_) value = std::abs(centered(value));

        auto [mad_lower, mad_upper]{ FindMiddleValues(scratch_) };
        const bool use_positive_deviations{ mad_upper == 0 };
        if (use_positive_deviations)
        {
            const auto positive_end{ std::remove(scratch_.begin(), scratch_.end(), 0.0L) };
            std::tie(mad_lower, mad_upper) = FindMiddleValues(
                std::span{ scratch_.begin(), positive_end });
        }
        const auto mad{ std::midpoint(mad_lower, mad_upper) };
        if (mad < std::numeric_limits<long double>::min())
            throw error_at(std::nullopt, "MAD scale underflows in scaled coordinates");
        const auto denominator{ kMadNormalConsistencyFactor * mad };

        constexpr auto max_double{ std::numeric_limits<double>::max() };
        const double min_gap{ std::sqrt(std::numeric_limits<double>::min()) };
        // Above this magnitude, even the smaller spacing at a power-of-two
        // boundary is at least min_gap.
        const double min_magnitude{ 2 * min_gap / std::numeric_limits<double>::epsilon() };
        double smallest{}, largest{};
        bool check_small_gaps{};
        for (std::size_t row = 0; row < rows_.size(); ++row)
        {
            const auto deviation{ centered(scale(input_at(row))) };
            const auto ratio{ deviation / denominator };
            if (!std::isfinite(ratio) || std::abs(ratio) > max_double)
                throw error_at(row, "standardized value exceeds the finite double range");
            const auto value{ static_cast<double>(ratio) };
            if (value == 0 && deviation != 0)
                throw error_at(row, "nonzero standardized value underflows to zero");
            output_at(row) = value;
            smallest = std::min(smallest, value);
            largest = std::max(largest, value);
            check_small_gaps |= value != 0 && std::abs(value) < min_magnitude;
        }

        const auto distance_error{ [&](std::string_view reason)
        {
            return std::runtime_error(std::string(stage_) + " Euclidean distance safety for "
                + column_location() + ": cannot certify squared distances (" + std::string(reason) + ").");
        } };
        // knncolle sums squared double differences in feature order.
        // Halve first to check the span without overflow, then round the
        // squared terms and running sum upward to bound both separate and fused operations.
        if (largest * 0.5 - smallest * 0.5 > std::sqrt(max_double) * 0.5)
            throw distance_error("bounding-box span is too large");
        constexpr auto infinity{ std::numeric_limits<double>::infinity() };
        const double span{ largest - smallest };
        const auto squared_span{ std::nextafter(span * span, infinity) };
        if (squared_span >= max_double - squared_distance_bound_)
            throw distance_error("bounding-box squared diagonal is too large");
        squared_distance_bound_ = std::nextafter(squared_distance_bound_ + squared_span, infinity);
        if (!std::isfinite(squared_distance_bound_))
            throw distance_error("bounding-box squared diagonal is too large");

        if (check_small_gaps)
        {
            // Only unusual magnitudes need a sort to check the closest distinct
            // coordinates. A subnormal square can lose precision or round to zero.
            for (std::size_t row = 0; row < rows_.size(); ++row) scratch_[row] = output_at(row);
            std::sort(scratch_.begin(), scratch_.end());
            for (std::size_t row = 1; row < scratch_.size(); ++row)
            {
                const auto delta{ scratch_[row] - scratch_[row - 1] };
                if (delta > 0 && delta < min_gap)
                    throw distance_error("a nonzero coordinate gap squares below the normal double range");
            }
        }

        if (use_positive_deviations)
        {
            matrix.zero_mad_fallback_features.push_back(name);
            Logger::Log(LogLevel::Warning, std::string(stage_) + " feature '" + std::string(name)
                + "' has zero MAD and varies; using nonzero-deviation median scaling.");
        }
        return true;
    }

private:
    std::span<const LocalFittingFeatureRow> rows_;
    std::string_view stage_;
    std::string_view model_key_tag_;
    std::vector<long double> scratch_;
    double squared_distance_bound_{};
};

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
                if (feature == 3 && !row.tail_peeling_ratio_unavailable_reason.empty())
                    reason += " (" + row.tail_peeling_ratio_unavailable_reason + ")";
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

    UmapFeatureMatrix features{
        std::vector<double>(feature_rows.size() * kSelectedFeatureCount), kSelectedFeatureCount, {} };
    std::vector<std::string_view> constant_features;
    try
    {
        MadRobustScaler scaler{ feature_rows, "First UMAP", model_key_tag };
        std::size_t column{};
        for (std::size_t feature = 0; feature < kInputFeatureCount; ++feature)
        {
            const auto & definition{ kFeatureDefinitions[feature] };
            if (!definition.include_in_umap) continue;
            if (!scaler.ScaleColumn(definition.name, column, features,
                [&](std::size_t row) { return feature_rows[row].features[feature]; }))
            {
                constant_features.push_back(definition.name);
            }
            ++column;
        }
    }
    catch (const std::runtime_error & error)
    {
        error_message = error.what();
        return std::nullopt;
    }
    if (constant_features.size() == kSelectedFeatureCount)
    {
        error_message = "All " + std::to_string(kSelectedFeatureCount)
            + " selected UMAP feature columns are constant; no embedding can be created.";
        return std::nullopt;
    }

    PreparedUmapInput prepared;
    prepared.feature_rows = std::move(feature_rows);
    prepared.configured_spot_indices = std::move(configured_spot_indices);
    prepared.elements = std::move(elements);
    prepared.features = std::move(features);
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
    const UmapBinning & first_binning,
    std::string_view model_key_tag)
{
    constexpr auto feature_count{ kSecondSelectedFeatureCount };
    constexpr auto missing{ std::numeric_limits<double>::quiet_NaN() };
    UmapFeatureMatrix features{ {}, feature_count, {} };
    if (rows.empty()) return features;
    features.values.resize(rows.size() * feature_count, missing);
    MadRobustScaler scaler{ rows, "Second UMAP", model_key_tag };
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
                    features.values[point * feature_count + column] = value;
            }
        }
        else
        {
            for (std::size_t point = 0; point < rows.size(); ++point)
                features.values[point * feature_count + column] = rows[point].features[definition.source_index];
        }

        const auto value_at{ [&](std::size_t point) -> double &
        {
            return features.values[point * feature_count + column];
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
        has_variable_feature |= scaler.ScaleColumn(definition.name, column, features, value_at);
        ++column;
    }
    if (!has_variable_feature) features.values.clear();
    return features;
}

} // namespace rhbm_gem::core::classifier
