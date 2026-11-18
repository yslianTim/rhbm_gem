#include "detail/CommandRunner.hpp"
#include "core/command/detail/LocalFittingFeatures.hpp"
#include <boost/json.hpp>
#include <map>

#include <rhbm_gem/data/io/DataRepository.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/data/object/ModelAnalysisView.hpp>
#include <rhbm_gem/utils/domain/ChemicalDataHelper.hpp>
#include <rhbm_gem/utils/domain/FilePathHelper.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>

#include <umappp/umappp.hpp>
#include <knncolle/knncolle.hpp>

#ifdef HAVE_ROOT
#include <rhbm_gem/utils/domain/ROOTHelper.hpp>
#include <rhbm_gem/utils/math/ArrayHelper.hpp>

#include <TCanvas.h>
#include <TColor.h>
#include <TGraphErrors.h>
#include <TH2.h>
#include <TLegend.h>
#endif

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <limits>
#include <memory>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <tuple>
#include <vector>

namespace rhbm_gem::core {

namespace {

constexpr std::size_t kOutputDimensionCount{ 2 };

struct UmapFeatureDefinition
{
    std::string_view name;
    bool include_in_umap;
};

// Keep these entries in local fitting feature order. Toggle include_in_umap to
// choose the features that are standardized and passed to UMAP, then rebuild.
constexpr std::array<UmapFeatureDefinition, detail::kLocalFittingFeatureCount>
    kFeatureDefinitions{
        UmapFeatureDefinition{ detail::kLocalFittingFeatureNames[0], true },
        UmapFeatureDefinition{ detail::kLocalFittingFeatureNames[1], true },
        UmapFeatureDefinition{ detail::kLocalFittingFeatureNames[2], false },
        UmapFeatureDefinition{ detail::kLocalFittingFeatureNames[3], true },
        UmapFeatureDefinition{ detail::kLocalFittingFeatureNames[4], false },
        UmapFeatureDefinition{ detail::kLocalFittingFeatureNames[5], false },
        UmapFeatureDefinition{ detail::kLocalFittingFeatureNames[6], false },
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

// Keep bin number first, then the local fitting features.
// Toggle include_in_umap to select features for the second UMAP, then rebuild.
constexpr std::array kSecondUmapFeatureDefinitions{
    UmapFeatureDefinition{ "bin number", true },
    UmapFeatureDefinition{ detail::kLocalFittingFeatureNames[0], true },
    UmapFeatureDefinition{ detail::kLocalFittingFeatureNames[1], false },
    UmapFeatureDefinition{ detail::kLocalFittingFeatureNames[2], true },
    UmapFeatureDefinition{ detail::kLocalFittingFeatureNames[3], false },
    UmapFeatureDefinition{ detail::kLocalFittingFeatureNames[4], false },
    UmapFeatureDefinition{ detail::kLocalFittingFeatureNames[5], true },
    UmapFeatureDefinition{ detail::kLocalFittingFeatureNames[6], true },
};
constexpr std::size_t kSecondSelectedFeatureCount{ static_cast<std::size_t>(std::count_if(
    kSecondUmapFeatureDefinitions.begin(),
    kSecondUmapFeatureDefinitions.end(),
    [](const UmapFeatureDefinition & feature) { return feature.include_in_umap; }))
};
static_assert(
    kSecondSelectedFeatureCount > 0,
    "At least one feature must be selected in kSecondUmapFeatureDefinitions for the second UMAP.");

struct ConfiguredUmapSpot
{
    std::string_view label;
    std::string_view spot;
    std::optional<std::string_view> residue;
};

constexpr std::array<ConfiguredUmapSpot, 5> kConfiguredUmapSpots{
    ConfiguredUmapSpot{ "C", "C", std::nullopt },
    ConfiguredUmapSpot{ "CA", "CA", std::nullopt },
    ConfiguredUmapSpot{ "N", "N", std::nullopt },
    ConfiguredUmapSpot{ "O", "O", std::nullopt },
    ConfiguredUmapSpot{ "O (Water)", "O", "HOH" },
};

std::optional<std::size_t> FindConfiguredUmapSpotIndex(
    std::string_view residue,
    std::string_view spot)
{
    // Match residue-specific categories before the generic spot categories.
    for (std::size_t index = 0; index < kConfiguredUmapSpots.size(); ++index)
    {
        const auto & configured_spot{ kConfiguredUmapSpots[index] };
        if (configured_spot.residue
            && residue == *configured_spot.residue
            && spot == configured_spot.spot)
        {
            return index;
        }
    }
    for (std::size_t index = 0; index < kConfiguredUmapSpots.size(); ++index)
    {
        const auto & configured_spot{ kConfiguredUmapSpots[index] };
        if (!configured_spot.residue && spot == configured_spot.spot) return index;
    }
    return std::nullopt;
}

std::string GetOtherUmapSpotLabel(Element element)
{
    return element == Element::UNK
        ? "Other Unknown"
        : "Side-Chain " + std::string(ChemicalDataHelper::GetLabel(element));
}

#ifdef HAVE_ROOT
struct UmapSpotPlotStyle
{
    std::string_view label;
    short color;
    short marker;
};

// Toggle this internal switch to omit all Other element graphs, then rebuild the project.
constexpr bool kDrawOtherUmapSpotGraph{ true };
constexpr std::array kUmapSpotPlotStyles{
    UmapSpotPlotStyle{ kConfiguredUmapSpots[0].label, kViolet + 1, 24 },
    UmapSpotPlotStyle{ kConfiguredUmapSpots[1].label, kRed + 1,    25 },
    UmapSpotPlotStyle{ kConfiguredUmapSpots[2].label, kGreen + 2,  26 },
    UmapSpotPlotStyle{ kConfiguredUmapSpots[3].label, kAzure + 2,  32 },
    UmapSpotPlotStyle{ kConfiguredUmapSpots[4].label, kOrange + 7, 28 },
};
static_assert(kUmapSpotPlotStyles.size() == kConfiguredUmapSpots.size());

int GetOtherElementColor(Element element)
{
    switch (element)
    {
    case Element::CARBON:     return kRed+1;
    case Element::NITROGEN:   return kGreen+2;
    case Element::OXYGEN:     return kAzure+2;
    case Element::PHOSPHORUS: return kOrange+7;
    case Element::SULFUR:     return kMagenta;
    case Element::UNK:        return kBlack;
    default: break;
    }
    const auto hue{ static_cast<float>(std::fmod(
        ChemicalDataHelper::GetAtomicNumber(element) * 137.508, 360.0)) };
    float red{}, green{}, blue{};
    TColor::HSV2RGB(hue, 0.65f, 0.75f, red, green, blue);
    return TColor::GetColor(red, green, blue);
}

short GetOtherElementMarker(Element element)
{
    switch (element)
    {
    case Element::CARBON:     return 2;
    case Element::NITROGEN:   return 5;
    case Element::OXYGEN:     return 27;
    case Element::SULFUR:     return 30;
    default: break;
    }
    return 5; // Default marker for other elements.
}
#endif

struct PreparedUmapInput
{
    std::vector<detail::LocalFittingFeatureRow> feature_rows;
    std::vector<std::optional<std::size_t>> configured_spot_indices;
    std::vector<Element> elements;
    std::vector<double> standardized_features;
    std::vector<std::string_view> constant_features;
    int effective_neighbors{ 0 };
    std::filesystem::path output_path;
    std::map<std::string, std::size_t> exclusions;
    boost::json::value map_normalization;
    std::size_t excluded_rows{};
};

std::string FormatFeatureLocation(
    std::string_view model_key_tag,
    int serial_id,
    std::string_view column_name)
{
    return "model key '" + std::string(model_key_tag)
        + "', atom serial " + std::to_string(serial_id)
        + ", column '" + std::string(column_name) + "'";
}

std::string JoinFeatureNames(const std::vector<std::string_view> & names)
{
    std::string result;
    for (std::size_t index = 0; index < names.size(); ++index)
    {
        if (index != 0) result += ", ";
        result += names[index];
    }
    return result;
}

std::optional<PreparedUmapInput> BuildAndStandardizeInput(
    const ModelObject & model_object,
    const UmapEmbeddingRequest & request,
    std::string & error_message)
{
    std::vector<detail::LocalFittingFeatureRow> source_rows;
    try
    {
        source_rows = detail::BuildLocalFittingFeatureRows(model_object);
    }
    catch (const std::exception & error)
    {
        error_message = "Failed to build local fitting features for model key '"
            + request.model_key_tag + "': " + error.what();
        return std::nullopt;
    }

    std::vector<detail::LocalFittingFeatureRow> feature_rows;
    std::vector<std::optional<std::size_t>> configured_spot_indices;
    std::vector<Element> elements;
    feature_rows.reserve(source_rows.size());
    configured_spot_indices.reserve(source_rows.size());
    elements.reserve(source_rows.size());
    std::size_t excluded_rows=0;
    std::map<std::string,std::size_t> exclusions;
    std::array<long double, kInputFeatureCount> means{};
    std::array<long double, kInputFeatureCount> sum_squared_differences{};

    for (auto & row : source_rows)
    {
        const auto configured_spot_index{
            FindConfiguredUmapSpotIndex(row.residue, row.spot)
        };
        bool valid=true;
        for (std::size_t feature = 0; feature < kInputFeatureCount; ++feature)
        {
            if (kFeatureDefinitions[feature].include_in_umap && !std::isfinite(row.features[feature]))
            {
                std::string reason="unavailable: " + std::string(kFeatureDefinitions[feature].name);
                if (feature == 3 && !row.tail_peeling_ratio_assessment.reason.empty())
                    reason += " (" + row.tail_peeling_ratio_assessment.reason + ")";
                ++exclusions[reason]; valid=false;
                Logger::Log(LogLevel::Warning, FormatFeatureLocation(request.model_key_tag,row.serial_id,kFeatureDefinitions[feature].name) + ": excluded; required feature unavailable.");
            }
        }
        if (!valid) { ++excluded_rows; continue; }

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
                    request.model_key_tag,
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
    prepared.standardized_features = std::move(standardized_features);
    prepared.constant_features = std::move(constant_features);
    prepared.effective_neighbors = std::min(
        request.num_neighbors,
        static_cast<int>(prepared.feature_rows.size() - 1));
    prepared.output_path = request.output_dir
        / ("umap_embedding_" + path_helper::EnsureSanitizedTag(request.model_key_tag) + ".csv");
    prepared.exclusions = std::move(exclusions);
    prepared.excluded_rows=excluded_rows;
    const auto & joint=model_object.GetAnalysisView().GetJointResult();
    if(joint && joint->metadata.map_normalization)
    {
        const auto & n=*joint->metadata.map_normalization;
        prepared.map_normalization=boost::json::object{{"requested",n.requested},{"applied",n.applied},{"divisor",n.divisor}};
    }
    return prepared;
}

void NormalizeAndValidateRequest(
    CommandRunner<UmapEmbeddingRequest> & runner,
    UmapEmbeddingRequest & request)
{
    runner.RequireExistingPath(request, &UmapEmbeddingRequest::database_path);
    runner.RequireNonEmptyList(request, &UmapEmbeddingRequest::model_key_tag);

    std::error_code file_error;
    if (!request.database_path.empty()
        && std::filesystem::exists(request.database_path, file_error)
        && !std::filesystem::is_regular_file(request.database_path, file_error))
    {
        runner.AddFieldValidationError(
            &UmapEmbeddingRequest::database_path,
            "database_path must name a regular file: " + request.database_path.string());
    }
    if (request.num_neighbors < 2)
    {
        runner.AddFieldValidationError(
            &UmapEmbeddingRequest::num_neighbors,
            "num_neighbors must be at least 2.");
    }
    if (!std::isfinite(request.min_dist)
        || request.min_dist < 0
        || request.min_dist > 1)
    {
        runner.AddFieldValidationError(
            &UmapEmbeddingRequest::min_dist,
            "min_dist must be a finite value within [0, 1].");
    }
    if (request.num_epochs < 0)
    {
        runner.AddFieldValidationError(
            &UmapEmbeddingRequest::num_epochs,
            "num_epochs must be non-negative.");
    }
}

std::vector<double> RunUmap(
    const std::vector<double> & features,
    std::size_t feature_count,
    int effective_neighbors,
    const UmapEmbeddingRequest & request)
{
    const auto point_count{ static_cast<int>(features.size() / feature_count) };
    std::vector<double> embedding(static_cast<std::size_t>(point_count) * kOutputDimensionCount);
    knncolle::VptreeBuilder<int, double, double> neighbor_builder(
        std::make_shared<knncolle::EuclideanDistance<double, double>>());
    umappp::Options options;
    options.num_neighbors = effective_neighbors;
    options.min_dist = request.min_dist;
    if (request.num_epochs > 0) options.num_epochs = request.num_epochs;
    const auto seed{ static_cast<umappp::RngEngine::result_type>(request.random_seed) };
    options.initialize_seed = seed;
    options.initialize_spectral_irlba_options.seed = seed;
    options.optimize_seed = seed;
    options.num_threads = request.job_count;
    options.num_threads_spectral = request.job_count;
    options.num_threads_optimize = 1;

    auto status{ umappp::initialize(feature_count, point_count, features.data(), neighbor_builder,
        kOutputDimensionCount, embedding.data(), options) };
    status.run(embedding.data());
    return embedding;
}

std::vector<double> RunUmapMeanShift(
    const std::vector<double> & embedding,
    const std::vector<detail::LocalFittingFeatureRow> & feature_rows,
    [[maybe_unused]] int thread_size)
{
    constexpr int iteration_count{ 100 };
    constexpr int requested_neighbors{ 200 };
    constexpr double alpha{ 2.0 };
    constexpr double sigma{ 1.5 };
    const int point_count{ static_cast<int>(embedding.size() / kOutputDimensionCount) };
    const int neighbors{ std::min(requested_neighbors, point_count) };
    if (neighbors != requested_neighbors)
    {
        Logger::Log(LogLevel::Warning,
            "UMAP mean shift neighbor count was limited to "
                + std::to_string(neighbors) + " for " + std::to_string(point_count)
                + " data rows.");
    }

    knncolle::VptreeBuilder<int, double, double> builder(
        std::make_shared<knncolle::EuclideanDistance<double, double>>());
    const auto index{ builder.build_unique(knncolle::SimpleMatrix<int, double>(
        kOutputDimensionCount, point_count, embedding.data())) };
    auto shifted{ embedding };
    std::atomic<bool> nonfinite_update{ false };
#ifdef USE_OPENMP
    #pragma omp parallel num_threads(thread_size)
#endif
    {
        auto searcher{ index->initialize() };
        std::vector<int> nearest;
        for (int iteration = 0; iteration < iteration_count; ++iteration)
        {
#ifdef USE_OPENMP
            #pragma omp for schedule(static)
#endif
            for (int point = 0; point < point_count; ++point)
            {
                const auto offset{ static_cast<std::size_t>(point) * kOutputDimensionCount };
                // Coordinate queries include coincident source points, including self.
                searcher->search(shifted.data() + offset, neighbors, &nearest, nullptr);
                std::array<double, kOutputDimensionCount> mean{};
                double weight_sum{};
                for (const int neighbor : nearest)
                {
                    const auto source{ static_cast<std::size_t>(neighbor) * kOutputDimensionCount };
                    const double dx{ shifted[offset] - embedding[source] };
                    const double dy{ shifted[offset + 1] - embedding[source + 1] };
                    const double distance_squared{ dx * dx + dy * dy };
                    const double weight{ std::exp(-alpha * distance_squared / (2.0 * sigma * sigma)) };
                    weight_sum += weight;
                    for (std::size_t dimension = 0; dimension < kOutputDimensionCount; ++dimension)
                    {
                        mean[dimension] += weight * embedding[source + dimension];
                    }
                }
                if (weight_sum == 0.0)
                {
                    Logger::Log(LogLevel::Warning,
                        "UMAP mean shift has zero total neighbor weight for atom serial "
                            + std::to_string(feature_rows[static_cast<std::size_t>(point)].serial_id)
                            + "; keeping its current position.");
                    continue;
                }
                for (auto & coordinate : mean) coordinate /= weight_sum;
                if (!std::isfinite(weight_sum)
                    || !std::all_of(mean.begin(), mean.end(),
                        [](double coordinate) { return std::isfinite(coordinate); }))
                {
                    nonfinite_update.store(true, std::memory_order_relaxed);
                    continue;
                }
                for (std::size_t dimension = 0; dimension < kOutputDimensionCount; ++dimension)
                {
                    shifted[offset + dimension] = mean[dimension];
                }
            }
        }
    }
    if (nonfinite_update.load(std::memory_order_relaxed))
    {
        throw std::runtime_error("Weighted mean shift produced a non-finite update.");
    }
    return shifted;
}

bool WriteEmbedding(
    const PreparedUmapInput & prepared,
    const std::vector<double> & embedding,
    const std::vector<double> & mean_shift_embedding,
    const std::vector<std::size_t> & bin_numbers,
    const std::vector<double> & second_embedding,
    std::string & error_message)
{
    std::ofstream output{ prepared.output_path, std::ios::out | std::ios::trunc };
    if (!output.is_open())
    {
        error_message = "Failed to open UMAP output file '"
            + prepared.output_path.string() + "'.";
        return false;
    }

    output << detail::BuildLocalFittingCsvHeader()
        << ",umap x,umap y,mean shift x,mean shift y,bin number,second umap x,second umap y\n";
    output << std::setprecision(std::numeric_limits<double>::max_digits10);
    for (std::size_t observation = 0;
        observation < prepared.feature_rows.size();
        ++observation)
    {
        const auto & row{ prepared.feature_rows[observation] };
        output << row.serial_id << ',' << row.residue << ',' << row.spot;
        for (std::size_t feature = 0;
            feature < detail::kLocalFittingFeatureCount;
            ++feature)
        {
            output << ',';
            if (!std::isfinite(row.features[feature])) continue;
            output << row.features[feature];
        }
        const auto offset{ observation * kOutputDimensionCount };
        output << ',' << embedding[offset] << ',' << embedding[offset + 1]
            << ',' << mean_shift_embedding[offset] << ',' << mean_shift_embedding[offset + 1]
            << ',' << bin_numbers[observation];
        if (second_embedding.empty()) output << ",,";
        else output << ',' << second_embedding[offset] << ',' << second_embedding[offset + 1];
        output << '\n';
    }
    output.close();
    boost::json::array features;
    for(const auto & feature:kFeatureDefinitions) if(feature.include_in_umap) features.push_back(boost::json::value(feature.name));
    boost::json::array second_umap_features;
    for (const auto & feature : kSecondUmapFeatureDefinitions)
    {
        if (feature.include_in_umap) second_umap_features.push_back(boost::json::value(feature.name));
    }
    boost::json::object exclusions;
    for(const auto & [reason,count]:prepared.exclusions) exclusions[reason]=count;
    std::ofstream metadata(prepared.output_path.string()+".metadata.json");
    metadata << boost::json::serialize(boost::json::object{{"estimator","joint-components"},{"features",features},
        {"second_umap_features",second_umap_features},
        {"peeling_mode","original-grid-v1"},{"map_normalization",prepared.map_normalization},
        {"column_semantics",boost::json::object{
            {detail::kLocalFittingFeatureNames[4],"Joint Gaussian amplitude A"},
            {detail::kLocalFittingFeatureNames[5],"Joint shared Gaussian/charge-kernel width tau"},
            {detail::kLocalFittingFeatureNames[6],"Joint charge-kernel coefficient C"}}},
        {"standardization","sample-z-score; constant columns zero"},{"included_rows",prepared.feature_rows.size()},{"excluded_rows",prepared.excluded_rows},{"exclusions",exclusions}});
    metadata.close();
    if(!metadata) { error_message="Failed to write UMAP metadata."; return false; }
    if (!output)
    {
        error_message = "Failed while writing UMAP output file '"
            + prepared.output_path.string() + "'.";
        return false;
    }
    return true;
}

struct UmapBin
{
    std::vector<std::size_t> points;
    std::size_t number{};
};

struct UmapBinning
{
    std::map<std::array<std::int64_t, kOutputDimensionCount>, UmapBin> bins;
    std::vector<std::size_t> point_bin_numbers;
};

UmapBinning BinMeanShiftEmbedding(
    const PreparedUmapInput & prepared,
    const std::vector<double> & mean_shift_embedding)
{
    constexpr double bin_width{ 0.01 };
    constexpr double min_bin_index{ static_cast<double>(std::numeric_limits<std::int64_t>::min()) };
    UmapBinning result;
    result.point_bin_numbers.resize(prepared.feature_rows.size());
    for (std::size_t point = 0; point < prepared.feature_rows.size(); ++point)
    {
        std::array<std::int64_t, kOutputDimensionCount> bin_index{};
        for (std::size_t axis = 0; axis < kOutputDimensionCount; ++axis)
        {
            const double index{ std::floor(mean_shift_embedding[point * kOutputDimensionCount + axis] / bin_width) };
            if (index < min_bin_index || index >= -min_bin_index)
            {
                throw std::runtime_error("Mean shift bin index is out of range for atom serial "
                    + std::to_string(prepared.feature_rows[point].serial_id) + ".");
            }
            bin_index[axis] = static_cast<std::int64_t>(index);
        }
        result.bins[bin_index].points.push_back(point);
    }
    std::size_t bin_number{};
    for (auto & [index, bin] : result.bins)
    {
        bin.number = ++bin_number;
        for (const auto point : bin.points) result.point_bin_numbers[point] = bin.number;
    }
    return result;
}

double MedianOfValues(std::vector<double> & values)
{
    std::sort(values.begin(), values.end());
    const auto middle{ values.size() / 2 };
    return std::midpoint(values[(values.size() - 1) / 2], values[middle]);
}

std::vector<double> BuildSecondUmapFeatures(
    const std::vector<detail::LocalFittingFeatureRow> & rows,
    const std::vector<std::size_t> & bin_numbers)
{
    constexpr auto feature_count{ kSecondSelectedFeatureCount };
    std::vector<double> features(rows.size() * feature_count);
    std::vector<double> finite_values;
    finite_values.reserve(rows.size());
    bool has_variable_feature{};
    std::size_t selected_feature{};
    for (std::size_t feature = 0; feature < kSecondUmapFeatureDefinitions.size(); ++feature)
    {
        const auto & definition{ kSecondUmapFeatureDefinitions[feature] };
        if (!definition.include_in_umap) continue;
        const auto column{ selected_feature++ };
        const std::string name{ definition.name };
        finite_values.clear();
        for (std::size_t point = 0; point < rows.size(); ++point)
        {
            const double value{ feature == 0
                ? static_cast<double>(bin_numbers[point])
                : rows[point].features[feature - 1] };
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
    return features;
}

bool WriteMeanShiftBins(
    const PreparedUmapInput & prepared,
    const UmapBinning & binning,
    std::string & error_message)
{
    const auto output_path{ prepared.output_path.parent_path()
        / (prepared.output_path.stem().string() + "_mean_shift_bins.csv") };
    std::ofstream output{ output_path };
    if (!output.is_open())
    {
        error_message = "Failed to open UMAP mean shift bin output file '" + output_path.string() + "'.";
        return false;
    }
    output << "bin number,bin x,bin y,point count,amplitude joint median,width joint median,offset joint median,"
        "spot class count,1st dominate spot,1st dominate spot count,1st dominate spot count ratio,"
        "2nd dominate spot,2nd dominate spot count,2nd dominate spot count ratio,"
        "3rd dominate spot,3rd dominate spot count,3rd dominate spot count ratio\n"
        << std::setprecision(std::numeric_limits<double>::max_digits10);
    for (const auto & [index, bin] : binning.bins)
    {
        std::array<std::vector<double>, 3> joint_values;
        std::map<std::string, std::size_t> spot_counts;
        for (const auto point : bin.points)
        {
            for (std::size_t feature = 0; feature < joint_values.size(); ++feature)
            {
                const double value{ prepared.feature_rows[point].features[4 + feature] };
                if (std::isfinite(value)) joint_values[feature].push_back(value);
            }
            const auto spot_index{ prepared.configured_spot_indices[point] };
            const auto label{ spot_index ? std::string(kConfiguredUmapSpots[*spot_index].label)
                : GetOtherUmapSpotLabel(prepared.elements[point]) };
            ++spot_counts[label];
        }
        output << bin.number << ',' << index[0] << ',' << index[1] << ',' << bin.points.size();
        for (auto & values : joint_values)
        {
            output << ',';
            if (!values.empty()) output << MedianOfValues(values);
        }
        std::vector<std::pair<std::string_view, std::size_t>> ranked_spots(
            spot_counts.begin(), spot_counts.end());
        std::sort(ranked_spots.begin(), ranked_spots.end(), [](const auto & lhs, const auto & rhs)
        {
            if (lhs.second != rhs.second) return lhs.second > rhs.second;
            return lhs.first < rhs.first;
        });
        output << ',' << ranked_spots.size();
        for (std::size_t rank = 0; rank < 3; ++rank)
        {
            if (rank >= ranked_spots.size())
            {
                output << ",,0,0";
                continue;
            }
            const auto & [label, count]{ ranked_spots[rank] };
            output << ',' << label << ',' << count << ','
                << static_cast<double>(count) / static_cast<double>(bin.points.size());
        }
        output << '\n';
    }
    output.close();
    if (!output)
    {
        error_message = "Failed while writing UMAP mean shift bin output file '" + output_path.string() + "'.";
        return false;
    }
    Logger::Log(LogLevel::Info, "UMAP mean shift bin count: " + std::to_string(binning.bins.size()));
    return true;
}

#ifdef HAVE_ROOT
struct UmapPlotRanges
{
    std::tuple<double, double> x;
    std::tuple<double, double> y;
};

std::optional<UmapPlotRanges> ComputeUmapPlotRanges(
    const PreparedUmapInput & prepared,
    const std::vector<double> & embedding)
{
    std::vector<double> plotted_x;
    std::vector<double> plotted_y;
    plotted_x.reserve(prepared.feature_rows.size());
    plotted_y.reserve(prepared.feature_rows.size());
    for (std::size_t point = 0; point < prepared.feature_rows.size(); ++point)
    {
        if (!prepared.configured_spot_indices[point] && !kDrawOtherUmapSpotGraph) continue;
        plotted_x.push_back(embedding[point * kOutputDimensionCount]);
        plotted_y.push_back(embedding[point * kOutputDimensionCount + 1]);
    }
    if (plotted_x.empty()) return std::nullopt;
    return UmapPlotRanges{
        array_helper::ComputeScalingRangeTuple(plotted_x, 0.1),
        array_helper::ComputeScalingRangeTuple(plotted_y, 0.1)
    };
}

bool PreparePlotOutputPath(
    const std::filesystem::path & plot_path,
    std::string & error_message)
{
    std::error_code file_error;
    const bool output_exists{ std::filesystem::exists(plot_path, file_error) };
    if (file_error)
    {
        error_message = "Failed to inspect UMAP plot output path '"
            + plot_path.string() + "': " + file_error.message() + ".";
        return false;
    }
    if (!output_exists) return true;

    if (!std::filesystem::is_regular_file(plot_path, file_error) || file_error)
    {
        error_message = "UMAP plot output path is not a regular file: "
            + plot_path.string() + ".";
        return false;
    }
    if (!std::filesystem::remove(plot_path, file_error) || file_error)
    {
        error_message = "Failed to replace UMAP plot output file '"
            + plot_path.string() + "': " + file_error.message() + ".";
        return false;
    }
    return true;
}

void RemoveStalePlotOutput(const std::filesystem::path & plot_path)
{
    std::error_code file_error;
    if (!std::filesystem::is_regular_file(plot_path, file_error)) return;
    std::filesystem::remove(plot_path, file_error);
    if (file_error)
    {
        Logger::Log(LogLevel::Warning,
            "Failed to remove stale UMAP plot output file '"
                + plot_path.string() + "': " + file_error.message() + ".");
    }
}

bool WriteEmbeddingPlot(
    const PreparedUmapInput & prepared,
    const std::vector<double> & embedding,
    const std::filesystem::path & plot_path,
    const UmapPlotRanges & ranges,
    std::string & error_message)
{
    try
    {
        std::vector<Element> other_elements;
        if (kDrawOtherUmapSpotGraph)
        {
            for (std::size_t observation = 0; observation < prepared.feature_rows.size(); ++observation)
            {
                if (!prepared.configured_spot_indices[observation])
                {
                    other_elements.push_back(prepared.elements[observation]);
                }
            }
            std::sort(other_elements.begin(), other_elements.end(), [](Element lhs, Element rhs)
            {
                if (lhs == Element::UNK) return false;
                if (rhs == Element::UNK) return true;
                return ChemicalDataHelper::GetAtomicNumber(lhs) < ChemicalDataHelper::GetAtomicNumber(rhs);
            });
            other_elements.erase(
                std::unique(other_elements.begin(), other_elements.end()),
                other_elements.end());
        }
        std::vector<std::unique_ptr<TGraphErrors>> graphs;
        std::vector<std::string> labels;

        for (const auto & style : kUmapSpotPlotStyles)
        {
            auto graph{ root_helper::CreateGraphErrors() };
            root_helper::SetMarkerAttribute(graph.get(), style.marker, 1.5f, style.color, 0.75f);
            graphs.emplace_back(std::move(graph));
            labels.emplace_back(style.label);
        }
        for (const auto element : other_elements)
        {
            auto graph{ root_helper::CreateGraphErrors() };
            root_helper::SetMarkerAttribute(
                graph.get(), GetOtherElementMarker(element), 1.5f,
                static_cast<short>(GetOtherElementColor(element)), 0.75f);
            graphs.emplace_back(std::move(graph));
            labels.emplace_back(GetOtherUmapSpotLabel(element));
        }

        for (std::size_t observation = 0;
            observation < prepared.feature_rows.size();
            ++observation)
        {
            const auto configured_spot_index{
                prepared.configured_spot_indices[observation]
            };
            if (!configured_spot_index && !kDrawOtherUmapSpotGraph) continue;

            const std::size_t category{ configured_spot_index ? *configured_spot_index
                : kUmapSpotPlotStyles.size() + static_cast<std::size_t>(std::distance(
                    other_elements.begin(),
                    std::find(other_elements.begin(), other_elements.end(), prepared.elements[observation]))) };
            const double x{ embedding[observation * kOutputDimensionCount] };
            const double y{ embedding[observation * kOutputDimensionCount + 1] };
            auto * const graph{ graphs[category].get() };
            const int point{ graph->GetN() };
            graph->SetPoint(point, x, y);
            graph->SetPointError(point, 0.0, 0.0);
        }
        if (!PreparePlotOutputPath(plot_path, error_message)) return false;

        const auto visible_categories{ static_cast<int>(std::count_if(
            graphs.begin(), graphs.end(), [](const auto & graph) { return graph->GetN() > 0; })) };
        const int legend_columns{ (visible_categories + 19) / 20 };
        const int canvas_width{ 900 + 230 * legend_columns };
        const float plot_right{ 900.0f / static_cast<float>(canvas_width) };
        auto canvas{ root_helper::CreateCanvas("umap_embedding_canvas", "", canvas_width, 700) };
        root_helper::SetCanvasDefaultStyle(canvas.get());
        root_helper::SetPadDefaultStyle(canvas.get());
        root_helper::SetPadMarginAttribute(
            canvas.get(), 108.0f / static_cast<float>(canvas_width),
            1.0f - plot_right, 0.16f, 0.08f);
        root_helper::SetPadLayout(canvas.get(), 1, 1, 0, 0, 1, 1);
        canvas->cd();

        auto frame{
            root_helper::CreateHist2D("umap_embedding_frame","",
                100, std::get<0>(ranges.x), std::get<1>(ranges.x),
                100, std::get<0>(ranges.y), std::get<1>(ranges.y))
        };
        frame->SetStats(0);
        frame->GetXaxis()->SetTitle("umap x");
        frame->GetYaxis()->SetTitle("umap y");
        frame->GetXaxis()->CenterTitle();
        frame->GetYaxis()->CenterTitle();
        root_helper::SetAxisTitleAttribute(frame->GetXaxis(), 45.0f, 1.1f);
        root_helper::SetAxisTitleAttribute(frame->GetYaxis(), 45.0f, 1.2f);
        root_helper::SetAxisLabelAttribute(frame->GetXaxis(), 40.0f, 0.01f);
        root_helper::SetAxisLabelAttribute(frame->GetYaxis(), 40.0f, 0.01f);
        frame->Draw();

        const int legend_rows{ (visible_categories + legend_columns - 1) / legend_columns };
        auto legend{ root_helper::CreateLegend(
            plot_right + 10.0 / canvas_width, 0.90 - 0.035 * legend_rows,
            1.0 - 10.0 / canvas_width, 0.90) };
        legend->SetNColumns(legend_columns);
        root_helper::SetLegendDefaultStyle(legend.get());
        root_helper::SetTextAttribute(legend.get(), 18.0f, 133, 12);
        for (std::size_t category = 0; category < graphs.size(); ++category)
        {
            if (graphs[category]->GetN() == 0) continue;
            graphs[category]->Draw("P SAME");
            legend->AddEntry(
                graphs[category].get(),
                labels[category].c_str(),
                "p");
        }
        legend->Draw();
        root_helper::PrintCanvasPad(canvas.get(), plot_path.string());

        std::error_code file_error;
        const bool is_regular_file{
            std::filesystem::is_regular_file(plot_path, file_error)
        };
        const auto file_size{
            is_regular_file ? std::filesystem::file_size(plot_path, file_error) : 0
        };
        if (file_error || !is_regular_file || file_size == 0)
        {
            error_message = "Failed to write UMAP plot output file '"
                + plot_path.string() + "'.";
            return false;
        }

        Logger::Log(LogLevel::Info,
            "UMAP scatter plot written to: " + plot_path.string());
        return true;
    }
    catch (const std::exception & error)
    {
        error_message = "Failed to create UMAP scatter plot '"
            + plot_path.string() + "': " + error.what();
        return false;
    }
}
#endif

bool WriteEmbeddingPlots(
    const PreparedUmapInput & prepared,
    const std::vector<double> & embedding,
    const std::vector<double> & mean_shift_embedding,
    const std::vector<double> & second_embedding,
    std::string & error_message)
{
#ifndef HAVE_ROOT
    (void)prepared;
    (void)embedding;
    (void)mean_shift_embedding;
    (void)second_embedding;
    (void)error_message;
    Logger::Log(LogLevel::Warning,
        "UMAP scatter plots skipped because ROOT support is unavailable.");
    return true;
#else
    auto plot_path{ prepared.output_path };
    plot_path.replace_extension(".pdf");
    const auto mean_shift_path{ prepared.output_path.parent_path()
        / (prepared.output_path.stem().string() + "_mean_shift.pdf") };
    const auto second_path{ prepared.output_path.parent_path()
        / (prepared.output_path.stem().string() + "_second_umap.pdf") };
    const auto ranges{ ComputeUmapPlotRanges(prepared, embedding) };
    if (!ranges)
    {
        RemoveStalePlotOutput(plot_path);
        RemoveStalePlotOutput(mean_shift_path);
        RemoveStalePlotOutput(second_path);
        Logger::Log(LogLevel::Warning,
            "UMAP scatter plots skipped because no configured spot categories were found.");
        return true;
    }
    if (second_embedding.empty()) RemoveStalePlotOutput(second_path);
    if (!WriteEmbeddingPlot(prepared, embedding, plot_path, *ranges, error_message)
        || !WriteEmbeddingPlot(prepared, mean_shift_embedding, mean_shift_path, *ranges, error_message))
    {
        return false;
    }
    if (second_embedding.empty()) return true;

    const auto second_ranges{ ComputeUmapPlotRanges(prepared, second_embedding) };
    return WriteEmbeddingPlot(prepared, second_embedding, second_path, *second_ranges, error_message);
#endif
}

bool ExecutePreparedRequest(
    CommandRunner<UmapEmbeddingRequest> & runner,
    const UmapEmbeddingRequest & request,
    const PreparedUmapInput & prepared)
{
    std::vector<double> embedding;
    try
    {
        embedding = RunUmap(prepared.standardized_features, kSelectedFeatureCount,
            prepared.effective_neighbors, request);
    }
    catch (const std::exception & error)
    {
        runner.RequirePrepareCondition(false,
            "UMAP embedding failed for model key '" + request.model_key_tag
                + "' in database '" + request.database_path.string()
                + "': " + error.what());
        return false;
    }

    for (const double coordinate : embedding)
    {
        if (std::isfinite(coordinate)) continue;
        runner.RequirePrepareCondition(false,
            "UMAP produced a non-finite embedding coordinate for model key '"
                + request.model_key_tag + "' in database '"
                + request.database_path.string() + "'.");
        return false;
    }

    std::vector<double> mean_shift_embedding;
    try
    {
        mean_shift_embedding = RunUmapMeanShift(embedding, prepared.feature_rows, request.job_count);
    }
    catch (const std::exception & error)
    {
        runner.RequirePrepareCondition(false,
            "UMAP mean shift failed for model key '" + request.model_key_tag
                + "' in database '" + request.database_path.string()
                + "': " + error.what());
        return false;
    }
    if (!std::all_of(mean_shift_embedding.begin(), mean_shift_embedding.end(),
        [](double coordinate) { return std::isfinite(coordinate); }))
    {
        runner.RequirePrepareCondition(false,
            "UMAP mean shift produced a non-finite coordinate for model key '"
                + request.model_key_tag + "' in database '"
                + request.database_path.string() + "'.");
        return false;
    }

    UmapBinning binning;
    std::vector<double> second_embedding;
    try
    {
        binning = BinMeanShiftEmbedding(prepared, mean_shift_embedding);
        const auto features{ BuildSecondUmapFeatures(
            prepared.feature_rows, binning.point_bin_numbers) };
        if (features.empty())
        {
            Logger::Log(LogLevel::Warning,
                "Second UMAP skipped because all selected feature columns are constant.");
        }
        else
        {
            second_embedding = RunUmap(features, kSecondSelectedFeatureCount,
                prepared.effective_neighbors, request);
            if (!std::all_of(second_embedding.begin(), second_embedding.end(),
                [](double coordinate) { return std::isfinite(coordinate); }))
            {
                throw std::runtime_error("Second UMAP produced a non-finite embedding coordinate.");
            }
        }
    }
    catch (const std::exception & error)
    {
        runner.RequirePrepareCondition(false,
            "UMAP binning or second embedding failed for model key '" + request.model_key_tag
                + "' in database '" + request.database_path.string() + "': " + error.what());
        return false;
    }

    std::string write_error;
    if (!WriteEmbedding(prepared, embedding, mean_shift_embedding, binning.point_bin_numbers,
        second_embedding, write_error))
    {
        runner.AddFieldValidationError(&CommandRequestBase::output_dir, write_error);
        return false;
    }
    Logger::Log(LogLevel::Info,
        "UMAP embedding written to: " + prepared.output_path.string());

    if (!WriteMeanShiftBins(prepared, binning, write_error))
    {
        runner.AddFieldValidationError(&CommandRequestBase::output_dir, write_error);
        return false;
    }

    std::string plot_error;
    if (!WriteEmbeddingPlots(prepared, embedding, mean_shift_embedding, second_embedding, plot_error))
    {
        runner.AddFieldValidationError(&CommandRequestBase::output_dir, plot_error);
        return false;
    }
    return true;
}

} // namespace

namespace command_internal {

CommandResult ExecuteUmapEmbeddingCommand(const UmapEmbeddingRequest & request)
{
    CommandRunner<UmapEmbeddingRequest> runner;
    std::optional<PreparedUmapInput> prepared_input;
    return runner.Run(
        request,
        NormalizeAndValidateRequest,
        [&prepared_input](
            CommandRunner<UmapEmbeddingRequest> & validation_runner,
            const UmapEmbeddingRequest & prepared_request)
        {
            std::unique_ptr<ModelObject> model_object;
            try
            {
                DataRepository repository{ prepared_request.database_path };
                try
                {
                    model_object = repository.LoadModel(prepared_request.model_key_tag);
                }
                catch (const std::exception & error)
                {
                    validation_runner.AddFieldValidationError(
                        &UmapEmbeddingRequest::model_key_tag,
                        "Failed to load model key '" + prepared_request.model_key_tag
                            + "': " + error.what());
                    return;
                }
            }
            catch (const std::exception & error)
            {
                validation_runner.AddFieldValidationError(
                    &UmapEmbeddingRequest::database_path,
                    "Failed to open or validate database '"
                        + prepared_request.database_path.string() + "': " + error.what());
                return;
            }

            std::string preparation_error;
            prepared_input = BuildAndStandardizeInput(
                *model_object,
                prepared_request,
                preparation_error);
            if (!prepared_input)
            {
                validation_runner.AddFieldValidationError(
                    &UmapEmbeddingRequest::model_key_tag,
                    preparation_error);
                return;
            }

            if (!prepared_input->constant_features.empty())
            {
                validation_runner.AddFieldNormalizationWarning(
                    &UmapEmbeddingRequest::model_key_tag,
                    "Constant feature columns were standardized to zero: "
                        + JoinFeatureNames(prepared_input->constant_features) + ".");
            }
            if (prepared_input->effective_neighbors != prepared_request.num_neighbors)
            {
                validation_runner.AddFieldNormalizationWarning(
                    &UmapEmbeddingRequest::num_neighbors,
                    "num_neighbors was limited to "
                        + std::to_string(prepared_input->effective_neighbors)
                        + " for " + std::to_string(prepared_input->feature_rows.size())
                        + " data rows.");
            }
        },
        [&runner, &prepared_input](const UmapEmbeddingRequest & prepared_request)
        {
            return ExecutePreparedRequest(runner, prepared_request, *prepared_input);
        });
}

} // namespace command_internal

} // namespace rhbm_gem::core
