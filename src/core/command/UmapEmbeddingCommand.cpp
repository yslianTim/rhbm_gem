#include "detail/CommandRunner.hpp"
#include "core/detail/LocalFittingFeatures.hpp"

#include <rhbm_gem/data/io/DataRepository.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/domain/FilePathHelper.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>

#include <umappp/umappp.hpp>

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
#include <cmath>
#include <cstddef>
#include <exception>
#include <fstream>
#include <iomanip>
#include <limits>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <system_error>
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
        UmapFeatureDefinition{ detail::kLocalFittingFeatureNames[2], true },
        UmapFeatureDefinition{ detail::kLocalFittingFeatureNames[3], true },
        UmapFeatureDefinition{ detail::kLocalFittingFeatureNames[4], true },
        UmapFeatureDefinition{ detail::kLocalFittingFeatureNames[5], true },
        UmapFeatureDefinition{ detail::kLocalFittingFeatureNames[6], false },
        UmapFeatureDefinition{ detail::kLocalFittingFeatureNames[7], true },
        UmapFeatureDefinition{ detail::kLocalFittingFeatureNames[8], true },
        UmapFeatureDefinition{ detail::kLocalFittingFeatureNames[9], false },
        UmapFeatureDefinition{ detail::kLocalFittingFeatureNames[10], false },
        UmapFeatureDefinition{ detail::kLocalFittingFeatureNames[11], true },
        UmapFeatureDefinition{ detail::kLocalFittingFeatureNames[12], false },
        UmapFeatureDefinition{ detail::kLocalFittingFeatureNames[13], false },
        UmapFeatureDefinition{ detail::kLocalFittingFeatureNames[14], true },
        UmapFeatureDefinition{ detail::kLocalFittingFeatureNames[15], false },
        UmapFeatureDefinition{ detail::kLocalFittingFeatureNames[16], false },
        UmapFeatureDefinition{ detail::kLocalFittingFeatureNames[17], true },
        UmapFeatureDefinition{ detail::kLocalFittingFeatureNames[18], false },
        UmapFeatureDefinition{ detail::kLocalFittingFeatureNames[19], false },
        UmapFeatureDefinition{ detail::kLocalFittingFeatureNames[20], true },
        UmapFeatureDefinition{ detail::kLocalFittingFeatureNames[21], false },
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

// Toggle this internal switch to restrict UMAP input to the configured spots,
// then rebuild the project.
constexpr bool kFilterUmapInputBySpot{ false };

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
    ConfiguredUmapSpot{ "HOH", "O", "HOH" },
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

#ifdef HAVE_ROOT
struct UmapSpotPlotStyle
{
    std::string_view label;
    short color;
};

// Toggle this internal switch to omit the Other graph, then rebuild the project.
constexpr bool kDrawOtherUmapSpotGraph{ true };
// Keep Other last: unmatched spot names fall back to that plot category.
constexpr std::array kUmapSpotPlotStyles{
    UmapSpotPlotStyle{ kConfiguredUmapSpots[0].label, kViolet + 1 },
    UmapSpotPlotStyle{ kConfiguredUmapSpots[1].label, kRed + 1 },
    UmapSpotPlotStyle{ kConfiguredUmapSpots[2].label, kGreen + 2 },
    UmapSpotPlotStyle{ kConfiguredUmapSpots[3].label, kAzure + 2 },
    UmapSpotPlotStyle{ kConfiguredUmapSpots[4].label, kOrange + 7 },
    UmapSpotPlotStyle{ "Other", kGray + 2 },
};
constexpr std::size_t kOtherUmapSpotPlotStyleIndex{
    kConfiguredUmapSpots.size()
};
static_assert(kUmapSpotPlotStyles.size() == kOtherUmapSpotPlotStyleIndex + 1);
#endif

struct PreparedUmapInput
{
    std::vector<detail::LocalFittingFeatureRow> feature_rows;
    std::vector<std::optional<std::size_t>> configured_spot_indices;
    std::vector<double> standardized_features;
    std::vector<std::string_view> constant_features;
    int effective_neighbors{ 0 };
    std::filesystem::path output_path;
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

std::filesystem::path BuildOutputPath(const UmapEmbeddingRequest & request)
{
    const std::string suffix{ path_helper::EnsureSanitizedTag(request.model_key_tag) };
    return request.output_dir / ("umap_embedding_" + suffix + ".csv");
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
        source_rows = detail::BuildLocalFittingFeatureRows(model_object, true);
    }
    catch (const std::exception & error)
    {
        error_message = "Failed to build local fitting features for model key '"
            + request.model_key_tag + "': " + error.what();
        return std::nullopt;
    }

    std::vector<detail::LocalFittingFeatureRow> feature_rows;
    std::vector<std::optional<std::size_t>> configured_spot_indices;
    feature_rows.reserve(source_rows.size());
    configured_spot_indices.reserve(source_rows.size());
    std::array<long double, kInputFeatureCount> means{};
    std::array<long double, kInputFeatureCount> sum_squared_differences{};

    for (auto & row : source_rows)
    {
        const auto configured_spot_index{
            FindConfiguredUmapSpotIndex(row.residue, row.spot)
        };
        if (kFilterUmapInputBySpot && !configured_spot_index) continue;

        for (std::size_t feature = 0; feature < kInputFeatureCount; ++feature)
        {
            if (!std::isfinite(row.features[feature]))
            {
                error_message = FormatFeatureLocation(
                    request.model_key_tag,
                    row.serial_id,
                    kFeatureDefinitions[feature].name)
                    + ": expected a finite value from the saved model.";
                return std::nullopt;
            }
        }

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
        if constexpr (kFilterUmapInputBySpot)
        {
            error_message = "UMAP embedding requires at least 3 data rows after spot "
                "filtering; found " + std::to_string(feature_rows.size()) + ".";
        }
        else
        {
            error_message = "UMAP embedding requires at least 3 data rows; found "
                + std::to_string(feature_rows.size()) + ".";
        }
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
    prepared.standardized_features = std::move(standardized_features);
    prepared.constant_features = std::move(constant_features);
    prepared.effective_neighbors = std::min(
        request.num_neighbors,
        static_cast<int>(prepared.feature_rows.size() - 1));
    prepared.output_path = BuildOutputPath(request);
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

bool WriteEmbedding(
    const PreparedUmapInput & prepared,
    const std::vector<double> & embedding,
    std::string & error_message)
{
    std::ofstream output{ prepared.output_path, std::ios::out | std::ios::trunc };
    if (!output.is_open())
    {
        error_message = "Failed to open UMAP output file '"
            + prepared.output_path.string() + "'.";
        return false;
    }

    output << detail::BuildLocalFittingCsvHeader() << ",umap x,umap y\n";
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
            if (detail::kLocalFittingFeatureIsIntegral[feature])
            {
                output << static_cast<long long>(row.features[feature]);
            }
            else
            {
                output << row.features[feature];
            }
        }
        output << ',' << embedding[observation * kOutputDimensionCount]
            << ',' << embedding[observation * kOutputDimensionCount + 1]
            << '\n';
    }
    output.close();
    if (!output)
    {
        error_message = "Failed while writing UMAP output file '"
            + prepared.output_path.string() + "'.";
        return false;
    }
    return true;
}

#ifdef HAVE_ROOT
std::filesystem::path BuildPlotOutputPath(const PreparedUmapInput & prepared)
{
    auto plot_path{ prepared.output_path };
    plot_path.replace_extension(".pdf");
    return plot_path;
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

#endif

bool WriteEmbeddingPlot(
    const PreparedUmapInput & prepared,
    const std::vector<double> & embedding,
    std::string & error_message)
{
#ifndef HAVE_ROOT
    (void)prepared;
    (void)embedding;
    (void)error_message;
    Logger::Log(LogLevel::Warning,
        "UMAP scatter plot skipped because ROOT support is unavailable.");
    return true;
#else
    const auto plot_path{ BuildPlotOutputPath(prepared) };
    try
    {
        std::array<std::unique_ptr<TGraphErrors>, kUmapSpotPlotStyles.size()> graphs;
        std::vector<double> plotted_x;
        std::vector<double> plotted_y;
        plotted_x.reserve(prepared.feature_rows.size());
        plotted_y.reserve(prepared.feature_rows.size());

        for (std::size_t category = 0; category < kUmapSpotPlotStyles.size(); ++category)
        {
            auto graph{ root_helper::CreateGraphErrors() };
            root_helper::SetMarkerAttribute(
                graph.get(), 20, 0.8f, kUmapSpotPlotStyles[category].color, 0.75f);
            graphs[category] = std::move(graph);
        }

        for (std::size_t observation = 0;
            observation < prepared.feature_rows.size();
            ++observation)
        {
            const auto configured_spot_index{
                prepared.configured_spot_indices[observation]
            };
            if (!configured_spot_index && !kDrawOtherUmapSpotGraph) continue;

            const std::size_t category{
                configured_spot_index.value_or(kOtherUmapSpotPlotStyleIndex)
            };
            const double x{ embedding[observation * kOutputDimensionCount] };
            const double y{ embedding[observation * kOutputDimensionCount + 1] };
            auto * const graph{ graphs[category].get() };
            const int point{ graph->GetN() };
            graph->SetPoint(point, x, y);
            graph->SetPointError(point, 0.0, 0.0);
            plotted_x.push_back(x);
            plotted_y.push_back(y);
        }
        if (plotted_x.empty())
        {
            RemoveStalePlotOutput(plot_path);
            Logger::Log(LogLevel::Warning,
                "UMAP scatter plot skipped because no configured spot categories were found.");
            return true;
        }
        if (!PreparePlotOutputPath(plot_path, error_message)) return false;

        const auto x_range{ array_helper::ComputeScalingRangeTuple(plotted_x, 0.1) };
        const auto y_range{ array_helper::ComputeScalingRangeTuple(plotted_y, 0.1) };
        auto canvas{ root_helper::CreateCanvas("umap_embedding_canvas", "", 900, 700) };
        root_helper::SetCanvasDefaultStyle(canvas.get());
        root_helper::SetPadDefaultStyle(canvas.get());
        root_helper::SetPadMarginAttribute(canvas.get(), 0.12f, 0.04f, 0.12f, 0.08f);
        root_helper::SetPadLayout(canvas.get(), 1, 1, 0, 0, 1, 1);
        canvas->cd();

        auto frame{ root_helper::CreateHist2D(
            "umap_embedding_frame",
            "UMAP embedding by spot",
            100,
            std::get<0>(x_range),
            std::get<1>(x_range),
            100,
            std::get<0>(y_range),
            std::get<1>(y_range)) };
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

        auto legend{ root_helper::CreateLegend(0.80, 0.72, 0.94, 0.90) };
        root_helper::SetLegendDefaultStyle(legend.get());
        root_helper::SetTextAttribute(legend.get(), 35.0f, 133, 12);
        for (std::size_t category = 0; category < kUmapSpotPlotStyles.size(); ++category)
        {
            if (graphs[category]->GetN() == 0) continue;
            graphs[category]->Draw("P SAME");
            legend->AddEntry(
                graphs[category].get(),
                kUmapSpotPlotStyles[category].label.data(),
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
#endif
}

bool ExecutePreparedRequest(
    CommandRunner<UmapEmbeddingRequest> & runner,
    const UmapEmbeddingRequest & request,
    const PreparedUmapInput & prepared)
{
    std::vector<double> embedding(
        prepared.feature_rows.size() * kOutputDimensionCount);
    try
    {
        knncolle::VptreeBuilder<int, double, double> neighbor_builder(
            std::make_shared<knncolle::EuclideanDistance<double, double>>());
        umappp::Options options;
        options.num_neighbors = prepared.effective_neighbors;
        options.min_dist = request.min_dist;
        if (request.num_epochs > 0) options.num_epochs = request.num_epochs;
        const auto seed{
            static_cast<umappp::RngEngine::result_type>(request.random_seed)
        };
        options.initialize_seed = seed;
        options.initialize_spectral_irlba_options.seed = seed;
        options.optimize_seed = seed;
        options.num_threads = request.job_count;
        options.num_threads_spectral = request.job_count;
        options.num_threads_optimize = 1;

        auto status{ umappp::initialize(
            static_cast<int>(kSelectedFeatureCount),
            static_cast<int>(prepared.feature_rows.size()),
            prepared.standardized_features.data(),
            neighbor_builder,
            kOutputDimensionCount,
            embedding.data(),
            options) };
        status.run(embedding.data());
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

    std::string write_error;
    if (!WriteEmbedding(prepared, embedding, write_error))
    {
        runner.AddFieldValidationError(&CommandRequestBase::output_dir, write_error);
        return false;
    }
    Logger::Log(LogLevel::Info,
        "UMAP embedding written to: " + prepared.output_path.string());

    std::string plot_error;
    if (!WriteEmbeddingPlot(prepared, embedding, plot_error))
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

            if (!model_object)
            {
                validation_runner.AddFieldValidationError(
                    &UmapEmbeddingRequest::model_key_tag,
                    "Loaded model key '" + prepared_request.model_key_tag
                        + "' did not produce a model object.");
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
            if (!prepared_input) return false;
            return ExecutePreparedRequest(runner, prepared_request, *prepared_input);
        });
}

} // namespace command_internal

} // namespace rhbm_gem::core
