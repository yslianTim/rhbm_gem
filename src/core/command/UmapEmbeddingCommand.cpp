#include "detail/CommandRunner.hpp"
#include "core/classifier/UmapEmbeddingAlgorithm.hpp"
#include "core/classifier/UmapEmbeddingFeatures.hpp"
#include "core/classifier/UmapSpotCategories.hpp"
#include <boost/json.hpp>
#include <map>

#include <rhbm_gem/data/io/DataRepository.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/data/object/ModelAnalysisView.hpp>
#include <rhbm_gem/utils/domain/ChemicalDataHelper.hpp>
#include <rhbm_gem/utils/domain/FilePathHelper.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>

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
#include <iterator>
#include <limits>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <tuple>
#include <vector>

namespace rhbm_gem::core {

namespace {

using namespace classifier;

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

    output << classifier::BuildLocalFittingCsvHeader()
        << ",umap x,umap y,mean shift x,mean shift y,bin number,second umap x,second umap y\n";
    output << std::setprecision(std::numeric_limits<double>::max_digits10);
    for (std::size_t observation = 0;
        observation < prepared.feature_rows.size();
        ++observation)
    {
        const auto & row{ prepared.feature_rows[observation] };
        output << row.serial_id << ',' << row.residue << ',' << row.spot;
        for (std::size_t feature = 0;
            feature < classifier::kLocalFittingFeatureCount;
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
    for (const auto name : GetUmapFeatureNames()) features.push_back(boost::json::value(name));
    boost::json::array second_umap_features;
    for (const auto name : GetSecondUmapFeatureNames()) second_umap_features.push_back(boost::json::value(name));
    boost::json::value map_normalization;
    if (const auto & normalization{ prepared.map_normalization }; normalization)
    {
        map_normalization = boost::json::object{
            { "requested", normalization->requested },
            { "applied", normalization->applied },
            { "divisor", normalization->divisor }
        };
    }
    boost::json::object exclusions;
    for(const auto & [reason,count]:prepared.exclusions) exclusions[reason]=count;
    std::ofstream metadata(prepared.output_path.string()+".metadata.json");
    metadata << boost::json::serialize(boost::json::object{{"estimator","joint-components"},{"features",features},
        {"second_umap_features",second_umap_features},
        {"peeling_mode","original-grid-v1"},{"map_normalization",map_normalization},
        {"column_semantics",boost::json::object{
            {classifier::kLocalFittingFeatureNames[4],"Joint Gaussian amplitude A"},
            {classifier::kLocalFittingFeatureNames[5],"Joint shared Gaussian/charge-kernel width tau"},
            {classifier::kLocalFittingFeatureNames[6],"Joint charge-kernel coefficient C"}}},
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

bool WriteMeanShiftBins(
    const PreparedUmapInput & prepared,
    const UmapBinning & binning,
    const std::filesystem::path & output_path,
    std::string & error_message)
{
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
    const std::vector<double> & second_mean_shift_embedding,
    std::string & error_message)
{
#ifndef HAVE_ROOT
    (void)prepared;
    (void)embedding;
    (void)mean_shift_embedding;
    (void)second_embedding;
    (void)second_mean_shift_embedding;
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
    const auto second_mean_shift_path{ prepared.output_path.parent_path()
        / (prepared.output_path.stem().string() + "_second_umap_mean_shift.pdf") };
    const auto ranges{ ComputeUmapPlotRanges(prepared, embedding) };
    if (!ranges)
    {
        RemoveStalePlotOutput(plot_path);
        RemoveStalePlotOutput(mean_shift_path);
        RemoveStalePlotOutput(second_path);
        RemoveStalePlotOutput(second_mean_shift_path);
        Logger::Log(LogLevel::Warning,
            "UMAP scatter plots skipped because no configured spot categories were found.");
        return true;
    }
    if (second_embedding.empty())
    {
        RemoveStalePlotOutput(second_path);
        RemoveStalePlotOutput(second_mean_shift_path);
    }
    if (!WriteEmbeddingPlot(prepared, embedding, plot_path, *ranges, error_message)
        || !WriteEmbeddingPlot(prepared, mean_shift_embedding, mean_shift_path, *ranges, error_message))
    {
        return false;
    }
    if (second_embedding.empty()) return true;

    const auto second_ranges{ ComputeUmapPlotRanges(prepared, second_embedding) };
    return WriteEmbeddingPlot(prepared, second_embedding, second_path, *second_ranges, error_message)
        && WriteEmbeddingPlot(prepared, second_mean_shift_embedding, second_mean_shift_path,
            *second_ranges, error_message);
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
        embedding = RunUmap(prepared.features, prepared.effective_neighbors,
            request.min_dist, request.num_epochs, request.random_seed, request.job_count);
    }
    catch (const std::exception & error)
    {
        runner.RequirePrepareCondition(false,
            "UMAP embedding failed for model key '" + request.model_key_tag
                + "' in database '" + request.database_path.string()
                + "': " + error.what());
        return false;
    }

    if (!HasFiniteCoordinates(embedding))
    {
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
    if (!HasFiniteCoordinates(mean_shift_embedding))
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
        binning = BinMeanShiftEmbedding(prepared.feature_rows, mean_shift_embedding);
        const auto second_features{ BuildSecondUmapFeatures(
            prepared.feature_rows, binning.point_bin_numbers) };
        if (second_features.values.empty())
        {
            Logger::Log(LogLevel::Warning,
                "Second UMAP skipped because all selected feature columns are constant.");
        }
        else
        {
            second_embedding = RunUmap(second_features, prepared.effective_neighbors,
                request.min_dist, request.num_epochs, request.random_seed, request.job_count);
            if (!HasFiniteCoordinates(second_embedding))
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

    std::vector<double> second_mean_shift_embedding;
    UmapBinning second_binning;
    if (!second_embedding.empty())
    {
        try
        {
            second_mean_shift_embedding = RunUmapMeanShift(
                second_embedding, prepared.feature_rows, request.job_count);
            if (!HasFiniteCoordinates(second_mean_shift_embedding))
            {
                throw std::runtime_error("Mean shift produced a non-finite coordinate.");
            }
            second_binning = BinMeanShiftEmbedding(prepared.feature_rows, second_mean_shift_embedding);
        }
        catch (const std::exception & error)
        {
            runner.RequirePrepareCondition(false,
                "Second UMAP mean shift or binning failed for model key '" + request.model_key_tag
                    + "' in database '" + request.database_path.string() + "': " + error.what());
            return false;
        }
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

    const auto output_stem{ prepared.output_path.parent_path() / prepared.output_path.stem() };
    if (!WriteMeanShiftBins(prepared, binning, output_stem.string() + "_mean_shift_bins.csv", write_error))
    {
        runner.AddFieldValidationError(&CommandRequestBase::output_dir, write_error);
        return false;
    }
    Logger::Log(LogLevel::Info, "UMAP mean shift bin count: " + std::to_string(binning.bins.size()));

    if (!WriteMeanShiftBins(prepared, second_binning,
        output_stem.string() + "_second_umap_mean_shift_bins.csv", write_error))
    {
        runner.AddFieldValidationError(&CommandRequestBase::output_dir, write_error);
        return false;
    }
    Logger::Log(LogLevel::Info,
        "Second UMAP mean shift bin count: " + std::to_string(second_binning.bins.size()));

    std::string plot_error;
    if (!WriteEmbeddingPlots(prepared, embedding, mean_shift_embedding, second_embedding,
        second_mean_shift_embedding, plot_error))
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

            std::vector<LocalFittingFeatureRow> source_rows;
            try
            {
                source_rows = BuildLocalFittingFeatureRows(*model_object);
            }
            catch (const std::exception & error)
            {
                validation_runner.AddFieldValidationError(
                    &UmapEmbeddingRequest::model_key_tag,
                    "Failed to build local fitting features for model key '"
                        + prepared_request.model_key_tag + "': " + error.what());
                return;
            }

            std::string preparation_error;
            prepared_input = BuildAndStandardizeInput(
                *model_object, std::move(source_rows), prepared_request.model_key_tag,
                prepared_request.num_neighbors, preparation_error);
            if (!prepared_input)
            {
                validation_runner.AddFieldValidationError(
                    &UmapEmbeddingRequest::model_key_tag,
                    preparation_error);
                return;
            }

            prepared_input->output_path = prepared_request.output_dir
                / ("umap_embedding_" + path_helper::EnsureSanitizedTag(prepared_request.model_key_tag) + ".csv");

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
