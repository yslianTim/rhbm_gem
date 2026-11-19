#include "core/classifier/UmapEmbeddingPlot.hpp"

#include <rhbm_gem/utils/domain/Logger.hpp>

#ifdef HAVE_ROOT
#include "core/classifier/UmapEmbeddingFeatures.hpp"
#include "core/classifier/UmapSpotCategories.hpp"

#include <rhbm_gem/utils/domain/ROOTHelper.hpp>
#include <rhbm_gem/utils/math/ArrayHelper.hpp>

#include <TCanvas.h>
#include <TColor.h>
#include <TGraphErrors.h>
#include <TH2.h>
#include <TLegend.h>

#include <algorithm>
#include <cmath>
#include <exception>
#include <iterator>
#include <memory>
#include <optional>
#include <system_error>
#include <tuple>
#include <utility>
#include <vector>
#endif

namespace rhbm_gem::core::classifier {

#ifdef HAVE_ROOT
namespace {

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

} // namespace
#endif

bool WriteEmbeddingPlots(
    const std::filesystem::path & output_path,
    const PreparedUmapInput & prepared,
    const UmapEmbeddingResult & first,
    const UmapEmbeddingResult & second,
    std::string & error_message)
{
#ifndef HAVE_ROOT
    (void)output_path;
    (void)prepared;
    (void)first;
    (void)second;
    (void)error_message;
    Logger::Log(LogLevel::Warning,
        "UMAP scatter plots skipped because ROOT support is unavailable.");
    return true;
#else
    auto plot_path{ output_path };
    plot_path.replace_extension(".pdf");
    const auto mean_shift_path{ output_path.parent_path()
        / (output_path.stem().string() + "_mean_shift.pdf") };
    const auto second_path{ output_path.parent_path()
        / (output_path.stem().string() + "_second_umap.pdf") };
    const auto second_mean_shift_path{ output_path.parent_path()
        / (output_path.stem().string() + "_second_umap_mean_shift.pdf") };
    const auto ranges{ ComputeUmapPlotRanges(prepared, first.coordinates) };
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
    if (second.coordinates.empty())
    {
        RemoveStalePlotOutput(second_path);
        RemoveStalePlotOutput(second_mean_shift_path);
    }
    if (!WriteEmbeddingPlot(prepared, first.coordinates, plot_path, *ranges, error_message)
        || !WriteEmbeddingPlot(prepared, first.mean_shift_coordinates, mean_shift_path, *ranges, error_message))
    {
        return false;
    }
    if (second.coordinates.empty()) return true;

    const auto second_ranges{ ComputeUmapPlotRanges(prepared, second.coordinates) };
    return WriteEmbeddingPlot(prepared, second.coordinates, second_path, *second_ranges, error_message)
        && WriteEmbeddingPlot(prepared, second.mean_shift_coordinates, second_mean_shift_path,
            *second_ranges, error_message);
#endif
}

} // namespace rhbm_gem::core::classifier
