#include <gtest/gtest.h>

#include <chrono>
#include <filesystem>
#include <string>

#include <rhbm_gem/utils/domain/LocalPainter.hpp>

#ifdef HAVE_ROOT
#include <TROOT.h>
#endif

namespace {

class ScopedTempDir
{
public:
    ScopedTempDir()
    {
        const auto now{
            std::chrono::steady_clock::now().time_since_epoch().count()
        };
        m_path =
            std::filesystem::temp_directory_path() /
            ("LocalPainterTest_" + std::to_string(now));
        std::filesystem::create_directories(m_path);
    }

    ~ScopedTempDir()
    {
        std::error_code ec;
        std::filesystem::remove_all(m_path, ec);
    }

    const std::filesystem::path & path() const
    {
        return m_path;
    }

private:
    std::filesystem::path m_path;
};

rhbm_gem::LinePlotRequest BuildValidLinePlotRequest(const std::filesystem::path & output_path)
{
    rhbm_gem::LinePlotRequest request;
    request.output_path = output_path;
    request.x_axis.title = "Alpha";
    request.panels = {
        rhbm_gem::LinePlotPanel{
            "Amplitude",
            rhbm_gem::AxisSpec{ "Bias", std::nullopt, true },
            {
                rhbm_gem::LineSeries{
                    "series",
                    { 0.0, 0.5, 1.0 },
                    { 1.0, 0.8, 0.6 },
                    std::nullopt
                }
            }
        }
    };
    return request;
}

rhbm_gem::HeatmapRequest BuildValidHeatmapRequest(const std::filesystem::path & output_path)
{
    rhbm_gem::HeatmapRequest request;
    request.output_path = output_path;
    request.values = (Eigen::MatrixXd(2, 2) << 1.0, 2.0, 3.0, 4.0).finished();
    request.x_axis.title = "X";
    request.y_axis.title = "Y";
    request.z_axis.title = "Z";
    return request;
}

} // namespace

TEST(LocalPainterTest, SaveLinePlotRejectsMissingPanels)
{
    rhbm_gem::LinePlotRequest request;
    request.output_path = "line.pdf";
    request.x_axis.title = "Alpha";

    const auto result{ rhbm_gem::local_painter::SaveLinePlot(request) };

    EXPECT_EQ(rhbm_gem::PlotStatus::INVALID_INPUT, result.status);
    EXPECT_FALSE(result.Succeeded());
}

TEST(LocalPainterTest, SaveLinePlotRejectsMismatchedSeriesLength)
{
    auto request{ BuildValidLinePlotRequest("line.pdf") };
    request.panels.front().series.front().y_values.pop_back();

    const auto result{ rhbm_gem::local_painter::SaveLinePlot(request) };

    EXPECT_EQ(rhbm_gem::PlotStatus::INVALID_INPUT, result.status);
    EXPECT_FALSE(result.Succeeded());
}

TEST(LocalPainterTest, SaveHeatmapRejectsEmptyMatrix)
{
    rhbm_gem::HeatmapRequest request;
    request.output_path = "heatmap.pdf";

    const auto result{ rhbm_gem::local_painter::SaveHeatmap(request) };

    EXPECT_EQ(rhbm_gem::PlotStatus::INVALID_INPUT, result.status);
    EXPECT_FALSE(result.Succeeded());
}

#ifndef HAVE_ROOT

TEST(LocalPainterTest, SaveLinePlotReportsRootUnavailableWithoutRoot)
{
    const auto result{
        rhbm_gem::local_painter::SaveLinePlot(BuildValidLinePlotRequest("line.pdf"))
    };

    EXPECT_EQ(rhbm_gem::PlotStatus::ROOT_UNAVAILABLE, result.status);
    EXPECT_FALSE(result.Succeeded());
}

TEST(LocalPainterTest, SaveHeatmapReportsRootUnavailableWithoutRoot)
{
    const auto result{
        rhbm_gem::local_painter::SaveHeatmap(BuildValidHeatmapRequest("heatmap.pdf"))
    };

    EXPECT_EQ(rhbm_gem::PlotStatus::ROOT_UNAVAILABLE, result.status);
    EXPECT_FALSE(result.Succeeded());
}

#else

TEST(LocalPainterTest, SaveLinePlotWritesOutputFile)
{
    gROOT->SetBatch(kTRUE);
    ScopedTempDir temp_dir;
    const auto output_path{ temp_dir.path() / "line_plot.pdf" };

    const auto result{
        rhbm_gem::local_painter::SaveLinePlot(BuildValidLinePlotRequest(output_path))
    };

    EXPECT_TRUE(result.Succeeded());
    EXPECT_TRUE(std::filesystem::exists(output_path));
}

TEST(LocalPainterTest, SaveHeatmapWritesOutputFile)
{
    gROOT->SetBatch(kTRUE);
    ScopedTempDir temp_dir;
    const auto output_path{ temp_dir.path() / "heatmap.pdf" };

    const auto result{
        rhbm_gem::local_painter::SaveHeatmap(BuildValidHeatmapRequest(output_path))
    };

    EXPECT_TRUE(result.Succeeded());
    EXPECT_TRUE(std::filesystem::exists(output_path));
}

#endif
