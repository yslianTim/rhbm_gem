#pragma once

#include "IPlotRenderBackend.hpp"

namespace rhbm_gem {

#ifdef HAVE_ROOT
class RootPlotRenderBackend final : public IPlotRenderBackend
{
public:
    std::unique_ptr<::TH1D> CreateHist1D(
        const std::string & name,
        const std::string & title,
        int bins,
        double x_min,
        double x_max) const override;
    std::unique_ptr<::TH2D> CreateHist2D(
        const std::string & name,
        const std::string & title,
        int x_bins,
        double x_min,
        double x_max,
        int y_bins,
        double y_min,
        double y_max) const override;
    std::unique_ptr<::TGraphErrors> CreateGraphErrors(int size = 0) const override;
    std::unique_ptr<::TF1> CreateLinearModelFunction(
        const std::string & name,
        double beta_0,
        double beta_1,
        double x_min = 0.0,
        double x_max = 1.5) const override;
    std::unique_ptr<::TF1> CreateGaus3DFunctionIn1D(
        const std::string & name,
        double amplitude,
        double width) const override;
    std::unique_ptr<::TF1> CreateGaus2DFunctionIn1D(
        const std::string & name,
        double amplitude,
        double width) const override;
};
#endif

} // namespace rhbm_gem
