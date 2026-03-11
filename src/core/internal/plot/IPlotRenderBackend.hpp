#pragma once

#include <memory>
#include <string>

#ifdef HAVE_ROOT
class TH1D;
class TH2D;
class TGraphErrors;
class TF1;
#endif

namespace rhbm_gem {

class IPlotRenderBackend
{
public:
    virtual ~IPlotRenderBackend() = default;

#ifdef HAVE_ROOT
    virtual std::unique_ptr<::TH1D> CreateHist1D(
        const std::string & name,
        const std::string & title,
        int bins,
        double x_min,
        double x_max) const = 0;
    virtual std::unique_ptr<::TH2D> CreateHist2D(
        const std::string & name,
        const std::string & title,
        int x_bins,
        double x_min,
        double x_max,
        int y_bins,
        double y_min,
        double y_max) const = 0;
    virtual std::unique_ptr<::TGraphErrors> CreateGraphErrors(int size = 0) const = 0;
    virtual std::unique_ptr<::TF1> CreateLinearModelFunction(
        const std::string & name,
        double beta_0,
        double beta_1,
        double x_min = 0.0,
        double x_max = 1.5) const = 0;
    virtual std::unique_ptr<::TF1> CreateGaus3DFunctionIn1D(
        const std::string & name,
        double amplitude,
        double width) const = 0;
    virtual std::unique_ptr<::TF1> CreateGaus2DFunctionIn1D(
        const std::string & name,
        double amplitude,
        double width) const = 0;
#endif
};

} // namespace rhbm_gem
