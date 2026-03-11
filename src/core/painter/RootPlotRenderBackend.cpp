#include "internal/plot/RootPlotRenderBackend.hpp"

#ifdef HAVE_ROOT

#include <rhbm_gem/utils/domain/ROOTHelper.hpp>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>

namespace rhbm_gem {

std::unique_ptr<::TH1D> RootPlotRenderBackend::CreateHist1D(
    const std::string & name,
    const std::string & title,
    int bins,
    double x_min,
    double x_max) const
{
    return ROOTHelper::CreateHist1D(name, title, bins, x_min, x_max);
}

std::unique_ptr<::TH2D> RootPlotRenderBackend::CreateHist2D(
    const std::string & name,
    const std::string & title,
    int x_bins,
    double x_min,
    double x_max,
    int y_bins,
    double y_min,
    double y_max) const
{
    return ROOTHelper::CreateHist2D(name, title, x_bins, x_min, x_max, y_bins, y_min, y_max);
}

std::unique_ptr<::TGraphErrors> RootPlotRenderBackend::CreateGraphErrors(int size) const
{
    return ROOTHelper::CreateGraphErrors(size);
}

std::unique_ptr<::TF1> RootPlotRenderBackend::CreateLinearModelFunction(
    const std::string & name,
    double beta_0,
    double beta_1,
    double x_min,
    double x_max) const
{
    return ROOTHelper::CreateLinearModelFunction(name, beta_0, beta_1, x_min, x_max);
}

std::unique_ptr<::TF1> RootPlotRenderBackend::CreateGaus3DFunctionIn1D(
    const std::string & name,
    double amplitude,
    double width) const
{
    return ROOTHelper::CreateGaus3DFunctionIn1D(name, amplitude, width);
}

std::unique_ptr<::TF1> RootPlotRenderBackend::CreateGaus2DFunctionIn1D(
    const std::string & name,
    double amplitude,
    double width) const
{
    return ROOTHelper::CreateGaus2DFunctionIn1D(name, amplitude, width);
}

} // namespace rhbm_gem

#endif
