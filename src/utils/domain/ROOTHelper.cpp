#include <rhbm_gem/utils/domain/ROOTHelper.hpp>

#ifdef HAVE_ROOT

#include <TAxis.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TF1.h>
#include <TF2.h>
#include <TFile.h>
#include <TFrame.h>
#include <TGaxis.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TGraphErrors.h>
#include <TGraph2DErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMarker.h>
#include <TMath.h>
#include <TPaletteAxis.h>
#include <TPaveStats.h>
#include <TPaveText.h>
#include <TProfile.h>
#include <TString.h>
#include <TStyle.h>
#include <TTree.h>
#include <TROOT.h>

#include <stdexcept>

using namespace ROOT;

std::unique_ptr<TCanvas> ROOTHelper::CreateCanvas(
    const std::string & name, const std::string & title, int width, int height)
{
    return std::make_unique<TCanvas>(name.data(), title.data(), width, height);
}

std::unique_ptr<TPad> ROOTHelper::CreatePad(
    const std::string & name, const std::string & title,
    double x_low, double y_low, double x_up, double y_up)
{
    return std::make_unique<TPad>(name.data(), title.data(), x_low, y_low, x_up, y_up);
}

std::unique_ptr<TH1D> ROOTHelper::CreateHist1D(
    const std::string & name, const std::string & title,
    int x_bin, double x_min, double x_max)
{
    return std::make_unique<TH1D>(name.data(), title.data(), x_bin, x_min, x_max);
}

std::unique_ptr<TH2D> ROOTHelper::CreateHist2D(
    const std::string & name, const std::string & title,
    int x_bin, double x_min, double x_max,
    int y_bin, double y_min, double y_max)
{
    return std::make_unique<TH2D>(name.data(), title.data(), x_bin, x_min, x_max, y_bin, y_min, y_max);
}

std::unique_ptr<TGraphErrors> ROOTHelper::CreateGraphErrors()
{
    return std::make_unique<TGraphErrors>();
}

std::unique_ptr<TGraphErrors> ROOTHelper::CreateGraphErrors(const int & point_size)
{
    return std::make_unique<TGraphErrors>(point_size);
}

std::unique_ptr<TGraph2DErrors> ROOTHelper::CreateGraph2DErrors()
{
    return std::make_unique<TGraph2DErrors>();
}

std::unique_ptr<TGraph2DErrors> ROOTHelper::CreateGraph2DErrors(const int & point_size)
{
    return std::make_unique<TGraph2DErrors>(point_size);
}

std::unique_ptr<TGraph2DErrors> ROOTHelper::CreateGraph2DErrors(
    const int & point_size,
    std::vector<double> & x, std::vector<double> & y, std::vector<double> & z)
{
    return std::make_unique<TGraph2DErrors>(point_size, &x[0], &y[0], &z[0]);
}

std::unique_ptr<TPaveText> ROOTHelper::CreatePaveText(
    double x1, double y1, double x2, double y2, std::string option, bool in_partition)
{
    if (in_partition == true)
    {
        x1 = GetXtoPadInCanvasPartition(x1);
        x2 = GetXtoPadInCanvasPartition(x2);
        y1 = GetYtoPadInCanvasPartition(y1);
        y2 = GetYtoPadInCanvasPartition(y2);
    }
    return std::make_unique<TPaveText>(x1, y1, x2, y2, option.data());
}

std::unique_ptr<TLegend> ROOTHelper::CreateLegend(
    double x1, double y1, double x2, double y2, bool in_partition)
{
    if (in_partition == true)
    {
        x1 = GetXtoPadInCanvasPartition(x1);
        x2 = GetXtoPadInCanvasPartition(x2);
        y1 = GetYtoPadInCanvasPartition(y1);
        y2 = GetYtoPadInCanvasPartition(y2);
    }
    return std::make_unique<TLegend>(x1, y1, x2, y2);
}

std::unique_ptr<TF1> ROOTHelper::CreateFunction1D(const std::string & name, const std::string & form)
{
    return std::make_unique<TF1>(name.data(), form.data());
}

std::unique_ptr<TF1> ROOTHelper::CreateLinearModelFunction(
    const std::string & name, double beta_0, double beta_1, double x_min, double x_max)
{
    auto function{ std::make_unique<TF1>(name.data(), LinearModelFunction, x_min, x_max, 2) };
    function->SetParameter(0, beta_0);
    function->SetParameter(1, beta_1);
    return function;
}

std::unique_ptr<TF1> ROOTHelper::CreateGaus2DFunctionIn1D(
    const std::string & name, double amplitude, double width, double x_min, double x_max)
{
    auto function{ std::make_unique<TF1>(name.data(), Gaus2DModelFunction, x_min, x_max, 2) };
    function->SetParameter(0, amplitude);
    function->SetParameter(1, width);
    return function;
}

std::unique_ptr<TF1> ROOTHelper::CreateGaus3DFunctionIn1D(
    const std::string & name, double amplitude, double width, double x_min, double x_max)
{
    auto function{ std::make_unique<TF1>(name.data(), Gaus3DModelFunction, x_min, x_max, 2) };
    function->SetParameter(0, amplitude);
    function->SetParameter(1, width);
    return function;
}

std::unique_ptr<TLine> ROOTHelper::CreateLine(double x1, double y1, double x2, double y2)
{
    return std::make_unique<TLine>(x1, y1, x2, y2);
}

void ROOTHelper::PrintCanvasOpen(TCanvas * canvas, const std::string & name)
{
    gErrorIgnoreLevel = kWarning;
    canvas->cd();
    canvas->Print(std::string(name).append("[").data());
}

void ROOTHelper::PrintCanvasPad(TCanvas * canvas, const std::string & name)
{
    canvas->Update();
    canvas->Print(name.data());
}

void ROOTHelper::PrintCanvasClose(TCanvas * canvas, const std::string & name)
{
    canvas->Print(std::string(name).append("]").data());
}

void ROOTHelper::SetPadFrameAttribute(
    TAttPad * pad, int mode, short size, short fill_style, short fill_color,
    short line_style, short line_width, short line_color)
{
    pad->SetFrameBorderMode(mode);
    pad->SetFrameBorderSize(size);
    pad->SetFrameFillStyle(fill_style);
    pad->SetFrameFillColor(fill_color);
    pad->SetFrameLineStyle(line_style);
    pad->SetFrameLineWidth(line_width);
    pad->SetFrameLineColor(line_color);
}

void ROOTHelper::SetPadMarginAttribute(
    TAttPad * pad, float left, float right, float bottom, float top)
{
    pad->SetLeftMargin(left);
    pad->SetRightMargin(right);
    pad->SetBottomMargin(bottom);
    pad->SetTopMargin(top);
}

void ROOTHelper::SetPadMarginInCanvas(
    TVirtualPad * pad, double left, double right, double bottom, double top)
{
    if (pad == nullptr)
    {
        throw std::invalid_argument("ROOTHelper::SetPadMarginInCanvas() requires a non-null pad.");
    }
    if (pad->GetCanvas() == nullptr)
    {
        throw std::runtime_error(
            "ROOTHelper::SetPadMarginInCanvas() requires the pad to be attached to a canvas.");
    }
    pad->GetCanvas()->Update();

    // Size of canvas's weight and height (unit: pixel)
    auto canvas_width_pixel { pad->GetCanvas()->GetWw() };
    auto canvas_height_pixel{ pad->GetCanvas()->GetWh() };

    auto pad_width_pixel { static_cast<float>(pad->GetAbsWNDC() * canvas_width_pixel) };
    auto pad_height_pixel{ static_cast<float>(pad->GetAbsHNDC() * canvas_height_pixel) };

    auto desired_left_margin_pixel  { static_cast<float>(left   * canvas_width_pixel) };
    auto desired_right_margin_pixel { static_cast<float>(right  * canvas_width_pixel) };
    auto desired_top_margin_pixel   { static_cast<float>(top    * canvas_height_pixel) };
    auto desired_bottom_margin_pixel{ static_cast<float>(bottom * canvas_height_pixel) };

    auto pad_left_margin  { desired_left_margin_pixel   / pad_width_pixel  };
    auto pad_right_margin { desired_right_margin_pixel  / pad_width_pixel  };
    auto pad_top_margin   { desired_top_margin_pixel    / pad_height_pixel };
    auto pad_bottom_margin{ desired_bottom_margin_pixel / pad_height_pixel };
    SetPadMarginAttribute(pad, pad_left_margin, pad_right_margin, pad_bottom_margin, pad_top_margin);
}

void ROOTHelper::SetPaveTextMarginInCanvas(
    TVirtualPad * pad, TPaveText * pave, double left, double right, double bottom, double top)
{
    if (!pave)
    {
        throw std::invalid_argument(
            "ROOTHelper::SetPaveTextMarginInCanvas() requires a non-null TPaveText.");
    }
    if (!pad)
    {
        throw std::invalid_argument(
            "ROOTHelper::SetPaveTextMarginInCanvas() requires a non-null TVirtualPad.");
    }
    pave->Draw();
    if (pad->GetCanvas()) pad->GetCanvas()->Update();

    double x1_canvas{ pad->GetXlowNDC() + left };          
    double x2_canvas{ pad->GetXlowNDC() + pad->GetAbsWNDC() - right };
    double y1_canvas{ pad->GetYlowNDC() + bottom };          
    double y2_canvas{ pad->GetYlowNDC() + pad->GetAbsHNDC() - top };

    pave->SetX1NDC((x1_canvas - pad->GetXlowNDC()) / pad->GetAbsWNDC());
    pave->SetX2NDC((x2_canvas - pad->GetXlowNDC()) / pad->GetAbsWNDC());
    pave->SetY1NDC((y1_canvas - pad->GetYlowNDC()) / pad->GetAbsHNDC());
    pave->SetY2NDC((y2_canvas - pad->GetYlowNDC()) / pad->GetAbsHNDC());
}

void ROOTHelper::SetLegendMarginInCanvas(
    TVirtualPad * pad, TLegend * legend, double left, double right, double bottom, double top)
{
    if (!legend)
    {
        throw std::invalid_argument(
            "ROOTHelper::SetLegendMarginInCanvas() requires a non-null TLegend.");
    }
    //legend->Draw();

    if (!pad)
    {
        throw std::invalid_argument(
            "ROOTHelper::SetLegendMarginInCanvas() requires a non-null TVirtualPad.");
    }
    if (pad->GetCanvas()) pad->GetCanvas()->Update();

    double x1_canvas{ pad->GetXlowNDC() + left };          
    double x2_canvas{ pad->GetXlowNDC() + pad->GetAbsWNDC() - right };
    double y1_canvas{ pad->GetYlowNDC() + bottom };          
    double y2_canvas{ pad->GetYlowNDC() + pad->GetAbsHNDC() - top };

    legend->SetX1NDC((x1_canvas - pad->GetXlowNDC()) / pad->GetAbsWNDC());
    legend->SetX2NDC((x2_canvas - pad->GetXlowNDC()) / pad->GetAbsWNDC());
    legend->SetY1NDC((y1_canvas - pad->GetYlowNDC()) / pad->GetAbsHNDC());
    legend->SetY2NDC((y2_canvas - pad->GetYlowNDC()) / pad->GetAbsHNDC());
}

void ROOTHelper::SetAxisTitleAttribute(TAttAxis * axis, float size, float offset, short font, short color)
{
    axis->SetTitleFont(font);
    axis->SetTitleSize(size);
    axis->SetTitleOffset(offset);
    axis->SetTitleColor(color);
}

void ROOTHelper::SetAxisLabelAttribute(TAttAxis * axis, float size, float offset, short font, short color)
{
    axis->SetLabelFont(font);
    axis->SetLabelSize(size);
    axis->SetLabelOffset(offset);
    axis->SetLabelColor(color);
}

void ROOTHelper::SetAxisTickAttribute(TAttAxis * axis, float length, int division)
{
    axis->SetTickLength(length);
    axis->SetNdivisions(division);
}

void ROOTHelper::SetTextAttribute(
    TAttText * text,
    float size, short font, short align, float angle, short color, float transparent)
{
    text->SetTextSize(size);
    text->SetTextFont(font);
    text->SetTextColorAlpha(color, transparent);
    text->SetTextAlign(align);
    text->SetTextAngle(angle);
}

void ROOTHelper::SetLineAttribute(
    TAttLine * line, short style, short width, short color, float transparent)
{
    line->SetLineStyle(style);
    line->SetLineWidth(width);
    line->SetLineColorAlpha(color, transparent);
}

void ROOTHelper::SetFillAttribute(TAttFill * fill, short style, short color, float transparent)
{
    fill->SetFillStyle(style);
    fill->SetFillColorAlpha(color, transparent);
}

void ROOTHelper::SetMarkerAttribute(TAttMarker * marker, short style, float size, short color, float transparent)
{
    marker->SetMarkerStyle(style);
    marker->SetMarkerSize(size);
    marker->SetMarkerColorAlpha(color, transparent);
}

void ROOTHelper::SetPaveAttribute(TPave * pave, int bordersize, double radius)
{
    pave->SetBorderSize(bordersize);
    pave->SetCornerRadius(radius);
}

void ROOTHelper::SetBoxAttribute(
    TBox * box, double left, double right, double bottom, double top)
{
    gPad->Update();
    box->SetBBoxX1(gPad->XtoPixel(left));
    box->SetBBoxX2(gPad->XtoPixel(right));
    box->SetBBoxY1(gPad->YtoPixel(bottom));
    box->SetBBoxY2(gPad->YtoPixel(top));
}

void ROOTHelper::SetPadInCanvas(
    TVirtualPad * pad, int width, int height,
    int division_x, int division_y, float margin_x, float margin_y)
{
    pad->SetCanvasSize(static_cast<unsigned int>(width), static_cast<unsigned int>(height));
    pad->Divide(division_x, division_y, margin_x, margin_y);
}

void ROOTHelper::SetPadDefaultStyle(TVirtualPad * pad)
{
    ROOTHelper::SetPadFrameAttribute(pad);
    ROOTHelper::SetFillAttribute(pad);
}

void ROOTHelper::SetPadRangeInCanvas(TVirtualPad * pad, double x1, double y1, double x2, double y2)
{
    auto x1_default = pad->GetAbsXlowNDC();
    auto y1_default = pad->GetAbsYlowNDC();
    auto x2_default = pad->GetAbsXlowNDC() + pad->GetAbsWNDC();
    auto y2_default = pad->GetAbsYlowNDC() + pad->GetAbsHNDC();
    if (x1==0.0 && y1 == 0.0 && x2==0.0 && y2 == 0.0)
    {
        pad->SetPad(x1_default, y1_default, x2_default, y2_default);
    }
    else
    {
        pad->SetPad(x1, y1, x2, y2);
    }
}

void ROOTHelper::SetPadLayout(
    TVirtualPad * pad, int grid_x, int grid_y, int log_x, int log_y, int tick_x, int tick_y)
{
    pad->SetGrid(grid_x, grid_y);
    pad->SetLogx(log_x);
    pad->SetLogy(log_y);
    pad->SetTicks(tick_x, tick_y);
}

void ROOTHelper::SetCanvasDefaultStyle(TCanvas * canvas)
{
    ROOTHelper::SetFillAttribute(canvas);
}

void ROOTHelper::SetPaveTextDefaultStyle(TPaveText * text)
{
    ROOTHelper::SetPaveAttribute(text);
    ROOTHelper::SetFillAttribute(text);
    ROOTHelper::SetTextAttribute(text);
}

void ROOTHelper::SetLegendDefaultStyle(TLegend * legend)
{
    ROOTHelper::SetPaveAttribute(legend);
    ROOTHelper::SetFillAttribute(legend);
    ROOTHelper::SetTextAttribute(legend);
}

void ROOTHelper::SetCanvasPartition(
    TCanvas * canvas, const int Nx, const int Ny,
    float lMargin, float rMargin, float bMargin, float tMargin,
    float vSpacing, float hSpacing)
{
    if (canvas == nullptr) return;

    auto x_dim{ static_cast<float>(Nx) };
    auto y_dim{ static_cast<float>(Ny) };
    auto vStep{ (1.0f - bMargin - tMargin - (y_dim - 1.0f) * vSpacing) / y_dim };
    auto hStep{ (1.0f - lMargin - rMargin - (x_dim - 1.0f) * hSpacing) / x_dim };
    float vposd, vposu, vfactor;
    float hposl, hposr, hfactor;
    for (int i = 0; i < Nx; i++)
    {
        if (i == 0)
        {
            hposl = 0.0;
            hposr = (Nx == 1) ? 1.0f : static_cast<float>(lMargin + hStep);
            hfactor = hposr - hposl;
        }
        else if (i == Nx - 1)
        {
            hposl = hposr + hSpacing;
            hposr = (hposl + hStep + rMargin > 1.0) ? 1.0f : static_cast<float>(hposl + hStep + rMargin);
            hfactor = hposr - hposl;
        }
        else
        {
            hposl = hposr + hSpacing;
            hposr = hposl + hStep;
            hfactor = hposr - hposl;
        }

        for (int j = 0; j < Ny; j++)
        {
            if (j == 0)
            {
                vposd = 0.0;
                vposu = (Ny == 1) ? 1.0f : static_cast<float>(bMargin + vStep);
                vfactor = vposu - vposd;
            }
            else if (j == Ny - 1)
            {
                vposd = vposu + vSpacing;
                vposu = (vposd + vStep + tMargin > 1.0) ? 1.0f : static_cast<float>(vposd + vStep + tMargin);
                vfactor = vposu - vposd;
            }
            else
            {
                vposd = vposu + vSpacing;
                vposu = vposd + vStep;
                vfactor = vposu - vposd;
            }

            canvas->cd(0);
            auto name{ TString::Format("pad_%d_%d", i, j) };
            auto pad{ dynamic_cast<TPad*>(canvas->FindObject(name.Data())) };
            if (pad) delete pad;
            pad = new TPad(name.Data(), "", hposl, vposd, hposr, vposu);
            auto left_margin{ (i == 0) ? static_cast<float>(lMargin / hfactor) : 0.0f };
            auto right_margin{ (i == Nx - 1) ? static_cast<float>(rMargin / hfactor) : 0.0f };
            auto bottom_margin{ (j == 0) ? static_cast<float>(bMargin / vfactor) : 0.0f };
            auto top_margin{ (j == Ny - 1) ? static_cast<float>(tMargin / vfactor) : 0.0f };
            SetPadMarginAttribute(pad, left_margin, right_margin, bottom_margin, top_margin);
            SetPadFrameAttribute(pad);
            SetFillAttribute(pad);
            pad->SetBorderMode(0);
            pad->SetBorderSize(0);
            pad->Draw();
        }
    }
}

void ROOTHelper::FindPadInCanvasPartition(TCanvas * canvas, int id_x, int id_y)
{
    canvas->cd(0);
    auto pad{ dynamic_cast<TPad*>(canvas->FindObject(
        TString::Format("pad_%d_%d", id_x, id_y).Data())) };
    if (pad == nullptr)
    {
        throw std::runtime_error("ROOTHelper::FindPadInCanvasPartition(): pad not found");
    }
    pad->Draw();
    pad->cd();
}

float ROOTHelper::GetPadXfactorInCanvasPartition(TCanvas * canvas, TVirtualPad * pad)
{
    auto pad_origin{ dynamic_cast<TPad*>(canvas->FindObject("pad_0_0")) };
    if (pad_origin == nullptr)
    {
        throw std::runtime_error("ROOTHelper::GetPadXfactorInCanvasPartition(): pad_0_0 not found");
    }
    return static_cast<float>(pad_origin->GetAbsWNDC() / pad->GetAbsWNDC());
}

float ROOTHelper::GetPadYfactorInCanvasPartition(TCanvas * canvas, TVirtualPad * pad)
{
    auto pad_origin{ dynamic_cast<TPad*>(canvas->FindObject("pad_0_0")) };
    if (pad_origin == nullptr)
    {
        throw std::runtime_error("ROOTHelper::GetPadYfactorInCanvasPartition(): pad_0_0 not found");
    }
    return static_cast<float>(pad_origin->GetAbsHNDC() / pad->GetAbsHNDC());
}

double ROOTHelper::GetXtoPadInCanvasPartition(double x)
{
    double xl, yl, xu, yu;
    gPad->GetPadPar(xl, yl, xu, yu);
    double pw{ xu - xl };
    double lm{ gPad->GetLeftMargin() };
    double rm{ gPad->GetRightMargin() };
    double fw{ pw - pw * lm - pw * rm };
    return (x * fw + pw * lm) / pw;
}

double ROOTHelper::GetYtoPadInCanvasPartition(double y)
{
    double xl, yl, xu, yu;
    gPad->GetPadPar(xl, yl, xu, yu);
    double ph{ yu - yl };
    double tm{ gPad->GetTopMargin() };
    double bm{ gPad->GetBottomMargin() };
    double fh{ ph - ph * bm - ph * tm };
    return (y * fh + bm * ph) / ph;
}

double ROOTHelper::ConvertGlobalTickLengthToPadTickLength(
    TVirtualPad * pad, double global_tick_length, bool use_width)
{
    auto fraction_pad_w{ 1.0 - pad->GetLeftMargin() - pad->GetRightMargin() };
    auto fraction_pad_h{ 1.0 - pad->GetBottomMargin() - pad->GetTopMargin() };
    auto tick_scale_x{ fraction_pad_w * pad->GetAbsHNDC() };
    auto tick_scale_y{ fraction_pad_h * pad->GetAbsWNDC() };
    return (use_width) ? global_tick_length / tick_scale_y : global_tick_length / tick_scale_x;
}

double ROOTHelper::GetLinearRegressionRSquare(const TGraphErrors * graph, const TF1 * function)
{
    auto data_size{ graph->GetN() };
    auto y_mean{ graph->GetMean(2) };
    auto residual_sum{ 0.0 };
    auto total_sum{ 0.0 };
    for (int i = 0; i < data_size; i++)
    {
        auto x_data{ graph->GetX()[i] };
        auto y_data{ graph->GetY()[i] };
        auto y_model{ function->Eval(x_data) };
        auto residual{ y_data - y_model };
        residual_sum += residual * residual;
        total_sum += std::pow(y_data - y_mean, 2);
    }
    auto r_square{ (total_sum == 0.0) ? 1.0 : 1.0 - residual_sum/total_sum };
    return r_square;
}

double ROOTHelper::PerformLinearRegression(
    TGraphErrors * graph, double & slope, double & intercept)
{
    auto fit_tmp{ CreateFunction1D("fit_linear","[1]*x+[0]") };
    graph->Fit(fit_tmp.get(), "0 WQ");
    double par[2];
    fit_tmp->GetParameters(par);
    intercept = par[0];
    slope = par[1];
    return GetLinearRegressionRSquare(graph, fit_tmp.get());
}

std::tuple<double, double> ROOTHelper::GetRangeInGraph(TGraphErrors * graph)
{
    double y_min{ std::numeric_limits<double>::max() };
    double y_max{ std::numeric_limits<double>::lowest() };
    auto data_size{ graph->GetN() };
    for (int i = 0; i < data_size; i++)
    {
        auto y_data{ graph->GetY()[i] };
        if (y_data < y_min) y_min = y_data;
        if (y_data > y_max) y_max = y_data;
    }
    return std::make_tuple(y_min, y_max);
}

double ROOTHelper::LinearModelFunction(double * x, double * par)
{
    double y{ par[0] + par[1]*x[0] };
    return y;
}

double ROOTHelper::Gaus2DModelFunction(double * x, double * par)
{
    double tau_square{ par[1]*par[1] };
    double y{ par[0] / (2.0*TMath::Pi()*tau_square) * std::exp(-x[0]*x[0]/(2.0*tau_square)) };
    return y;
}

double ROOTHelper::Gaus3DModelFunction(double * x, double * par)
{
    double tau_square{ par[1]*par[1] };
    double y{ par[0] * std::pow(2.0*TMath::Pi()*tau_square, -1.5) * std::exp(-x[0]*x[0]/(2.0*tau_square)) };
    return y;
}

#endif
