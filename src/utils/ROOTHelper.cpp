#include "ROOTHelper.hpp"

#include <TAxis.h>
#include <TCanvas.h>
#include <TExec.h>
#include <TF1.h>
#include <TF2.h>
#include <TFile.h>
#include <TFrame.h>
#include <TGaxis.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TGraphErrors.h>
#include <TGraph2DErrors.h>
#include <TGraphMultiErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMarker.h>
#include <TMath.h>
#include <TPad.h>
#include <TPaletteAxis.h>
#include <TPaveStats.h>
#include <TPaveText.h>
#include <TProfile.h>
#include <TString.h>
#include <TStyle.h>
#include <TTree.h>
#include <TROOT.h>

using namespace ROOT;

std::unique_ptr<TCanvas> ROOTHelper::CreateCanvas(
    const std::string & name, const std::string & title, int width, int height)
{
    return std::make_unique<TCanvas>(name.data(), title.data(), width, height);
}

std::unique_ptr<TGraphErrors> CreateGraphErrors(void)
{
    return std::make_unique<TGraphErrors>();
}

std::unique_ptr<TGraphErrors> CreateGraphErrors(const int & point_size)
{
    return std::make_unique<TGraphErrors>(point_size);
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
    TAttText * _text,
    float size, short font, short align, float angle, short color, float transparent)
{
    _text->SetTextSize(size);
    _text->SetTextFont(font);
    _text->SetTextColorAlpha(color, transparent);
    _text->SetTextAlign(align);
    _text->SetTextAngle(angle);
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
    TPad * pad, int width, int height,
    int division_x, int division_y, float margin_x, float margin_y)
{
    pad->SetCanvasSize(width, height);
    pad->Divide(division_x, division_y, margin_x, margin_y);
}

void ROOTHelper::SetPadStyle(TVirtualPad * pad, double x1, double y1, double x2, double y2)
{
    ROOTHelper::SetPadRangeInCanvas(pad, x1, y1, x2, y2);
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

void ROOTHelper::SetCanvasStyle(TCanvas * canvas)
{
    ROOTHelper::SetFillAttribute(canvas);
}

void ROOTHelper::SetPaveTextStyle(TPaveText * text)
{
    ROOTHelper::SetPaveAttribute(text);
    ROOTHelper::SetFillAttribute(text);
    ROOTHelper::SetTextAttribute(text);
}

void ROOTHelper::SetLegendStyle(TLegend * legend)
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

    auto vStep{ (1. - bMargin - tMargin - (Ny - 1) * vSpacing) / Ny };
    auto hStep{ (1. - lMargin - rMargin - (Nx - 1) * hSpacing) / Nx };
    float vposd, vposu, vfactor;
    float hposl, hposr, hfactor;
    for (int i = 0; i < Nx; i++)
    {
        if (i == 0)
        {
            hposl = 0.0;
            hposr = (Nx == 1) ? 1.0 : lMargin + hStep;
            hfactor = hposr - hposl;
        }
        else if (i == Nx - 1)
        {
            hposl = hposr + hSpacing;
            hposr = (hposl + hStep + rMargin > 1.0) ? 1.0 : hposl + hStep + rMargin;
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
                vposu = (Ny == 1) ? 1.0 : bMargin + vStep;
                vfactor = vposu - vposd;
            }
            else if (j == Ny - 1)
            {
                vposd = vposu + vSpacing;
                vposu = (vposd + vStep + tMargin > 1.0) ? 1.0 : vposd + vStep + tMargin;
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
            auto leftMargin{ (i == 0) ? lMargin / hfactor : 0.0 };
            auto rightMargin{ (i == Nx - 1) ? rMargin / hfactor : 0.0 };
            auto bottomMargin{ (j == 0) ? bMargin / vfactor : 0.0 };
            auto topMargin{ (j == Ny - 1) ? tMargin / vfactor : 0.0 };
            SetPadMarginAttribute(pad, leftMargin, rightMargin, bottomMargin, topMargin);
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
    auto pad{ dynamic_cast<TPad*>(canvas->FindObject(TString::Format("pad_%d_%d", id_x, id_y).Data())) };
    pad->Draw();
    pad->cd();
}

float ROOTHelper::GetPadXfactorInCanvasPartition(TCanvas * canvas, TVirtualPad * pad)
{
    auto pad_origin{ dynamic_cast<TPad*>(canvas->FindObject("pad_0_0")) };
    return pad_origin->GetAbsWNDC() / pad->GetAbsWNDC();
}

float ROOTHelper::GetPadYfactorInCanvasPartition(TCanvas * canvas, TVirtualPad * pad)
{
    auto pad_origin{ dynamic_cast<TPad*>(canvas->FindObject("pad_0_0")) };
    return pad_origin->GetAbsHNDC() / pad->GetAbsHNDC();
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