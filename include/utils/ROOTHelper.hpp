#pragma once

#include <iostream>
#include <memory>
#include <string>

#ifdef HAVE_ROOT
class TCanvas;
class TPad;
class TH2D;
class TAttPad;
class TAttAxis;
class TAttText;
class TAttLine;
class TAttFill;
class TAttMarker;
class TPave;
class TBox;
class TPaveText;
class TLegend;
class TVirtualPad;
class TGraphErrors;
#endif

class ROOTHelper
{

public:
    ROOTHelper(void) = default;
    ~ROOTHelper() = default;
    #ifdef HAVE_ROOT
    static std::unique_ptr<TCanvas> CreateCanvas(const std::string & name, const std::string & title, int width=600, int height=600);
    static std::unique_ptr<TPad> CreatePad(const std::string & name, const std::string & title, double x_low, double y_low, double x_up, double y_up);
    static std::unique_ptr<TH2D> CreateHist2D(const std::string & name, const std::string & title, int x_bin, double x_min, double x_max, int y_bin, double y_min, double y_max);
    static std::unique_ptr<TGraphErrors> CreateGraphErrors(void);
    static std::unique_ptr<TGraphErrors> CreateGraphErrors(const int & point_size);
    static std::unique_ptr<TPaveText> CreatePaveText(double x1, double y1, double x2, double y2, std::string option="nbNDC", bool in_partition=false);
    static void PrintCanvasOpen(TCanvas * canvas, const std::string & name);
    static void PrintCanvasPad(TCanvas * canvas, const std::string & name);
    static void PrintCanvasClose(TCanvas * canvas, const std::string & name);
    static void SetPadFrameAttribute(
        TAttPad * pad, int mode=0, short size=0, short fill_style=4000, short fill_color=0,
        short line_style=0, short line_width=1, short line_color=1);
    static void SetPadMarginAttribute(
        TAttPad * pad, float left, float right, float bottom, float top);
    static void SetPadMarginInCanvas(TPad * pad, double left, double right, double bottom, double top);
    static void SetPaveTextMarginInCanvas(TPad * pad, TPaveText * pave, double left, double right, double bottom, double top);
    static void SetAxisTitleAttribute(
        TAttAxis * axis, float size, float offset=1.0f, short font=133, short color=1);
    static void SetAxisLabelAttribute(
        TAttAxis * axis, float size, float offset=0.005f, short font=133, short color=1);
    static void SetAxisTickAttribute(
        TAttAxis * axis, float length=0.03, int division=510);
    static void SetTextAttribute(
        TAttText * _text,
        float size=1.0f, short font=133, short align=12, float angle=0.0f, short color=1, float transparent=1.0f);
    static void SetLineAttribute(
        TAttLine * line, short style=1, short width=1, short color=1, float transparent=1.0f);
    static void SetFillAttribute(
        TAttFill * fill, short style=4000, short color=0, float transparent=1.0f);
    static void SetMarkerAttribute(
        TAttMarker * marker, short style, float size=1.0f, short color=1, float transparent=1.0f);
    static void SetPaveAttribute(
        TPave * pave, int border_size=0, double radius=0.2);
    static void SetBoxAttribute(
        TBox * box, double left, double right, double bottom, double top);
    static void SetPadInCanvas(TPad * pad, int width, int height, int division_x=1, int divition_y=1, float margin_x=0.0f, float margin_y=0.0f);
    static void SetPadDefaultStyle(TVirtualPad * pad);
    static void SetPadRangeInCanvas(TVirtualPad * pad, double x1, double y1, double x2, double y2);
    static void SetPadLayout(
        TVirtualPad * pad, int grid_x=0, int grid_y=0, int log_x=0, int log_y=0, int tick_x=0, int tick_y=1);
    static void SetCanvasDefaultStyle(TCanvas * canvas);
    static void SetPaveTextDefaultStyle(TPaveText * text);
    static void SetLegendDefaultStyle(TLegend * legend);
    static void SetCanvasPartition(
        TCanvas * canvas, const int Nx, const int Ny,
        float lMargin, float rMargin, float bMargin, float tMargin,
        float vSpacing = 0.0, float hSpacing = 0.0);
    static void FindPadInCanvasPartition(TCanvas * canvas, int id_x, int id_y);
    static float GetPadXfactorInCanvasPartition(TCanvas * canvas, TVirtualPad * pad);
    static float GetPadYfactorInCanvasPartition(TCanvas * canvas, TVirtualPad * pad);
    static double GetXtoPadInCanvasPartition(double x);
    static double GetYtoPadInCanvasPartition(double y);
    static double ConvertGlobalTickLengthToPadTickLength(TPad * pad, double global_tick_length, bool use_width=true);
    #endif

};