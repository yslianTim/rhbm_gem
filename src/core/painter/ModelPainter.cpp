#include "PainterFunctions.hpp"

#include <rhbm_gem/data/object/ModelAnalysisView.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/utils/domain/ChemicalDataHelper.hpp>
#include <rhbm_gem/utils/domain/ComponentHelper.hpp>
#include <rhbm_gem/utils/domain/GlobalEnumClass.hpp>
#include <rhbm_gem/utils/domain/AtomKeySystem.hpp>
#include <rhbm_gem/utils/domain/StringHelper.hpp>
#include <rhbm_gem/utils/domain/FilePathHelper.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/math/ArrayHelper.hpp>
#include <rhbm_gem/utils/hrl/LocalPotentialSeries.hpp>

#include "PotentialPlotBuilder.hpp"
#include "core/painter/AtomStyleCatalog.hpp"
#include "data/detail/AtomClassifier.hpp"
#include "detail/PainterModelValidation.hpp"
#include "detail/PainterOutputLabel.hpp"

#ifdef HAVE_ROOT
#include <rhbm_gem/utils/domain/ROOTHelper.hpp>
#include <TStyle.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TGraphErrors.h>
#include <TGraph2DErrors.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TColor.h>
#include <TMarker.h>
#include <TAxis.h>
#include <TH2.h>
#include <TF1.h>
#include <TLine.h>
#include <TH3F.h>
#include <TH1.h>
#endif

#include <vector>

namespace rhbm_gem {

namespace {

double ComputeAtomGausEstimateMinimum(
    const ModelObject & model_object,
    int par_id,
    Element element)
{
    std::vector<double> gaus_estimate_list;
    gaus_estimate_list.reserve(model_object.GetSelectedAtomCount());
    for (const auto * atom : model_object.GetSelectedAtoms())
    {
        if (atom->GetElement() != element) continue;
        gaus_estimate_list.emplace_back(
            AtomLocalPotentialView::RequireFor(*atom).GetEstimateMDPDE().GetDisplayParameter(par_id));
    }
    return array_helper::ComputeMin(gaus_estimate_list.data(), gaus_estimate_list.size());
}

} // namespace

class ModelPainter
{
    std::vector<ModelObject *> m_model_object_list;
    std::string m_folder_path{ "./" };

public:
    ModelPainter(
        const painter_internal::ModelObjectList & model_objects,
        const std::string & output_folder);
    void Run();

private:
    void AddModel(ModelObject & data_object);
    void PaintAtomGroupGausMainChain(ModelObject * model_object, const std::string & name);
    void PaintAtomGroupGausNucleotideMainChain(ModelObject * model_object, const std::string & name);
    void PaintAtomMapValueMainChain(ModelObject * model_object, const std::string & name);
    void PaintGroupWidthScatterPlot(
        ModelObject * model_object,
        const std::string & name,
        int par_id=0,
        bool draw_box_plot=false);
    void PaintAtomXYPosition(ModelObject * model_object, const std::string & name);
    void PaintAtomGausScatterPlot(
        ModelObject * model_object,
        const std::string & name,
        bool do_normalize=false);
    void PaintAtomGausMainChain(ModelObject * model_object, const std::string & name);
    void PaintAtomRankMainChain(ModelObject * model_object, const std::string & name);

#ifdef HAVE_ROOT
    void PrintAmplitudePad(::TPad * pad, ::TH2 * hist);
    void PrintWidthPad(::TPad * pad, ::TH2 * hist);
    void PrintAmplitudeSummaryPad(::TPad * pad, ::TH2 * hist);
    void PrintAtomWidthSummaryPad(::TPad * pad, ::TH2 * hist);
    void PrintGausSummaryPad(::TPad * pad, ::TH2 * hist);

    void ModifyAxisLabelSideChain(
        ::TPad * pad,
        ::TH2 * hist,
        Residue residue,
        const std::vector<std::string> & label_list);
#endif
};

ModelPainter::ModelPainter(
    const painter_internal::ModelObjectList & model_objects,
    const std::string & output_folder) :
    m_folder_path{ path_helper::EnsureTrailingSlash(output_folder) }
{
    for (auto * model_object : model_objects)
    {
        AddModel(*model_object);
    }
}

void ModelPainter::AddModel(ModelObject & data_object)
{
    painter_internal::RequireGroupedAnalyzedModel(data_object, "ModelPainter");
    m_model_object_list.push_back(&data_object);
}

void ModelPainter::Run()
{
    Logger::Log(LogLevel::Info, "ModelPainter::Run() called.");
    Logger::Log(LogLevel::Info, "Folder path: " + m_folder_path);
    Logger::Log(LogLevel::Info, "Number of atom objects to be painted: "
                + std::to_string(m_model_object_list.size()));

    for (auto model_object : m_model_object_list)
    {
        auto label{ painter_internal::BuildPainterOutputLabel(*model_object) };
        PaintAtomXYPosition(model_object, "atom_position_"+ label);
        PaintAtomGausScatterPlot(model_object, "atom_gaus_scatter_"+ label, false);
        PaintAtomMapValueMainChain(model_object, "atom_map_value_main_chain_"+ label);
        PaintAtomRankMainChain(model_object, "atom_rank_main_chain_"+ label);
    }
}

void ModelPainter::PaintAtomMapValueMainChain(ModelObject * model_object, const std::string & name)
{
    auto file_path{ m_folder_path + name };
    Logger::Log(LogLevel::Info, " ModelPainter::PaintAtomMapValueMainChain");
    (void)model_object;

    #ifdef HAVE_ROOT
    const auto & class_key{ ChemicalDataHelper::GetSimpleAtomClassKey() };
    auto entry_iter{ std::make_unique<ModelAnalysisView>(*model_object) };
    auto plot_builder{ std::make_unique<PotentialPlotBuilder>(model_object) };
    const int col_size{ 4 };
    const int row_size{ 1 };
    gStyle->SetLineScalePS(1.5);
    gStyle->SetGridColor(kGray);

    auto canvas{ root_helper::CreateCanvas("test","", 1000, 350) };
    root_helper::SetCanvasDefaultStyle(canvas.get());
    root_helper::SetCanvasPartition(canvas.get(), col_size, row_size, 0.10f, 0.02f, 0.20f, 0.20f, 0.01f, 0.01f);
    root_helper::PrintCanvasOpen(canvas.get(), file_path);

    const size_t main_chain_element_size{ 4 };

    std::vector<std::unique_ptr<TGraphErrors>> map_value_graph_list[main_chain_element_size];
    std::unique_ptr<TF1> gaus_function[main_chain_element_size];
    double amplitude_prior[main_chain_element_size];
    double width_prior[main_chain_element_size];
    std::vector<double> y_array;
    y_array.reserve(model_object->GetSelectedAtomCount());
        for (size_t k = 0; k < main_chain_element_size; k++)
    {
        auto group_key{ AtomClassifier::GetMainChainSimpleAtomClassGroupKey(k) };
        if (!entry_iter->HasAtomGroup(group_key, class_key)) continue;
        for (auto atom : entry_iter->GetAtomObjectList(group_key, class_key))
        {
            auto atom_plot_builder{ std::make_unique<PotentialPlotBuilder>(atom) };
            auto graph{ atom_plot_builder->CreateBinnedDistanceToMapValueGraph() };
            root_helper::SetLineAttribute(graph.get(), 1, 2, static_cast<short>(kAzure-7), 0.3f);
            map_value_graph_list[k].emplace_back(std::move(graph));
            const auto atom_view{ AtomLocalPotentialView::RequireFor(*atom) };
            auto map_value_range{
                local_potential_series::ComputeResponseRange(atom_view.GetSamplingEntries(), 0.0)
            };
            y_array.emplace_back(std::get<0>(map_value_range));
            y_array.emplace_back(std::get<1>(map_value_range));
        }
        gaus_function[k] = plot_builder->CreateAtomGroupGausFunctionPrior(group_key, class_key);
        amplitude_prior[k] = entry_iter->GetAtomGroupPrior(group_key, class_key).GetDisplayParameter(0);
        width_prior[k] = entry_iter->GetAtomGroupPrior(group_key, class_key).GetDisplayParameter(1);
    }

    auto y_range{ array_helper::ComputeScalingRangeTuple(y_array, 0.15) };
    auto x_min{ 0.01 };
    auto x_max{ 1.49 };
    auto y_min{ std::get<0>(y_range) };
    auto y_max{ std::get<1>(y_range) };
    auto line_ref{ root_helper::CreateLine(x_min, 0.0, x_max, 0.0) };

    std::unique_ptr<TH2> frame[col_size][row_size];
    std::unique_ptr<TPaveText> element_text[col_size];
    std::unique_ptr<TPaveText> result_text[col_size];
    for (int i = 0; i < col_size; i++)
    {
        for (int j = 0; j < row_size; j++)
        {
            root_helper::FindPadInCanvasPartition(canvas.get(), i, j);
            root_helper::SetPadLayout(gPad, 1, 1, 0, 0, 0, 0);
            root_helper::SetPadFrameAttribute(gPad, 0, 0, 4000, 0, 0, 0);
            auto x_factor{ root_helper::GetPadXfactorInCanvasPartition(canvas.get(), gPad) };
            auto y_factor{ root_helper::GetPadYfactorInCanvasPartition(canvas.get(), gPad) };
            frame[i][j] = root_helper::CreateHist2D(Form("hist_%d_%d", i, j),"", 100, x_min, x_max, 100, y_min, y_max);
            root_helper::SetAxisTitleAttribute(frame[i][j]->GetXaxis(), 0.0f);
            root_helper::SetAxisLabelAttribute(frame[i][j]->GetXaxis(), 35.0f, 0.005f, 133);
            root_helper::SetAxisTickAttribute(frame[i][j]->GetXaxis(), static_cast<float>(y_factor*0.05/x_factor), 505);
            root_helper::SetAxisTitleAttribute(frame[i][j]->GetYaxis(), 40.0f, 1.2f);
            root_helper::SetAxisLabelAttribute(frame[i][j]->GetYaxis(), 35.0f, 0.01f, 133);
            root_helper::SetAxisTickAttribute(frame[i][j]->GetYaxis(), static_cast<float>(x_factor*0.05/y_factor), 506);
            root_helper::SetLineAttribute(frame[i][j].get(), 1, 0);
            frame[i][j]->SetStats(0);
            frame[i][j]->GetXaxis()->SetTitle("");
            frame[i][j]->GetYaxis()->SetTitle("Map Value");
            frame[i][j]->Draw();

            for (auto & graph : map_value_graph_list[i])
            {
                graph->Draw("L X0");
            }

            root_helper::SetLineAttribute(line_ref.get(), 2, 1, kBlack);
            line_ref->Draw();

            root_helper::SetLineAttribute(gaus_function[i].get(), 2, 2, kRed);
            gaus_function[i]->Draw("SAME");

            element_text[i] = root_helper::CreatePaveText(0.70, 0.75, 0.99, 0.99, "nbNDC ARC", true);
            root_helper::SetPaveTextDefaultStyle(element_text[i].get());
            root_helper::SetPaveAttribute(element_text[i].get(), 0, 0.2);
            root_helper::SetTextAttribute(element_text[i].get(), 40.0f, 103, 22);
            root_helper::SetFillAttribute(element_text[i].get(), 1001, AtomStyleCatalog::GetMainChainElementColor(static_cast<size_t>(i)), 0.5f);
            element_text[i]->AddText(AtomStyleCatalog::GetMainChainElementLabel(static_cast<size_t>(i)).data());
            element_text[i]->Draw();

            result_text[i] = root_helper::CreatePaveText(0.05, 0.08, 0.95, 0.25, "nbNDC", true);
            root_helper::SetPaveTextDefaultStyle(result_text[i].get());
            root_helper::SetTextAttribute(result_text[i].get(), 20.0f, 133, 22, 0.0, kRed);
            root_helper::SetFillAttribute(result_text[i].get(), 4000);
            result_text[i]->AddText(Form("#font[2]{A} = %.2f  ;  #tau = %.2f", amplitude_prior[i], width_prior[i]));
            result_text[i]->Draw();
        }
    }

    canvas->cd();
    auto pad_title{ root_helper::CreatePad("pad_title","", 0.00, 0.80, 1.00, 1.00) };
    root_helper::SetPadDefaultStyle(pad_title.get());
    pad_title->Draw();
    pad_title->cd();

    auto subtitle1_text{ root_helper::CreatePaveText(0.10, 0.10, 0.25, 0.90, "nbNDC ARC", false) };
    root_helper::SetPaveTextDefaultStyle(subtitle1_text.get());
    root_helper::SetPaveAttribute(subtitle1_text.get(), 0, 0.2);
    root_helper::SetFillAttribute(subtitle1_text.get(), 1001, kAzure-7);
    root_helper::SetTextAttribute(subtitle1_text.get(), 40.0f, 133, 22, 0.0, kYellow-10);
    subtitle1_text->AddText(Form("%.2f #AA", model_object->GetResolution()));
    subtitle1_text->Draw();

    auto subtitle2_text{ root_helper::CreatePaveText(0.26, 0.10, 0.45, 0.90, "nbNDC ARC", false) };
    root_helper::SetPaveTextDefaultStyle(subtitle2_text.get());
    root_helper::SetPaveAttribute(subtitle2_text.get(), 0, 0.2);
    root_helper::SetFillAttribute(subtitle2_text.get(), 1001, kAzure-7, 0.5f);
    root_helper::SetTextAttribute(subtitle2_text.get(), 35.0f, 103, 22);
    subtitle2_text->AddText(Form("PDB-%s", model_object->GetPdbID().data()));
    subtitle2_text->Draw();

    auto subtitle3_text{ root_helper::CreatePaveText(0.46, 0.10, 0.68, 0.90, "nbNDC ARC", false) };
    root_helper::SetPaveTextDefaultStyle(subtitle3_text.get());
    root_helper::SetPaveAttribute(subtitle3_text.get(), 0, 0.2);
    root_helper::SetFillAttribute(subtitle3_text.get(), 1001, kAzure-7, 0.5f);
    root_helper::SetTextAttribute(subtitle3_text.get(), 35.0f, 103, 22);
    subtitle3_text->AddText(model_object->GetEmdID().data());
    subtitle3_text->Draw();

    auto legend{ root_helper::CreateLegend(0.69, 0.10, 0.99, 0.90, false) };
    root_helper::SetLegendDefaultStyle(legend.get());
    root_helper::SetFillAttribute(legend.get(), 4000);
    root_helper::SetTextAttribute(legend.get(), 25.0f, 133, 12, 0.0);
    legend->SetMargin(0.15f);
    legend->AddEntry(gaus_function[0].get(),
        "Gaussian Model #color[633]{#phi (#font[1]{A},#font[1]{#tau})}", "l");
    legend->AddEntry(map_value_graph_list[0].at(0).get(),
        "Map Value", "l");
    legend->Draw();

    canvas->cd();
    auto pad_bottom{ root_helper::CreatePad("pad_bottom","", 0.10, 0.00, 0.98, 0.11) };
    root_helper::SetPadDefaultStyle(pad_bottom.get());
    pad_bottom->Draw();
    pad_bottom->cd();
    auto x_title_text{ root_helper::CreatePaveText(0.00, 0.00, 1.00, 1.00, "nbNDC", false) };
    root_helper::SetPaveTextDefaultStyle(x_title_text.get());
    root_helper::SetFillAttribute(x_title_text.get(), 4000);
    root_helper::SetTextAttribute(x_title_text.get(), 35.0f, 133, 22);
    x_title_text->AddText("Radial Distance #[]{#AA}");
    x_title_text->Draw();

    root_helper::PrintCanvasPad(canvas.get(), file_path);
    root_helper::PrintCanvasClose(canvas.get(), file_path);
    Logger::Log(LogLevel::Info, " Output file: " + file_path);
    #endif
}

void ModelPainter::PaintAtomXYPosition(
    ModelObject * model_object, const std::string & name)
{
    auto file_path{ m_folder_path + name };
    Logger::Log(LogLevel::Info, " ModelPainter::PaintAtomXYPosition");
    auto entry_iter{ std::make_unique<ModelAnalysisView>(*model_object) };
    auto plot_builder{ std::make_unique<PotentialPlotBuilder>(model_object) };

    #ifdef HAVE_ROOT

    gStyle->SetLineScalePS(1.5);
    gStyle->SetGridColor(kGray);

    auto canvas{ root_helper::CreateCanvas("test","", 1000, 1000) };
    root_helper::SetCanvasDefaultStyle(canvas.get());
    root_helper::PrintCanvasOpen(canvas.get(), file_path);

    auto normalized_z_position{ 0.5 };
    auto z_ratio_window{ 0.1 };
    auto graph{ plot_builder->CreateAtomXYPositionTomographyGraph(
        normalized_z_position, z_ratio_window, true)
    };
    auto z_position{ model_object->GetModelPosition(2, normalized_z_position) };
    auto z_window{ model_object->GetModelLength(2) * z_ratio_window };
    auto com_pos{ model_object->GetCenterOfMassPosition() };
    std::vector<double> x_array, y_array;
    x_array.reserve(static_cast<size_t>(graph->GetN()));
    y_array.reserve(static_cast<size_t>(graph->GetN()));
    for (int p = 0; p < graph->GetN(); p++)
    {
        x_array.emplace_back(graph->GetPointX(p));
        y_array.emplace_back(graph->GetPointY(p));
    }

    auto x_range{ array_helper::ComputeScalingRangeTuple(x_array, 0.1) };
    double x_min{ std::get<0>(x_range) };
    double x_max{ std::get<1>(x_range) };

    auto y_range{ array_helper::ComputeScalingRangeTuple(y_array, 0.1) };
    double y_min{ std::get<0>(y_range) };
    double y_max{ std::get<1>(y_range) };

    canvas->cd();
    root_helper::SetPadMarginInCanvas(gPad, 0.15, 0.05, 0.12, 0.12);
    root_helper::SetPadLayout(gPad, 1, 1, 0, 0);

    auto frame{ root_helper::CreateHist2D("frame","", 500, x_min, x_max, 500, y_min, y_max) };
    root_helper::SetAxisTitleAttribute(frame->GetXaxis(), 50.0f, 1.1f);
    root_helper::SetAxisTitleAttribute(frame->GetYaxis(), 50.0f, 1.4f);
    root_helper::SetAxisLabelAttribute(frame->GetXaxis(), 45.0f, 0.01f);
    root_helper::SetAxisLabelAttribute(frame->GetYaxis(), 45.0f, 0.01f);
    frame->SetStats(0);
    frame->GetXaxis()->SetTitle("Atom Position #font[2]{x} #[]{#AA}");
    frame->GetYaxis()->SetTitle("Atom Position #font[2]{y} #[]{#AA}");
    frame->Draw();
    root_helper::SetMarkerAttribute(graph.get(), 6, 1.0f, kAzure-7);
    graph->Draw("P X0");

    auto pos_text{ root_helper::CreatePaveText(0.0, 0.0, 1.0, 1.0, "nbNDC", false) };
    root_helper::SetPaveTextMarginInCanvas(gPad, pos_text.get(), 0.45, 0.06, 0.83, 0.12);
    root_helper::SetPaveTextDefaultStyle(pos_text.get());
    root_helper::SetFillAttribute(pos_text.get(), 4000);
    root_helper::SetTextAttribute(pos_text.get(), 35, 133, 32);
    pos_text->AddText(
        Form("Atom Position #font[2]{z} = %.1f #pm %.1f #[]{#AA}",
        z_position - com_pos.at(2), z_window*0.5));
    pos_text->Draw();

    auto title_text{ root_helper::CreatePaveText(0.02, 0.89, 0.67, 0.98, "nbNDC ARC", false) };
    root_helper::SetPaveTextDefaultStyle(title_text.get());
    root_helper::SetPaveAttribute(title_text.get(), 0, 0.2);
    root_helper::SetLineAttribute(title_text.get(), 1, 0);
    root_helper::SetFillAttribute(title_text.get(), 1001, kAzure-7);
    root_helper::SetTextAttribute(title_text.get(), 47, 23, 22, 0.0, kYellow-10);
    title_text->AddText("Tomography of Atomic Model");
    title_text->Draw();

    auto model_text{ root_helper::CreatePaveText(0.69, 0.89, 0.95, 0.98, "nbNDC ARC", false) };
    root_helper::SetPaveTextDefaultStyle(model_text.get());
    root_helper::SetPaveAttribute(model_text.get(), 0, 0.2);
    root_helper::SetFillAttribute(model_text.get(), 1001, kAzure-7, 0.5f);
    root_helper::SetTextAttribute(model_text.get(), 50, 103, 22);
    model_text->AddText(Form("PDB-%s", model_object->GetPdbID().data()));
    model_text->Draw();

    root_helper::PrintCanvasPad(canvas.get(), file_path);
    root_helper::PrintCanvasClose(canvas.get(), file_path);
    Logger::Log(LogLevel::Info, " Output file: " + file_path);
    #endif
}

void ModelPainter::PaintAtomGausScatterPlot(
    ModelObject * model_object, const std::string & name, bool do_normalize)
{
    auto file_path{ m_folder_path + name };
    Logger::Log(LogLevel::Info, " ModelPainter::PaintAtomGausScatterPlot");
    (void)model_object;
    (void)do_normalize;

    #ifdef HAVE_ROOT
    auto entry_iter{ std::make_unique<ModelAnalysisView>(*model_object) };
    auto plot_builder{ std::make_unique<PotentialPlotBuilder>(model_object) };
    auto amplitude_min{ ComputeAtomGausEstimateMinimum(*model_object, 0, Element::OXYGEN) };

    gStyle->SetLineScalePS(1.5);
    gStyle->SetGridColor(kGray);

    auto canvas{ root_helper::CreateCanvas("test","", 1300, 1000) };
    root_helper::SetCanvasDefaultStyle(canvas.get());
    root_helper::PrintCanvasOpen(canvas.get(), file_path);

    std::unique_ptr<TH2> frame;
    std::unordered_map<Element, std::unique_ptr<TGraphErrors>> graph_map;
    std::vector<double> x_array, y_array;
    x_array.reserve(model_object->GetSelectedAtomCount());
    y_array.reserve(model_object->GetSelectedAtomCount());
    for (auto & [element_type, element_name] : ChemicalDataHelper::GetElementLabelMap())
    {
        auto graph
        {
            (do_normalize == true) ?
            plot_builder->CreateNormalizedAtomGausEstimateScatterGraph(element_type, amplitude_min) :
            plot_builder->CreateAtomGausEstimateScatterGraph(element_type)
        };
        auto atomic_number{ ChemicalDataHelper::GetAtomicNumber(element_type) };
        auto marker_size{ (atomic_number <= 8) ? 1.2f : 2.0f };
        short marker_color{ (ChemicalDataHelper::IsStandardElement(element_type)) ?
            ChemicalDataHelper::GetDisplayColor(element_type) : static_cast<short>(kRed)
        };
        auto transparency{ (atomic_number <= 8) ? 0.2f : 1.0f };
        root_helper::SetMarkerAttribute(graph.get(),
            ChemicalDataHelper::GetDisplayMarker(element_type), marker_size,
            marker_color, transparency);
        for (int p = 0; p < graph->GetN(); p++)
        {
            auto x{ graph->GetPointX(p) };
            auto y{ graph->GetPointY(p) };
            x_array.emplace_back(x);
            y_array.emplace_back(y);
        }
        if (graph->GetN() == 0) continue;
        graph_map.emplace(element_type, std::move(graph));
    }

    auto x_range{ array_helper::ComputeScalingPercentileRangeTuple(x_array, 0.2, 0.005, 0.995) };
    double x_min{ std::get<0>(x_range) };
    double x_max{ std::get<1>(x_range) };

    auto y_range{ array_helper::ComputeScalingPercentileRangeTuple(y_array, 0.2, 0.005, 0.995) };
    double y_min{ std::get<0>(y_range) };
    double y_max{ std::get<1>(y_range) };

    canvas->cd();
    root_helper::SetPadMarginInCanvas(gPad, 0.12, 0.20, 0.12, 0.10);
    root_helper::SetPadLayout(gPad, 1, 1, 0, 0);

    frame = root_helper::CreateHist2D("frame","", 500, x_min, x_max, 500, y_min, y_max);
    root_helper::SetAxisTitleAttribute(frame->GetXaxis(), 50.0f, 1.1f);
    root_helper::SetAxisTitleAttribute(frame->GetYaxis(), 50.0f, 1.4f);
    root_helper::SetAxisLabelAttribute(frame->GetXaxis(), 45.0f, 0.01f);
    root_helper::SetAxisLabelAttribute(frame->GetYaxis(), 45.0f, 0.01f);
    frame->SetStats(0);
    frame->GetXaxis()->SetTitle("Amplitude #font[2]{A}");
    frame->GetYaxis()->SetTitle("Width #font[2]{#tau}");
    frame->GetXaxis()->CenterTitle();
    frame->GetYaxis()->CenterTitle();
    frame->Draw();
    auto legend{ root_helper::CreateLegend(0.81, 0.12, 0.98, 0.90, false) };
    root_helper::SetLegendDefaultStyle(legend.get());
    root_helper::SetTextAttribute(legend.get(), 40.0f, 133, 12);
    root_helper::SetFillAttribute(legend.get(), 4000);
    for (int i = 1; i <= 90; i++)
    {
        auto element{ ChemicalDataHelper::GetElementFromAtomicNumber(i) };
        if (graph_map.find(element) == graph_map.end()) continue;
        auto & graph{ graph_map.at(element) };
        graph->Draw("P X0");
        legend->AddEntry(graph.get(),
                         Form("#font[102]{%s} (%d)",
                            ChemicalDataHelper::GetElementLabelMap().at(element).data(),
                            graph->GetN()), "p");
    }
    legend->Draw();

    auto title_text{ root_helper::CreatePaveText(0.12, 0.91, 0.98, 0.99, "nbNDC ARC", false) };
    root_helper::SetPaveTextDefaultStyle(title_text.get());
    root_helper::SetPaveAttribute(title_text.get(), 0, 0.2);
    root_helper::SetLineAttribute(title_text.get(), 1, 0);
    root_helper::SetFillAttribute(title_text.get(), 1001, kAzure-7);
    root_helper::SetTextAttribute(title_text.get(), 50.0f, 23, 22, 0.0, kYellow-10);
    title_text->AddText("Atom-Based Main/Side-Chain #font[2]{A}#minus#font[2]{#tau} Scatter Plot");
    title_text->Draw();

    root_helper::PrintCanvasPad(canvas.get(), file_path);
    root_helper::PrintCanvasClose(canvas.get(), file_path);
    Logger::Log(LogLevel::Info, " Output file: " + file_path);
    #endif
}

void ModelPainter::PaintAtomRankMainChain(
    ModelObject * model_object, const std::string & name)
{
    auto file_path{ m_folder_path + name };
    Logger::Log(LogLevel::Info, " ModelPainter::PaintAtomRankMainChain");
    auto entry_iter{ std::make_unique<ModelAnalysisView>(*model_object) };
    auto plot_builder{ std::make_unique<PotentialPlotBuilder>(model_object) };

    #ifdef HAVE_ROOT

    gStyle->SetLineScalePS(1.5);
    gStyle->SetGridColor(kGray);
    const int col_size{ 4 };
    const int row_size{ 3 };
    auto canvas{ root_helper::CreateCanvas("test","", 1000, 800) };
    root_helper::SetCanvasDefaultStyle(canvas.get());
    root_helper::SetCanvasPartition(
        canvas.get(), col_size, row_size, 0.08f, 0.05f, 0.10f, 0.10f, 0.02f, 0.01f);
    root_helper::PrintCanvasOpen(canvas.get(), file_path);

    const std::string y_title[row_size]
    {
        "Intensity",
        "Width",
        "Amplitude"
    };

    int chain_size;
    std::unique_ptr<TH2> frame[col_size][row_size];
    std::vector<std::unique_ptr<TH1D>> rank_hist_list[row_size][3];
    std::unique_ptr<TH1D> rank_reference[row_size];
    std::vector<Residue> veto_residues_list{  };
    rank_hist_list[2][0] = plot_builder->CreateMainChainAtomGausRankHistogram(0, chain_size, Residue::UNK, 2, veto_residues_list);
    rank_hist_list[1][0] = plot_builder->CreateMainChainAtomGausRankHistogram(1, chain_size, Residue::UNK, 1, veto_residues_list);
    rank_hist_list[0][0] = plot_builder->CreateMainChainAtomGausRankHistogram(2, chain_size, Residue::UNK, 0, veto_residues_list);
    rank_hist_list[2][1] = plot_builder->CreateMainChainAtomGausRankHistogram(0, chain_size, Residue::GLY);
    rank_hist_list[1][1] = plot_builder->CreateMainChainAtomGausRankHistogram(1, chain_size, Residue::GLY);
    rank_hist_list[0][1] = plot_builder->CreateMainChainAtomGausRankHistogram(2, chain_size, Residue::GLY);
    rank_hist_list[2][2] = plot_builder->CreateMainChainAtomGausRankHistogram(0, chain_size, Residue::PRO);
    rank_hist_list[1][2] = plot_builder->CreateMainChainAtomGausRankHistogram(1, chain_size, Residue::PRO);
    rank_hist_list[0][2] = plot_builder->CreateMainChainAtomGausRankHistogram(2, chain_size, Residue::PRO);
    
    for (int j = 0; j < row_size; j++)
    {
        for (auto & hist : rank_hist_list[j][0]) hist->Scale(100.0 / hist->GetEntries());
        for (auto & hist : rank_hist_list[j][1]) hist->Scale(100.0 / hist->GetEntries());
        for (auto & hist : rank_hist_list[j][2]) hist->Scale(100.0 / hist->GetEntries());
        rank_reference[j] = root_helper::CreateHist1D(Form("rank_reference_%d", j),"", 4, 0.5, 4.5);
        for (int p = 1; p <= rank_reference[j]->GetNbinsX(); p++)
        {
            rank_reference[j]->SetBinContent(p, 150.0);
        }
    }

    gStyle->SetTextFont(22);
    gStyle->SetPaintTextFormat(".1f");
    std::unique_ptr<TPaveText> element_text[col_size];
    std::unique_ptr<TPaveText> par_text[row_size];
    std::unique_ptr<TLegend> legend;
    for (int i = 0; i < col_size; i++)
    {
        auto element_color{ AtomStyleCatalog::GetMainChainElementColor(static_cast<size_t>(i)) };
        for (int j = 0; j < row_size; j++)
        {
            root_helper::FindPadInCanvasPartition(canvas.get(), i, j);
            root_helper::SetPadLayout(gPad, 1, 0, 0, 0, 0, 0);
            root_helper::SetPadFrameAttribute(gPad, 0, 0, 4000, 0, 0, 0, 0);
            auto x_factor{ root_helper::GetPadXfactorInCanvasPartition(canvas.get(), gPad) };
            auto y_factor{ root_helper::GetPadYfactorInCanvasPartition(canvas.get(), gPad) };
            frame[i][j] = root_helper::CreateHist2D(Form("frame_%d_%d", i, j),"", 100, 0.1, 4.9, 100, 0.01, 110.0);
            root_helper::SetAxisTitleAttribute(frame[i][j]->GetXaxis(), 45.0f, 0.8f, 133);
            root_helper::SetAxisLabelAttribute(frame[i][j]->GetXaxis(), 40.0f, 0.01f, 103, kCyan+3);
            root_helper::SetAxisTickAttribute(frame[i][j]->GetXaxis(), static_cast<float>(y_factor*0.05/x_factor), 5);
            root_helper::SetAxisTitleAttribute(frame[i][j]->GetYaxis(), 28.0f, 1.3f, 133);
            root_helper::SetAxisLabelAttribute(frame[i][j]->GetYaxis(), 25.0f, 0.005f, 133);
            root_helper::SetAxisTickAttribute(frame[i][j]->GetYaxis(), static_cast<float>(x_factor*0.05/y_factor), 505);
            root_helper::SetLineAttribute(frame[i][j].get(), 1, 0);
            frame[i][j]->GetXaxis()->SetTitle("Rank");
            frame[i][j]->GetYaxis()->SetTitle("Count Ratio (%)");
            frame[i][j]->GetXaxis()->CenterTitle();
            frame[i][j]->GetYaxis()->CenterTitle();
            frame[i][j]->SetStats(0);
            frame[i][j]->Draw();

            root_helper::SetFillAttribute(rank_reference[j].get(), 1001, kGray, 0.2f);
            rank_reference[j]->SetBarWidth(0.9f);
            rank_reference[j]->SetBarOffset(0.05f);
            rank_reference[j]->Draw("BAR SAME");

            auto text_size{ 0.2f * 4.0f/static_cast<float>(gPad->GetAbsHNDC()) };
            root_helper::SetFillAttribute(rank_hist_list[j][0].at(static_cast<size_t>(i)).get(), 1001, kAzure-7, 0.5f);
            root_helper::SetLineAttribute(rank_hist_list[j][0].at(static_cast<size_t>(i)).get(), 1, 0);
            root_helper::SetMarkerAttribute(rank_hist_list[j][0].at(static_cast<size_t>(i)).get(), 1, text_size, kAzure-7);
            rank_hist_list[j][0].at(static_cast<size_t>(i))->SetBarWidth(0.25f);
            rank_hist_list[j][0].at(static_cast<size_t>(i))->SetBarOffset(0.125f);
            rank_hist_list[j][0].at(static_cast<size_t>(i))->Draw("BAR TEXT0 SAME");

            root_helper::SetFillAttribute(rank_hist_list[j][1].at(static_cast<size_t>(i)).get(), 1001, kOrange-6, 0.5f);
            root_helper::SetLineAttribute(rank_hist_list[j][1].at(static_cast<size_t>(i)).get(), 1, 0);
            root_helper::SetMarkerAttribute(rank_hist_list[j][1].at(static_cast<size_t>(i)).get(), 1, text_size, kOrange-6);
            rank_hist_list[j][1].at(static_cast<size_t>(i))->SetBarWidth(0.25f);
            rank_hist_list[j][1].at(static_cast<size_t>(i))->SetBarOffset(0.375f);
            rank_hist_list[j][1].at(static_cast<size_t>(i))->Draw("BAR TEXT0 SAME");

            root_helper::SetFillAttribute(rank_hist_list[j][2].at(static_cast<size_t>(i)).get(), 1001, kMagenta-2, 0.5f);
            root_helper::SetLineAttribute(rank_hist_list[j][2].at(static_cast<size_t>(i)).get(), 1, 0);
            root_helper::SetMarkerAttribute(rank_hist_list[j][2].at(static_cast<size_t>(i)).get(), 1, text_size, kMagenta-2);
            rank_hist_list[j][2].at(static_cast<size_t>(i))->SetBarWidth(0.25f);
            rank_hist_list[j][2].at(static_cast<size_t>(i))->SetBarOffset(0.625f);
            rank_hist_list[j][2].at(static_cast<size_t>(i))->Draw("BAR TEXT0 SAME");

            if (j == row_size - 1)
            {
                element_text[i] = root_helper::CreatePaveText(0.01, 1.01, 0.99, 1.20, "nbNDC ARC", true);
                root_helper::SetPaveTextDefaultStyle(element_text[i].get());
                root_helper::SetPaveAttribute(element_text[i].get(), 0, 0.2);
                root_helper::SetTextAttribute(element_text[i].get(), 30.0f, 133, 22);
                root_helper::SetFillAttribute(element_text[i].get(), 1001, element_color, 0.5f);
                element_text[i]->AddText(AtomStyleCatalog::GetMainChainElementLabel(static_cast<size_t>(i)).data());
                element_text[i]->Draw();
            }
            if (i == col_size - 1)
            {
                par_text[j] = root_helper::CreatePaveText(1.01, 0.01, 1.22, 0.99, "nbNDC ARC", true);
                root_helper::SetPaveTextDefaultStyle(par_text[j].get());
                root_helper::SetPaveAttribute(par_text[j].get(), 0, 0.2);
                root_helper::SetTextAttribute(par_text[j].get(), 30.0f, 133, 22, 0.0f, kYellow-10);
                root_helper::SetLineAttribute(par_text[j].get(), 1, 0);
                root_helper::SetFillAttribute(par_text[j].get(), 1001, kAzure-7);
                par_text[j]->AddText(y_title[j].data());
                auto text{ par_text[j]->GetLineWith(y_title[j].data()) };
                text->SetTextAngle(90.0f);
                par_text[j]->Draw();
            }
        }
    }
    root_helper::PrintCanvasPad(canvas.get(), file_path);
    root_helper::PrintCanvasClose(canvas.get(), file_path);
    Logger::Log(LogLevel::Info, " Output file: " + file_path);
    #endif
}


} // namespace rhbm_gem

namespace rhbm_gem {

#ifdef HAVE_ROOT

void ModelPainter::PrintAmplitudePad(TPad * pad, TH2 * hist)
{
    pad->cd();
    root_helper::SetPadMarginInCanvas(pad, 0.07, 0.005, 0.01, 0.01);
    root_helper::SetPadLayout(pad, 1, 1, 0, 0);
    root_helper::SetAxisTitleAttribute(hist->GetXaxis(), 0.0f);
    root_helper::SetAxisTitleAttribute(hist->GetYaxis(), 35.0f, 1.4f);
    root_helper::SetAxisLabelAttribute(hist->GetXaxis(), 0.0f);
    root_helper::SetAxisLabelAttribute(hist->GetYaxis(), 30.0f, 0.01f);

    auto x_tick_length{ root_helper::ConvertGlobalTickLengthToPadTickLength(pad, 0.0, 0) };
    auto y_tick_length{ root_helper::ConvertGlobalTickLengthToPadTickLength(pad, 0.015, 1) };
    root_helper::SetAxisTickAttribute(hist->GetXaxis(), static_cast<float>(x_tick_length), 21);
    root_helper::SetAxisTickAttribute(hist->GetYaxis(), static_cast<float>(y_tick_length), 506);
    hist->GetXaxis()->SetLimits(-1.0, 20.0);
    hist->SetStats(0);
    hist->GetYaxis()->SetTitle("Amplitude");
    hist->GetYaxis()->CenterTitle();
    hist->Draw();
}

void ModelPainter::PrintWidthPad(TPad * pad, TH2 * hist)
{
    pad->cd();
    root_helper::SetPadMarginInCanvas(pad, 0.07, 0.005, 0.11, 0.01);
    root_helper::SetPadLayout(pad, 1, 1, 0, 0);
    root_helper::SetAxisTitleAttribute(hist->GetXaxis(), 0.0);
    root_helper::SetAxisTitleAttribute(hist->GetYaxis(), 35.0, 1.4f);
    root_helper::SetAxisLabelAttribute(hist->GetXaxis(), 32.0, 0.11f, 103, kCyan+3);
    root_helper::SetAxisLabelAttribute(hist->GetYaxis(), 30.0, 0.01f);

    auto x_tick_length{ root_helper::ConvertGlobalTickLengthToPadTickLength(pad, 0.0, 0) };
    auto y_tick_length{ root_helper::ConvertGlobalTickLengthToPadTickLength(pad, 0.015, 1) };
    root_helper::SetAxisTickAttribute(hist->GetXaxis(), static_cast<float>(x_tick_length), 21);
    root_helper::SetAxisTickAttribute(hist->GetYaxis(), static_cast<float>(y_tick_length), 505);
    hist->GetXaxis()->SetLimits(-1.0, 20.0);
    hist->GetXaxis()->ChangeLabel(1, -1.0, 0.0);
    hist->GetXaxis()->ChangeLabel(-1, -1.0, 0.0);
    for (size_t i = 0; i < ChemicalDataHelper::GetStandardAminoAcidCount(); i++)
    {
        auto residue{ ChemicalDataHelper::GetStandardAminoAcidList().at(i) };
        auto label{ ChemicalDataHelper::GetLabel(residue) };
        auto label_index{ static_cast<int>(i) + 2 };
        hist->GetXaxis()->ChangeLabel(label_index, 90.0, -1, 12, -1, -1, label.data());
    }

    hist->SetStats(0);
    hist->GetYaxis()->SetTitle("Width");
    hist->GetYaxis()->CenterTitle();
    hist->Draw();
}

void ModelPainter::PrintAmplitudeSummaryPad(TPad * pad, TH2 * hist)
{
    pad->cd();
    root_helper::SetPadMarginInCanvas(pad, 0.005, 0.005, 0.01, 0.01);
    root_helper::SetPadLayout(pad, 1, 1, 0, 0);
    root_helper::SetAxisTitleAttribute(hist->GetXaxis(), 0.0);
    root_helper::SetAxisTitleAttribute(hist->GetYaxis(), 0.0);
    root_helper::SetAxisLabelAttribute(hist->GetXaxis(), 0.0);
    root_helper::SetAxisLabelAttribute(hist->GetYaxis(), 0.0);

    auto x_tick_length{ root_helper::ConvertGlobalTickLengthToPadTickLength(pad, 0.0, 0) };
    auto y_tick_length{ root_helper::ConvertGlobalTickLengthToPadTickLength(pad, 0.015, 1) };
    root_helper::SetAxisTickAttribute(hist->GetXaxis(), static_cast<float>(x_tick_length), 5);
    root_helper::SetAxisTickAttribute(hist->GetYaxis(), static_cast<float>(y_tick_length), 506);
    hist->GetXaxis()->SetLimits(-1.0, 4.0);
    hist->SetStats(0);
    hist->Draw();
}

void ModelPainter::PrintAtomWidthSummaryPad(TPad * pad, TH2 * hist)
{
    pad->cd();
    root_helper::SetPadMarginInCanvas(pad, 0.005, 0.005, 0.11, 0.01);
    root_helper::SetPadLayout(pad, 1, 1, 0, 0);
    root_helper::SetAxisTitleAttribute(hist->GetXaxis(), 35.0f, 1.0f);
    root_helper::SetAxisTitleAttribute(hist->GetYaxis(), 0.0f);
    root_helper::SetAxisLabelAttribute(hist->GetXaxis(), 35.0f, 0.005f, 103);
    root_helper::SetAxisLabelAttribute(hist->GetYaxis(), 0.0f);

    auto x_tick_length{ root_helper::ConvertGlobalTickLengthToPadTickLength(pad, 0.0, 0) };
    auto y_tick_length{ root_helper::ConvertGlobalTickLengthToPadTickLength(pad, 0.015, 1) };
    root_helper::SetAxisTickAttribute(hist->GetXaxis(), static_cast<float>(x_tick_length), 5);
    root_helper::SetAxisTickAttribute(hist->GetYaxis(), static_cast<float>(y_tick_length), 505);
    hist->GetXaxis()->SetLimits(-1.0, 4.0);
    hist->GetXaxis()->ChangeLabel(1, -1.0, 0.0);
    hist->GetXaxis()->ChangeLabel(-1, -1.0, 0.0);

    for (size_t i = 0; i < 4; i++)
    {
        auto index{ static_cast<int>(i) + 2 };
        auto color{ AtomStyleCatalog::GetMainChainElementColor(i) };
        auto label{ AtomStyleCatalog::GetMainChainElementLabel(i) };
        hist->GetXaxis()->ChangeLabel(index, 0.0, -1, -1, color, -1, label.data());
    }
    hist->SetStats(0);
    hist->GetXaxis()->SetTitle("Element");
    hist->GetXaxis()->CenterTitle();
    hist->Draw();
}

void ModelPainter::PrintGausSummaryPad(TPad * pad, TH2 * hist)
{
    pad->cd();
    root_helper::SetPadMarginInCanvas(pad, 0.005, 0.07, 0.11, 0.01);
    root_helper::SetPadLayout(pad, 1, 1, 0, 0);
    root_helper::SetAxisTitleAttribute(hist->GetXaxis(), 35.0f, 1.0f);
    root_helper::SetAxisTitleAttribute(hist->GetYaxis(), 35.0f, 1.5f);
    root_helper::SetAxisLabelAttribute(hist->GetXaxis(), 35.0f);
    root_helper::SetAxisLabelAttribute(hist->GetYaxis(), 35.0f, 0.01f);

    auto x_tick_length{ root_helper::ConvertGlobalTickLengthToPadTickLength(pad, 0.03, 0) };
    auto y_tick_length{ root_helper::ConvertGlobalTickLengthToPadTickLength(pad, 0.015, 1) };
    root_helper::SetAxisTickAttribute(hist->GetXaxis(), static_cast<float>(x_tick_length), 505);
    root_helper::SetAxisTickAttribute(hist->GetYaxis(), static_cast<float>(y_tick_length), 505);
    hist->SetStats(0);
    hist->GetXaxis()->SetTitle("Amplitude Estimate");
    hist->GetYaxis()->SetTitle("Width Estimate");
    hist->GetXaxis()->CenterTitle();
    hist->GetYaxis()->CenterTitle();
    hist->Draw("Y+");
}

void ModelPainter::ModifyAxisLabelSideChain(
    TPad * pad, TH2 * hist, Residue residue, const std::vector<std::string> & label_list)
{
    if (ChemicalDataHelper::IsStandardResidue(residue) == false) return;

    auto x_tick_length{ root_helper::ConvertGlobalTickLengthToPadTickLength(pad, 0.0, 0) };
    auto y_tick_length{ root_helper::ConvertGlobalTickLengthToPadTickLength(pad, 0.015, 1) };
    auto label_size{ static_cast<int>(label_list.size()) };
    root_helper::SetAxisTickAttribute(hist->GetXaxis(), static_cast<float>(x_tick_length), label_size+1);
    root_helper::SetAxisTickAttribute(hist->GetYaxis(), static_cast<float>(y_tick_length), 506);

    hist->GetXaxis()->SetLimits(-1.0, static_cast<double>(label_list.size()));
    hist->GetXaxis()->ChangeLabel(1, -1.0, 0.0);
    hist->GetXaxis()->ChangeLabel(-1, -1.0, 0.0);

    for (int i = 0; i < label_size; i++)
    {
        hist->GetXaxis()->ChangeLabel(i+2, 0.0, -1, -1, -1, -1, label_list.at(static_cast<size_t>(i)).data());
    }
}

#endif

namespace painter_internal {

void PaintModel(const ModelObjectList & model_objects, const std::string & output_folder)
{
    ModelPainter painter{ model_objects, output_folder };
    painter.Run();
}

} // namespace painter_internal

} // namespace rhbm_gem
