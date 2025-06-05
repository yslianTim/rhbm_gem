#include "ChargeAnalysisVisitor.hpp"
#include "AtomObject.hpp"
#include "ModelObject.hpp"
#include "MapObject.hpp"
#include "DataObjectManager.hpp"
#include "MapInterpolationVisitor.hpp"
#include "HRLModelHelper.hpp"
#include "ScopeTimer.hpp"
#include "AtomicChargeEntry.hpp"
#include "GroupChargeEntry.hpp"
#include "AtomicPotentialEntry.hpp"
#include "GroupPotentialEntry.hpp"
#include "AtomicInfoHelper.hpp"
#include "KeyPacker.hpp"
#include "ArrayStats.hpp"
#include "AtomClassifier.hpp"
#include "GlobalEnumClass.hpp"
#include "AminoAcidInfoHelper.hpp"
#include "ElectricPotential.hpp"

#include <iostream>
#include <tuple>
#include <fstream>
#include <unordered_set>

#ifdef HAVE_ROOT
#include "ROOTHelper.hpp"
#include <TStyle.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLine.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TColor.h>
#include <TMarker.h>
#include <TAxis.h>
#include <TH2.h>
#include <TF1.h>
#endif

ChargeAnalysisVisitor::ChargeAnalysisVisitor(void) :
    m_thread_size{ 1 },
    m_alpha_r{ 0.0 }, m_alpha_g{ 0.0 }, m_x_min{ 0.0 }, m_x_max{ 0.0 }
{

}

ChargeAnalysisVisitor::~ChargeAnalysisVisitor()
{

}

void ChargeAnalysisVisitor::VisitAtomObject(AtomObject * data_object)
{
    if (data_object == nullptr) return;
    if (data_object->GetAtomicPotentialEntry() == nullptr)
    {
        data_object->SetSelectedFlag(false);
        return;
    }
    data_object->SetSelectedFlag(true);
}

void ChargeAnalysisVisitor::VisitModelObject(ModelObject * data_object)
{
    ScopeTimer timer("ChargeAnalysisVisitor::VisitModelObject");
    if (data_object == nullptr) return;
    data_object->Update();
    m_selected_atom_list = data_object->GetSelectedAtomList();
    for (auto & atom : m_selected_atom_list)
    {
        auto atom_charge_entry{ std::make_unique<AtomicChargeEntry>() };
        atom->AddAtomicChargeEntry(std::move(atom_charge_entry));
    }
    std::cout <<" Number of selected atom = "<< m_selected_atom_list.size() << std::endl;
}

void ChargeAnalysisVisitor::VisitMapObject(MapObject * data_object)
{
    if (data_object == nullptr) return;
}

void ChargeAnalysisVisitor::Analysis(DataObjectManager * data_manager)
{
    std::cout <<"- Analysis..." << std::endl;
    try
    {
        const auto & model_object{ data_manager->GetDataObjectRef(m_model_key_tag) };
        model_object->Accept(this);

        const auto & neutral_model_object{ data_manager->GetDataObjectRef(m_neutral_model_key_tag) };
        const auto & positive_model_object{ data_manager->GetDataObjectRef(m_positive_model_key_tag) };
        const auto & negative_model_object{ data_manager->GetDataObjectRef(m_negative_model_key_tag) };

        RunChargeFitting(dynamic_cast<ModelObject*>(neutral_model_object.get()),
                         dynamic_cast<ModelObject*>(positive_model_object.get()),
                         dynamic_cast<ModelObject*>(negative_model_object.get()));
    }
    catch(const std::exception & e)
    {
        std::cerr << e.what() << std::endl;
    }
}

void ChargeAnalysisVisitor::RunChargeFitting(
    ModelObject * model_neutral, ModelObject * model_pos, ModelObject * model_neg)
{
    ScopeTimer timer("ChargeAnalysisVisitor::RunChargeFitting");
    auto electric_potential{ std::make_unique<ElectricPotential>() };
    electric_potential->SetModelChoice(1);
    electric_potential->SetBlurringWidth(0.50);

    auto gaus_potential{ std::make_unique<ElectricPotential>() };
    gaus_potential->SetModelChoice(2);

    auto neutral_atom_map{ BuildSerialIDAtomObjectMap(model_neutral) };
    auto positive_atom_map{ BuildSerialIDAtomObjectMap(model_pos) };
    auto negative_atom_map{ BuildSerialIDAtomObjectMap(model_neg) };

    // Atom Classification
    std::unordered_map<size_t, std::vector<AtomObject *>> atom_list_map;
    size_t type_id{ 999 };
    for (auto atom : m_selected_atom_list)
    {
        if (atom->GetSpecialAtomFlag() == true) continue;
        if (AtomClassifier::IsMainChainMember(atom->GetElement(), atom->GetRemoteness(), type_id) == false) continue;
        atom_list_map[type_id].emplace_back(atom);
    }

    std::unordered_map<size_t, std::vector<std::tuple<float, float, float>>> regression_data_list_map;

    

    // Group Charge Fitting
    for (const auto & [id, atom_list] : atom_list_map)
    {
        auto group_size{ atom_list.size() };
        Eigen::ArrayXf distance_array{ Eigen::ArrayXf::Zero(static_cast<int>(group_size)) };
        distance_array.setRandom();

        std::vector<std::tuple<std::vector<Eigen::VectorXd>, std::string>> data_array;
        std::vector<std::tuple<float, float, float>> regression_data_list;
        std::vector<Eigen::VectorXd> sampling_entry_list; // same as regression_data_list
        data_array.reserve(group_size);
        regression_data_list.reserve(group_size);
        sampling_entry_list.reserve(group_size);

        auto count{ 0 };
        for (auto atom : atom_list)
        {
            auto potential_entry{ atom->GetAtomicPotentialEntry() };
            if (potential_entry == nullptr)
            {
                std::cout <<"Warning: Atomic potential entry is nullptr for atom "
                        << atom->GetInfo() << ". This atom will be skipped." << std::endl;
                continue;
            }
            auto serial_id{ atom->GetSerialID() };
            auto atom_neutral{ neutral_atom_map.at(serial_id) };
            auto atom_positive{ positive_atom_map.at(serial_id) };
            auto atom_negative{ negative_atom_map.at(serial_id) };

            auto potential_entry_neutral{ atom_neutral->GetAtomicPotentialEntry() };
            auto potential_entry_positive{ atom_positive->GetAtomicPotentialEntry() };
            auto potential_entry_negative{ atom_negative->GetAtomicPotentialEntry() };

            auto charge{ AminoAcidInfoHelper::GetPartialCharge(
                atom->GetResidue(),
                atom->GetElement(),
                atom->GetRemoteness(),
                atom->GetBranch(),
                atom->GetStructure()) };
            auto is_negative_charged{ charge < 0.0 };
            
            auto distance{ distance_array(count) };
            auto func_phi_0{ electric_potential->GetPotentialValue(atom->GetElement(), distance, 0.0) };
            auto func_phi_pos{ electric_potential->GetPotentialValue(atom->GetElement(), distance, 1.0) };
            auto func_phi_neg{ electric_potential->GetPotentialValue(atom->GetElement(), distance, -1.0) };
            //std::cout << func_phi_0 << "\t" << func_phi_pos << "\t" << func_phi_neg << std::endl;

            //auto x0{ potential_entry_neutral->GetAmplitudeEstimateMDPDE() / std::pow(potential_entry_neutral->GetWidthEstimateMDPDE(), 3) };
            //auto x1_pos{ potential_entry_positive->GetAmplitudeEstimateMDPDE() / std::pow(potential_entry_positive->GetWidthEstimateMDPDE(), 3) };
            //auto x1_neg{ potential_entry_negative->GetAmplitudeEstimateMDPDE() / std::pow(potential_entry_negative->GetWidthEstimateMDPDE(), 3) };
            auto x0{ func_phi_0 };
            //auto x1_pos{ func_phi_pos - func_phi_0 };
            //auto x1_neg{ func_phi_0 - func_phi_neg };
            auto x1_pos{ func_phi_pos };
            auto x1_neg{ func_phi_neg };
            auto x1{ (is_negative_charged) ? -x1_neg : x1_pos };
            auto amplitude{ potential_entry->GetAmplitudeEstimateMDPDE() };
            auto width{ potential_entry->GetWidthEstimateMDPDE() };
            auto y{ gaus_potential->GetPotentialValue(atom->GetElement(), distance, 0.0, amplitude, width) };
            //auto y{ amplitude / std::pow(width, 3) };

            Eigen::VectorXd sampling_entry(4);
            sampling_entry(0) = 1.0;
            sampling_entry(1) = x0;
            sampling_entry(2) = x1;
            sampling_entry(3) = y;
            sampling_entry_list.emplace_back(sampling_entry);
            regression_data_list.emplace_back(std::make_tuple(x0, x1, y));
            data_array.emplace_back(std::make_tuple(sampling_entry_list, ""));

            count++;
        }
        regression_data_list_map[id] = regression_data_list;

        auto model_estimator{ std::make_unique<HRLModelHelper>(3, static_cast<int>(group_size)) };
        model_estimator->SetDataArray(data_array);
        model_estimator->RunEstimation(m_alpha_r, m_alpha_g);

        auto model_group_mean{ model_estimator->GetMuVectorMean() };
        auto model_group_mdpde{ model_estimator->GetMuVectorMDPDE() };
        auto prior_estimate{ model_estimator->GetMuVectorPrior() };
        auto prior_variance{ model_estimator->GetCapitalLambdaMatrix() };
        //std::cout <<"[ID-"<< id <<"] "<< model_group_mdpde.transpose() <<"\t"<< model_group_mdpde(2)/(model_group_mdpde(1)+model_group_mdpde(2)) << std::endl;
        std::cout <<"[ID-"<< id <<"] "<< regression_data_list.size() <<"\t"<< prior_estimate.transpose() <<"\t"<< prior_estimate(2)/(prior_estimate(1)+prior_estimate(2)) << std::endl;
        
    }
    PrintRegressionResult(regression_data_list_map);
}

void ChargeAnalysisVisitor::SetFitRange(double x_min, double x_max)
{
    m_x_min = x_min;
    m_x_max = x_max;
}

std::unordered_map<int, AtomObject *>
ChargeAnalysisVisitor::BuildSerialIDAtomObjectMap(ModelObject * model_object)
{
    std::unordered_map<int, AtomObject *> map;
    for (auto & atom : model_object->GetSelectedAtomList())
    {
        if (atom->GetAtomicPotentialEntry() == nullptr) continue;
        auto serial_id{ atom->GetSerialID() };
        if (map.find(serial_id) == map.end())
        {
            map[serial_id] = atom;
        }
        else
        {
            std::cerr << "Warning: Duplicate serial ID found: " << serial_id << std::endl;
        }
    }
    //std::cout << "Number of atoms in map: " << map.size() << std::endl;
    return map;
}

void ChargeAnalysisVisitor::PrintRegressionResult(
    const std::unordered_map<size_t, std::vector<std::tuple<float, float, float>>> & data_map)
{
    auto file_path{ "/Users/yslian/Downloads/regression_data.pdf" };

    const int col_size{ 2 };
    const std::string col_label[col_size]
    {
        "#font[2]{x}_{1}",
        "#font[2]{x}_{2}"
    };

    #ifdef HAVE_ROOT
    gStyle->SetGridColor(kGray);
    auto canvas{ ROOTHelper::CreateCanvas("test","", 700, 400) };
    ROOTHelper::SetCanvasDefaultStyle(canvas.get());
    ROOTHelper::SetCanvasPartition(canvas.get(), col_size, 1, 0.12f, 0.02f, 0.20f, 0.02f, 0.02f, 0.02f);
    ROOTHelper::PrintCanvasOpen(canvas.get(), file_path);

    std::unique_ptr<TH2> frame[col_size];
    for (int i = 0; i < col_size; i++)
    {
        ROOTHelper::FindPadInCanvasPartition(canvas.get(), i, 0);
        ROOTHelper::SetPadLayout(gPad, 1, 1, 0, 0, 0, 0);
        ROOTHelper::SetPadFrameAttribute(gPad, 0, 0, 4000, 0, 0, 0);
        auto x_factor{ ROOTHelper::GetPadXfactorInCanvasPartition(canvas.get(), gPad) };
        auto y_factor{ ROOTHelper::GetPadYfactorInCanvasPartition(canvas.get(), gPad) };
        frame[i] = ROOTHelper::CreateHist2D(Form("hist_%d_%d", i, 0),"", 100, 0.0, 1.0, 100, 0.0, 1.0);
        ROOTHelper::SetAxisTitleAttribute(frame[i]->GetXaxis(), 35.0f, 0.9f, 133);
        ROOTHelper::SetAxisLabelAttribute(frame[i]->GetXaxis(), 35.0f, 0.005f, 133);
        ROOTHelper::SetAxisTickAttribute(frame[i]->GetXaxis(), static_cast<float>(y_factor*0.05/x_factor), 505);
        ROOTHelper::SetAxisTitleAttribute(frame[i]->GetYaxis(), 40.0f, 1.0f, 133);
        ROOTHelper::SetAxisLabelAttribute(frame[i]->GetYaxis(), 35.0f, 0.01f, 133);
        ROOTHelper::SetAxisTickAttribute(frame[i]->GetYaxis(), static_cast<float>(x_factor*0.05/y_factor), 505);
        ROOTHelper::SetLineAttribute(frame[i].get(), 1, 0);
        frame[i]->SetStats(0);
        frame[i]->GetXaxis()->SetTitle(col_label[i].data());
        frame[i]->GetXaxis()->CenterTitle();
        frame[i]->GetYaxis()->CenterTitle();
    }

    std::vector<std::unique_ptr<TGraphErrors>> graph_list;
    for (auto & [id, data_array] : data_map)
    {
        auto label{ AtomClassifier::GetMainChainElementLabel(id) };
        

        std::unique_ptr<TGraphErrors> graph[col_size];
        graph[0] = ROOTHelper::CreateGraphErrors();
        graph[1] = ROOTHelper::CreateGraphErrors();
        auto count{ 0 };
        for (auto & [x1, x2, y] : data_array)
        {
            graph[0]->SetPoint(count, x1, y);
            graph[1]->SetPoint(count, x2, y);
            count++;
        }

        for (int i = 0; i < col_size; i++)
        {
            ROOTHelper::FindPadInCanvasPartition(canvas.get(), i, 0);
            double x_min, x_max, y_min, y_max;
            graph[i]->ComputeRange(x_min, y_min, x_max, y_max);
            auto x_range{ x_max - x_min };

            frame[i]->GetYaxis()->SetTitle(Form("#font[2]{y} (#font[102]{%s})", label.data()));
            frame[i]->GetXaxis()->SetLimits(x_min - 0.1*x_range, x_max + 0.1*x_range);
            frame[i]->GetYaxis()->SetLimits(y_min, y_max);
            frame[i]->Draw();

            ROOTHelper::SetMarkerAttribute(graph[i].get(), 20, 1.0f, kAzure-7);
            graph[i]->Draw("P");
            graph_list.emplace_back(std::move(graph[i]));
        }

        ROOTHelper::PrintCanvasPad(canvas.get(), file_path);
    }
    ROOTHelper::PrintCanvasClose(canvas.get(), file_path);
    #endif
}