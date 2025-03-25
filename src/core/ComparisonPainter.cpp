#include "ComparisonPainter.hpp"
#include "ModelObject.hpp"
#include "DataObjectBase.hpp"
#include "PotentialEntryIterator.hpp"
#include "FilePathHelper.hpp"
#include "AtomicInfoHelper.hpp"
#include "ArrayStats.hpp"

#ifdef HAVE_ROOT
#include "ROOTHelper.hpp"
#include <TStyle.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TColor.h>
#include <TMarker.h>
#include <TAxis.h>
#include <TH2.h>
#endif

#include <vector>
#include <tuple>

using ElementKeyType = GroupKeyMapping<ElementGroupClassifierTag>::type;
using ResidueKeyType = GroupKeyMapping<ResidueGroupClassifierTag>::type;

ComparisonPainter::ComparisonPainter(void) :
    m_folder_path{ "./" }
{

}

ComparisonPainter::~ComparisonPainter()
{

}

void ComparisonPainter::SetFolder(const std::string & folder_path)
{
    m_folder_path = FilePathHelper::EnsureTrailingSlash(folder_path);
}

void ComparisonPainter::AddDataObject(DataObjectBase * data_object)
{
    m_model_object_list.push_back(dynamic_cast<ModelObject *>(data_object));
}

void ComparisonPainter::AddReferenceDataObject(DataObjectBase * data_object, const std::string & label)
{
    m_ref_model_object_list_map[label].emplace_back(dynamic_cast<ModelObject *>(data_object));
}

void ComparisonPainter::Painting(void)
{
    std::cout <<"- ComparisonPainter::Painting"<<std::endl;
    std::cout <<"  Folder path: "<< m_folder_path << std::endl;
    std::cout <<"  Number of model objects to be painted: "<< m_model_object_list.size() << std::endl;
    std::cout <<"  Number of reference types to be painted: "<< m_ref_model_object_list_map.size() << std::endl;
}

void ComparisonPainter::PaintSimulationGaus(const std::string & name)
{
    auto file_path{ m_folder_path + name };
    std::cout <<"- ComparisonPainter::PaintSimulationGaus"<< std::endl;

    #ifdef HAVE_ROOT
    const int primary_element_size{ 4 };
    const char * element_label[primary_element_size]
    {
        "Alpha Carbon",
        "Carbonyl Carbon",
        "Peptide Nitrogen",
        "Carbonyl Oxygen"
    };
    int element_index[primary_element_size]   { 6, 6, 7, 8 };
    int remoteness_index[primary_element_size]{ 1, 0, 0, 0 };
    int color_element[primary_element_size]   { kRed+1, kOrange+1, kGreen+2, kAzure+2 };
    int marker_element[primary_element_size]  { 54, 53, 55, 59 };

    std::unique_ptr<TGraphErrors> amplitude_graph[primary_element_size];
    std::unique_ptr<TGraphErrors> width_graph[primary_element_size];
    std::vector<double> amplitude_array, width_array;
    amplitude_array.reserve(primary_element_size * AtomicInfoHelper::GetStandardResidueCount());
    width_array.reserve(primary_element_size * AtomicInfoHelper::GetStandardResidueCount());
    for (int i = 0; i < primary_element_size; i++)
    {
        auto element{ element_index[i] };
        auto remoteness{ remoteness_index[i] };
        amplitude_graph[i] = ROOTHelper::CreateGraphErrors();
        width_graph[i] = ROOTHelper::CreateGraphErrors();
        auto group_key{ std::make_tuple(element, remoteness, false) };
        auto count{ 0 };
        for (auto & model_object : m_model_object_list)
        {
            auto entry_iter{ std::make_unique<PotentialEntryIterator>(model_object) };
            if (entry_iter->IsAvailableGroupKey(group_key) == false) continue;
            auto gaus_estimate{ entry_iter->GetGausEstimatePrior(group_key) };
            auto gaus_variance{ entry_iter->GetGausVariancePrior(group_key) };
            amplitude_array.push_back(std::get<0>(gaus_estimate));
            width_array.push_back(std::get<1>(gaus_estimate));
            amplitude_graph[i]->SetPoint(count, count, std::get<0>(gaus_estimate));
            amplitude_graph[i]->SetPointError(count, 0.0, std::get<0>(gaus_variance));
            width_graph[i]->SetPoint(count, count, std::get<1>(gaus_estimate));
            width_graph[i]->SetPointError(count, 0.0, std::get<1>(gaus_variance));
            count++;
        }
        ROOTHelper::SetMarkerAttribute(amplitude_graph[i].get(), marker_element[i], 1.3, color_element[i]);
        ROOTHelper::SetMarkerAttribute(width_graph[i].get(), marker_element[i], 1.3, color_element[i]);
        ROOTHelper::SetLineAttribute(amplitude_graph[i].get(), 1, 1, color_element[i]);
        ROOTHelper::SetLineAttribute(width_graph[i].get(), 1, 1, color_element[i]);
    }
    #endif
}


#ifdef HAVE_ROOT


#endif