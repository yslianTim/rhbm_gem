#include "AtomPainter.hpp"
#include "AtomObject.hpp"
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

AtomPainter::AtomPainter(void) :
    m_folder_path{ "./" }
{

}

AtomPainter::~AtomPainter()
{

}

void AtomPainter::SetFolder(const std::string & folder_path)
{
    m_folder_path = FilePathHelper::EnsureTrailingSlash(folder_path);
}

void AtomPainter::AddDataObject(DataObjectBase * data_object)
{
    m_atom_object_list.push_back(dynamic_cast<AtomObject *>(data_object));
}

void AtomPainter::AddReferenceDataObject(DataObjectBase * data_object, const std::string & label)
{
    m_ref_atom_object_map[label] = dynamic_cast<AtomObject *>(data_object);
}

void AtomPainter::Painting(void)
{
    std::cout <<"- AtomPainter::Painting"<<std::endl;
    std::cout <<"  Folder path: "<< m_folder_path << std::endl;
    std::cout <<"  Number of atom objects to be painted: "<< m_atom_object_list.size() << std::endl;
    
    auto plot_name{ "sampling_map_value.pdf" };
    PaintSamplingMapValue(plot_name);
}

void AtomPainter::PaintSamplingMapValue(const std::string & name)
{
    auto file_path{ m_folder_path + name };

    #ifdef HAVE_ROOT

    auto canvas{ ROOTHelper::CreateCanvas("test","", 1000, 1000) };
    ROOTHelper::SetCanvasDefaultStyle(canvas.get());
    ROOTHelper::PrintCanvasOpen(canvas.get(), file_path);
    for (auto atom_object : m_atom_object_list)
    {
        auto atom_iter{ std::make_unique<PotentialEntryIterator>(atom_object) };
        auto data_array{ atom_iter->GetDistanceAndMapValueList() };

        std::vector<float> distance_array, map_value_array;
        distance_array.reserve(data_array.size());
        map_value_array.reserve(data_array.size());
        for (auto & data : data_array)
        {
            distance_array.emplace_back(std::get<0>(data));
            map_value_array.emplace_back(std::get<1>(data));
        }
        auto scaling{ 0.0 };
        auto map_value_range{ ArrayStats<float>::ComputeScalingRangeTuple(map_value_array, scaling) };
        auto distance_range{ ArrayStats<float>::ComputeScalingRangeTuple(distance_array, scaling) };
        std::cout << "  Map value range: " << std::get<0>(map_value_range) << " - " << std::get<1>(map_value_range) << std::endl;
        std::cout << "  Distance range: " << std::get<0>(distance_range) << " - " << std::get<1>(distance_range) << std::endl;

        auto map_value_hist{ ROOTHelper::CreateHist2D("hist_0","", 15, std::get<0>(distance_range), std::get<1>(distance_range), 100, std::get<0>(map_value_range), std::get<1>(map_value_range)) };
        for (auto & [distance, map_value] : data_array)
        {
            map_value_hist->Fill(distance, map_value);
        }

        map_value_hist->Draw("COLZ");

        ROOTHelper::PrintCanvasPad(canvas.get(), file_path);
    }
    ROOTHelper::PrintCanvasClose(canvas.get(), file_path);
    #endif
}

#ifdef HAVE_ROOT

#endif