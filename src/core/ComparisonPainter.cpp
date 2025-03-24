#include "ComparisonPainter.hpp"
#include "ModelObject.hpp"
#include "DataObjectBase.hpp"
#include "GroupPotentialEntry.hpp"
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




#ifdef HAVE_ROOT


#endif