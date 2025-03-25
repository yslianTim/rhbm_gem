#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include "PainterBase.hpp"

class ModelObject;

#ifdef HAVE_ROOT
class TPad;
class TH2;
class TPaveText;
#endif

class ComparisonPainter : public PainterBase
{
    std::vector<ModelObject *> m_model_object_list;
    std::unordered_map<std::string, std::vector<ModelObject *>> m_ref_model_object_list_map;
    std::string m_folder_path;

public:
    ComparisonPainter(void);
    ~ComparisonPainter();
    void SetFolder(const std::string & folder_path) override;
    void AddDataObject(DataObjectBase * data_object) override;
    void AddReferenceDataObject(DataObjectBase * data_object, const std::string & label) override;
    void Painting(void) override;

private:
    void PaintSimulationGaus(const std::string & name);

    #ifdef HAVE_ROOT


    #endif

};