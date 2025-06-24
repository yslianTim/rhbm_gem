#pragma once

#include <string>
#include <vector>
#include <unordered_map>
#include "PainterBase.hpp"

class AtomObject;

#ifdef HAVE_ROOT
class TPad;
class TH2;
class TPaveText;
#endif

class AtomPainter : public PainterBase
{
    std::vector<AtomObject *> m_atom_object_list;
    std::unordered_map<std::string, AtomObject *> m_ref_atom_object_map;
    std::string m_folder_path;

public:
    AtomPainter(void);
    ~AtomPainter();
    void SetFolder(const std::string & folder_path) override;
    void AddDataObject(DataObjectBase * data_object) override;
    void AddReferenceDataObject(DataObjectBase * data_object, const std::string & label) override;
    void Painting(void) override;

private:
    void PaintDemoPlot(const std::string & name);
    void PaintRegressionCheckPlot(const std::string & name);

};
