#pragma once

#include <iostream>
#include "DataObjectVisitorBase.hpp"

class PotentialDisplayVisitor : public DataObjectVisitorBase
{
    std::string m_model_key_tag;

public:
    PotentialDisplayVisitor(void);
    ~PotentialDisplayVisitor();

    void VisitAtomObject(AtomObject * data_object) override;
    void VisitModelObject(ModelObject * data_object) override;
    void VisitMapObject(MapObject * data_object) override;
    void Analysis(DataObjectManager * data_manager) override;

    void SetModelObjectKeyTag(const std::string & value) { m_model_key_tag = value; }

};