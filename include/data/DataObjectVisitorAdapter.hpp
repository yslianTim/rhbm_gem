#pragma once

#include "DataObjectVisitorBase.hpp"

class DataObjectVisitorAdapter : public DataObjectVisitorBase
{
public:
    ~DataObjectVisitorAdapter() override = default;
    void VisitAtomObject(AtomObject * data_object) override { (void)data_object; }
    void VisitModelObject(ModelObject * data_object) override { (void)data_object; }
    void VisitMapObject(MapObject * data_object) override { (void)data_object; }
    void VisitDataObjectManager(DataObjectManager * data_manager) override { (void)data_manager; }
};
