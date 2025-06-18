#pragma once

class AtomObject;
class ModelObject;
class MapObject;
class DataObjectManager;

class DataObjectVisitorBase
{
public:
    virtual ~DataObjectVisitorBase() = default;
    virtual void VisitAtomObject(AtomObject * data_object) = 0;
    virtual void VisitModelObject(ModelObject * data_object) = 0;
    virtual void VisitMapObject(MapObject * data_object) = 0;
    virtual void Analysis(DataObjectManager * data_manager) = 0;
};
