#pragma once

class AtomObject;
class ModelObject;
class MapObject;
class DataObjectManager;

class DataObjectVisitorBase
{
public:
    virtual ~DataObjectVisitorBase() = default;
    virtual void VisitAtomObject(AtomObject * data_object) { (void)data_object; }
    virtual void VisitModelObject(ModelObject * data_object){ (void)data_object; }
    virtual void VisitMapObject(MapObject * data_object) { (void)data_object; }
};
