#pragma once

namespace rhbm_gem {

class AtomObject;
class BondObject;
class ModelObject;
class MapObject;

class DataObjectVisitorBase
{
public:
    virtual ~DataObjectVisitorBase() = default;
    virtual void VisitAtomObject(AtomObject * data_object) { (void)data_object; }
    virtual void VisitBondObject(BondObject * data_object) { (void)data_object; }
    virtual void VisitModelObject(ModelObject * data_object){ (void)data_object; }
    virtual void VisitMapObject(MapObject * data_object) { (void)data_object; }
};

} // namespace rhbm_gem
