#pragma once

namespace rhbm_gem {

class AtomObject;
class BondObject;
class ModelObject;
class MapObject;

class DataObjectVisitor
{
public:
    virtual ~DataObjectVisitor() = default;
    virtual void VisitAtomObject(AtomObject & data_object) = 0;
    virtual void VisitBondObject(BondObject & data_object) = 0;
    virtual void VisitModelObject(ModelObject & data_object) = 0;
    virtual void VisitMapObject(MapObject & data_object) = 0;
};

class ConstDataObjectVisitor
{
public:
    virtual ~ConstDataObjectVisitor() = default;
    virtual void VisitAtomObject(const AtomObject & data_object) = 0;
    virtual void VisitBondObject(const BondObject & data_object) = 0;
    virtual void VisitModelObject(const ModelObject & data_object) = 0;
    virtual void VisitMapObject(const MapObject & data_object) = 0;
};

} // namespace rhbm_gem
