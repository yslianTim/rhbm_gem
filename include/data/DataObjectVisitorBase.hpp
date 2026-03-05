#pragma once

#include <stdexcept>
#include <string>

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

class StrictDataObjectVisitorBase : public DataObjectVisitorBase
{
public:
    ~StrictDataObjectVisitorBase() override = default;

    void VisitAtomObject(AtomObject * data_object) override
    {
        (void)data_object;
        throw std::logic_error(
            "StrictDataObjectVisitorBase::VisitAtomObject() is not implemented.");
    }

    void VisitBondObject(BondObject * data_object) override
    {
        (void)data_object;
        throw std::logic_error(
            "StrictDataObjectVisitorBase::VisitBondObject() is not implemented.");
    }

    void VisitModelObject(ModelObject * data_object) override
    {
        (void)data_object;
        throw std::logic_error(
            "StrictDataObjectVisitorBase::VisitModelObject() is not implemented.");
    }

    void VisitMapObject(MapObject * data_object) override
    {
        (void)data_object;
        throw std::logic_error(
            "StrictDataObjectVisitorBase::VisitMapObject() is not implemented.");
    }
};

} // namespace rhbm_gem
