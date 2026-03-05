#include "DataObjectDispatch.hpp"

#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>

#include "AtomObject.hpp"
#include "BondObject.hpp"
#include "DataObjectBase.hpp"
#include "DataObjectVisitor.hpp"
#include "MapObject.hpp"
#include "ModelObject.hpp"

namespace rhbm_gem {
namespace {

[[noreturn]] void ThrowTypeMismatch(
    std::string_view context,
    const char * expected_type,
    const char * actual_type)
{
    throw std::runtime_error(
        std::string(context) + ": expected " + expected_type + " but got " + actual_type + ".");
}

class ExpectModelObjectVisitor : public ConstDataObjectVisitor
{
    std::string m_context;
    const ModelObject * m_model_object{ nullptr };

public:
    explicit ExpectModelObjectVisitor(std::string_view context) :
        m_context{ context }
    {
    }

    void VisitAtomObject(const AtomObject & data_object) override
    {
        (void)data_object;
        ThrowTypeMismatch(m_context, "ModelObject", "AtomObject");
    }

    void VisitBondObject(const BondObject & data_object) override
    {
        (void)data_object;
        ThrowTypeMismatch(m_context, "ModelObject", "BondObject");
    }

    void VisitModelObject(const ModelObject & data_object) override
    {
        m_model_object = &data_object;
    }

    void VisitMapObject(const MapObject & data_object) override
    {
        (void)data_object;
        ThrowTypeMismatch(m_context, "ModelObject", "MapObject");
    }

    const ModelObject & GetModelObject() const
    {
        if (m_model_object == nullptr)
        {
            throw std::runtime_error(
                m_context + ": expected ModelObject but got unsupported DataObject type.");
        }
        return *m_model_object;
    }
};

class ExpectMapObjectVisitor : public ConstDataObjectVisitor
{
    std::string m_context;
    const MapObject * m_map_object{ nullptr };

public:
    explicit ExpectMapObjectVisitor(std::string_view context) :
        m_context{ context }
    {
    }

    void VisitAtomObject(const AtomObject & data_object) override
    {
        (void)data_object;
        ThrowTypeMismatch(m_context, "MapObject", "AtomObject");
    }

    void VisitBondObject(const BondObject & data_object) override
    {
        (void)data_object;
        ThrowTypeMismatch(m_context, "MapObject", "BondObject");
    }

    void VisitModelObject(const ModelObject & data_object) override
    {
        (void)data_object;
        ThrowTypeMismatch(m_context, "MapObject", "ModelObject");
    }

    void VisitMapObject(const MapObject & data_object) override
    {
        m_map_object = &data_object;
    }

    const MapObject & GetMapObject() const
    {
        if (m_map_object == nullptr)
        {
            throw std::runtime_error(
                m_context + ": expected MapObject but got unsupported DataObject type.");
        }
        return *m_map_object;
    }
};

class CatalogTypeNameVisitor : public ConstDataObjectVisitor
{
    std::string m_type_name;

public:
    void VisitAtomObject(const AtomObject & data_object) override
    {
        (void)data_object;
        throw std::runtime_error("GetCatalogTypeName(): AtomObject is not a top-level catalog type.");
    }

    void VisitBondObject(const BondObject & data_object) override
    {
        (void)data_object;
        throw std::runtime_error("GetCatalogTypeName(): BondObject is not a top-level catalog type.");
    }

    void VisitModelObject(const ModelObject & data_object) override
    {
        (void)data_object;
        m_type_name = "model";
    }

    void VisitMapObject(const MapObject & data_object) override
    {
        (void)data_object;
        m_type_name = "map";
    }

    std::string GetTypeName() const
    {
        if (m_type_name.empty())
        {
            throw std::runtime_error("GetCatalogTypeName(): no catalog type resolved.");
        }
        return m_type_name;
    }
};

} // namespace

const ModelObject & ExpectModelObject(
    const DataObjectBase & data_object,
    std::string_view context)
{
    ExpectModelObjectVisitor visitor{ context };
    data_object.Accept(visitor, ModelVisitMode::SelfOnly);
    return visitor.GetModelObject();
}

ModelObject & ExpectModelObject(
    DataObjectBase & data_object,
    std::string_view context)
{
    const auto & const_ref{
        ExpectModelObject(static_cast<const DataObjectBase &>(data_object), context) };
    return const_cast<ModelObject &>(const_ref);
}

const MapObject & ExpectMapObject(
    const DataObjectBase & data_object,
    std::string_view context)
{
    ExpectMapObjectVisitor visitor{ context };
    data_object.Accept(visitor, ModelVisitMode::SelfOnly);
    return visitor.GetMapObject();
}

MapObject & ExpectMapObject(
    DataObjectBase & data_object,
    std::string_view context)
{
    const auto & const_ref{
        ExpectMapObject(static_cast<const DataObjectBase &>(data_object), context) };
    return const_cast<MapObject &>(const_ref);
}

std::string GetCatalogTypeName(const DataObjectBase & data_object)
{
    CatalogTypeNameVisitor visitor;
    data_object.Accept(visitor, ModelVisitMode::SelfOnly);
    return visitor.GetTypeName();
}

} // namespace rhbm_gem
