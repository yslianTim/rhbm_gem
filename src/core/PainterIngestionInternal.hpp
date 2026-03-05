#pragma once

#include <stdexcept>
#include <string>
#include <string_view>

#include "DataObjectBase.hpp"
#include "DataObjectVisitor.hpp"
#include "ModelVisitMode.hpp"

namespace rhbm_gem::painter_internal {

template <typename IngestMode>
[[noreturn]] void ThrowInvalidTypeForMode(
    IngestMode ingest_mode,
    IngestMode reference_mode,
    std::string_view painter_name)
{
    const std::string method_name{
        ingest_mode == reference_mode ? "AddReferenceDataObject" : "AddDataObject"
    };
    throw std::runtime_error(
        std::string(painter_name) + "::" + method_name + "(): invalid data_object type");
}

template <typename IngestMode>
void AddDataObject(
    DataObjectBase * data_object,
    DataObjectVisitor & visitor,
    IngestMode & ingest_mode,
    IngestMode data_mode,
    IngestMode reference_mode,
    std::string & ingest_label,
    std::string_view painter_name)
{
    ingest_mode = data_mode;
    ingest_label.clear();
    if (data_object == nullptr)
    {
        ThrowInvalidTypeForMode(ingest_mode, reference_mode, painter_name);
    }
    data_object->Accept(visitor, ModelVisitMode::SelfOnly);
}

template <typename IngestMode>
void AddReferenceDataObject(
    DataObjectBase * data_object,
    const std::string & label,
    DataObjectVisitor & visitor,
    IngestMode & ingest_mode,
    IngestMode reference_mode,
    std::string & ingest_label,
    std::string_view painter_name)
{
    ingest_mode = reference_mode;
    ingest_label = label;
    if (data_object == nullptr)
    {
        ThrowInvalidTypeForMode(ingest_mode, reference_mode, painter_name);
    }
    data_object->Accept(visitor, ModelVisitMode::SelfOnly);
}

} // namespace rhbm_gem::painter_internal
