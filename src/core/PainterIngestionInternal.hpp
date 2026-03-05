#pragma once

#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>

#include "DataObjectBase.hpp"

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

template <typename ExpectedType, typename IngestMode, typename OnTypedDataObject>
void AddDataObject(
    DataObjectBase * data_object,
    IngestMode & ingest_mode,
    IngestMode data_mode,
    IngestMode reference_mode,
    std::string & ingest_label,
    std::string_view painter_name,
    OnTypedDataObject && on_typed_data_object)
{
    ingest_mode = data_mode;
    ingest_label.clear();
    if (data_object == nullptr)
    {
        ThrowInvalidTypeForMode(ingest_mode, reference_mode, painter_name);
    }
    auto typed_data_object{ dynamic_cast<ExpectedType *>(data_object) };
    if (typed_data_object == nullptr)
    {
        ThrowInvalidTypeForMode(ingest_mode, reference_mode, painter_name);
    }
    std::forward<OnTypedDataObject>(on_typed_data_object)(*typed_data_object);
}

template <typename ExpectedType, typename IngestMode, typename OnTypedDataObject>
void AddReferenceDataObject(
    DataObjectBase * data_object,
    const std::string & label,
    IngestMode & ingest_mode,
    IngestMode reference_mode,
    std::string & ingest_label,
    std::string_view painter_name,
    OnTypedDataObject && on_typed_data_object)
{
    ingest_mode = reference_mode;
    ingest_label = label;
    if (data_object == nullptr)
    {
        ThrowInvalidTypeForMode(ingest_mode, reference_mode, painter_name);
    }
    auto typed_data_object{ dynamic_cast<ExpectedType *>(data_object) };
    if (typed_data_object == nullptr)
    {
        ThrowInvalidTypeForMode(ingest_mode, reference_mode, painter_name);
    }
    std::forward<OnTypedDataObject>(on_typed_data_object)(*typed_data_object);
}

} // namespace rhbm_gem::painter_internal
