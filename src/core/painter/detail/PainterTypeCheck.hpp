#pragma once

#include <stdexcept>
#include <string>
#include <string_view>

#include <rhbm_gem/data/object/DataObjectBase.hpp>

namespace rhbm_gem::painter_internal {

[[noreturn]] inline void ThrowInvalidPainterObjectType(
    std::string_view painter_name,
    std::string_view method_name)
{
    throw std::runtime_error(
        std::string(painter_name) + "::" + std::string(method_name)
        + "(): invalid data_object type");
}

template <typename ExpectedType>
ExpectedType & RequirePainterObject(
    DataObjectBase * data_object,
    std::string_view painter_name,
    std::string_view method_name)
{
    if (data_object == nullptr)
    {
        ThrowInvalidPainterObjectType(painter_name, method_name);
    }
    auto typed_data_object{ dynamic_cast<ExpectedType *>(data_object) };
    if (typed_data_object == nullptr)
    {
        ThrowInvalidPainterObjectType(painter_name, method_name);
    }
    return *typed_data_object;
}

} // namespace rhbm_gem::painter_internal
