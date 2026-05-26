#pragma once

#include <string>

#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/domain/StringHelper.hpp>

namespace rhbm_gem::painter_internal {

inline std::string BuildPainterOutputLabel(const ModelObject & model_object)
{
    auto is_simulation{ model_object.GetEmdID().find("Simulation") != std::string::npos };
    auto label{ model_object.GetPdbID() };
    if (is_simulation)
    {
        label += "_" + model_object.GetEmdID()
                 + "_bw" + string_helper::ToStringWithPrecision(model_object.GetResolution(), 2);
    }
    label += ".pdf";
    return label;
}

} // namespace rhbm_gem::painter_internal
