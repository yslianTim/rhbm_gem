#pragma once

#include <string>
#include <unordered_map>
#include <vector>

namespace rhbm_gem {

class ModelObject;

namespace core {

using ModelObjectList = std::vector<ModelObject *>;
using ReferenceModelGroupMap = std::unordered_map<std::string, std::vector<ModelObject *>>;

void PaintAtom(const ModelObjectList & model_objects, const std::string & output_folder);
void PaintGaus(const ModelObjectList & model_objects, const std::string & output_folder);
void PaintComparison(
    const ModelObjectList & model_objects,
    const ReferenceModelGroupMap & reference_model_groups,
    const std::string & output_folder);
void PaintDemo(
    const ModelObjectList & model_objects,
    const ReferenceModelGroupMap & reference_model_groups,
    const std::string & output_folder);

} // namespace core
} // namespace rhbm_gem
