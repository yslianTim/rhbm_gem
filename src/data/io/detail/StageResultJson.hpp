#pragma once
#include <string>
namespace rhbm_gem { class ModelObject;
namespace stage_result_io {
std::string Encode(const ModelObject &);
void Decode(ModelObject &, const std::string &);
void AdaptLegacy(ModelObject &);
void MapSnapshot(ModelObject &);
void Validate(const ModelObject &);
}}
