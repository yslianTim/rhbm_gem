#pragma once
#include <string>
#include <vector>
#include <rhbm_gem/utils/domain/GlobalEnumClass.hpp>
namespace rhbm_gem { class ModelObject;
namespace stage_result_io {
std::string Encode(const ModelObject &);
void Decode(ModelObject &, const std::string &, int expected_version, const std::vector<GroupKey> & legacy_groups);
void AdaptLegacy(ModelObject &);
void MapSnapshot(ModelObject &);
void Validate(const ModelObject &);
}}
