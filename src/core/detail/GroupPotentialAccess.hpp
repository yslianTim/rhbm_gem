#pragma once

#include <string>

#include "data/detail/ModelAnalysisAccess.hpp"
#include "data/detail/ModelAnalysisData.hpp"
#include "data/detail/GroupPotentialEntry.hpp"
#include <rhbm_gem/data/object/ModelObject.hpp>

namespace rhbm_gem {

inline GroupPotentialEntry & EnsureAtomGroupPotentialEntry(
    ModelObject & model_object,
    const std::string & class_key)
{
    return ModelAnalysisAccess::Mutable(model_object).Atoms().EnsureGroupEntry(class_key);
}

inline GroupPotentialEntry & EnsureBondGroupPotentialEntry(
    ModelObject & model_object,
    const std::string & class_key)
{
    return ModelAnalysisAccess::Mutable(model_object).Bonds().EnsureGroupEntry(class_key);
}

inline GroupPotentialEntry * FindAtomGroupPotentialEntry(
    ModelObject & model_object,
    const std::string & class_key)
{
    return ModelAnalysisAccess::Mutable(model_object).Atoms().FindGroupEntry(class_key);
}

inline const GroupPotentialEntry * FindAtomGroupPotentialEntry(
    const ModelObject & model_object,
    const std::string & class_key)
{
    return ModelAnalysisAccess::Read(model_object).Atoms().FindGroupEntry(class_key);
}

inline GroupPotentialEntry * FindBondGroupPotentialEntry(
    ModelObject & model_object,
    const std::string & class_key)
{
    return ModelAnalysisAccess::Mutable(model_object).Bonds().FindGroupEntry(class_key);
}

inline const GroupPotentialEntry * FindBondGroupPotentialEntry(
    const ModelObject & model_object,
    const std::string & class_key)
{
    return ModelAnalysisAccess::Read(model_object).Bonds().FindGroupEntry(class_key);
}

} // namespace rhbm_gem
