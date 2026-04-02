#include "data/detail/ModelAnalysisAccess.hpp"

#include "data/detail/ModelAnalysisData.hpp"
#include "data/detail/ModelSpatialCache.hpp"
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/BondObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/math/KDTreeAlgorithm.hpp>

namespace rhbm_gem {

ModelAnalysisData & ModelAnalysisAccess::Mutable(ModelObject & model_object)
{
    return *model_object.m_analysis_data;
}

const ModelAnalysisData & ModelAnalysisAccess::Read(const ModelObject & model_object)
{
    return *model_object.m_analysis_data;
}

void ModelAnalysisAccess::ClearFitStates(ModelObject & model_object)
{
    Mutable(model_object).ClearFitStates();
}

ModelObject * ModelAnalysisAccess::OwnerOf(const AtomObject & atom_object)
{
    return atom_object.m_owner_model;
}

ModelObject * ModelAnalysisAccess::OwnerOf(const BondObject & bond_object)
{
    return bond_object.m_owner_model;
}

std::vector<AtomObject *> ModelAnalysisAccess::FindAtomsInRange(
    ModelObject & model_object,
    const AtomObject & center_atom,
    double range)
{
    model_object.EnsureKDTreeRoot();
    if (model_object.m_spatial_cache == nullptr ||
        model_object.m_spatial_cache->kd_tree_root == nullptr)
    {
        return {};
    }

    return KDTreeAlgorithm<AtomObject>::RangeSearch(
        model_object.m_spatial_cache->kd_tree_root.get(),
        const_cast<AtomObject *>(&center_atom),
        range);
}

} // namespace rhbm_gem
