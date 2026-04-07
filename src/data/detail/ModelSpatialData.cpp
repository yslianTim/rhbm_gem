#include "data/detail/ModelSpatialData.hpp"

#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>

namespace rhbm_gem {

ModelSpatialData::ModelSpatialData() = default;

ModelSpatialData::~ModelSpatialData() = default;

ModelSpatialData & ModelSpatialData::Of(ModelObject & model_object)
{
    if (model_object.m_spatial_data == nullptr)
    {
        model_object.m_spatial_data = std::make_unique<ModelSpatialData>();
    }
    return *model_object.m_spatial_data;
}

const ModelSpatialData & ModelSpatialData::Of(const ModelObject & model_object)
{
    return Of(const_cast<ModelObject &>(model_object));
}

std::vector<AtomObject *> ModelSpatialData::FindAtomsInRange(
    ModelObject & model_object,
    const AtomObject & center_atom,
    double range)
{
    auto & spatial_data{ Of(model_object) };
    spatial_data.EnsureKDTreeRoot(model_object);
    if (spatial_data.m_kd_tree_root == nullptr)
    {
        return {};
    }

    return KDTreeAlgorithm<AtomObject>::RangeSearch(
        spatial_data.m_kd_tree_root.get(),
        const_cast<AtomObject *>(&center_atom),
        range);
}

void ModelSpatialData::Clear()
{
    m_kd_tree_root.reset();
}

void ModelSpatialData::EnsureKDTreeRoot(ModelObject & model_object)
{
    if (m_kd_tree_root != nullptr)
    {
        return;
    }

    std::vector<AtomObject *> atom_ptr_list;
    atom_ptr_list.reserve(model_object.m_atom_list.size());
    for (auto & atom : model_object.m_atom_list)
    {
        atom_ptr_list.emplace_back(atom.get());
    }
    m_kd_tree_root = KDTreeAlgorithm<AtomObject>::BuildKDTree(atom_ptr_list, 0);
}

} // namespace rhbm_gem
