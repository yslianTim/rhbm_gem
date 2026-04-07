#include "data/detail/ModelDerivedState.hpp"

#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/math/ArrayStats.hpp>

namespace rhbm_gem {

ModelDerivedState::ModelDerivedState() = default;

ModelDerivedState::~ModelDerivedState() = default;

ModelDerivedState & ModelDerivedState::Of(ModelObject & model_object)
{
    if (model_object.m_derived_state == nullptr)
    {
        model_object.m_derived_state = std::make_unique<ModelDerivedState>();
    }
    return *model_object.m_derived_state;
}

const ModelDerivedState & ModelDerivedState::Of(const ModelObject & model_object)
{
    return Of(const_cast<ModelObject &>(model_object));
}

std::array<float, 3> ModelDerivedState::GetCenterOfMassPosition(ModelObject & model_object)
{
    if (m_center_of_mass_position.has_value())
    {
        return *m_center_of_mass_position;
    }

    std::vector<float> pos_x, pos_y, pos_z;
    pos_x.reserve(model_object.m_atom_list.size());
    pos_y.reserve(model_object.m_atom_list.size());
    pos_z.reserve(model_object.m_atom_list.size());
    for (auto & atom : model_object.m_atom_list)
    {
        const auto & pos{ atom->GetPositionRef() };
        pos_x.emplace_back(pos.at(0));
        pos_y.emplace_back(pos.at(1));
        pos_z.emplace_back(pos.at(2));
    }

    m_center_of_mass_position = std::array<float, 3>{
        ArrayStats<float>::ComputeMean(pos_x.data(), pos_x.size()),
        ArrayStats<float>::ComputeMean(pos_y.data(), pos_y.size()),
        ArrayStats<float>::ComputeMean(pos_z.data(), pos_z.size())
    };
    return *m_center_of_mass_position;
}

std::tuple<double, double> ModelDerivedState::GetModelPositionRange(
    ModelObject & model_object,
    int axis)
{
    if (m_model_position_range[static_cast<size_t>(axis)].has_value())
    {
        return *m_model_position_range[static_cast<size_t>(axis)];
    }

    std::vector<double> position_list;
    position_list.reserve(model_object.m_atom_list.size());
    for (auto & atom : model_object.m_atom_list)
    {
        position_list.emplace_back(atom->GetPosition().at(static_cast<size_t>(axis)));
    }
    m_model_position_range[static_cast<size_t>(axis)] =
        ArrayStats<double>::ComputeScalingRangeTuple(position_list, 0.0);
    return *m_model_position_range[static_cast<size_t>(axis)];
}

std::vector<AtomObject *> ModelDerivedState::FindAtomsInRange(
    ModelObject & model_object,
    const AtomObject & center_atom,
    double range)
{
    EnsureKDTreeRoot(model_object);
    if (m_kd_tree_root == nullptr)
    {
        return {};
    }

    return KDTreeAlgorithm<AtomObject>::RangeSearch(
        m_kd_tree_root.get(),
        const_cast<AtomObject *>(&center_atom),
        range);
}

void ModelDerivedState::Clear()
{
    m_kd_tree_root.reset();
    m_center_of_mass_position.reset();
    for (auto & axis_range : m_model_position_range)
    {
        axis_range.reset();
    }
}

void ModelDerivedState::EnsureKDTreeRoot(ModelObject & model_object)
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
