#pragma once

#include <cstddef>
#include <vector>

#include <rhbm_gem/data/object/ModelObject.hpp>

namespace rhbm_gem {

class AtomObject;
class BondObject;

class ModelSelectionView
{
public:
    static const std::vector<AtomObject *> & SelectedAtoms(const ModelObject & model_object)
    {
        return model_object.m_selected_atom_list;
    }

    static const std::vector<BondObject *> & SelectedBonds(const ModelObject & model_object)
    {
        return model_object.m_selected_bond_list;
    }

    static size_t SelectedAtomCount(const ModelObject & model_object)
    {
        return model_object.m_selected_atom_list.size();
    }

    static size_t SelectedBondCount(const ModelObject & model_object)
    {
        return model_object.m_selected_bond_list.size();
    }
};

} // namespace rhbm_gem
