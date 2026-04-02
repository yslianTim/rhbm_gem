#pragma once

#include <set>
#include <unordered_set>
#include <utility>

#include <rhbm_gem/data/object/ModelObject.hpp>

namespace rhbm_gem {

class ModelSelectionAccess
{
public:
    static void ApplyLoadedAtomSelection(
        ModelObject & model_object,
        const std::unordered_set<int> & selected_serial_ids)
    {
        model_object.SetAtomSelectionBulk(selected_serial_ids);
    }

    static void ApplyLoadedBondSelection(
        ModelObject & model_object,
        const std::set<std::pair<int, int>> & selected_serial_pairs)
    {
        model_object.SetBondSelectionBulk(selected_serial_pairs);
    }
};

} // namespace rhbm_gem
