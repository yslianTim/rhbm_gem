#pragma once

#include <vector>

#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>

namespace rhbm_gem {

std::vector<AtomObject *> FindAtomsInRange(
    ModelObject & model_object,
    const AtomObject & center_atom,
    double range);

} // namespace rhbm_gem
