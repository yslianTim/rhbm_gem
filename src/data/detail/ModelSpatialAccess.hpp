#pragma once

#include <vector>

namespace rhbm_gem {

class AtomObject;
class ModelObject;

class ModelSpatialAccess
{
public:
    static std::vector<AtomObject *> FindAtomsInRange(
        ModelObject & model_object,
        const AtomObject & center_atom,
        double range);
};

} // namespace rhbm_gem
