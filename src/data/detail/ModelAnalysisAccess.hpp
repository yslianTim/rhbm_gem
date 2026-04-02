#pragma once

#include <vector>

namespace rhbm_gem {

class AtomObject;
class BondObject;
class ModelAnalysisData;
class ModelObject;

class ModelAnalysisAccess
{
public:
    static ModelAnalysisData & Mutable(ModelObject & model_object);
    static const ModelAnalysisData & Read(const ModelObject & model_object);
    static void ClearFitStates(ModelObject & model_object);
    static ModelObject * OwnerOf(const AtomObject & atom_object);
    static ModelObject * OwnerOf(const BondObject & bond_object);
    static std::vector<AtomObject *> FindAtomsInRange(
        ModelObject & model_object,
        const AtomObject & center_atom,
        double range);
};

} // namespace rhbm_gem
