#pragma once

#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/domain/AtomKeySystem.hpp>
#include <rhbm_gem/utils/domain/BondKeySystem.hpp>
#include <rhbm_gem/utils/domain/ComponentKeySystem.hpp>

namespace rhbm_gem {

class ModelAnalysisState;
template <typename T> struct KDNode;

class ModelObjectAccess
{
public:
    static ModelAnalysisState & AnalysisState(ModelObject & model_object)
    {
        return *model_object.m_analysis_state;
    }

    static const ModelAnalysisState & AnalysisState(const ModelObject & model_object)
    {
        return *model_object.m_analysis_state;
    }

    static ::KDNode<AtomObject> * KDTreeRoot(const ModelObject & model_object)
    {
        return model_object.m_kd_tree_root.get();
    }

    static ComponentKeySystem & ComponentKeySystemRef(ModelObject & model_object)
    {
        return *model_object.m_component_key_system;
    }

    static AtomKeySystem & AtomKeySystemRef(ModelObject & model_object)
    {
        return *model_object.m_atom_key_system;
    }

    static BondKeySystem & BondKeySystemRef(ModelObject & model_object)
    {
        return *model_object.m_bond_key_system;
    }

    static void ClearAnalysisFitStates(ModelObject & model_object);
};

} // namespace rhbm_gem
