#pragma once

#include <memory>

#include "data/detail/GroupPotentialEntry.hpp"
#include "data/detail/LocalPotentialEntry.hpp"
#include "data/detail/LocalPotentialFitState.hpp"
#include "data/detail/ModelAnalysisState.hpp"
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/BondObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>

namespace rhbm_gem {

class ModelAnalysisAccess
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

    static LocalPotentialEntry & EnsureLocalEntry(ModelObject & model_object, const AtomObject & atom_object)
    {
        return AnalysisState(model_object).Atoms().EnsureLocalEntry(atom_object);
    }

    static LocalPotentialEntry & EnsureLocalEntry(ModelObject & model_object, const BondObject & bond_object)
    {
        return AnalysisState(model_object).Bonds().EnsureLocalEntry(bond_object);
    }

    static void SetLocalEntry(
        ModelObject & model_object,
        const AtomObject & atom_object,
        std::unique_ptr<LocalPotentialEntry> entry)
    {
        AnalysisState(model_object).Atoms().SetLocalEntry(atom_object, std::move(entry));
    }

    static void SetLocalEntry(
        ModelObject & model_object,
        const BondObject & bond_object,
        std::unique_ptr<LocalPotentialEntry> entry)
    {
        AnalysisState(model_object).Bonds().SetLocalEntry(bond_object, std::move(entry));
    }

    static LocalPotentialEntry * FindLocalEntry(ModelObject & model_object, const AtomObject & atom_object)
    {
        return AnalysisState(model_object).Atoms().FindLocalEntry(atom_object);
    }

    static const LocalPotentialEntry * FindLocalEntry(
        const ModelObject & model_object,
        const AtomObject & atom_object)
    {
        return AnalysisState(model_object).Atoms().FindLocalEntry(atom_object);
    }

    static LocalPotentialEntry * FindLocalEntry(ModelObject & model_object, const BondObject & bond_object)
    {
        return AnalysisState(model_object).Bonds().FindLocalEntry(bond_object);
    }

    static const LocalPotentialEntry * FindLocalEntry(
        const ModelObject & model_object,
        const BondObject & bond_object)
    {
        return AnalysisState(model_object).Bonds().FindLocalEntry(bond_object);
    }

    static LocalPotentialEntry * FindLocalEntry(const AtomObject & atom_object)
    {
        auto * owner{ atom_object.GetOwnerModel() };
        return owner == nullptr ? nullptr : FindLocalEntry(*owner, atom_object);
    }

    static const LocalPotentialEntry * FindLocalEntryConst(const AtomObject & atom_object)
    {
        const auto * owner{ atom_object.GetOwnerModel() };
        return owner == nullptr ? nullptr : FindLocalEntry(*owner, atom_object);
    }

    static const LocalPotentialEntry * FindLocalEntry(const BondObject & bond_object)
    {
        auto * owner{ bond_object.GetOwnerModel() };
        return owner == nullptr ? nullptr : FindLocalEntry(*owner, bond_object);
    }

    static const LocalPotentialEntry * FindLocalEntryConst(const BondObject & bond_object)
    {
        const auto * owner{ bond_object.GetOwnerModel() };
        return owner == nullptr ? nullptr : FindLocalEntry(*owner, bond_object);
    }

    static LocalPotentialFitState & EnsureFitState(ModelObject & model_object, const AtomObject & atom_object)
    {
        return AnalysisState(model_object).Atoms().EnsureFitState(atom_object);
    }

    static LocalPotentialFitState & EnsureFitState(ModelObject & model_object, const BondObject & bond_object)
    {
        return AnalysisState(model_object).Bonds().EnsureFitState(bond_object);
    }

    static void ClearAnalysisFitStates(ModelObject & model_object)
    {
        model_object.m_analysis_state->ClearFitStates();
    }
};

} // namespace rhbm_gem
