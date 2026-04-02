#pragma once

#include <memory>
#include <string>
#include <vector>

namespace rhbm_gem {

class AtomObject;
class BondObject;
class GroupPotentialEntry;
class LocalPotentialEntry;
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

    static LocalPotentialEntry & EnsureLocalEntry(
        ModelObject & model_object,
        const AtomObject & atom_object);
    static LocalPotentialEntry & EnsureLocalEntry(
        ModelObject & model_object,
        const BondObject & bond_object);
    static void SetLocalEntry(
        ModelObject & model_object,
        const AtomObject & atom_object,
        std::unique_ptr<LocalPotentialEntry> entry);
    static void SetLocalEntry(
        ModelObject & model_object,
        const BondObject & bond_object,
        std::unique_ptr<LocalPotentialEntry> entry);
    static LocalPotentialEntry * FindLocalEntry(
        ModelObject & model_object,
        const AtomObject & atom_object);
    static const LocalPotentialEntry * FindLocalEntry(
        const ModelObject & model_object,
        const AtomObject & atom_object);
    static LocalPotentialEntry * FindLocalEntry(
        ModelObject & model_object,
        const BondObject & bond_object);
    static const LocalPotentialEntry * FindLocalEntry(
        const ModelObject & model_object,
        const BondObject & bond_object);
    static const LocalPotentialEntry * FindLocalEntry(const AtomObject & atom_object);
    static const LocalPotentialEntry * FindLocalEntry(const BondObject & bond_object);
    static const LocalPotentialEntry & RequireLocalEntry(
        const LocalPotentialEntry * entry,
        const char * context);
    static const LocalPotentialEntry & RequireLocalEntry(const AtomObject & atom_object);
    static const LocalPotentialEntry & RequireLocalEntry(const BondObject & bond_object);

    static GroupPotentialEntry & EnsureAtomGroupEntry(
        ModelObject & model_object,
        const std::string & class_key);
    static GroupPotentialEntry & EnsureBondGroupEntry(
        ModelObject & model_object,
        const std::string & class_key);
    static GroupPotentialEntry * FindAtomGroupEntry(
        ModelObject & model_object,
        const std::string & class_key);
    static const GroupPotentialEntry * FindAtomGroupEntry(
        const ModelObject & model_object,
        const std::string & class_key);
    static GroupPotentialEntry * FindBondGroupEntry(
        ModelObject & model_object,
        const std::string & class_key);
    static const GroupPotentialEntry * FindBondGroupEntry(
        const ModelObject & model_object,
        const std::string & class_key);

    static std::vector<AtomObject *> FindAtomsInRange(
        ModelObject & model_object,
        const AtomObject & center_atom,
        double range);
};

} // namespace rhbm_gem
