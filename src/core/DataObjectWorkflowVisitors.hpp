#pragma once

#include <array>
#include <tuple>
#include <unordered_map>
#include <vector>

#include "DataObjectVisitor.hpp"
#include "MapInterpolationVisitor.hpp"
#include "OptionEnumClass.hpp"

class SamplerBase;
class AtomSelector;

namespace rhbm_gem {

class MapObject;
class ModelObject;
class AtomObject;
class BondObject;

class MapNormalizeVisitor : public DataObjectVisitor
{
public:
    void VisitAtomObject(AtomObject & data_object) override;
    void VisitBondObject(BondObject & data_object) override;
    void VisitModelObject(ModelObject & data_object) override;
    void VisitMapObject(MapObject & data_object) override;
};

struct ModelPreparationOptions
{
    bool select_all_atoms{ false };
    bool select_all_bonds{ false };
    bool apply_atom_symmetry_filter{ false };
    bool apply_bond_symmetry_filter{ false };
    bool asymmetry_flag{ false };
    bool update_model{ false };
    bool initialize_atom_local_entries{ false };
    bool initialize_bond_local_entries{ false };
};

class ModelPreparationVisitor : public DataObjectVisitor
{
    ModelPreparationOptions m_options;

public:
    explicit ModelPreparationVisitor(ModelPreparationOptions options);

    void VisitAtomObject(AtomObject & data_object) override;
    void VisitBondObject(BondObject & data_object) override;
    void VisitModelObject(ModelObject & data_object) override;
    void VisitMapObject(MapObject & data_object) override;
};

class ModelSelectionVisitor : public DataObjectVisitor
{
    ::AtomSelector * m_selector;

public:
    explicit ModelSelectionVisitor(::AtomSelector * selector);

    void VisitAtomObject(AtomObject & data_object) override;
    void VisitBondObject(BondObject & data_object) override;
    void VisitModelObject(ModelObject & data_object) override;
    void VisitMapObject(MapObject & data_object) override;
};

struct ModelAtomCollectorOptions
{
    bool selected_only{ true };
    bool require_local_potential_entry{ false };
};

class ModelSelectedAtomCollectorVisitor : public DataObjectVisitor
{
    ModelAtomCollectorOptions m_options;
    std::vector<AtomObject *> m_atom_list;

public:
    explicit ModelSelectedAtomCollectorVisitor(ModelAtomCollectorOptions options = {});

    void VisitAtomObject(AtomObject & data_object) override;
    void VisitBondObject(BondObject & data_object) override;
    void VisitModelObject(ModelObject & data_object) override;
    void VisitMapObject(MapObject & data_object) override;

    const std::vector<AtomObject *> & GetAtomList() const { return m_atom_list; }
};

struct SimulationAtomPreparationOptions
{
    PartialCharge partial_charge_choice{ PartialCharge::PARTIAL };
    bool include_unknown_atoms{ true };
};

class SimulationAtomPreparationVisitor : public DataObjectVisitor
{
    SimulationAtomPreparationOptions m_options;
    std::vector<AtomObject *> m_atom_list;
    std::unordered_map<int, double> m_atom_charge_map;
    std::array<float, 3> m_range_minimum;
    std::array<float, 3> m_range_maximum;
    bool m_has_atom{ false };

public:
    explicit SimulationAtomPreparationVisitor(
        SimulationAtomPreparationOptions options = {});

    void VisitAtomObject(AtomObject & data_object) override;
    void VisitBondObject(BondObject & data_object) override;
    void VisitModelObject(ModelObject & data_object) override;
    void VisitMapObject(MapObject & data_object) override;

    const std::vector<AtomObject *> & GetAtomList() const { return m_atom_list; }
    const std::unordered_map<int, double> & GetAtomChargeMap() const { return m_atom_charge_map; }
    bool HasAtom() const { return m_has_atom; }
    std::array<float, 3> GetRangeMinimum() const { return m_range_minimum; }
    std::array<float, 3> GetRangeMaximum() const { return m_range_maximum; }
};

class ModelAtomBondContextVisitor : public DataObjectVisitor
{
    std::unordered_map<int, AtomObject *> m_atom_map;
    std::unordered_map<int, std::vector<BondObject *>> m_bond_map;

public:
    void VisitAtomObject(AtomObject & data_object) override;
    void VisitBondObject(BondObject & data_object) override;
    void VisitModelObject(ModelObject & data_object) override;
    void VisitMapObject(MapObject & data_object) override;

    const std::unordered_map<int, AtomObject *> & GetAtomMap() const { return m_atom_map; }
    const std::unordered_map<int, std::vector<BondObject *>> & GetBondMap() const { return m_bond_map; }
};

class MapSamplingWorkflow
{
    MapInterpolationVisitor m_interpolation_visitor;

public:
    explicit MapSamplingWorkflow(::SamplerBase * sampler);

    std::vector<std::tuple<float, float>> Sample(
        const MapObject & map_object,
        const std::array<float, 3> & position,
        const std::array<float, 3> & axis_vector={ 0.0f, 0.0f, 0.0f });
};

void NormalizeMapObject(MapObject & map_object);
void PrepareModelObject(ModelObject & model_object, const ModelPreparationOptions & options);
void ApplyModelSelection(ModelObject & model_object, ::AtomSelector & selector);

} // namespace rhbm_gem
