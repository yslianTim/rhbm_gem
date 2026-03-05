#pragma once

#include <array>
#include <tuple>
#include <vector>

#include "DataObjectVisitor.hpp"
#include "MapInterpolationVisitor.hpp"

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
