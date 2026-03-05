#include "DataObjectWorkflowVisitors.hpp"

#include <memory>
#include <stdexcept>

#include "AtomObject.hpp"
#include "BondObject.hpp"
#include "LocalPotentialEntry.hpp"
#include "MapObject.hpp"
#include "ModelObject.hpp"

namespace rhbm_gem {

namespace {

[[noreturn]] void ThrowUnsupportedType(const char * visitor_name, const char * supported_type)
{
    throw std::logic_error(
        std::string(visitor_name) + " supports " + supported_type + " only.");
}

} // namespace

void MapNormalizeVisitor::VisitAtomObject(AtomObject & data_object)
{
    (void)data_object;
    ThrowUnsupportedType("MapNormalizeVisitor", "MapObject");
}

void MapNormalizeVisitor::VisitBondObject(BondObject & data_object)
{
    (void)data_object;
    ThrowUnsupportedType("MapNormalizeVisitor", "MapObject");
}

void MapNormalizeVisitor::VisitModelObject(ModelObject & data_object)
{
    (void)data_object;
    ThrowUnsupportedType("MapNormalizeVisitor", "MapObject");
}

void MapNormalizeVisitor::VisitMapObject(MapObject & data_object)
{
    data_object.MapValueArrayNormalization();
}

ModelPreparationVisitor::ModelPreparationVisitor(ModelPreparationOptions options) :
    m_options{ options }
{
}

void ModelPreparationVisitor::VisitAtomObject(AtomObject & data_object)
{
    (void)data_object;
    ThrowUnsupportedType("ModelPreparationVisitor", "ModelObject");
}

void ModelPreparationVisitor::VisitBondObject(BondObject & data_object)
{
    (void)data_object;
    ThrowUnsupportedType("ModelPreparationVisitor", "ModelObject");
}

void ModelPreparationVisitor::VisitModelObject(ModelObject & data_object)
{
    if (m_options.select_all_atoms)
    {
        for (auto & atom : data_object.GetAtomList())
        {
            atom->SetSelectedFlag(true);
        }
    }
    if (m_options.select_all_bonds)
    {
        for (auto & bond : data_object.GetBondList())
        {
            bond->SetSelectedFlag(true);
        }
    }

    if (m_options.apply_atom_symmetry_filter)
    {
        data_object.FilterAtomFromSymmetry(m_options.asymmetry_flag);
    }
    if (m_options.apply_bond_symmetry_filter)
    {
        data_object.FilterBondFromSymmetry(m_options.asymmetry_flag);
    }

    bool has_updated{ false };
    if (m_options.update_model)
    {
        data_object.Update();
        has_updated = true;
    }

    if ((m_options.initialize_atom_local_entries || m_options.initialize_bond_local_entries)
        && !has_updated)
    {
        data_object.Update();
    }

    if (m_options.initialize_atom_local_entries)
    {
        for (auto * atom : data_object.GetSelectedAtomList())
        {
            atom->AddLocalPotentialEntry(std::make_unique<LocalPotentialEntry>());
        }
    }
    if (m_options.initialize_bond_local_entries)
    {
        for (auto * bond : data_object.GetSelectedBondList())
        {
            bond->AddLocalPotentialEntry(std::make_unique<LocalPotentialEntry>());
        }
    }
}

void ModelPreparationVisitor::VisitMapObject(MapObject & data_object)
{
    (void)data_object;
    ThrowUnsupportedType("ModelPreparationVisitor", "ModelObject");
}

MapSamplingWorkflow::MapSamplingWorkflow(::SamplerBase * sampler) :
    m_interpolation_visitor{ sampler }
{
}

std::vector<std::tuple<float, float>> MapSamplingWorkflow::Sample(
    const MapObject & map_object,
    const std::array<float, 3> & position,
    const std::array<float, 3> & axis_vector)
{
    m_interpolation_visitor.SetPosition(position);
    m_interpolation_visitor.SetAxisVector(axis_vector);
    map_object.Accept(m_interpolation_visitor);
    return m_interpolation_visitor.ConsumeSamplingDataList();
}

} // namespace rhbm_gem
