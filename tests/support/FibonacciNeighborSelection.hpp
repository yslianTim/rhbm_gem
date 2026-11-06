#pragma once
// Fixed-geometry research controls; deliberately outside the Production policy API.
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/domain/Constants.hpp>
#include <rhbm_gem/utils/domain/SampleFilter.hpp>
#include <rhbm_gem/utils/math/ArrayHelper.hpp>
#include <rhbm_gem/utils/math/NumericValidation.hpp>
#include <map>
#include <algorithm>

namespace second_stage_test {
struct NeighborSelectionRow
{
    bool f0{}, ownership{}, cone{}, qscore{}, shared_qscore{};
};
struct NeighborSelectionShell
{
    std::vector<int> fixed_ids, qscore_ids;
};
struct NeighborSelectionMasks
{
    std::vector<NeighborSelectionRow> rows;
    std::map<double, NeighborSelectionShell> shells;
};
inline NeighborSelectionMasks FibonacciNeighborMasks(
    const rhbm_gem::AtomObject & atom, const rhbm_gem::ModelObject & model,
    const LocalPotentialSampleList & raw)
{
    using namespace rhbm_gem;
    const auto & center=atom.GetPositionRef();
    numeric_validation::RequireAllFinite(center,"target position");
    auto fixed=atom.FindNeighborAtoms();
    fixed.erase(std::remove_if(fixed.begin(),fixed.end(),[](const auto * a) {
        return a->GetElement()==Element::HYDROGEN;
    }),fixed.end());
    NeighborSelectionMasks masks;
    SamplingPointList f0_points;
    std::vector<std::array<double,3>> fixed_positions;
    for (const auto & sample:raw) f0_points.push_back(sample.point);
    for (const auto * a:fixed) fixed_positions.push_back(a->GetPosition());
    sample_filter::FilterSamplingPointList(f0_points,center,fixed_positions);
    std::map<double,std::vector<AtomObject *>> shell_neighbors;
    for (const auto & sample:raw)
    {
        const auto & point=sample.point;
        const double r=point.distance, exclusion=.9*r;
        numeric_validation::RequireFiniteNonNegative(r,"sample radius");
        numeric_validation::RequireAllFinite(point.position,"sample position");
        if (!shell_neighbors.contains(r))
        {
            auto candidates=model.FindNeighborAtoms(atom,2.*r);
            candidates.erase(std::remove_if(candidates.begin(),candidates.end(),[](const auto * a) {
                return a->GetElement()==Element::HYDROGEN;
            }),candidates.end());
            auto & shell=masks.shells[r];
            for (const auto * a:fixed) shell.fixed_ids.push_back(a->GetSerialID());
            for (const auto * a:candidates) shell.qscore_ids.push_back(a->GetSerialID());
            shell_neighbors[r]=std::move(candidates);
        }
        NeighborSelectionRow row; row.f0=f0_points[masks.rows.size()].is_selected;
        const double distance=array_helper::ComputeNorm(point.position,center);
        const auto direction=array_helper::ComputeVector(center,point.position,true);
        for (const auto * neighbor:fixed)
        {
            const auto & position=neighbor->GetPositionRef();
            // Match the F0 coincident-position contract for ownership/cone only.
            if (numeric_validation::IsNonEqual(center,position))
            {
                row.ownership |= distance>array_helper::ComputeNorm(point.position,position);
                row.cone |= array_helper::ComputeDotProduct(direction,
                    array_helper::ComputeVector(center,position,true))>std::cos(30.*Constants::pi/180.);
            }
            row.shared_qscore |= sample_filter::IsInsideNeighborExclusion(point.position,position,exclusion*exclusion);
        }
        for (const auto * neighbor:shell_neighbors.at(r))
            row.qscore |= sample_filter::IsInsideNeighborExclusion(
                point.position,neighbor->GetPositionRef(),exclusion*exclusion);
        masks.rows.push_back(row);
    }
    return masks;
}
inline void ApplyFibonacciNeighborSelection(
    LocalPotentialSampleList & raw,const NeighborSelectionMasks & masks,bool qscore)
{
    if (raw.size()!=masks.rows.size()) throw std::invalid_argument("Selection mask size mismatch.");
    for (std::size_t i=0;i<raw.size();++i)
        raw[i].point.is_selected=!(qscore ? masks.rows[i].qscore : masks.rows[i].ownership);
}
} // namespace second_stage_test
