#include "PostJointQScore.hpp"
#include "data/detail/JointStageAdapter.hpp"
#include <rhbm_gem/core/QScoreHelper.hpp>
#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <map>
#include <set>
#include <stdexcept>

namespace rhbm_gem::core::detail {
namespace {
double SquaredDistance(const std::array<double,3> & a, const std::array<double,3> & b)
{
    double square{};
    for (std::size_t axis=0;axis<3;++axis) { const double delta=a[axis]-b[axis]; square+=delta*delta; }
    return square;
}
bool SameEndpoint(const LocalStageEstimate & a, const LocalStageEstimate & b)
{
    return a.source.method==b.source.method && a.source.role==b.source.role &&
        a.source.atom_id==b.source.atom_id && a.source.component_id==b.source.component_id &&
        a.source.run_id==b.source.run_id && !a.source.fixed_from_initial_stage &&
        a.convergence==JointCheckStatus::Passed && b.convergence==JointCheckStatus::Passed &&
        a.point && b.point && a.point->ToVector()==b.point->ToVector();
}
struct ContributorModel
{
    std::optional<GaussianModel3D> point;
    bool halo{};
    std::vector<std::string> extra_donors;
};
}

std::vector<std::pair<int, PostJointQScoreResult>> BuildPostJointQScores(
    const MapObject & map, const ModelObject & model, const JointProblem & problem,
    const JointAnalysisResult & snapshot, const HaloMedianBuildResult * halo_build,
    JointPostprocessingMetrics & metrics)
{
    const auto start=std::chrono::steady_clock::now();
    if (!IsQualifiedJointAnalysisResult(snapshot) || snapshot.atom_ids!=problem.OriginalInput().atom_ids ||
        !snapshot.selection_domain || snapshot.selection_domain->target_indices.empty())
        throw std::invalid_argument("QScore requires a qualified matching Joint snapshot.");
    const auto & original=problem.OriginalInput();
    if (!original.selection_domain ||
        snapshot.selection_domain->target_indices!=original.selection_domain->target_indices ||
        snapshot.selection_domain->contract!=original.selection_domain->contract ||
        snapshot.selection_domain->contributor_policy!=original.selection_domain->contributor_policy ||
        snapshot.selection_domain->observation_radius!=original.selection_domain->observation_radius ||
        snapshot.selection_domain->support_radius!=original.selection_domain->support_radius ||
        snapshot.row_ids!=original.row_ids || snapshot.metadata.grid_size!=map.GetGridSize() ||
        snapshot.metadata.grid_spacing!=map.GetGridSpacing() || snapshot.metadata.origin!=map.GetOrigin())
        throw std::invalid_argument("QScore map/Original-domain identity mismatch.");
    for (std::size_t row=0;row<original.row_ids.size();++row)
    {
        const auto grid=std::stoull(original.row_ids[row]);
        if (grid>=map.GetMapValueArraySize() || !std::isfinite(map.GetMapValue(grid)) ||
            map.GetMapValue(grid)!=original.observations.at(row))
            throw std::invalid_argument("QScore Original observations differ from map.");
    }
    const auto & target_indices=snapshot.selection_domain->target_indices;
    const auto * first=model.FindAtomPtr(std::stoi(snapshot.atom_ids.at(target_indices.front())));
    const auto run=AtomLocalPotentialView::For(*first).GetJointStageEstimate().source.run_id;
    if (run.empty()) throw std::invalid_argument("QScore requires Joint run identity.");
    const auto expected=data_internal::BuildJointStageEstimates(snapshot,run);
    const auto rebuilt_donors=halo_build ? std::vector<HaloMedianTargetDonor>{} :
        QualifiedJointTargetDonors(model,problem,snapshot);
    const auto & donors=halo_build ? halo_build->target_donors : rebuilt_donors;
    std::set<std::size_t> qualified_donors;
    if (snapshot.target_evidence)
        for (const auto & point:snapshot.target_evidence->point_atoms)
            if (point.component_convergence==JointCheckStatus::Passed &&
                point.parameter_validity==JointCheckStatus::Passed && point.local_conditioning=="Adequate")
                qualified_donors.insert(point.atom);
    bool donors_valid=donors.size()==qualified_donors.size();
    std::set<std::size_t> seen_donors;
    for (const auto & donor:donors)
    {
        const auto found=expected.find(std::stoi(donor.atom_id));
        const auto * atom=model.FindAtomPtr(std::stoi(donor.atom_id));
        const auto local=AtomLocalPotentialView::For(*atom);
        if (!qualified_donors.contains(donor.contributor_index) ||
            !seen_donors.insert(donor.contributor_index).second ||
            snapshot.atom_ids.at(donor.contributor_index)!=donor.atom_id || found==expected.end() ||
            !found->second.point || donor.position!=atom->GetPositionRef() ||
            donor.point!=std::array{found->second.point->GetAmplitude(),found->second.point->GetWidth(),
                found->second.point->GetOffset()} || !local.IsAvailable() ||
            !SameEndpoint(local.GetJointStageEstimate(),found->second)) donors_valid=false;
    }
    // One donor index for all requested Halos, including extra-domain coverage.
    const GeometryLocalTargetDonorIndex donor_index(donors);
    std::map<int,const JointHaloMedianReconstruction *> saved_halos;
    if (snapshot.halo_median_reconstruction)
    {
        const auto & saved=*snapshot.halo_median_reconstruction;
        if (saved.joint_run_id!=run || saved.contract!=HaloMedianReconstructionContractVersion ||
            saved.donor_policy!=GeometryLocalTargetMedianDonorPolicy ||
            saved.neighbor_count!=GeometryLocalTargetMedianNeighborCount || saved.valid_target_count!=donors.size())
            throw std::invalid_argument("QScore Halo provenance mismatch.");
        for (const auto & halo:saved.halos) saved_halos.emplace(std::stoi(halo.atom_id),&halo);
    }
    std::map<int,ContributorModel> cache;
    const auto contributor_model=[&](const AtomObject & atom) -> const ContributorModel & {
        const int id=atom.GetSerialID();
        if (const auto found=cache.find(id);found!=cache.end()) return found->second;
        auto & result=cache[id];
        const auto endpoint=expected.find(id);
        if (endpoint!=expected.end() && endpoint->second.source.role==FittingRole::Target)
        {
            const auto view=AtomLocalPotentialView::For(atom);
            if (view.IsAvailable() && SameEndpoint(view.GetJointStageEstimate(),endpoint->second))
                result.point=endpoint->second.point;
            return result;
        }
        result.halo=true;
        if (!donors_valid) return result;
        const auto nearest=donor_index.Query(atom.GetPositionRef());
        if (nearest.empty()) return result;
        const auto median=GeometryLocalTargetMedian(nearest);
        std::vector<std::string> ids;
        for (const auto * donor:nearest) ids.push_back(donor->atom_id);
        if (const auto saved=saved_halos.find(id);saved!=saved_halos.end())
        {
            const auto & halo=*saved->second;
            if (halo.joint_run_id!=run || halo.contract!=HaloMedianReconstructionContractVersion ||
                halo.source!=HaloMedianSource::GeometryLocalTargetMedian || halo.donor_atom_ids!=ids ||
                halo.donor_count!=ids.size() ||
                std::array{halo.amplitude,halo.width,halo.charge}!=median)
                throw std::invalid_argument("QScore saved Halo representative mismatch.");
            result.point=GaussianModel3D{halo.amplitude,halo.width,halo.charge};
        }
        else if (endpoint==expected.end())
        {
            result.point=GaussianModel3D{median[0],median[1],median[2]};
            result.extra_donors=std::move(ids);
        }
        // A recorded Halo without its saved reconstruction is unavailable.
        return result;
    };
    std::vector<std::pair<int,PostJointQScoreResult>> results;
    results.reserve(target_indices.size());
    for (const auto index:target_indices)
    {
        const int id=std::stoi(snapshot.atom_ids.at(index));
        const auto & atom=*model.FindAtomPtr(id);
        const auto & stage=AtomLocalPotentialView::For(atom).GetJointStageEstimate();
        auto & result=results.emplace_back(id,PostJointQScoreResult{}).second;
        result.source=stage.source;
        result.reference_height=model.GetReferenceHeight();
        result.reference_offset=model.GetReferenceOffset();
        if (!SameEndpoint(stage,expected.at(id))) { result.reason="target-endpoint-unavailable"; continue; }
        LocalPotentialSampleList samples;
        try { samples=PrepareQScoreSamples(atom,map,model); }
        catch (const std::out_of_range &) { result.reason="sample-outside-map"; continue; }
        result.sample_count=samples.size();
        auto neighbors=model.FindNeighborAtoms(atom,4.5+1.0e-12);
        std::erase_if(neighbors,[](const auto * neighbor){return neighbor->GetElement()==Element::HYDROGEN;});
        std::sort(neighbors.begin(),neighbors.end(),[](const auto * a,const auto * b){return a->GetSerialID()<b->GetSerialID();});
        std::set<int> used, missing, halos;
        std::vector<double> observations, joint_reference;
        observations.reserve(samples.size()); joint_reference.reserve(samples.size());
        for (const auto & sample:samples)
        {
            long double prediction{};
            if (!std::isfinite(sample.response)) result.reason="nonfinite-response";
            for (const auto * neighbor:neighbors)
            {
                if (neighbor==&atom) continue;
                const double square=SquaredDistance(sample.point.position,neighbor->GetPositionRef());
                if (square>6.25) continue;
                const int neighbor_id=neighbor->GetSerialID();
                const auto & contributor=contributor_model(*neighbor);
                if (used.insert(neighbor_id).second && contributor.halo)
                {
                    halos.insert(neighbor_id);
                    if (!contributor.extra_donors.empty())
                        result.extra_halos.push_back({neighbor_id,contributor.extra_donors});
                }
                if (!contributor.point)
                {
                    missing.insert(neighbor_id);
                    if (result.reason.empty()) result.reason="contributor-unavailable";
                    continue;
                }
                ++result.evaluated_contributions;
                const double response=contributor.point->ResponseAtDistance(std::sqrt(square));
                if (!std::isfinite(response)) result.reason="nonfinite-contribution";
                else prediction+=response;
            }
            const double peeled=static_cast<double>(static_cast<long double>(sample.response)-prediction);
            if (!std::isfinite(peeled)) result.reason="nonfinite-response";
            observations.push_back(peeled);
            joint_reference.push_back(stage.point->ResponseAtDistance(sample.point.distance));
        }
        result.contributor_count=used.size();
        result.missing_contributor_count=missing.size();
        result.halo_count=halos.size();
        std::sort(result.extra_halos.begin(),result.extra_halos.end(),[](const auto & a,const auto & b){return a.atom_id<b.atom_id;});
        if (result.reason.empty())
        {
            const auto green=CalculateQScoreCorrelation(observations,
                StandardQScoreReference(samples,result.reference_height,result.reference_offset));
            const auto red=CalculateQScoreCorrelation(observations,joint_reference);
            if (!green.value) result.reason="standard:"+green.reason;
            else if (!red.value) result.reason="joint:"+red.reason;
            else
            {
                result.scores=PostJointQScoreValues{*green.value,*red.value};
                result.status=EvidenceStatus::Available;
            }
        }
        metrics.estimated_qscore_peak_bytes=std::max(metrics.estimated_qscore_peak_bytes,
            samples.capacity()*sizeof(LocalPotentialSample)+samples.size()*3*sizeof(double)+
            neighbors.capacity()*sizeof(const AtomObject *));
    }
    for (const auto & [id,result]:results)
    {
        (void)id;
        ++metrics.qscore_targets;
        if (result.scores) ++metrics.qscore_available; else ++metrics.qscore_unavailable;
        metrics.qscore_samples+=result.sample_count;
        metrics.qscore_contributions+=result.evaluated_contributions;
    }
    std::size_t dynamic_payload{};
    for (const auto & [id,entry]:cache)
    {
        (void)id;
        dynamic_payload+=entry.extra_donors.capacity()*sizeof(std::string);
        for (const auto & donor:entry.extra_donors) dynamic_payload+=donor.capacity();
        if (!entry.halo) continue;
        if (!entry.point) { ++metrics.qscore_unavailable_halos; continue; }
        ++metrics.qscore_halos;
        if (!entry.extra_donors.empty()) ++metrics.qscore_extra_halos;
    }
    for (const auto & [id,result]:results)
    {
        (void)id;
        dynamic_payload+=result.source.atom_id.capacity()+result.source.component_id.capacity()+
            result.source.run_id.capacity()+result.reason.capacity()+result.contract.capacity()+
            result.extra_halos.capacity()*sizeof(QScoreExtraHaloSource);
        for (const auto & halo:result.extra_halos)
        {
            dynamic_payload+=halo.donor_atom_ids.capacity()*sizeof(std::string);
            for (const auto & donor:halo.donor_atom_ids) dynamic_payload+=donor.capacity();
        }
    }
    // Container/node payload estimate; excludes allocator overhead and borrowed map/problem.
    metrics.estimated_qscore_peak_bytes+=dynamic_payload+
        cache.size()*(sizeof(decltype(cache)::value_type)+3*sizeof(void *))+
        expected.size()*(sizeof(decltype(expected)::value_type)+3*sizeof(void *))+
        donors.size()*(sizeof(HaloMedianTargetDonor)+8*sizeof(void *))+
        results.capacity()*sizeof(decltype(results)::value_type);
    metrics.qscore_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    return results;
}
} // namespace rhbm_gem::core::detail
