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
    const JointAnalysisResult & snapshot, const ContributorReconstructionBuildResult * reconstruction,
    JointPostprocessingMetrics & metrics)
{
    const auto start=std::chrono::steady_clock::now();
    if (!IsQualifiedJointAnalysisResult(snapshot) || snapshot.atom_ids!=problem.OriginalInput().atom_ids ||
        !snapshot.selection_domain || snapshot.selection_domain->output_indices.empty())
        throw std::invalid_argument("QScore requires a qualified matching Joint snapshot.");
    const auto & original=problem.OriginalInput();
    if (snapshot.selection_domain!=original.selection_domain ||
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
    const auto & output_indices=snapshot.selection_domain->output_indices;
    const auto * first=model.FindAtomPtr(std::stoi(snapshot.atom_ids.at(output_indices.front())));
    const auto run=AtomLocalPotentialView::For(*first).GetJointStageEstimate().source.run_id;
    if (run.empty()) throw std::invalid_argument("QScore requires Joint run identity.");
    const auto expected=data_internal::BuildJointStageEstimates(snapshot,run);
    std::optional<ContributorReconstructionBuildResult> rebuilt;
    if (!reconstruction || reconstruction->published_snapshot!=&snapshot)
        rebuilt=BuildContributorReconstruction(model,problem,snapshot,run);
    const auto & workspace=rebuilt ? *rebuilt : *reconstruction;
    const auto & donor_index=*workspace.donor_index;
    const auto & donors=donor_index.Donors();
    bool donors_valid=true;
    for (const auto & donor:donors)
    {
        const auto * atom=model.FindAtomPtr(std::stoi(donor.atom_id));
        const auto view=AtomLocalPotentialView::For(*atom);
        const auto & endpoint=expected.at(atom->GetSerialID());
        if (!endpoint.point || donor.position!=atom->GetPositionRef() ||
            donor.point!=std::array{endpoint.point->GetAmplitude(),endpoint.point->GetWidth(),endpoint.point->GetOffset()} ||
            !view.IsAvailable() || !SameEndpoint(view.GetJointStageEstimate(),endpoint))
            donors_valid=false;
    }
    std::map<int,const JointContributorReplacement *> replacements;
    if (snapshot.contributor_reconstruction)
    {
        const auto & saved=*snapshot.contributor_reconstruction;
        if (saved.joint_run_id!=run)
            throw std::invalid_argument("QScore contributor reconstruction run mismatch.");
        for (const auto & replacement:saved.replacements)
            replacements.emplace(std::stoi(replacement.atom_id),&replacement);
    }
    else if (workspace.audit.direct_contributor_count!=snapshot.atom_ids.size())
        throw std::invalid_argument("QScore is missing Joint contributor reconstruction.");
    std::map<int,ContributorModel> cache;
    const auto contributor_model=[&](const AtomObject & atom) -> const ContributorModel & {
        const int id=atom.GetSerialID();
        if (const auto found=cache.find(id);found!=cache.end()) return found->second;
        auto & result=cache[id];
        const auto endpoint=expected.find(id);
        result.halo=endpoint==expected.end() || endpoint->second.source.role==FittingRole::ExternalHalo;
        if (!result.halo)
        {
            const auto view=AtomLocalPotentialView::For(atom);
            if (!view.IsAvailable() || !SameEndpoint(view.GetJointStageEstimate(),endpoint->second))
                return result;
        }
        if (const auto saved=replacements.find(id);saved!=replacements.end())
        {
            const auto & replacement=*saved->second;
            if (donors_valid && replacement.source==ReconstructionSource::GeometryLocalFittingMedian)
                result.point=GaussianModel3D{replacement.amplitude,replacement.width,replacement.charge};
        }
        else if (!result.halo)
            result.point=endpoint->second.point;
        else if (endpoint==expected.end() && donors_valid)
        {
            const auto nearest=donor_index.Query(atom.GetPositionRef());
            if (!nearest.empty())
            {
                const auto median=GeometryLocalFittingMedian(nearest);
                result.point=GaussianModel3D{median[0],median[1],median[2]};
                for (const auto * donor:nearest) result.extra_donors.push_back(donor->atom_id);
            }
        }
        return result;
    };
    std::vector<std::pair<int,PostJointQScoreResult>> results;
    results.reserve(output_indices.size());
    for (const auto index:output_indices)
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
        donors.size()*(sizeof(ReconstructionDonor)+8*sizeof(void *))+
        results.capacity()*sizeof(decltype(results)::value_type);
    metrics.qscore_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    return results;
}
} // namespace rhbm_gem::core::detail
