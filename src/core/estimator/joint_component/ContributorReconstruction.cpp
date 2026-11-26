#include "ContributorReconstruction.hpp"

#include "WidthBounds.hpp"
#include "TargetEvidence.hpp"
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/algorithm/KDTreeAlgorithm.hpp>
#include <rhbm_gem/utils/math/ArrayHelper.hpp>
#include <algorithm>
#include <array>
#include <cmath>
#include <chrono>
#include <set>
#include <stdexcept>

namespace rhbm_gem::core::detail {
namespace {
struct TargetNode
{
    const ReconstructionDonor * donor{};
    std::array<double,3> position{};
    const std::array<double,3> & GetPosition() const { return position; }
};

double SquaredDistance(const std::array<double,3> & left,const std::array<double,3> & right)
{
    double square{};
    for (std::size_t axis=0;axis<3;++axis)
    {
        const auto delta=left[axis]-right[axis];
        square+=delta*delta;
    }
    return square;
}

struct ContributorAssessment
{
    std::vector<std::optional<std::array<double,3>>> joint_points;
    std::vector<unsigned> reasons;
};

using joint_component::WidthBoundaryStatus;
using joint_component::WidthBounds;
using Status=JointCheckStatus;

WidthBounds AcceptedWidthBounds(const JointAnalysisResult & result)
{
    WidthBounds bounds;
    bool seen_upper{}, seen_lower{};
    for (const auto & check:result.evidence)
    {
        const bool upper=check.name=="global-width-calibration";
        if (!upper && check.name!="global-width-lower-bound") continue;
        auto & seen=upper ? seen_upper : seen_lower;
        if (seen || check.scope!=JointEvidenceScope::AssembledGlobal)
            throw std::invalid_argument("Invalid accepted Joint width bound provenance.");
        seen=true;
        if (!upper && check.status==Status::NotRun && !check.threshold && !check.value && !check.reason.empty())
            continue;
        if (check.status!=Status::Passed || !check.threshold || !std::isfinite(*check.threshold) ||
            *check.threshold<=0)
            throw std::invalid_argument("Invalid accepted Joint width bound provenance.");
        (upper ? bounds.upper : bounds.lower)=std::log(*check.threshold);
    }
    if (!bounds.IsValid()) throw std::invalid_argument("Invalid accepted Joint width interval.");
    return bounds;
}

ContributorAssessment AssessContributors(const JointAnalysisResult & result)
{
    const auto count=result.atom_ids.size();
    if (!IsQualifiedJointAnalysisResult(result) || !result.selection_domain ||
        !result.selection_domain->IsValid(count) || !result.target_evidence ||
        result.target_evidence->contract!=JointTargetPointValidationContractVersion)
        throw std::invalid_argument("Contributor reconstruction requires a qualified Joint endpoint.");
    const auto & domain=*result.selection_domain;
    const auto bounds=AcceptedWidthBounds(result);
    ContributorAssessment assessment;
    assessment.joint_points.resize(count);
    assessment.reasons.resize(count,ExternalHaloReconstruction);
    std::vector<const JointAtomPointValidation *> evidence(count);
    for (const auto & point:result.target_evidence->point_atoms)
    {
        if (point.atom>=count || !IsFittingTarget(domain.RoleOf(point.atom)) || evidence[point.atom])
            throw std::invalid_argument("Invalid Joint contributor evidence mapping.");
        evidence[point.atom]=&point;
    }
    for (const auto & component:result.components)
    {
        const auto & atoms=component.layout ? component.layout->full_atoms : component.atoms;
        if (!component.state || component.runtime_convergence!=Status::Passed ||
            component.target_runtime_convergence!=Status::Passed ||
            component.state->ac.size()!=2*atoms.size() || component.state->tau.size()!=atoms.size() ||
            component.state->log_tau.size()!=atoms.size())
            throw std::invalid_argument("Invalid certified Joint contributor state.");
        const auto & state=*component.state;
        const auto * component_points=component.target_evidence ? &component.target_evidence->point_atoms : nullptr;
        if (component_points && component_points->size()!=atoms.size())
            throw std::invalid_argument("Invalid certified Joint component point coverage.");
        for (std::size_t local=0;local<atoms.size();++local)
        {
            const auto atom=atoms[local];
            if (atom>=count || !evidence[atom] || assessment.joint_points[atom])
                throw std::invalid_argument("Invalid Joint contributor state mapping.");
            if (component_points)
            {
                const auto & local_point=(*component_points)[local];
                if (local_point.atom!=atom || local_point.component_convergence!=Status::Passed ||
                    local_point.parameter_validity!=Status::Passed)
                    throw std::invalid_argument("Invalid certified Joint component point evidence.");
            }
            const auto & quality=*evidence[atom];
            const auto a=state.ac[2*local], tau=state.tau[local], c=state.ac[2*local+1];
            if (quality.component_convergence!=Status::Passed || quality.parameter_validity!=Status::Passed ||
                !joint_component::ValidJointPointParameters(a,tau,c,state.log_tau[local]) ||
                !bounds.ContainsWithRoundoff(std::log(tau)))
                throw std::invalid_argument("Invalid certified Joint contributor endpoint.");
            unsigned reasons=0;
            switch (bounds.Classify(state.log_tau[local]))
            {
            case WidthBoundaryStatus::LowerBound: reasons=WidthLowerBoundReconstruction; break;
            case WidthBoundaryStatus::UpperBound: reasons=WidthUpperBoundReconstruction; break;
            case WidthBoundaryStatus::Unavailable:
                throw std::invalid_argument("Invalid accepted Joint contributor width.");
            case WidthBoundaryStatus::Interior: break;
            }
            if (quality.local_conditioning!="Adequate")
            {
                if (quality.local_conditioning!="Warning" && quality.local_conditioning!="Unavailable" &&
                    quality.local_conditioning!="NotAssessed")
                    throw std::invalid_argument("Invalid Joint contributor conditioning evidence.");
                reasons|=ConditioningNotAdequateReconstruction;
            }
            assessment.joint_points[atom]=std::array{a,tau,c};
            assessment.reasons[atom]=reasons;
        }
    }
    for (std::size_t atom=0;atom<count;++atom)
        if (IsFittingTarget(domain.RoleOf(atom)) && !assessment.joint_points[atom])
            throw std::invalid_argument("Missing certified Joint contributor endpoint.");
    return assessment;
}

std::array<double,3> ContributorPosition(const ModelObject & model,const std::string & id)
{
    const auto * atom=model.FindAtomPtr(std::stoi(id));
    if (!atom || std::to_string(atom->GetSerialID())!=id || atom->GetElement()==Element::HYDROGEN)
        throw std::invalid_argument("Invalid Joint contributor identity.");
    const auto position=atom->GetPosition();
    if (!std::all_of(position.begin(),position.end(),[](double value){return std::isfinite(value);}))
        throw std::invalid_argument("Nonfinite Joint contributor position.");
    return position;
}

std::vector<ReconstructionDonor> CollectDonors(const JointAnalysisResult & result,
    const ContributorAssessment & assessment,const ModelObject * model)
{
    std::vector<ReconstructionDonor> out;
    for (std::size_t atom=0;atom<result.atom_ids.size();++atom)
    {
        if (assessment.reasons[atom]) continue;
        const auto & id=result.atom_ids[atom];
        out.push_back({atom,id,model ? ContributorPosition(*model,id) : std::array<double,3>{},
            *assessment.joint_points[atom]});
    }
    std::sort(out.begin(),out.end(),[](const auto & left,const auto & right) { return left.atom_id<right.atom_id; });
    return out;
}

JointContributorReplacement MakeReplacement(const JointAnalysisResult & result,std::size_t atom,
    unsigned reasons,std::string_view run,const std::vector<const ReconstructionDonor *> & donors)
{
    JointContributorReplacement replacement;
    replacement.contributor_index=atom;
    replacement.atom_id=result.atom_ids[atom];
    replacement.role=result.selection_domain->RoleOf(atom);
    replacement.rejection_reasons=reasons;
    replacement.joint_run_id=run;
    replacement.donor_count=donors.size();
    if (donors.empty())
        replacement.source=ReconstructionSource::Unavailable;
    else
    {
        const auto median=GeometryLocalFittingMedian(donors);
        replacement.amplitude=median[0];
        replacement.width=median[1];
        replacement.charge=median[2];
        for (const auto * donor:donors) replacement.donor_atom_ids.push_back(donor->atom_id);
    }
    return replacement;
}
}

struct ReconstructionDonorIndex::Impl
{
    std::vector<ReconstructionDonor> donors;
    std::vector<TargetNode> nodes;
    std::unique_ptr<KDNode<TargetNode>> tree;

    explicit Impl(std::vector<ReconstructionDonor> source):donors(std::move(source))
    {
        nodes.reserve(donors.size());
        for (const auto & donor:donors) nodes.push_back({&donor,donor.position});
        tree=KDTreeAlgorithm<TargetNode>::BuildKDTree(nodes);
    }
};

ReconstructionDonorIndex::ReconstructionDonorIndex(
    std::vector<ReconstructionDonor> donors):m_impl(std::make_unique<Impl>(std::move(donors)))
{
}

ReconstructionDonorIndex::~ReconstructionDonorIndex()=default;

const std::vector<ReconstructionDonor> & ReconstructionDonorIndex::Donors() const
{
    return m_impl->donors;
}

std::vector<const ReconstructionDonor *> ReconstructionDonorIndex::Query(
    const std::array<double,3> & position) const
{
    const auto & nodes=m_impl->nodes;
    if (nodes.empty()) return {};
    TargetNode query{nullptr,position};
    const auto count=std::min(GeometryLocalFittingMedianNeighborCount,nodes.size());
    const auto nearest=KDTreeAlgorithm<TargetNode>::KNearestNeighbors(m_impl->tree.get(),&query,count);
    if (nearest.empty()) return {};
    double radius{};
    for (const auto * donor:nearest)
        radius=std::max(radius,std::sqrt(SquaredDistance(position,donor->position)));
    auto candidates=KDTreeAlgorithm<TargetNode>::RangeSearch(m_impl->tree.get(),&query,radius+1.0e-12);
    std::sort(candidates.begin(),candidates.end(),[&](const auto * left,const auto * right) {
        const auto left_distance=SquaredDistance(position,left->position);
        const auto right_distance=SquaredDistance(position,right->position);
        if (left_distance!=right_distance) return left_distance<right_distance;
        return left->donor->atom_id<right->donor->atom_id;
    });
    if (candidates.size()>count) candidates.resize(count);
    std::vector<const ReconstructionDonor *> selected;
    selected.reserve(candidates.size());
    for (const auto * candidate:candidates) selected.push_back(candidate->donor);
    return selected;
}

std::array<double,3> GeometryLocalFittingMedian(const std::vector<const ReconstructionDonor *> & donors)
{
    if (donors.empty()) throw std::invalid_argument("Cannot compute contributor median from empty donors.");
    std::array<double,3> result{};
    for (std::size_t parameter=0;parameter<3;++parameter)
    {
        std::vector<double> values;
        values.reserve(donors.size());
        for (const auto * donor:donors)
        {
            const auto value=donor->point[parameter];
            if (!std::isfinite(value)) throw std::invalid_argument("Nonfinite contributor median donor.");
            values.push_back(value);
        }
        result[parameter]=array_helper::ComputeMedian(values);
    }
    return result;
}

JointContributorReconstructionSnapshot ContributorReconstructionBuildResult::TakeSnapshot()
{
    if (!snapshot) throw std::logic_error("contributor median snapshot was already transferred.");
    auto result = std::move(*snapshot);
    snapshot.reset();
    return result;
}

std::vector<ReconstructionDonor> CollectQualifiedReconstructionDonors(
    const ModelObject & model,const JointAnalysisResult & result)
{
    return CollectDonors(result,AssessContributors(result),&model);
}

ContributorReconstructionBuildResult BuildContributorReconstruction(
    const ModelObject & model,const JointProblem & problem,const JointAnalysisResult & result,
    std::string_view joint_run_id)
{
    const auto & input=problem.OriginalInput();
    if (joint_run_id.empty() || result.atom_ids!=input.atom_ids || result.selection_domain!=input.selection_domain)
        throw std::invalid_argument("Contributor reconstruction Joint run identity mismatch.");

    using Clock=std::chrono::steady_clock;
    auto phase=Clock::now();
    auto assessment=AssessContributors(result);
    auto donors=CollectDonors(result,assessment,&model);
    ContributorReconstructionBuildResult out;
    out.audit.collection_seconds=std::chrono::duration<double>(Clock::now()-phase).count();
    out.snapshot.emplace();
    auto & snapshot=*out.snapshot;
    snapshot.joint_run_id=joint_run_id;
    snapshot.qualified_donor_count=out.audit.qualified_donor_count=donors.size();
    phase=Clock::now();
    out.donor_index=std::make_unique<ReconstructionDonorIndex>(std::move(donors));
    out.audit.index_seconds=std::chrono::duration<double>(Clock::now()-phase).count();

    phase=Clock::now();
    for (std::size_t atom=0;atom<result.atom_ids.size();++atom)
    {
        const auto reasons=assessment.reasons[atom];
        if (!reasons) continue;
        const auto nearest=out.donor_index->Query(ContributorPosition(model,result.atom_ids[atom]));
        auto replacement=MakeReplacement(result,atom,reasons,joint_run_id,nearest);
        const bool lower=(reasons&WidthLowerBoundReconstruction)!=0;
        const bool upper=(reasons&WidthUpperBoundReconstruction)!=0;
        const bool conditioning=(reasons&ConditioningNotAdequateReconstruction)!=0;
        out.audit.lower_bound+=lower;
        out.audit.upper_bound+=upper;
        out.audit.conditioning_rejected+=conditioning;
        out.audit.overlapping_rejections+=(lower || upper) && conditioning;
        if (nearest.empty()) ++out.audit.unavailable;
        else if (replacement.role==FittingRole::OutputTarget) ++out.audit.replaced_outputs;
        else if (replacement.role==FittingRole::BufferTarget) ++out.audit.replaced_buffers;
        else ++out.audit.reconstructed_halos;
        snapshot.replacements.push_back(std::move(replacement));
    }
    out.audit.query_seconds=std::chrono::duration<double>(Clock::now()-phase).count();
    if (result.contributor_reconstruction &&
        (*result.contributor_reconstruction!=snapshot ||
            HasUnavailableJointContributorReconstruction(result)!=(out.audit.unavailable!=0)))
        throw std::invalid_argument("Invalid Joint contributor reconstruction provenance.");
    out.joint_points=std::move(assessment.joint_points);
    return out;
}

void ValidateSavedContributorReconstruction(const JointAnalysisResult & result,const ModelObject * model)
{
    const auto & saved=result.contributor_reconstruction;
    bool unavailable=false;
    for (const auto & check:result.evidence)
    {
        if (check.name!=JointContributorReconstructionEvidenceName) continue;
        if (unavailable || check.status!=Status::Unavailable || check.scope!=JointEvidenceScope::AssembledGlobal ||
            check.value || check.threshold || check.reason!=NoQualifiedReconstructionDonors)
            throw std::invalid_argument("Invalid contributor reconstruction unavailable evidence.");
        unavailable=true;
    }
    if (!IsQualifiedJointAnalysisResult(result))
    {
        if (saved || unavailable)
            throw std::invalid_argument("Unqualified Joint endpoint has contributor reconstruction.");
        return;
    }
    const auto assessment=AssessContributors(result);
    const auto needed=static_cast<std::size_t>(std::count_if(assessment.reasons.begin(),assessment.reasons.end(),
        [](unsigned reasons){return reasons!=0;}));
    if (!saved)
    {
        if (needed || unavailable) throw std::invalid_argument("Missing Joint contributor reconstruction.");
        return;
    }
    auto donors=CollectDonors(result,assessment,model);
    const auto donor_count=std::min(GeometryLocalFittingMedianNeighborCount,donors.size());
    if (saved->joint_run_id.empty() || saved->contract!=ContributorReconstructionContractVersion ||
        saved->donor_policy!=ContributorReconstructionDonorPolicy ||
        saved->neighbor_count!=GeometryLocalFittingMedianNeighborCount ||
        saved->qualified_donor_count!=donors.size() || saved->replacements.size()!=needed ||
        unavailable!=(needed!=0 && donors.empty()))
        throw std::invalid_argument("Invalid Joint contributor reconstruction coverage or contract.");

    std::optional<ReconstructionDonorIndex> index;
    if (model && needed && !donors.empty()) index.emplace(std::move(donors));
    auto record=saved->replacements.begin();
    for (std::size_t atom=0;atom<assessment.reasons.size();++atom)
    {
        const auto reasons=assessment.reasons[atom];
        if (!reasons) continue;
        if (record->donor_atom_ids.size()!=donor_count)
            throw std::invalid_argument("Invalid Joint contributor donor count.");
        std::vector<const ReconstructionDonor *> selected;
        if (index) selected=index->Query(ContributorPosition(*model,result.atom_ids[atom]));
        else
        {
            std::set<std::string> seen;
            for (const auto & id:record->donor_atom_ids)
            {
                const auto donor=std::lower_bound(donors.begin(),donors.end(),id,
                    [](const auto & candidate,const auto & key){return candidate.atom_id<key;});
                if (donor==donors.end() || donor->atom_id!=id || !seen.insert(id).second)
                    throw std::invalid_argument("Invalid Joint contributor donor identity.");
                selected.push_back(&*donor);
            }
        }
        if (*record!=MakeReplacement(result,atom,reasons,saved->joint_run_id,selected))
            throw std::invalid_argument("Invalid Joint contributor reconstruction provenance.");
        ++record;
    }
}
}
