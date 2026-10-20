#include "HaloMedianReconstruction.hpp"

#include <rhbm_gem/core/JointComponentEstimator.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/algorithm/KDTreeAlgorithm.hpp>
#include <rhbm_gem/utils/math/ArrayHelper.hpp>
#include <rhbm_gem/utils/math/GaussianModel3D.hpp>
#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <set>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

namespace rhbm_gem::core::detail {
namespace {
struct TargetNode
{
    const HaloMedianTargetDonor * donor{};
    std::array<double,3> position{};
    const std::array<double,3> & GetPosition() const { return position; }
};

std::array<double,3> Point(const GaussianModel3D & point)
{
    return {point.GetAmplitude(),point.GetWidth(),point.GetOffset()};
}

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

std::array<double,3> Median(const std::vector<const HaloMedianTargetDonor *> & donors)
{
    if (donors.empty()) throw std::invalid_argument("Cannot compute Halo median from empty donors.");
    std::array<double,3> result{};
    for (std::size_t parameter=0;parameter<3;++parameter)
    {
        std::vector<double> values;
        values.reserve(donors.size());
        for (const auto * donor:donors)
        {
            const auto value=donor->point[parameter];
            if (!std::isfinite(value)) throw std::invalid_argument("Nonfinite Halo median donor.");
            values.push_back(value);
        }
        result[parameter]=array_helper::ComputeMedian(values);
    }
    return result;
}

std::array<double,3> Dispersion(
    const std::vector<const HaloMedianTargetDonor *> & donors,const std::array<double,3> & median)
{
    std::array<double,3> result{};
    for (std::size_t parameter=0;parameter<3;++parameter)
        for (const auto * donor:donors)
            result[parameter]=std::max(result[parameter],std::abs(donor->point[parameter]-median[parameter]));
    return result;
}

std::map<std::size_t,GaussianModel3D> TargetPoints(const JointAnalysisResult & result)
{
    std::map<std::size_t,GaussianModel3D> points;
    for (const auto & component:result.components)
    {
        if (!component.state) continue;
        const auto & full=component.layout ? component.layout->full_atoms : component.atoms;
        if (component.state->ac.size()!=2*full.size() || component.state->b.size()!=full.size())
            throw std::invalid_argument("Invalid Joint target state for Halo median reconstruction.");
        for (std::size_t k=0;k<full.size();++k)
            if (!points.emplace(full[k],GaussianModel3D{component.state->ac[2*k],component.state->b[k],
                component.state->ac[2*k+1]}).second)
                throw std::invalid_argument("Duplicate Joint target state for Halo median reconstruction.");
    }
    return points;
}

bool QualifiedTarget(const JointAnalysisResult & result,std::size_t atom)
{
    if (result.runtime_convergence!=JointCheckStatus::Passed ||
        result.target_runtime_convergence!=JointCheckStatus::Passed || !result.target_evidence)
        return false;
    const auto found=std::find_if(result.target_evidence->point_atoms.begin(),result.target_evidence->point_atoms.end(),
        [atom](const auto & point){return point.atom==atom;});
    return found!=result.target_evidence->point_atoms.end() &&
        found->component_convergence==JointCheckStatus::Passed &&
        found->parameter_validity==JointCheckStatus::Passed &&
        found->local_conditioning=="Adequate";
}

bool SameSelectionDomain(const JointSelectionDomain & left,const JointSelectionDomain & right)
{
    const bool same_pure=left.pure_halo_support_exclusion.has_value()==right.pure_halo_support_exclusion.has_value() &&
        (!left.pure_halo_support_exclusion || *left.pure_halo_support_exclusion==*right.pure_halo_support_exclusion);
    const bool same_fixed=left.fixed_halo_selection.has_value()==right.fixed_halo_selection.has_value() &&
        (!left.fixed_halo_selection || *left.fixed_halo_selection==*right.fixed_halo_selection);
    return left.contract==right.contract && left.target_indices==right.target_indices &&
        left.observation_radius==right.observation_radius && left.support_radius==right.support_radius &&
        left.contributor_policy==right.contributor_policy &&
        left.halo_eligibility.has_value()==right.halo_eligibility.has_value() && same_fixed && same_pure;
}

bool ValidPoint(const GaussianModel3D & point)
{
    return std::isfinite(point.GetAmplitude()) && point.GetAmplitude()>=0 &&
        std::isfinite(point.GetWidth()) && point.GetWidth()>0 && std::isfinite(point.GetOffset());
}
}

struct GeometryLocalTargetDonorIndex::Impl
{
    const std::vector<HaloMedianTargetDonor> & donors;
    std::vector<TargetNode> nodes;
    std::unique_ptr<KDNode<TargetNode>> tree;

    explicit Impl(const std::vector<HaloMedianTargetDonor> & source):donors(source)
    {
        nodes.reserve(donors.size());
        for (const auto & donor:donors) nodes.push_back({&donor,donor.position});
        tree=KDTreeAlgorithm<TargetNode>::BuildKDTree(nodes);
    }
};

GeometryLocalTargetDonorIndex::GeometryLocalTargetDonorIndex(
    const std::vector<HaloMedianTargetDonor> & donors):m_impl(std::make_unique<Impl>(donors))
{
}

GeometryLocalTargetDonorIndex::~GeometryLocalTargetDonorIndex()=default;

std::vector<const HaloMedianTargetDonor *> GeometryLocalTargetDonorIndex::Query(
    const std::array<double,3> & position) const
{
    const auto & nodes=m_impl->nodes;
    if (nodes.empty()) return {};
    TargetNode query{nullptr,position};
    const auto count=std::min(GeometryLocalTargetMedianNeighborCount,nodes.size());
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
    std::vector<const HaloMedianTargetDonor *> selected;
    selected.reserve(candidates.size());
    for (const auto * candidate:candidates) selected.push_back(candidate->donor);
    return selected;
}

std::array<double,3> GeometryLocalTargetMedian(
    const std::vector<const HaloMedianTargetDonor *> & donors)
{
    return Median(donors);
}

std::vector<std::size_t> DistinctOriginalSupportRows(const std::vector<JointSupport> & support)
{
    std::vector<std::size_t> rows;
    rows.reserve(support.size());
    for(const auto & entry:support) rows.push_back(entry.row);
    std::sort(rows.begin(),rows.end());
    rows.erase(std::unique(rows.begin(),rows.end()),rows.end());
    return rows;
}

std::vector<HaloMedianTargetDonor> QualifiedJointTargetDonors(
    const ModelObject & model,const JointProblem & problem,const JointAnalysisResult & result)
{
    const auto & input=problem.OriginalInput();
    if (!input.selection_domain || input.selection_domain->target_indices.empty())
        throw std::invalid_argument("halo-median-no-selection");
    if (result.atom_ids!=input.atom_ids || !result.selection_domain ||
        !SameSelectionDomain(*result.selection_domain,*input.selection_domain))
        throw std::invalid_argument("halo-median-joint-run-identity-mismatch");
    const auto points=TargetPoints(result);
    std::vector<HaloMedianTargetDonor> out;
    for (const auto target:input.selection_domain->target_indices)
    {
        if (!QualifiedTarget(result,target)) continue;
        const auto found=points.find(target);
        if (found==points.end() || !ValidPoint(found->second)) continue;
        const auto atom_id=input.atom_ids[target];
        const auto * atom=model.FindAtomPtr(std::stoi(atom_id));
        if (!atom) throw std::invalid_argument("halo-median-unknown-target");
        const auto position=atom->GetPosition();
        if (!std::all_of(position.begin(),position.end(),[](double value){return std::isfinite(value);}))
            throw std::invalid_argument("halo-median-nonfinite-target-position");
        out.push_back({target,atom_id,position,Point(found->second)});
    }
    std::sort(out.begin(),out.end(),[](const auto & left,const auto & right) {
        return left.atom_id<right.atom_id;
    });
    return out;
}

HaloFullABCPreclassification BuildHaloFullABCPreclassification(
    const ModelObject & model,const JointProblem & problem,const JointAnalysisResult & result)
{
    const auto & input=problem.OriginalInput();
    if(!input.selection_domain || input.selection_domain->target_indices.empty())
        throw std::invalid_argument("halo-median-no-selection");
    const auto donors=QualifiedJointTargetDonors(model,problem,result);
    HaloFullABCPreclassification out;
    if(donors.empty()) out.unavailable_reason="no-qualified-joint-target-donors";
    else
    {
        std::vector<const HaloMedianTargetDonor *> donor_refs;
        donor_refs.reserve(donors.size());
        for(const auto & donor:donors) donor_refs.push_back(&donor);
        const auto median=GeometryLocalTargetMedian(donor_refs);
        if(!std::isfinite(median[0]) || median[0]<0 || !std::isfinite(median[1]) || !(median[1]>0) ||
            !std::isfinite(median[2]))
            out.unavailable_reason="invalid-qualified-joint-target-median";
        else
        {
            out.available=true;
            out.global_target_median=median;
        }
    }

    std::unordered_set<int> target_ids;
    for(const auto target:input.selection_domain->target_indices)
        target_ids.insert(std::stoi(input.atom_ids.at(target)));
    std::unordered_map<int,std::size_t> contributor_indices;
    contributor_indices.reserve(input.atom_ids.size());
    for(std::size_t index=0;index<input.atom_ids.size();++index)
        contributor_indices.emplace(std::stoi(input.atom_ids[index]),index);
    std::vector<const AtomObject *> candidates;
    for(const auto & atom:model.GetAtomList())
        if(atom->GetElement()!=Element::HYDROGEN && !target_ids.contains(atom->GetSerialID()))
            candidates.push_back(atom.get());
    std::sort(candidates.begin(),candidates.end(),[](const auto * left,const auto * right) {
        return left->GetSerialID()<right->GetSerialID();
    });
    out.halos.reserve(candidates.size());
    for(const auto * atom:candidates)
    {
        HaloFullABCClassificationRecord record;
        record.atom_id=std::to_string(atom->GetSerialID());
        if(const auto found=contributor_indices.find(atom->GetSerialID());found!=contributor_indices.end())
        {
            record.original_contributor_index=found->second;
            record.original_support_rows=DistinctOriginalSupportRows(input.support.at(found->second));
            for(const auto row:record.original_support_rows)
            {
                if(row>=input.row_ids.size())
                    throw std::invalid_argument("Invalid Halo observation support row.");
            }
        }
        record.distinct_original_observations=record.original_support_rows.size();
        if(record.distinct_original_observations<3)
        {
            record.role=HaloFullABCRole::FixedGlobalTargetMedian;
            ++out.fixed_halo_count;
        }
        else
        {
            record.role=HaloFullABCRole::FreeFullABC;
            ++out.free_halo_count;
        }
        if(out.available) record.parameter_profile=*out.global_target_median;
        out.halos.push_back(std::move(record));
    }
    return out;
}

HaloMedianBuildResult BuildHaloMedianReconstruction(
    const ModelObject & model,const JointProblem & problem,const JointAnalysisResult & result,
    std::string_view joint_run_id)
{
    const auto & input=problem.OriginalInput();
    if (!input.selection_domain || input.selection_domain->target_indices.empty())
        throw std::invalid_argument("halo-median-no-selection");
    HaloMedianBuildResult out;
    out.target_donors=QualifiedJointTargetDonors(model,problem,result);
    if(out.target_donors.empty()) throw std::invalid_argument("halo-median-no-valid-targets");
    out.snapshot.joint_run_id=std::string(joint_run_id);
    out.snapshot.neighbor_count=GeometryLocalTargetMedianNeighborCount;
    out.snapshot.donor_policy=std::string(GeometryLocalTargetMedianDonorPolicy);
    std::vector<const HaloMedianTargetDonor *> all_donors;
    all_donors.reserve(out.target_donors.size());
    for(const auto & donor:out.target_donors) all_donors.push_back(&donor);
    out.snapshot.valid_target_count=all_donors.size();
    out.snapshot.global_median=Median(all_donors);
    out.audit.valid_target_count=all_donors.size();
    out.audit.global_median=out.snapshot.global_median;
    const GeometryLocalTargetDonorIndex donor_index(out.target_donors);

    const auto & targets=input.selection_domain->target_indices;
    for (std::size_t atom=0;atom<input.atom_ids.size();++atom)
    {
        if (std::binary_search(targets.begin(),targets.end(),atom)) continue;
        const auto * object=model.FindAtomPtr(std::stoi(input.atom_ids[atom]));
        if (!object) throw std::invalid_argument("halo-median-unknown-halo");
        const auto donors=donor_index.Query(object->GetPosition());
        if (donors.empty()) throw std::invalid_argument("halo-median-no-valid-targets");
        const auto median=GeometryLocalTargetMedian(donors);
        JointHaloMedianReconstruction halo;
        halo.halo_contributor_index=atom;
        halo.atom_id=input.atom_ids[atom];
        halo.amplitude=median[0];
        halo.width=median[1];
        halo.charge=median[2];
        halo.source=HaloMedianSource::GeometryLocalTargetMedian;
        halo.donor_count=donors.size();
        halo.joint_run_id=std::string(joint_run_id);
        halo.contract=std::string(HaloMedianReconstructionContractVersion);
        for (const auto * donor:donors) halo.donor_atom_ids.push_back(donor->atom_id);
        out.snapshot.halos.push_back(std::move(halo));

        std::set<std::size_t> rows;
        for (const auto & support:input.support[atom]) rows.insert(support.row);
        out.audit.halos.push_back({input.atom_ids[atom],donors.size(),rows.size(),median,
            Dispersion(donors,median)});
    }
    out.audit.geometry_local_halos=out.snapshot.halos.size();
    return out;
}
}
