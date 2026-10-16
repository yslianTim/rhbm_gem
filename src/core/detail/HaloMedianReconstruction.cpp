#include "HaloMedianReconstruction.hpp"

#include "data/detail/AtomClassifier.hpp"
#include <rhbm_gem/core/JointComponentEstimator.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/math/ArrayHelper.hpp>
#include <rhbm_gem/utils/math/GaussianModel3D.hpp>
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>

namespace rhbm_gem::core::detail {
namespace {
struct Donor
{
    std::size_t atom{};
    GroupKey group_key{};
    GaussianModel3D point{};
};

std::array<double,3> Median(const std::vector<Donor> & donors, std::size_t offset)
{
    std::vector<double> values;
    values.reserve(donors.size());
    for (const auto & donor:donors)
    {
        const auto value=offset==0 ? donor.point.GetAmplitude() : offset==1 ? donor.point.GetWidth() : donor.point.GetOffset();
        if (!std::isfinite(value)) throw std::invalid_argument("Nonfinite Halo median donor.");
        values.push_back(value);
    }
    if (values.empty()) throw std::invalid_argument("Cannot compute Halo median from empty donors.");
    return {array_helper::ComputeMedian(values),0,0};
}

std::array<double,3> Median3(const std::vector<Donor> & donors)
{
    std::array<double,3> result{};
    for (std::size_t k=0;k<3;++k) result[k]=Median(donors,k)[0];
    return result;
}

std::array<double,3> Dispersion(const std::vector<Donor> & donors, const std::array<double,3> & median)
{
    std::array<double,3> result{};
    for (std::size_t k=0;k<3;++k)
    {
        double maximum{};
        for (const auto & donor:donors)
        {
            const double value=k==0 ? donor.point.GetAmplitude() : k==1 ? donor.point.GetWidth() : donor.point.GetOffset();
            maximum=std::max(maximum,std::abs(value-median[k]));
        }
        result[k]=maximum;
    }
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
        {
            if (!points.emplace(full[k],GaussianModel3D{component.state->ac[2*k],component.state->b[k],
                component.state->ac[2*k+1]}).second)
                throw std::invalid_argument("Duplicate Joint target state for Halo median reconstruction.");
        }
    }
    return points;
}

bool QualifiedTarget(const JointAnalysisResult & result, std::size_t atom)
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

bool SameSelectionDomain(const JointSelectionDomain & left, const JointSelectionDomain & right)
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

void ValidatePoint(const GaussianModel3D & point)
{
    if (!std::isfinite(point.GetAmplitude()) || point.GetAmplitude()<0 ||
        !std::isfinite(point.GetWidth()) || point.GetWidth()<=0 || !std::isfinite(point.GetOffset()))
        throw std::invalid_argument("Invalid Joint target donor parameters.");
}
}

HaloMedianBuildResult BuildHaloMedianReconstruction(
    const ModelObject & model, const JointProblem & problem, const JointAnalysisResult & result,
    std::string_view joint_run_id)
{
    const auto & input=problem.OriginalInput();
    if (!input.selection_domain || input.selection_domain->target_indices.empty())
        throw std::invalid_argument("halo-median-no-selection");
    if (result.atom_ids!=input.atom_ids || !result.selection_domain || !SameSelectionDomain(*result.selection_domain,*input.selection_domain))
        throw std::invalid_argument("halo-median-joint-run-identity-mismatch");
    const auto points=TargetPoints(result);
    std::map<GroupKey,std::vector<Donor>> donors_by_group;
    std::vector<Donor> global_donors;
    for (const auto target:input.selection_domain->target_indices)
    {
        if (!QualifiedTarget(result,target)) continue;
        const auto found=points.find(target);
        if (found==points.end()) continue;
        ValidatePoint(found->second);
        const auto * atom=model.FindAtomPtr(std::stoi(input.atom_ids[target]));
        if (!atom) throw std::invalid_argument("halo-median-unknown-target");
        const auto key=data_internal::GetGroupKey(atom);
        Donor donor{target,key,found->second};
        donors_by_group[key].push_back(donor);
        global_donors.push_back(std::move(donor));
    }
    if (global_donors.empty()) throw std::invalid_argument("halo-median-no-valid-targets");

    HaloMedianBuildResult out;
    out.snapshot.joint_run_id=std::string(joint_run_id);
    out.snapshot.valid_target_count=global_donors.size();
    out.snapshot.global_median=Median3(global_donors);
    out.audit.valid_target_count=global_donors.size();
    out.audit.global_median=out.snapshot.global_median;
    out.audit.group_count=donors_by_group.size();

    std::map<GroupKey,std::array<double,3>> group_medians;
    for (const auto & [key, donors]:donors_by_group)
    {
        const auto median=Median3(donors);
        group_medians.emplace(key,median);
        out.audit.groups.push_back({key,donors.size(),0,0,median,Dispersion(donors,median),false});
    }
    std::map<GroupKey,std::size_t> group_audit_index;
    for (std::size_t i=0;i<out.audit.groups.size();++i) group_audit_index.emplace(out.audit.groups[i].group_key,i);

    const auto & targets=input.selection_domain->target_indices;
    for (std::size_t atom=0;atom<input.atom_ids.size();++atom)
    {
        if (std::binary_search(targets.begin(),targets.end(),atom)) continue;
        const auto * object=model.FindAtomPtr(std::stoi(input.atom_ids[atom]));
        if (!object) throw std::invalid_argument("halo-median-unknown-halo");
        const auto key=data_internal::GetGroupKey(object);
        const auto found=group_medians.find(key);
        const bool fallback=found==group_medians.end();
        const auto median=fallback ? out.snapshot.global_median : found->second;
        const auto donor_count=fallback ? global_donors.size() : donors_by_group.at(key).size();
        auto & audit=out.audit.groups;
        if (fallback)
        {
            out.audit.fallback_halos++;
            std::set<std::size_t> rows;
            for (const auto & support:input.support[atom]) rows.insert(support.row);
            out.audit.fallback_support_rows+=rows.size();
        }
        else
        {
            out.audit.same_group_halos++;
            auto & record=audit.at(group_audit_index.at(key));
            ++record.halo_count;
            std::set<std::size_t> rows;
            for (const auto & support:input.support[atom]) rows.insert(support.row);
            record.affected_support_rows+=rows.size();
        }
        out.snapshot.halos.push_back({atom,input.atom_ids[atom],key,median[0],median[1],median[2],
            fallback ? HaloMedianSource::GlobalTargetMedian : HaloMedianSource::GroupKeyMedian,
            donor_count,std::string(joint_run_id),std::string(HaloMedianReconstructionContractVersion)});
    }
    return out;
}
}
