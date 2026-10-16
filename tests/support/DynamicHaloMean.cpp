#include "DynamicHaloMean.hpp"
#include "core/detail/joint_component/Numerics.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <set>
#include <stdexcept>
#include <unordered_map>
#include <utility>

namespace joint_test::dynamic_halo {
namespace {
bool IncludesCharge(const Mode mode) {return mode!=Mode::AB;}
const HaloOwner & FindHalo(const Ownership & ownership,std::size_t atom)
{
    const auto found=std::lower_bound(ownership.halos.begin(),ownership.halos.end(),atom,
        [](const HaloOwner & item,std::size_t value) {return item.halo<value;});
    if(found==ownership.halos.end() || found->halo!=atom)
        throw std::logic_error("Owned Dynamic Halo metadata is missing.");
    return *found;
}
std::vector<std::size_t> UniqueRows(const rhbm_gem::core::JointProblemInput & input,
    std::size_t atom)
{
    std::vector<std::size_t> rows;
    rows.reserve(input.support.at(atom).size());
    for(const auto & support:input.support[atom]) rows.push_back(support.row);
    std::sort(rows.begin(),rows.end());
    rows.erase(std::unique(rows.begin(),rows.end()),rows.end());
    return rows;
}
}

Ownership BuildOwnership(const rhbm_gem::core::JointProblemInput & input,std::vector<Block> blocks,
    OwnerSelection selection)
{
    if(!input.selection_domain || input.atom_ids.size()!=input.support.size() ||
        input.row_ids.size()!=input.observations.size())
        throw std::invalid_argument("Dynamic Halo ownership requires a complete selected problem.");

    std::vector<bool> is_target(input.atom_ids.size());
    for(const auto target:input.selection_domain->target_indices)
    {
        if(target>=input.atom_ids.size() || is_target[target])
            throw std::invalid_argument("Dynamic Halo ownership has an invalid target identity.");
        is_target[target]=true;
    }
    std::vector<std::size_t> target_memberships(input.atom_ids.size());
    std::set<std::string> block_ids;
    for(const auto & block:blocks)
    {
        if(block.id.empty() || !block_ids.insert(block.id).second)
            throw std::invalid_argument("Dynamic Halo blocks require unique stable identities.");
        for(const auto target:block.targets)
        {
            if(target>=is_target.size() || !is_target[target])
                throw std::invalid_argument("Dynamic Halo block contains a non-target atom.");
            ++target_memberships[target];
        }
    }
    for(std::size_t atom=0;atom<is_target.size();++atom)
        if(is_target[atom] && target_memberships[atom]!=1)
            throw std::invalid_argument("Every Dynamic Halo target must belong to exactly one block.");
    std::vector<std::size_t> halos;
    for(std::size_t atom=0;atom<is_target.size();++atom)
        if(!is_target[atom] && !input.support[atom].empty()) halos.push_back(atom);

    std::vector<std::vector<std::size_t>> block_rows(blocks.size());
    std::vector<std::vector<std::size_t>> row_blocks(input.observations.size());
    for(std::size_t block=0;block<blocks.size();++block)
    {
        auto & rows=block_rows[block];
        for(const auto target:blocks[block].targets)
        {
            if(target>=is_target.size() || !is_target[target])
                throw std::invalid_argument("Dynamic Halo block contains a non-target atom.");
            for(const auto & support:input.support[target])
            {
                if(support.row>=input.observations.size())
                    throw std::invalid_argument("Dynamic Halo target support row is out of range.");
                rows.push_back(support.row);
            }
        }
        std::sort(rows.begin(),rows.end());
        rows.erase(std::unique(rows.begin(),rows.end()),rows.end());
        for(const auto row:rows) row_blocks[row].push_back(block);
    }

    Ownership result;
    result.blocks=std::move(blocks);
    result.block_halos.resize(result.blocks.size());
    result.block_dynamic_rows=std::move(block_rows);
    for(const auto & memberships:row_blocks) result.row_block_entries+=memberships.size();

    result.halos.reserve(halos.size());
    for(const auto halo:halos)
    {
        const auto halo_rows=UniqueRows(input,halo);
        std::vector<std::size_t> overlap(result.blocks.size());
        std::vector<std::size_t> intersecting;
        for(const auto row:halo_rows)
        {
            if(row>=row_blocks.size())
                throw std::invalid_argument("Dynamic Halo support row is out of range.");
            for(const auto block:row_blocks[row])
            {
                ++overlap[block];
                intersecting.push_back(block);
            }
        }
        std::sort(intersecting.begin(),intersecting.end());
        intersecting.erase(std::unique(intersecting.begin(),intersecting.end()),intersecting.end());
        if(intersecting.empty())
            throw std::invalid_argument("Supported Dynamic Halo has no target block owner.");

        auto candidates=intersecting;
        std::sort(candidates.begin(),candidates.end(),[&](std::size_t left,std::size_t right) {
            if(overlap[left]!=overlap[right]) return overlap[left]>overlap[right];
            return result.blocks[left].id<result.blocks[right].id;
        });
        const auto owner_rank=selection==OwnerSelection::SecondBestOverlap && candidates.size()>1 ? 1u : 0u;
        const auto owner=candidates[owner_rank];

        std::vector<std::size_t> components;
        for(const auto block:intersecting) components.push_back(result.blocks[block].component);
        std::sort(components.begin(),components.end());
        components.erase(std::unique(components.begin(),components.end()),components.end());
        if(components.size()>1) ++result.cross_component_halos;

        std::vector<std::size_t> overlapping_targets;
        for(const auto block:intersecting)
            for(const auto target:result.blocks[block].targets)
            {
                const auto target_rows=UniqueRows(input,target);
                std::size_t i{},j{};
                while(i<target_rows.size() && j<halo_rows.size())
                {
                    if(target_rows[i]==halo_rows[j])
                    {overlapping_targets.push_back(target); break;}
                    if(target_rows[i]<halo_rows[j]) ++i; else ++j;
                }
            }
        std::sort(overlapping_targets.begin(),overlapping_targets.end());
        overlapping_targets.erase(std::unique(overlapping_targets.begin(),overlapping_targets.end()),
            overlapping_targets.end());

        std::vector<std::size_t> owner_mean_targets;
        for(const auto target:result.blocks[owner].targets)
            if(std::binary_search(overlapping_targets.begin(),overlapping_targets.end(),target))
                owner_mean_targets.push_back(target);
        std::sort(owner_mean_targets.begin(),owner_mean_targets.end());
        if(owner_mean_targets.empty())
            throw std::logic_error("Halo owner has no overlapping Target for an overlap mean.");

        result.halos.push_back({halo,owner,overlap[owner],std::move(intersecting),
            std::move(overlapping_targets),std::move(owner_mean_targets)});
        result.block_halos[owner].push_back(halo);
        auto & dynamic_rows=result.block_dynamic_rows[owner];
        dynamic_rows.insert(dynamic_rows.end(),halo_rows.begin(),halo_rows.end());
    }
    for(auto & rows:result.block_dynamic_rows)
    {
        std::sort(rows.begin(),rows.end());
        rows.erase(std::unique(rows.begin(),rows.end()),rows.end());
    }
    for(auto & owned:result.block_halos) std::sort(owned.begin(),owned.end());
    std::sort(result.halos.begin(),result.halos.end(),[](const auto & a,const auto & b) {
        return a.halo<b.halo;
    });
    return result;
}

AugmentedDesign BuildAugmentedDesign(const rhbm_gem::core::JointProblemInput & input,
    const Ownership & ownership,std::size_t block,std::span<const Parameters> target_parameters,Mode mode)
{
    if(block>=ownership.blocks.size() || target_parameters.size()!=input.atom_ids.size())
        throw std::invalid_argument("Invalid Dynamic Halo block design input.");
    const auto & targets=ownership.blocks[block].targets;
    const auto & halos=ownership.block_halos[block];
    AugmentedDesign result;
    result.rows=ownership.block_dynamic_rows.at(block);
    result.values=Eigen::MatrixXd::Zero(static_cast<Eigen::Index>(result.rows.size()),
        static_cast<Eigen::Index>(2*targets.size()));
    std::unordered_map<std::size_t,Eigen::Index> local_rows;
    local_rows.reserve(result.rows.size());
    for(std::size_t row=0;row<result.rows.size();++row)
        local_rows.emplace(result.rows[row],static_cast<Eigen::Index>(row));
    const auto add=[&](std::size_t atom,std::size_t target_slot,double width,double amplitude_share) {
        for(const auto & support:input.support.at(atom))
        {
            const auto row=local_rows.find(support.row);
            if(row==local_rows.end()) throw std::logic_error("Dynamic Halo block omitted a supported row.");
            const auto basis=rhbm_gem::core::joint_component::EvaluateKernel(support.squared_distance,width,2.5);
            result.values(row->second,static_cast<Eigen::Index>(2*target_slot))+=amplitude_share*basis.gaussian;
            if(mode==Mode::ABC || atom==targets[target_slot])
                result.values(row->second,static_cast<Eigen::Index>(2*target_slot+1))+=amplitude_share*basis.charge;
        }
    };
    for(std::size_t slot=0;slot<targets.size();++slot)
        add(targets[slot],slot,target_parameters[targets[slot]].b,1.0);
    if(mode!=Mode::StaticABC)
    {
        for(const auto halo:halos)
        {
            const auto & owner=FindHalo(ownership,halo);
            const auto & mean_targets=MeanTargets(ownership,owner);
            const auto mean=HaloMean(ownership,owner,target_parameters,IncludesCharge(mode));
            const double share=1.0/static_cast<double>(mean_targets.size());
            for(std::size_t slot=0;slot<targets.size();++slot)
                if(std::find(mean_targets.begin(),mean_targets.end(),targets[slot])!=mean_targets.end())
                    add(halo,slot,mean.b,share);
        }
    }
    return result;
}

Eigen::MatrixXd CoupledWidthJacobian(const rhbm_gem::core::JointProblemInput & input,
    const Ownership & ownership,std::size_t block,std::span<const Parameters> target_parameters,
    const Eigen::VectorXd & beta,Mode mode)
{
    if(block>=ownership.blocks.size() || target_parameters.size()!=input.atom_ids.size())
        throw std::invalid_argument("Invalid Dynamic Halo width Jacobian input.");
    const auto & targets=ownership.blocks[block].targets;
    if(beta.size()!=static_cast<Eigen::Index>(2*targets.size()))
        throw std::invalid_argument("Dynamic Halo width Jacobian coefficient dimension mismatch.");
    const auto & halos=ownership.block_halos[block];
    const auto & rows=ownership.block_dynamic_rows[block];
    Eigen::MatrixXd result=Eigen::MatrixXd::Zero(static_cast<Eigen::Index>(rows.size()),
        static_cast<Eigen::Index>(targets.size()));
    std::unordered_map<std::size_t,Eigen::Index> local_rows;
    local_rows.reserve(rows.size());
    for(std::size_t row=0;row<rows.size();++row)
        local_rows.emplace(rows[row],static_cast<Eigen::Index>(row));
    for(std::size_t slot=0;slot<targets.size();++slot)
    {
        const auto target=targets[slot];
        for(const auto & support:input.support.at(target))
        {
            const auto row=local_rows.at(support.row);
            const auto basis=rhbm_gem::core::joint_component::EvaluateKernel(support.squared_distance,
                target_parameters[target].b,2.5);
            result(row,static_cast<Eigen::Index>(slot))+=beta(static_cast<Eigen::Index>(2*slot))*basis.gaussian_log_width+
                beta(static_cast<Eigen::Index>(2*slot+1))*basis.charge_log_width;
        }
    }
    if(!halos.empty() && mode!=Mode::StaticABC)
    {
        for(const auto halo:halos)
        {
            const auto & owner=FindHalo(ownership,halo);
            const auto & mean_targets=MeanTargets(ownership,owner);
            const auto mean=HaloMean(ownership,owner,target_parameters,IncludesCharge(mode));
            const double share=1.0/static_cast<double>(mean_targets.size());
            for(std::size_t slot=0;slot<targets.size();++slot)
            {
                if(std::find(mean_targets.begin(),mean_targets.end(),targets[slot])==mean_targets.end()) continue;
                const double chain=target_parameters[targets[slot]].b*share/mean.b;
                const double halo_a=mean.a;
                const double halo_c=mode==Mode::AB ? 0.0 : mean.c;
                for(const auto & support:input.support.at(halo))
                {
                    const auto row=local_rows.at(support.row);
                    const auto basis=rhbm_gem::core::joint_component::EvaluateKernel(support.squared_distance,mean.b,2.5);
                    result(row,static_cast<Eigen::Index>(slot))+=chain*(halo_a*basis.gaussian_log_width+
                        halo_c*basis.charge_log_width);
                }
            }
        }
    }
    return result;
}

std::vector<Parameters> DeriveHalos(const rhbm_gem::core::JointProblemInput & input,
    const Ownership & ownership,std::span<const Parameters> target_parameters,Mode mode)
{
    if(target_parameters.size()!=input.atom_ids.size())
        throw std::invalid_argument("Dynamic Halo state does not cover every contributor.");
    auto result=std::vector<Parameters>(target_parameters.begin(),target_parameters.end());
    for(const auto & owner:ownership.halos)
    {
        result[owner.halo]=HaloMean(ownership,owner,target_parameters,IncludesCharge(mode));
    }
    return result;
}

std::vector<double> ReplayPrediction(const rhbm_gem::core::JointProblemInput & input,
    std::span<const Parameters> parameters)
{
    if(parameters.size()!=input.atom_ids.size())
        throw std::invalid_argument("Dynamic Halo replay state does not cover every contributor.");
    std::vector<double> result(input.observations.size());
    for(std::size_t atom=0;atom<parameters.size();++atom)
    {
        const auto & parameter=parameters[atom];
        for(const auto & support:input.support[atom])
        {
            const auto basis=rhbm_gem::core::joint_component::EvaluateKernel(support.squared_distance,parameter.b,2.5);
            result[support.row]+=parameter.a*basis.gaussian+parameter.c*basis.charge;
        }
    }
    return result;
}

FullState BuildState(const rhbm_gem::core::JointProblemInput & input,const Ownership & ownership,
    std::span<const Parameters> targets,Mode mode,double scale)
{
    if(!(scale>0.0) || !std::isfinite(scale)) throw std::invalid_argument("Invalid Dynamic Halo observation scale.");
    FullState result;
    result.contributors=DeriveHalos(input,ownership,targets,mode);
    result.prediction=ReplayPrediction(input,result.contributors);
    long double squared{};
    for(std::size_t row=0;row<input.observations.size();++row)
    {
        const long double residual=static_cast<long double>(result.prediction[row])-input.observations[row];
        squared+=residual*residual;
    }
    result.objective=static_cast<double>(0.5L*squared/(static_cast<long double>(scale)*scale));
    if(!std::isfinite(result.objective)) throw std::runtime_error("Dynamic Halo replay objective is nonfinite.");
    return result;
}

bool CommitMonotone(FullState & current,FullState && candidate)
{
    if(candidate.contributors.size()!=current.contributors.size() ||
        candidate.prediction.size()!=current.prediction.size() || !std::isfinite(candidate.objective))
        throw std::invalid_argument("Dynamic Halo transaction state dimensions do not match.");
    const double enclosure=1e-12+2e-12*std::abs(current.objective);
    if(candidate.objective>current.objective+enclosure) return false;
    current=std::move(candidate);
    return true;
}

std::optional<std::size_t> DenseBytes(std::size_t rows,std::size_t columns)
{
    constexpr auto bytes_per_value=sizeof(double);
    if(columns!=0 && rows>std::numeric_limits<std::size_t>::max()/columns) return std::nullopt;
    const auto values=rows*columns;
    if(values>std::numeric_limits<std::size_t>::max()/bytes_per_value) return std::nullopt;
    return values*bytes_per_value;
}

Parameters Mean(std::span<const Parameters> values, bool include_charge)
{
    if(values.empty()) throw std::invalid_argument("Cannot compute an empty Target mean.");
    Parameters result{0.0,0.0,0.0};
    for(const auto & value:values)
    {
        if(!std::isfinite(value.a) || value.a<0 || !std::isfinite(value.b) || !(value.b>0) ||
            !std::isfinite(value.c))
            throw std::invalid_argument("Target mean contains an invalid parameter.");
        result.a+=value.a;
        result.b+=value.b;
        if(include_charge) result.c+=value.c;
    }
    const double count=static_cast<double>(values.size());
    result.a/=count; result.b/=count; result.c=include_charge ? result.c/count : 0.0;
    return result;
}

const std::vector<std::size_t> & MeanTargets(const Ownership & ownership,const HaloOwner & halo)
{
    return ownership.mean_scope==MeanScope::WholeCore ?
        ownership.blocks.at(halo.owner_block).targets : halo.owner_mean_targets;
}

Parameters HaloMean(const Ownership & ownership,const HaloOwner & halo,
    std::span<const Parameters> target_parameters,bool include_charge)
{
    if(target_parameters.empty()) throw std::invalid_argument("Dynamic Halo mean has no Target state.");
    const auto & targets=MeanTargets(ownership,halo);
    std::vector<Parameters> members;
    members.reserve(targets.size());
    for(const auto target:targets)
    {
        if(target>=target_parameters.size()) throw std::invalid_argument("Dynamic Halo mean Target is out of range.");
        members.push_back(target_parameters[target]);
    }
    return Mean(members,include_charge);
}

double MeanWidthLogDerivative(double target_width, std::size_t core_size)
{
    if(!std::isfinite(target_width) || !(target_width>0) || core_size==0)
        throw std::invalid_argument("Invalid Dynamic Halo width derivative input.");
    return target_width/static_cast<double>(core_size);
}
}
