#include "DynamicHaloMean.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <set>
#include <stdexcept>
#include <utility>

namespace joint_test::dynamic_halo {
namespace {
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

Ownership BuildOwnership(const rhbm_gem::core::JointProblemInput & input, std::vector<Block> blocks)
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

        auto owner=intersecting.front();
        for(const auto block:intersecting)
            if(overlap[block]>overlap[owner] || (overlap[block]==overlap[owner] &&
                result.blocks[block].id<result.blocks[owner].id)) owner=block;

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

        result.halos.push_back({halo,owner,overlap[owner],std::move(intersecting),
            std::move(overlapping_targets)});
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

double MeanWidthLogDerivative(double target_width, std::size_t core_size)
{
    if(!std::isfinite(target_width) || !(target_width>0) || core_size==0)
        throw std::invalid_argument("Invalid Dynamic Halo width derivative input.");
    return target_width/static_cast<double>(core_size);
}
}
