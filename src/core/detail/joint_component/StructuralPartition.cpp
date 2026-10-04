#include "StructuralPartition.hpp"
#include <algorithm>
#include <deque>
#include <numeric>

namespace rhbm_gem::core::joint_component {
StructuralTopology::StructuralTopology(const JointProblemInput & input,const JointParameterLayout & layout)
    :input_(input),informative_rows_(input.observations.size(),false),row_offsets_(input.observations.size()+1),
     atom_order_(layout.full_atoms),neighbor_marks_(input.atom_ids.size())
{
    for(auto row:layout.informative_rows) informative_rows_.at(row)=true;
    std::sort(atom_order_.begin(),atom_order_.end(),[&](auto a,auto b) {
        return input_.atom_ids.at(a)<input_.atom_ids.at(b);
    });
    for(auto atom:atom_order_) for(const auto & support:input_.support.at(atom))
        if(informative_rows_.at(support.row)) ++row_offsets_[support.row+1];
    std::partial_sum(row_offsets_.begin(),row_offsets_.end(),row_offsets_.begin());
    row_atoms_.resize(row_offsets_.back()); auto cursor=row_offsets_;
    for(auto atom:atom_order_) for(const auto & support:input_.support.at(atom))
        if(informative_rows_[support.row]) row_atoms_[cursor[support.row]++]=static_cast<Eigen::Index>(atom);
}
bool StructuralTopology::LessAtom(Eigen::Index a,Eigen::Index b) const
{
    return input_.atom_ids.at(static_cast<std::size_t>(a))<input_.atom_ids.at(static_cast<std::size_t>(b));
}
Indices StructuralTopology::Neighbors(std::size_t atom)
{
    Indices out; ++neighbor_generation_;
    for(const auto & support:input_.support.at(atom)) if(informative_rows_[support.row])
        for(auto k=row_offsets_[support.row];k<row_offsets_[support.row+1];++k)
        {
            const auto other=row_atoms_[k]; const auto index=static_cast<std::size_t>(other);
            if(neighbor_marks_[index]!=neighbor_generation_)
            {neighbor_marks_[index]=neighbor_generation_; out.push_back(other);}
        }
    std::sort(out.begin(),out.end(),[&](auto a,auto b){return LessAtom(a,b);});
    return out;
}
StructuralBlockPartition StructuralTopology::BuildCores(std::size_t core_atoms,bool include_metadata)
{
    if(core_atoms==0) throw std::invalid_argument("Invalid structural partition policy");
    const auto count=input_.atom_ids.size(); std::vector<bool> assigned(count,false);
    StructuralBlockPartition result;
    for(auto seed:atom_order_) if(!assigned[seed])
    {
        StructuralCore core; core.id=input_.atom_ids[seed];
        std::deque<Eigen::Index> queue{static_cast<Eigen::Index>(seed)};
        std::vector<bool> queued(count,false); queued[seed]=true;
        while(!queue.empty() && core.atoms.size()<core_atoms)
        {
            const auto atom=queue.front(); queue.pop_front(); const auto index=static_cast<std::size_t>(atom);
            if(assigned[index]) continue;
            assigned[index]=true; core.atoms.push_back(atom);
            for(auto other:Neighbors(index))
            {
                const auto other_index=static_cast<std::size_t>(other);
                if(!assigned[other_index] && !queued[other_index])
                {queue.push_back(other); queued[other_index]=true;}
            }
        }
        std::sort(core.atoms.begin(),core.atoms.end(),[&](auto a,auto b){return LessAtom(a,b);});
        result.cores.push_back(std::move(core));
    }
    if(include_metadata) for(auto & core:result.cores)
    {
        std::vector<bool> in_core(count,false);
        for(auto atom:core.atoms)
        {
            in_core[static_cast<std::size_t>(atom)]=true;
            for(const auto & support:input_.support[static_cast<std::size_t>(atom)])
                if(informative_rows_[support.row]) core.affected_rows.push_back(static_cast<Eigen::Index>(support.row));
        }
        std::sort(core.affected_rows.begin(),core.affected_rows.end());
        core.affected_rows.erase(std::unique(core.affected_rows.begin(),core.affected_rows.end()),core.affected_rows.end());
        for(auto atom:core.atoms) for(auto other:Neighbors(static_cast<std::size_t>(atom)))
            if(!in_core[static_cast<std::size_t>(other)]) core.neighbor_atoms.push_back(other);
        std::sort(core.neighbor_atoms.begin(),core.neighbor_atoms.end(),[&](auto a,auto b){return LessAtom(a,b);});
        core.neighbor_atoms.erase(std::unique(core.neighbor_atoms.begin(),core.neighbor_atoms.end()),core.neighbor_atoms.end());
    }
    return result;
}
Indices StructuralTopology::ContextAtoms(const Indices & core_atoms,std::size_t context_hops,std::size_t max_block_atoms)
{
    const auto count=input_.atom_ids.size(); std::vector<bool> included(count,false);
    for(auto atom:core_atoms) included.at(static_cast<std::size_t>(atom))=true;
    Indices context,frontier=core_atoms;
    for(std::size_t hop=0;hop<context_hops && !frontier.empty();++hop)
    {
        Indices next;
        for(auto atom:frontier) for(auto other:Neighbors(static_cast<std::size_t>(atom)))
        {
            const auto index=static_cast<std::size_t>(other);
            if(included[index]) continue;
            included[index]=true; context.push_back(other); next.push_back(other);
            if(core_atoms.size()+context.size()>max_block_atoms) throw std::runtime_error("structural-context-limit");
        }
        frontier=std::move(next);
    }
    std::sort(context.begin(),context.end(),[&](auto a,auto b){return LessAtom(a,b);});
    return context;
}
StructuralBlockPartition BuildStructuralBlockPartition(const JointProblemInput & input,const JointParameterLayout & layout,
    std::size_t core_atoms,std::size_t context_hops,std::size_t max_block_atoms)
{
    if(core_atoms>max_block_atoms) throw std::invalid_argument("Invalid structural partition policy");
    StructuralTopology topology(input,layout); auto partition=topology.BuildCores(core_atoms);
    for(auto & core:partition.cores) core.context_atoms=topology.ContextAtoms(core.atoms,context_hops,max_block_atoms);
    return partition;
}
}
