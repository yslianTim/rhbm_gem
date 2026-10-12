#include "StructuralPartition.hpp"
#include <algorithm>
#include <deque>
#include <numeric>

namespace rhbm_gem::core::joint_component {
StructuralTopology::StructuralTopology(const Domain & domain,const Identities & identities,
    const JointParameterLayout & layout)
    :domain_(domain),identities_(identities),informative_rows_(static_cast<std::size_t>(domain.rows),false),
     row_offsets_(static_cast<std::size_t>(domain.rows)+1),atom_order_(layout.full_atoms),
     neighbor_marks_(domain.atoms.size())
{
    if(domain.rows<0 || identities.size()!=domain.atoms.size())
        throw std::invalid_argument("Invalid structural partition view.");
    for(auto row:layout.informative_rows) informative_rows_.at(row)=true;
    std::sort(atom_order_.begin(),atom_order_.end(),[&](auto a,auto b) {
        return identities_.at(a)<identities_.at(b);
    });
    for(auto atom:atom_order_) for(const auto & support:domain_.atoms[atom])
    {
        const auto row=static_cast<std::size_t>(support.row);
        if(informative_rows_.at(row)) ++row_offsets_[row+1];
    }
    std::partial_sum(row_offsets_.begin(),row_offsets_.end(),row_offsets_.begin());
    row_atoms_.resize(row_offsets_.back()); auto cursor=row_offsets_;
    for(auto atom:atom_order_) for(const auto & support:domain_.atoms[atom])
    {
        const auto row=static_cast<std::size_t>(support.row);
        if(informative_rows_[row]) row_atoms_[cursor[row]++]=static_cast<Eigen::Index>(atom);
    }
}
bool StructuralTopology::LessAtom(Eigen::Index a,Eigen::Index b) const
{
    return identities_.at(static_cast<std::size_t>(a))<identities_.at(static_cast<std::size_t>(b));
}
Indices StructuralTopology::Neighbors(std::size_t atom)
{
    Indices out; ++neighbor_generation_;
    for(const auto & support:domain_.atoms[atom])
    {
        const auto row=static_cast<std::size_t>(support.row);
        if(!informative_rows_[row]) continue;
        for(auto k=row_offsets_[row];k<row_offsets_[row+1];++k)
        {
            const auto other=row_atoms_[k]; const auto index=static_cast<std::size_t>(other);
            if(neighbor_marks_[index]!=neighbor_generation_)
            {neighbor_marks_[index]=neighbor_generation_; out.push_back(other);}
        }
    }
    std::sort(out.begin(),out.end(),[&](auto a,auto b){return LessAtom(a,b);});
    return out;
}
StructuralBlockPartition StructuralTopology::BuildCores(std::size_t core_atoms,bool include_metadata)
{
    if(core_atoms==0) throw std::invalid_argument("Invalid structural partition policy");
    const auto count=domain_.atoms.size(); std::vector<bool> assigned(count,false);
    StructuralBlockPartition result;
    for(auto seed:atom_order_) if(!assigned[seed])
    {
        StructuralCore core; core.id=identities_[seed];
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
    if(include_metadata) PopulateMetadata(result);
    return result;
}
StructuralBlockCensus StructuralTopology::Census(const StructuralBlockPartition & partition) const
{
    StructuralBlockCensus census;
    std::vector<std::size_t> row_marks(static_cast<std::size_t>(domain_.rows));
    std::vector<std::size_t> neighbor_marks(domain_.atoms.size());
    std::vector<bool> in_core(domain_.atoms.size(),false);
    std::size_t generation{};
    for(const auto & core:partition.cores)
    {
        if(generation==std::numeric_limits<std::size_t>::max()) throw std::overflow_error("Structural census generation overflow.");
        ++generation;
        std::size_t block_rows{},block_support{},block_neighbors{};
        for(auto atom:core.atoms)
        {
            const auto index=static_cast<std::size_t>(atom);
            in_core[index]=true;
            if(!SafeAdd(block_support,domain_.atoms[index].size(),block_support))
                throw std::overflow_error("Structural support census overflow.");
            for(const auto & support:domain_.atoms[static_cast<std::size_t>(atom)])
            {
                const auto row=static_cast<std::size_t>(support.row);
                if(informative_rows_[row] && row_marks[row]!=generation)
                {
                    row_marks[row]=generation;
                    if(!SafeAdd(block_rows,1,block_rows)) throw std::overflow_error("Structural row census overflow.");
                }
            }
        }
        for(auto atom:core.atoms)
        {
            const auto index=static_cast<std::size_t>(atom);
            for(const auto & support:domain_.atoms[index])
            {
                const auto row=static_cast<std::size_t>(support.row);
                if(!informative_rows_[row]) continue;
                for(auto k=row_offsets_[row];k<row_offsets_[row+1];++k)
                {
                    const auto other=static_cast<std::size_t>(row_atoms_[k]);
                    if(!in_core[other] && neighbor_marks[other]!=generation)
                    {
                        neighbor_marks[other]=generation;
                        if(!SafeAdd(block_neighbors,1,block_neighbors))
                            throw std::overflow_error("Structural neighbor census overflow.");
                    }
                }
            }
            in_core[index]=false;
        }
        if(!SafeAdd(census.affected_row_entries,block_rows,census.affected_row_entries))
            throw std::overflow_error("Structural row census overflow.");
        if(!SafeAdd(census.neighbor_atom_entries,block_neighbors,census.neighbor_atom_entries) ||
            !SafeAdd(census.total_block_support_entries,block_support,census.total_block_support_entries))
            throw std::overflow_error("Structural resource census overflow.");
        census.maximum_block_rows=std::max(census.maximum_block_rows,block_rows);
        census.maximum_block_support_entries=std::max(census.maximum_block_support_entries,block_support);
        census.maximum_neighbor_atoms=std::max(census.maximum_neighbor_atoms,block_neighbors);
    }
    return census;
}
std::size_t StructuralTopology::StorageBytes() const
{
    std::size_t total{};
    const auto add=[&](std::size_t count,std::size_t width) {
        std::size_t bytes{};
        return SafeMultiply(count,width,bytes) && SafeAdd(total,bytes,total);
    };
    if(!add(informative_rows_.capacity()/8+(informative_rows_.capacity()%8!=0),sizeof(std::uint8_t)) ||
        !add(row_offsets_.capacity(),sizeof(std::size_t)) || !add(row_atoms_.capacity(),sizeof(Eigen::Index)) ||
        !add(atom_order_.capacity(),sizeof(std::size_t)) || !add(neighbor_marks_.capacity(),sizeof(std::size_t)))
        return std::numeric_limits<std::size_t>::max();
    return total;
}
void StructuralTopology::PopulateMetadata(StructuralBlockPartition & partition)
{
    const auto count=domain_.atoms.size();
    for(auto & core:partition.cores)
    {
        std::vector<bool> in_core(count,false);
        for(auto atom:core.atoms)
        {
            in_core[static_cast<std::size_t>(atom)]=true;
            for(const auto & support:domain_.atoms[static_cast<std::size_t>(atom)])
                if(informative_rows_[static_cast<std::size_t>(support.row)]) core.affected_rows.push_back(support.row);
        }
        std::sort(core.affected_rows.begin(),core.affected_rows.end());
        core.affected_rows.erase(std::unique(core.affected_rows.begin(),core.affected_rows.end()),core.affected_rows.end());
        for(auto atom:core.atoms) for(auto other:Neighbors(static_cast<std::size_t>(atom)))
            if(!in_core[static_cast<std::size_t>(other)]) core.neighbor_atoms.push_back(other);
        std::sort(core.neighbor_atoms.begin(),core.neighbor_atoms.end(),[&](auto a,auto b){return LessAtom(a,b);});
        core.neighbor_atoms.erase(std::unique(core.neighbor_atoms.begin(),core.neighbor_atoms.end()),core.neighbor_atoms.end());
    }
}
Indices StructuralTopology::ContextAtoms(const Indices & core_atoms,std::size_t context_hops,std::size_t max_block_atoms)
{
    const auto count=domain_.atoms.size(); std::vector<bool> included(count,false);
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
    const auto owner=std::shared_ptr<const JointProblemInput>(&input,[](const JointProblemInput *){});
    const Domain domain(owner);
    const Identities identities(std::shared_ptr<const std::vector<std::string>>(owner,&input.atom_ids));
    return BuildStructuralBlockPartition(domain,identities,layout,core_atoms,context_hops,max_block_atoms);
}
StructuralBlockPartition BuildStructuralBlockPartition(const Domain & domain,const Identities & identities,
    const JointParameterLayout & layout,std::size_t core_atoms,std::size_t context_hops,std::size_t max_block_atoms)
{
    if(core_atoms>max_block_atoms) throw std::invalid_argument("Invalid structural partition policy");
    StructuralTopology topology(domain,identities,layout); auto partition=topology.BuildCores(core_atoms);
    for(auto & core:partition.cores) core.context_atoms=topology.ContextAtoms(core.atoms,context_hops,max_block_atoms);
    return partition;
}
}
