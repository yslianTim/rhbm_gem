#include "Problem.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>

namespace rhbm_gem::core::joint_component {
namespace {
std::size_t Index(Eigen::Index i) {return static_cast<std::size_t>(i);}
}
ComponentPartition Partition(const Domain & domain,const Identities & ids)
{
    const std::vector<std::string> identities=ids;
    const auto count=static_cast<Eigen::Index>(domain.atoms.size());
    if(domain.rows<0 || ids.size()!=domain.atoms.size() || std::set<std::string>(identities.begin(),identities.end()).size()!=ids.size())
        throw std::invalid_argument("Invalid structural partition identities.");
    std::vector<Eigen::Index> parent(static_cast<std::size_t>(count)); std::iota(parent.begin(),parent.end(),0);
    auto root=[&](Eigen::Index a) {while(parent[Index(a)]!=a) {parent[Index(a)]=parent[Index(parent[Index(a)])]; a=parent[Index(a)];} return a;};
    std::vector<std::vector<Eigen::Index>> contributors(static_cast<std::size_t>(domain.rows));
    ComponentPartition out; auto mappings=std::make_shared<PartitionMappings>(); out.mappings=mappings;
    mappings->atom_component.assign(ids.size(),-1); mappings->row_component.assign(static_cast<std::size_t>(domain.rows),-1);
    mappings->atom_to_local.assign(ids.size(),-1); mappings->row_to_local.assign(static_cast<std::size_t>(domain.rows),-1);
    for(Eigen::Index a=0;a<count;++a)
    {
        std::set<Eigen::Index> seen;
        for(const auto & s:domain.atoms[Index(a)])
        {
            if(s.row<0 || s.row>=domain.rows || !std::isfinite(s.square) || s.square<0 || s.square>6.25 || !seen.insert(s.row).second)
                throw std::invalid_argument("Invalid structural membership.");
            contributors[Index(s.row)].push_back(a);
        }
        if(domain.atoms[Index(a)].empty()) out.unobserved_atoms.push_back(a);
    }
    for(const auto & row:contributors) for(std::size_t k=1;k<row.size();++k) parent[Index(root(row[k]))]=root(row[0]);
    std::map<Eigen::Index,std::vector<Eigen::Index>> groups;
    for(Eigen::Index a=0;a<count;++a) groups[root(a)].push_back(a);
    for(const auto & [key,atoms]:groups)
    {
        (void)key; ComponentView v; v.atoms=atoms;
        v.id=ids[Index(*std::min_element(atoms.begin(),atoms.end(),[&](auto a,auto b){return ids[Index(a)]<ids[Index(b)];}))];
        out.components.push_back(std::move(v));
    }
    std::sort(out.components.begin(),out.components.end(),[](const auto & a,const auto & b){return a.id<b.id;});
    for(std::size_t c=0;c<out.components.size();++c) for(auto a:out.components[c].atoms) mappings->atom_component[Index(a)]=static_cast<Eigen::Index>(c);
    for(Eigen::Index r=0;r<domain.rows;++r)
    {
        if(contributors[Index(r)].empty()) {out.constant_rows.push_back(r); continue;}
        const auto c=mappings->atom_component[Index(contributors[Index(r)][0])]; mappings->row_component[Index(r)]=c; out.components[Index(c)].rows.push_back(r);
    }
    for(std::size_t c=0;c<out.components.size();++c)
    {
        auto & v=out.components[c]; v.mappings=mappings; v.component_index=static_cast<Eigen::Index>(c);
        for(std::size_t r=0;r<v.rows.size();++r) mappings->row_to_local[Index(v.rows[r])]=static_cast<Eigen::Index>(r);
        for(std::size_t a=0;a<v.atoms.size();++a) mappings->atom_to_local[Index(v.atoms[a])]=static_cast<Eigen::Index>(a);
        const auto rows=CompactRowMapping::FromParentRows(v.rows);
        v.domain=domain.Select(v.atoms,static_cast<Eigen::Index>(v.rows.size()),rows);
    }
    return out;
}
ComponentProblemView MakeComponentProblemView(const JointParameterLayout & parent_layout,const ComponentView & view)
{
    const auto parent=view.domain.atoms.Snapshot();
    ComponentProblemView out; out.parent_row_count=static_cast<Eigen::Index>(view.rows.size());
    Indices view_to_local_rows;
    Indices local_to_view_atoms;
    Indices identity_parent_atoms;
    out.local_to_parent_atoms.reserve(parent_layout.full_atoms.size());
    out.local_to_parent_rows.reserve(parent_layout.informative_rows.size());
    local_to_view_atoms.reserve(parent_layout.full_atoms.size());
    identity_parent_atoms.reserve(parent_layout.full_atoms.size());
    for(std::size_t k=0;k<parent_layout.informative_rows.size();++k)
    {
        const auto row=parent_layout.informative_rows[k];
        const auto component_row=view.LocalRow(static_cast<Eigen::Index>(row));
        if(component_row<0) throw std::invalid_argument("Component row is outside its view.");
        view_to_local_rows.push_back(component_row);
        out.local_to_parent_rows.push_back(static_cast<Eigen::Index>(row));
        out.local_layout.informative_rows.push_back(k);
    }
    for(std::size_t k=0;k<parent_layout.full_atoms.size();++k)
    {
        const auto atom=parent_layout.full_atoms[k];
        const auto component_atom=view.LocalAtom(static_cast<Eigen::Index>(atom));
        if(component_atom<0) throw std::invalid_argument("Component atom is outside its view.");
        out.local_to_parent_atoms.push_back(static_cast<Eigen::Index>(atom));
        local_to_view_atoms.push_back(component_atom);
        identity_parent_atoms.push_back(static_cast<Eigen::Index>(view.domain.atoms.ParentAtom(static_cast<std::size_t>(component_atom))));
        out.local_layout.full_atoms.push_back(k);
    }
    const auto row_mapping=CompactRowMapping::FromParentRows(view_to_local_rows);
    out.domain=view.domain.Select(local_to_view_atoms,static_cast<Eigen::Index>(out.local_layout.informative_rows.size()),row_mapping);
    const Identities parent_atoms(std::shared_ptr<const std::vector<std::string>>(parent,&parent->atom_ids));
    const Identities parent_rows(std::shared_ptr<const std::vector<std::string>>(parent,&parent->row_ids));
    out.atom_ids=parent_atoms.Select(identity_parent_atoms);
    out.row_ids=parent_rows.Select(out.local_to_parent_rows);
    return out;
}
}
