#pragma once
#include <rhbm_gem/core/JointComponentEstimator.hpp>
#include <Eigen/Core>
#include <algorithm>
#include <limits>
#include <numeric>
#include <span>
#include <stdexcept>

namespace rhbm_gem::core::joint_component {
using Vector=Eigen::VectorXd;
using VectorRef=Eigen::Ref<const Vector>;
using VectorMap=Eigen::Map<const Vector>;
using Indices=std::vector<Eigen::Index>;
inline bool SafeAdd(std::size_t a,std::size_t b,std::size_t & out)
{if(b>std::numeric_limits<std::size_t>::max()-a) return false; out=a+b; return true;}
inline bool SafeMultiply(std::size_t a,std::size_t b,std::size_t & out)
{if(a && b>std::numeric_limits<std::size_t>::max()/a) return false; out=a*b; return true;}
struct CompactRowMapping
{
    Indices parent_rows,local_rows;
    static std::shared_ptr<const CompactRowMapping> FromParentRows(const Indices & rows)
    {
        auto mapping=std::make_shared<CompactRowMapping>();
        std::vector<std::pair<Eigen::Index,Eigen::Index>> pairs; pairs.reserve(rows.size());
        for(std::size_t k=0;k<rows.size();++k) pairs.emplace_back(rows[k],static_cast<Eigen::Index>(k));
        std::sort(pairs.begin(),pairs.end(),[](const auto & a,const auto & b){return a.first<b.first;});
        for(std::size_t k=0;k<pairs.size();++k)
        {
            if(k && pairs[k-1].first==pairs[k].first) throw std::invalid_argument("Duplicate row mapping key.");
            mapping->parent_rows.push_back(pairs[k].first); mapping->local_rows.push_back(pairs[k].second);
        }
        return mapping;
    }
    Eigen::Index Find(Eigen::Index row) const
    {
        const auto found=std::lower_bound(parent_rows.begin(),parent_rows.end(),row);
        return found==parent_rows.end() || *found!=row ? -1 : local_rows[static_cast<std::size_t>(found-parent_rows.begin())];
    }
    std::size_t StorageBytes() const
    {
        std::size_t parent{},local{},total{};
        if(!SafeMultiply(parent_rows.capacity(),sizeof(Eigen::Index),parent) ||
            !SafeMultiply(local_rows.capacity(),sizeof(Eigen::Index),local) || !SafeAdd(parent,local,total))
            return std::numeric_limits<std::size_t>::max();
        return total;
    }
};
struct LocalSupportIndex {std::size_t support{}; Eigen::Index row{};};

// A context may outlive its problem handle. The map's deleter retains its owner.
inline std::shared_ptr<const VectorMap> Observe(std::shared_ptr<const JointProblemInput> input)
{
    auto * map=new VectorMap(input->observations.data(),static_cast<Eigen::Index>(input->observations.size()));
    return {map,[input=std::move(input)](const VectorMap * p){delete p;}};
}
class Identities
{
    std::shared_ptr<const std::vector<std::string>> values;
    std::shared_ptr<const Indices> indices;
public:
    Identities():values(std::make_shared<const std::vector<std::string>>()) {}
    Identities(std::vector<std::string> v):values(std::make_shared<const std::vector<std::string>>(std::move(v))) {}
    Identities(std::initializer_list<std::string> v):Identities(std::vector<std::string>(v)) {}
    explicit Identities(std::shared_ptr<const std::vector<std::string>> v):values(std::move(v)) {}
    std::size_t size() const {return indices ? indices->size() : values->size();}
    const std::string & at(std::size_t i) const {return values->at(indices ? static_cast<std::size_t>(indices->at(i)) : i);}
    const std::string & operator[](std::size_t i) const {return at(i);}
    Identities Select(const Indices & selected) const
    {
        Identities out(*this); auto map=std::make_shared<Indices>(); map->reserve(selected.size());
        for(auto i:selected) map->push_back(indices ? indices->at(static_cast<std::size_t>(i)) : i);
        out.indices=std::move(map); return out;
    }
    struct Iterator
    {
        const Identities * owner; std::size_t index;
        const std::string & operator*() const {return (*owner)[index];}
        Iterator & operator++() {++index; return *this;}
        bool operator!=(const Iterator & other) const {return index!=other.index;}
    };
    Iterator begin() const {return {this,0};} Iterator end() const {return {this,size()};}
    operator std::vector<std::string>() const {std::vector<std::string> out; for(const auto & id:*this) out.push_back(id); return out;}
    bool operator==(const Identities & other) const
    {if(size()!=other.size()) return false; for(std::size_t i=0;i<size();++i) if(at(i)!=other.at(i)) return false; return true;}
};
struct Support {Eigen::Index row; double square;};
// Component support is an indexed view of the original memberships, including
// zero-valued numerical basis entries. No local squared distances are copied.
class DomainAtoms
{
    std::shared_ptr<const JointProblemInput> input;
    std::shared_ptr<const DomainAtoms> parent;
    Indices selected;
    std::shared_ptr<const CompactRowMapping> row_to_local;
    std::shared_ptr<const std::vector<std::vector<LocalSupportIndex>>> local_supports;
    std::size_t RootAtom(std::size_t a) const {return parent ? parent->RootAtom(static_cast<std::size_t>(selected.at(a))) : a;}
    Eigen::Index LocalRow(Eigen::Index row) const
    {
        if(!parent) return row;
        const auto previous=parent->LocalRow(row);
        return previous<0 ? -1 : row_to_local->Find(previous);
    }
public:
    explicit DomainAtoms(std::shared_ptr<const JointProblemInput> snapshot):input(std::move(snapshot)) {}
    DomainAtoms Select(const Indices & atoms,std::shared_ptr<const CompactRowMapping> rows) const
    {
        DomainAtoms out(input); out.parent=std::make_shared<const DomainAtoms>(*this); out.selected=atoms; out.row_to_local=std::move(rows);
        auto indices=std::make_shared<std::vector<std::vector<LocalSupportIndex>>>(atoms.size());
        for(std::size_t a=0;a<atoms.size();++a)
        {
            const auto & support=input->support.at(out.RootAtom(a));
            for(std::size_t k=0;k<support.size();++k)
            {
                const auto local=out.LocalRow(static_cast<Eigen::Index>(support[k].row));
                if(local>=0) (*indices)[a].push_back({k,local});
            }
        }
        out.local_supports=std::move(indices);
        return out;
    }
    std::size_t size() const {return parent ? selected.size() : input->support.size();}
    struct Rows
    {
        const DomainAtoms * owner; std::span<const JointSupport> entries;
        const std::vector<LocalSupportIndex> * local{};
        std::size_t size() const {return local ? local->size() : entries.size();} bool empty() const {return size()==0;}
        Support operator[](std::size_t i) const
        {const auto & s=entries[local ? local->at(i).support : i]; return {local ? local->at(i).row : owner->LocalRow(static_cast<Eigen::Index>(s.row)),s.squared_distance};}
        Eigen::Index ParentRow(std::size_t i) const
        {return static_cast<Eigen::Index>(entries[local ? local->at(i).support : i].row);}
        struct Iterator
        {
            const Rows * rows; std::size_t index;
            Support operator*() const {return (*rows)[index];}
            Iterator & operator++() {++index; return *this;}
            bool operator!=(const Iterator & other) const {return index!=other.index;}
        };
        Iterator begin() const {return {this,0};} Iterator end() const {return {this,size()};}
    };
    Rows operator[](std::size_t a) const {return {this,input->support.at(RootAtom(a)),local_supports ? &local_supports->at(a) : nullptr};}
    struct Iterator
    {
        const DomainAtoms * owner; std::size_t index;
        Rows operator*() const {return (*owner)[index];}
        Iterator & operator++() {++index; return *this;}
        bool operator!=(const Iterator & other) const {return index!=other.index;}
    };
    Iterator begin() const {return {this,0};} Iterator end() const {return {this,size()};}
    std::shared_ptr<const JointProblemInput> Snapshot() const {return input;}
    std::size_t ParentAtom(std::size_t a) const {return RootAtom(a);}
    const JointProblemInput * Storage() const {return input.get();}
    std::size_t LocalSupportIndexEntries() const
    {
        std::size_t count{};
        if(local_supports) for(const auto & atom:*local_supports)
            if(!SafeAdd(count,atom.capacity(),count)) return std::numeric_limits<std::size_t>::max();
        return count;
    }
    std::size_t LocalSupportIndexBytes() const
    {
        std::size_t bytes{};
        if(local_supports && !SafeMultiply(local_supports->capacity(),sizeof(std::vector<LocalSupportIndex>),bytes))
            return std::numeric_limits<std::size_t>::max();
        if(local_supports) for(const auto & atom:*local_supports)
        {
            std::size_t item{};
            if(!SafeMultiply(atom.capacity(),sizeof(LocalSupportIndex),item) || !SafeAdd(bytes,item,bytes))
                return std::numeric_limits<std::size_t>::max();
        }
        return bytes;
    }
};
struct Domain
{
    Eigen::Index rows;
    DomainAtoms atoms;
    explicit Domain(std::shared_ptr<const JointProblemInput> input):rows(static_cast<Eigen::Index>(input->observations.size())),atoms(std::move(input)) {}
    Domain(Eigen::Index count,std::vector<std::vector<Support>> support):Domain([&] {
        auto input=std::make_shared<JointProblemInput>(); input->support.resize(support.size());
        for(std::size_t a=0;a<support.size();++a) for(const auto & s:support[a])
            input->support[a].push_back({static_cast<std::size_t>(s.row),s.square});
        return input;
    }()) {rows=count;}
    Domain Select(const Indices & selected,Eigen::Index count,std::shared_ptr<const CompactRowMapping> mapping) const
    {auto out=*this; out.rows=count; out.atoms=atoms.Select(selected,std::move(mapping)); return out;}
    std::vector<std::vector<Support>> CopySupport() const
    {
        std::vector<std::vector<Support>> out(atoms.size());
        for(std::size_t a=0;a<atoms.size();++a) for(const auto & s:atoms[a]) out[a].push_back(s);
        return out;
    }
};
struct PartitionMappings {Indices atom_component,row_component,atom_to_local,row_to_local;};
}
