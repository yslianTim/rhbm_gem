#pragma once
#include <rhbm_gem/core/JointComponentEstimator.hpp>
#include <Eigen/Core>
#include <algorithm>
#include <numeric>
#include <span>
#include <stdexcept>

namespace rhbm_gem::core::joint_component {
using Vector=Eigen::VectorXd;
using VectorRef=Eigen::Ref<const Vector>;
using VectorMap=Eigen::Map<const Vector>;
using Indices=std::vector<Eigen::Index>;

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
    std::shared_ptr<const Indices> row_to_local;
    std::size_t RootAtom(std::size_t a) const {return parent ? parent->RootAtom(static_cast<std::size_t>(selected.at(a))) : a;}
    Eigen::Index LocalRow(Eigen::Index row) const
    {return parent ? row_to_local->at(static_cast<std::size_t>(parent->LocalRow(row))) : row;}
public:
    explicit DomainAtoms(std::shared_ptr<const JointProblemInput> snapshot):input(std::move(snapshot)) {}
    DomainAtoms Select(const Indices & atoms,std::shared_ptr<const Indices> rows) const
    {DomainAtoms out(input); out.parent=std::make_shared<const DomainAtoms>(*this); out.selected=atoms; out.row_to_local=std::move(rows); return out;}
    std::size_t size() const {return parent ? selected.size() : input->support.size();}
    struct Rows
    {
        const DomainAtoms * owner; std::span<const JointSupport> entries;
        std::size_t size() const {return entries.size();} bool empty() const {return entries.empty();}
        Support operator[](std::size_t i) const
        {const auto & s=entries[i]; return {owner->LocalRow(static_cast<Eigen::Index>(s.row)),s.squared_distance};}
        struct Iterator
        {
            const Rows * rows; std::size_t index;
            Support operator*() const {return (*rows)[index];}
            Iterator & operator++() {++index; return *this;}
            bool operator!=(const Iterator & other) const {return index!=other.index;}
        };
        Iterator begin() const {return {this,0};} Iterator end() const {return {this,size()};}
    };
    Rows operator[](std::size_t a) const {return {this,input->support.at(RootAtom(a))};}
    struct Iterator
    {
        const DomainAtoms * owner; std::size_t index;
        Rows operator*() const {return (*owner)[index];}
        Iterator & operator++() {++index; return *this;}
        bool operator!=(const Iterator & other) const {return index!=other.index;}
    };
    Iterator begin() const {return {this,0};} Iterator end() const {return {this,size()};}
    const JointProblemInput * Storage() const {return input.get();}
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
    Domain Select(const Indices & selected,Eigen::Index count,std::shared_ptr<const Indices> mapping) const
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
