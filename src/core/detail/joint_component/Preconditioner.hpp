#pragma once
#include "ProfileJacobianOperator.hpp"
#include <cmath>
#include <set>

namespace rhbm_gem::core::joint_component {
enum class PreconditionerSpace {Width,FreeAC};
struct PreconditionerBlock
{
    std::string id;
    Indices core_atoms,overlap_atoms,informative_rows; // Parent snapshot indices.
};
// Structural data only. Keep through shared_ptr<const ...>; numerical factors
// and active-face mappings belong to a solver-specific linearization instead.
struct PreconditionerPartition
{
    const std::shared_ptr<const JointProblemInput> problem;
    const JointParameterLayout layout;
    const std::vector<PreconditionerBlock> blocks;
    std::vector<std::vector<std::size_t>> atom_blocks; // parent atom -> memberships
    PreconditionerPartition(std::shared_ptr<const JointProblemInput> input,JointParameterLayout l,
        std::vector<PreconditionerBlock> b):problem(std::move(input)),layout(std::move(l)),blocks(std::move(b))
    {
        if(!problem) throw std::invalid_argument("Missing preconditioner snapshot");
        atom_blocks.resize(problem->atom_ids.size());
        const std::set<std::size_t> full(layout.full_atoms.begin(),layout.full_atoms.end());
        const std::set<std::size_t> rows(layout.informative_rows.begin(),layout.informative_rows.end());
        std::set<std::size_t> cores; std::set<std::string> names;
        for(std::size_t k=0;k<blocks.size();++k)
        {
            const auto & block=blocks[k]; std::set<Eigen::Index> members;
            if(block.id.empty() || !names.insert(block.id).second || block.core_atoms.empty())
                throw std::invalid_argument("Invalid preconditioner block identity");
            for(const auto * atoms:{&block.core_atoms,&block.overlap_atoms}) for(auto a:*atoms)
            {
                if(a<0 || static_cast<std::size_t>(a)>=atom_blocks.size() || !full.contains(static_cast<std::size_t>(a)) || !members.insert(a).second)
                    throw std::invalid_argument("Invalid preconditioner atom mapping");
                atom_blocks[static_cast<std::size_t>(a)].push_back(k);
                if(atoms==&block.core_atoms && !cores.insert(static_cast<std::size_t>(a)).second)
                    throw std::invalid_argument("Overlapping preconditioner cores");
            }
            std::set<Eigen::Index> unique_rows;
            for(auto r:block.informative_rows)
                if(r<0 || static_cast<std::size_t>(r)>=problem->observations.size() || !rows.contains(static_cast<std::size_t>(r)) || !unique_rows.insert(r).second)
                    throw std::invalid_argument("Invalid preconditioner row mapping");
        }
        if(cores!=full) throw std::invalid_argument("Uncovered preconditioner core atom");
    }
};
struct BlockCoordinates
{
    Indices global; // local -> solver coordinate; never parent atom indices.
    Vector weights;
    Eigen::Index dimension{};
    Vector Restrict(VectorRef v) const
    {
        if(v.size()!=dimension) throw std::invalid_argument("Invalid global RHS size");
        Vector out(static_cast<Eigen::Index>(global.size()));
        for(std::size_t k=0;k<global.size();++k) out(static_cast<Eigen::Index>(k))=weights(static_cast<Eigen::Index>(k))*v(global[k]);
        return out;
    }
    void Scatter(VectorRef local,Vector & global_vector) const
    {
        if(local.size()!=static_cast<Eigen::Index>(global.size()) || global_vector.size()!=dimension) throw std::invalid_argument("Invalid block RHS size");
        for(std::size_t k=0;k<global.size();++k) global_vector(global[k])+=weights(static_cast<Eigen::Index>(k))*local(static_cast<Eigen::Index>(k));
    }
};
struct SolverBlockMapping
{
    PreconditionerSpace space;
    Eigen::Index dimension;
    std::vector<BlockCoordinates> blocks;
    // Construction checks coverage and creates symmetric partition-of-unity
    // weights. It neither partitions the graph nor builds Schwarz factors.
    SolverBlockMapping(PreconditionerSpace s,Eigen::Index size,const std::vector<Indices> & coordinates):space(s),dimension(size)
    {
        if(size<=0) throw std::invalid_argument("Invalid preconditioner dimension");
        Vector membership=Vector::Zero(size);
        for(const auto & block:coordinates)
        {
            std::set<Eigen::Index> unique;
            if(block.empty()) throw std::invalid_argument("Empty preconditioner block");
            for(auto k:block)
            {
                if(k<0 || k>=size || !unique.insert(k).second) throw std::invalid_argument("Invalid block coordinate");
                membership(k)+=1;
            }
        }
        if((membership.array()==0).any()) throw std::invalid_argument("Uncovered preconditioner coordinate");
        for(const auto & block:coordinates)
        {
            BlockCoordinates out{block,Vector(static_cast<Eigen::Index>(block.size())),size};
            for(std::size_t k=0;k<block.size();++k) out.weights(static_cast<Eigen::Index>(k))=1/std::sqrt(membership(block[k]));
            blocks.push_back(std::move(out));
        }
    }
};
struct PreconditionerContext
{
    std::shared_ptr<const LinearizationIdentity> linearization;
    PreconditionerSpace space{PreconditionerSpace::Width};
    Vector metric; // Positive trust-region diagonal, in solver coordinates.
    double damping{};
    bool Valid() const
    {return linearization && metric.size()>0 && metric.allFinite() && (metric.array()>0).all() && std::isfinite(damping) && damping>=0;}
    bool Matches(const PreconditionerContext & other) const
    {
        return linearization==other.linearization && space==other.space && damping==other.damping &&
            metric.size()==other.metric.size() && (metric.array()==other.metric.array()).all();
    }
};
// Solver owns the operator and preconditioner separately. Numeric actions are
// immutable during a Krylov solve; a changed context requires a new instance.
class IdentityPreconditioner
{
    const PreconditionerContext context_;
public:
    explicit IdentityPreconditioner(PreconditionerContext context):context_(std::move(context))
    {if(!context_.Valid()) throw std::invalid_argument("Invalid preconditioner context");}
    Vector ApplyInverse(VectorRef rhs,const PreconditionerContext & context) const
    {
        if(!context_.Matches(context)) throw std::logic_error("Stale preconditioner context");
        if(rhs.size()!=context_.metric.size() || !rhs.allFinite()) throw std::invalid_argument("Invalid preconditioner RHS");
        return rhs;
    }
};
}
