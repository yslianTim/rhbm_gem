#pragma once
#include "SnapshotViews.hpp"
#include <rhbm_gem/core/JointComponentEstimator.hpp>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace rhbm_gem::core::joint_component {
struct StructuralCore
{
    std::string id;
    Indices atoms,affected_rows,neighbor_atoms,context_atoms; // Parent snapshot indices.
};
struct StructuralBlockPartition
{
    std::vector<StructuralCore> cores;
};

class StructuralTopology
{
    const JointProblemInput & input_;
    std::vector<bool> informative_rows_;
    std::vector<std::size_t> row_offsets_;
    Indices row_atoms_;
    std::vector<std::size_t> atom_order_,neighbor_marks_;
    std::size_t neighbor_generation_{};
    bool LessAtom(Eigen::Index,Eigen::Index) const;
    Indices Neighbors(std::size_t);
public:
    StructuralTopology(const JointProblemInput &,const JointParameterLayout &);
    StructuralBlockPartition BuildCores(std::size_t core_atoms,bool include_metadata=true);
    Indices ContextAtoms(const Indices & core_atoms,std::size_t context_hops,std::size_t max_block_atoms);
};

// Builds deterministic cores from informative-row incidence. Context atoms
// are the requested graph-hop expansion and remain structural data only.
StructuralBlockPartition BuildStructuralBlockPartition(const JointProblemInput &,const JointParameterLayout &,
    std::size_t core_atoms,std::size_t context_hops=0,
    std::size_t max_block_atoms=std::numeric_limits<std::size_t>::max());
}
