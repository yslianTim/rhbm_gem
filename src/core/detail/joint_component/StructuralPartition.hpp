#pragma once
#include "SnapshotViews.hpp"
#include "core/detail/joint_component/JointComponentEstimator.hpp"
#include <string>
#include <vector>

namespace rhbm_gem::core::joint_component {
struct StructuralCore
{
    std::string id;
    Indices atoms,affected_rows,neighbor_atoms; // Indices in the supplied problem view.
};
struct StructuralBlockPartition
{
    std::vector<StructuralCore> cores;
};
struct StructuralBlockCensus
{
    std::size_t affected_row_entries{},neighbor_atom_entries{},total_block_support_entries{},maximum_block_rows{},
        maximum_block_support_entries{},maximum_neighbor_atoms{};
};

class StructuralTopology
{
    const Domain & domain_;
    const Identities & identities_;
    std::vector<bool> informative_rows_;
    std::vector<std::size_t> row_offsets_;
    Indices row_atoms_;
    std::vector<std::size_t> atom_order_,neighbor_marks_;
    std::size_t neighbor_generation_{};
    bool LessAtom(Eigen::Index,Eigen::Index) const;
    Indices Neighbors(std::size_t);
public:
    StructuralTopology(const Domain &,const Identities &,const JointParameterLayout &);
    StructuralBlockPartition BuildCores(std::size_t core_atoms,bool include_metadata=true);
    StructuralBlockCensus Census(const StructuralBlockPartition &) const;
    std::size_t StorageBytes() const;
    void PopulateMetadata(StructuralBlockPartition &);
};
}
