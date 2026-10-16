#pragma once

#include <rhbm_gem/core/JointComponentEstimator.hpp>
#include <Eigen/Dense>
#include <array>
#include <cstddef>
#include <optional>
#include <span>
#include <string>
#include <vector>

namespace joint_test::dynamic_halo {

struct Block
{
    std::string id;
    std::size_t component{};
    std::vector<std::size_t> targets;
};

struct HaloOwner
{
    std::size_t halo{};
    std::size_t owner_block{};
    std::size_t owner_overlap{};
    std::vector<std::size_t> intersecting_blocks;
    std::vector<std::size_t> overlapping_targets;
};

struct Ownership
{
    std::vector<Block> blocks;
    std::vector<HaloOwner> halos;
    std::vector<std::vector<std::size_t>> block_halos;
    std::vector<std::vector<std::size_t>> block_dynamic_rows;
    std::size_t cross_component_halos{};
    std::size_t row_block_entries{};
};

Ownership BuildOwnership(const rhbm_gem::core::JointProblemInput &, std::vector<Block>);

struct Parameters
{
    double a{}, b{1}, c{};
};

enum class Mode {AB,ABC,StaticABC};
struct AugmentedDesign
{
    std::vector<std::size_t> rows;
    Eigen::MatrixXd values;
};
struct FullState
{
    std::vector<Parameters> contributors;
    std::vector<double> prediction;
    double objective{};
};

AugmentedDesign BuildAugmentedDesign(const rhbm_gem::core::JointProblemInput &,const Ownership &,
    std::size_t,std::span<const Parameters>,Mode);
Eigen::MatrixXd CoupledWidthJacobian(const rhbm_gem::core::JointProblemInput &,const Ownership &,
    std::size_t,std::span<const Parameters>,const Eigen::VectorXd &,Mode);
std::vector<Parameters> DeriveHalos(const rhbm_gem::core::JointProblemInput &,const Ownership &,
    std::span<const Parameters>,Mode);
std::vector<double> ReplayPrediction(const rhbm_gem::core::JointProblemInput &,std::span<const Parameters>);
FullState BuildState(const rhbm_gem::core::JointProblemInput &,const Ownership &,
    std::span<const Parameters>,Mode,double);
bool CommitMonotone(FullState &,FullState &&);
std::optional<std::size_t> DenseBytes(std::size_t rows,std::size_t columns);

Parameters Mean(std::span<const Parameters>, bool include_charge);
double MeanWidthLogDerivative(double target_width, std::size_t core_size);

}
