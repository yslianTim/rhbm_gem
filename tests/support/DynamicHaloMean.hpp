#pragma once

#include <rhbm_gem/core/JointComponentEstimator.hpp>
#include <array>
#include <cstddef>
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

Parameters Mean(std::span<const Parameters>, bool include_charge);
double MeanWidthLogDerivative(double target_width, std::size_t core_size);

}
