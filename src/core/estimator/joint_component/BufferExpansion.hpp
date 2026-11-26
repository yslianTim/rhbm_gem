#pragma once

#include <cstddef>
#include <optional>
#include <span>
#include <vector>

namespace rhbm_gem::core {
struct JointProblemInput;
namespace joint_component {

struct BufferSelection
{
    std::optional<std::size_t> buffer;
    std::vector<bool> retained_rows;
    std::size_t recovered_rows{};
};

// Uses the Original snapshot's native support and its canonical atom order.
BufferSelection SelectBuffer(const JointProblemInput &,std::span<const std::size_t> output_indices);
std::optional<std::size_t> EstimateBufferSelectionScratchBytes(std::size_t rows,std::size_t contributors);

}
}
