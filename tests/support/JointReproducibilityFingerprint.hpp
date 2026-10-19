#pragma once

#include "core/detail/joint_component/Problem.hpp"

#include <boost/json/object.hpp>
#include <filesystem>
#include <string>
#include <vector>

namespace rhbm_gem::joint_reproducibility {

struct FingerprintOptions
{
    std::string execution_role;
    bool instrumentation_enabled{};
    std::string fixed_neighbor_local_search{"profile-lm"};
};

boost::json::object WriteFingerprint(
    const std::filesystem::path & output,
    const std::filesystem::path & model_path,
    const std::filesystem::path & map_path,
    double normalization_divisor,
    const MapObject & map,
    const ModelObject & model,
    const core::JointProblem & problem,
    const std::vector<double> & initial_widths,
    const FingerprintOptions & options);

} // namespace rhbm_gem::joint_reproducibility
