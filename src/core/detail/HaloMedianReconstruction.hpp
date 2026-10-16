#pragma once

#include <rhbm_gem/data/object/JointAnalysisResult.hpp>
#include <cstddef>
#include <map>
#include <string>
#include <string_view>
#include <vector>

namespace rhbm_gem {
class ModelObject;
namespace core { class JointProblem; }
namespace core::detail {

struct HaloMedianCoverageRecord
{
    GroupKey group_key{};
    std::size_t target_donor_count{};
    std::size_t halo_count{};
    std::size_t affected_support_rows{};
    std::array<double,3> median{};
    std::array<double,3> dispersion{};
    bool uses_global_fallback{};
};

struct HaloMedianCoverageAudit
{
    std::size_t valid_target_count{};
    std::size_t group_count{};
    std::size_t same_group_halos{};
    std::size_t fallback_halos{};
    std::size_t fallback_support_rows{};
    std::array<double,3> global_median{};
    std::vector<HaloMedianCoverageRecord> groups;
};

struct HaloMedianBuildResult
{
    JointHaloMedianSnapshot snapshot;
    HaloMedianCoverageAudit audit;
};

// Builds the immutable post-Joint Halo reconstruction snapshot. Only target
// point evidence from this JointAnalysisResult may enter the donor pool.
HaloMedianBuildResult BuildHaloMedianReconstruction(
    const ModelObject &, const JointProblem &, const JointAnalysisResult &,
    std::string_view joint_run_id);

} // namespace core::detail
} // namespace rhbm_gem
