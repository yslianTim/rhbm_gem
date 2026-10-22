#pragma once

#include "core/detail/joint_component/JointComponentEstimator.hpp"
#include <rhbm_gem/data/object/JointAnalysisResult.hpp>
#include <cstddef>
#include <array>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace rhbm_gem {
class ModelObject;
namespace core { class JointProblem; }
namespace core::detail {

struct HaloMedianCoverageAudit
{
    std::size_t valid_target_count{};
    std::size_t geometry_local_halos{};
};

struct HaloMedianTargetDonor
{
    std::size_t contributor_index{};
    std::string atom_id;
    std::array<double,3> position{};
    std::array<double,3> point{};
};

struct HaloMedianBuildResult
{
    // TakeSnapshot() transfers this value to the persisted Joint result; the
    // donor workspace remains alive for postprocessing.
    std::optional<JointHaloMedianSnapshot> snapshot;
    HaloMedianCoverageAudit audit;
    std::vector<HaloMedianTargetDonor> target_donors;

    JointHaloMedianSnapshot TakeSnapshot();
};

class GeometryLocalTargetDonorIndex
{
public:
    explicit GeometryLocalTargetDonorIndex(const std::vector<HaloMedianTargetDonor> &);
    ~GeometryLocalTargetDonorIndex();
    GeometryLocalTargetDonorIndex(const GeometryLocalTargetDonorIndex &) = delete;
    GeometryLocalTargetDonorIndex & operator=(const GeometryLocalTargetDonorIndex &) = delete;
    std::vector<const HaloMedianTargetDonor *> Query(const std::array<double,3> &) const;

private:
    struct Impl;
    std::unique_ptr<Impl> m_impl;
};
std::array<double,3> GeometryLocalTargetMedian(
    const std::vector<const HaloMedianTargetDonor *> &);
std::vector<HaloMedianTargetDonor> QualifiedJointTargetDonors(
    const ModelObject &,const JointProblem &,const JointAnalysisResult &);
// Builds the immutable post-Joint Halo reconstruction snapshot. Only target
// point evidence from this JointAnalysisResult may enter the donor pool.
HaloMedianBuildResult BuildHaloMedianReconstruction(
    const ModelObject &, const JointProblem &, const JointAnalysisResult &,
    std::string_view joint_run_id);

} // namespace core::detail
} // namespace rhbm_gem
