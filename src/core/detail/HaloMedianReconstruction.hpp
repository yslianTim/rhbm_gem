#pragma once

#include <rhbm_gem/core/JointComponentEstimator.hpp>
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

struct HaloMedianCoverageRecord
{
    std::string halo_atom_id;
    std::size_t target_donor_count{};
    std::size_t affected_support_rows{};
    std::array<double,3> median{};
    std::array<double,3> dispersion{};
};

struct HaloMedianCoverageAudit
{
    std::size_t valid_target_count{};
    std::size_t geometry_local_halos{};
    std::array<double,3> global_median{};
    std::vector<HaloMedianCoverageRecord> halos;
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
    JointHaloMedianSnapshot snapshot;
    HaloMedianCoverageAudit audit;
    std::vector<HaloMedianTargetDonor> target_donors;
};

enum class HaloFullABCRole {FixedGlobalTargetMedian,FreeFullABC};
struct HaloFullABCClassificationRecord
{
    std::string atom_id;
    std::optional<std::size_t> original_contributor_index;
    std::size_t distinct_original_observations{};
    HaloFullABCRole role{HaloFullABCRole::FixedGlobalTargetMedian};
    std::vector<std::size_t> original_support_rows;
    std::optional<std::array<double,3>> parameter_profile;
    bool operator==(const HaloFullABCClassificationRecord &) const = default;
};
struct HaloFullABCPreclassification
{
    bool available{};
    std::string unavailable_reason;
    std::optional<std::array<double,3>> global_target_median;
    std::size_t fixed_halo_count{},free_halo_count{};
    std::vector<HaloFullABCClassificationRecord> halos;
    bool operator==(const HaloFullABCPreclassification &) const = default;
};
std::vector<std::size_t> DistinctOriginalSupportRows(const std::vector<JointSupport> &);

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
HaloFullABCPreclassification BuildHaloFullABCPreclassification(
    const ModelObject &,const JointProblem &,const JointAnalysisResult &);

// Builds the immutable post-Joint Halo reconstruction snapshot. Only target
// point evidence from this JointAnalysisResult may enter the donor pool.
HaloMedianBuildResult BuildHaloMedianReconstruction(
    const ModelObject &, const JointProblem &, const JointAnalysisResult &,
    std::string_view joint_run_id);

} // namespace core::detail
} // namespace rhbm_gem
