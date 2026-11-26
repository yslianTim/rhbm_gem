#pragma once

#include "core/estimator/joint_component/JointComponentEstimator.hpp"
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

struct ReconstructionCoverageAudit
{
    std::size_t direct_contributor_count{}, qualified_donor_count{};
    std::size_t reconstructed_halos{}, unavailable{};
    std::size_t lower_bound{}, upper_bound{}, conditioning_rejected{}, overlapping_rejections{};
    double assessment_seconds{}, collection_seconds{}, index_seconds{}, query_seconds{};
};

struct ReconstructionDonor
{
    std::size_t contributor_index{};
    std::string atom_id;
    std::array<double,3> position{};
    std::array<double,3> point{};
};

class ReconstructionDonorIndex
{
public:
    explicit ReconstructionDonorIndex(std::vector<ReconstructionDonor>);
    ~ReconstructionDonorIndex();
    ReconstructionDonorIndex(const ReconstructionDonorIndex &) = delete;
    ReconstructionDonorIndex & operator=(const ReconstructionDonorIndex &) = delete;
    const std::vector<ReconstructionDonor> & Donors() const;
    std::vector<const ReconstructionDonor *> Query(const std::array<double,3> &) const;

private:
    struct Impl;
    std::unique_ptr<Impl> m_impl;
};
struct ContributorReconstructionBuildResult
{
    std::optional<JointContributorReconstructionSnapshot> snapshot;
    ReconstructionCoverageAudit audit;
    std::vector<std::optional<std::array<double,3>>> joint_points;
    std::unique_ptr<ReconstructionDonorIndex> donor_index;
    // Set by the workflow after publishing the internally constructed snapshot.
    const JointAnalysisResult * published_snapshot{};

    JointContributorReconstructionSnapshot TakeSnapshot();
};

std::array<double,3> GeometryLocalFittingMedian(
    const std::vector<const ReconstructionDonor *> &);
std::vector<ReconstructionDonor> CollectQualifiedReconstructionDonors(
    const ModelObject &,const JointAnalysisResult &);
ContributorReconstructionBuildResult BuildContributorReconstruction(
    const ModelObject &, const JointProblem &, const JointAnalysisResult &, std::string_view joint_run_id);
void ValidateSavedContributorReconstruction(const JointAnalysisResult &,const ModelObject * = nullptr);

} // namespace core::detail
} // namespace rhbm_gem
