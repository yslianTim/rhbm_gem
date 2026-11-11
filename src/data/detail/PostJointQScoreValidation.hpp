#pragma once

#include <rhbm_gem/data/object/PostJointQScoreResult.hpp>
#include "core/estimator/joint_component/HaloMedianReconstruction.hpp"
#include <map>
#include <set>

namespace rhbm_gem {
class ModelObject;
class AtomObject;
struct JointAnalysisResult;
namespace data_internal {
// One immutable snapshot lookup per batch/serialization, not per atom.
class PostJointQScoreValidation
{
    const ModelObject & m_model;
    const JointAnalysisResult * m_snapshot;
    std::map<int,LocalStageEstimate> m_expected;
    std::set<std::string> m_donors, m_contributors;
    std::vector<core::detail::HaloMedianTargetDonor> m_donor_records;
    std::optional<core::detail::GeometryLocalTargetDonorIndex> m_donor_index;
    mutable std::map<int,std::vector<std::string>> m_extra_donor_ids;
public:
    explicit PostJointQScoreValidation(const ModelObject &);
    std::string Validate(const PostJointQScoreResult &, const AtomObject &, const LocalStageEstimate &) const;
};
}
}
