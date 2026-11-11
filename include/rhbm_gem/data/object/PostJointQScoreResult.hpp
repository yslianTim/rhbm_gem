#pragma once

#include <rhbm_gem/utils/hrl/GaussianEstimationTypes.hpp>

namespace rhbm_gem {

inline constexpr char PostJointQScoreContract[] = "joint-peeled-qscore-v1";

struct PostJointQScoreValues
{
    double standard{};
    double joint{};
};

// Only extra-domain Halos need additional provenance; recorded Halos are owned
// by the immutable Joint snapshot. Representatives are never optimizer endpoints.
struct QScoreExtraHaloSource
{
    int atom_id{};
    std::vector<std::string> donor_atom_ids;
};

struct PostJointQScoreResult
{
    EstimateSource source;
    std::string contract{PostJointQScoreContract};
    EvidenceStatus status{EvidenceStatus::Unavailable};
    std::optional<PostJointQScoreValues> scores;
    std::string reason;
    double reference_height{}, reference_offset{};
    std::size_t sample_count{}, contributor_count{}, missing_contributor_count{};
    std::size_t evaluated_contributions{}, halo_count{};
    std::vector<QScoreExtraHaloSource> extra_halos;
};

} // namespace rhbm_gem
