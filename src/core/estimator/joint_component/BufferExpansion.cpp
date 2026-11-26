#include "BufferExpansion.hpp"

#include "JointComponentEstimator.hpp"
#include "ResourceWork.hpp"
#include <algorithm>
#include <cmath>
#include <limits>

namespace rhbm_gem::core::joint_component {
namespace {
constexpr auto NoBlocker=std::numeric_limits<std::size_t>::max();
constexpr auto MultipleBlockers=NoBlocker-1;
constexpr double CoreRadiusSquared=0.5*0.5;

struct Candidate
{
    long double score{};
    std::size_t retained_support{};
    std::size_t core_hits{};
};
}

std::optional<std::size_t> EstimateBufferSelectionScratchBytes(std::size_t rows,std::size_t contributors)
{
    std::size_t blocker_bytes{},candidate_bytes{},total{};
    if(!CheckedMultiplyBytes(rows,sizeof(std::size_t),blocker_bytes) ||
        !CheckedMultiplyBytes(contributors,sizeof(Candidate)+sizeof(std::size_t),candidate_bytes) ||
        !CheckedAddBytes(blocker_bytes,candidate_bytes,total)) return {};
    return total;
}

BufferSelection SelectBuffer(const JointProblemInput & input,std::span<const std::size_t> output_indices)
{
    BufferSelection result;
    result.retained_rows.assign(input.observations.size(),true);
    if(output_indices.size()==input.support.size()) return result;

    std::vector<std::size_t> blockers(input.observations.size(),NoBlocker);
    for(std::size_t atom=0;atom<input.support.size();++atom)
    {
        if(std::binary_search(output_indices.begin(),output_indices.end(),atom)) continue;
        for(const auto & entry:input.support[atom])
        {
            auto & blocker=blockers[entry.row];
            if(blocker==NoBlocker) blocker=atom;
            else if(blocker!=atom) blocker=MultipleBlockers;
        }
    }

    std::vector<Candidate> candidates(input.support.size());
    std::vector<std::size_t> touched;
    touched.reserve(candidates.size());
    for(const auto blocker:blockers)
        if(blocker<MultipleBlockers) ++candidates[blocker].retained_support;
    for(const auto atom:output_indices)
    {
        std::size_t core_count{};
        for(const auto & entry:input.support[atom])
        {
            if(entry.squared_distance>CoreRadiusSquared) continue;
            ++core_count;
            const auto blocker=blockers[entry.row];
            if(blocker<MultipleBlockers && candidates[blocker].core_hits++==0)
                touched.push_back(blocker);
        }
        for(const auto blocker:touched)
        {
            auto & candidate=candidates[blocker];
            candidate.score+=static_cast<long double>(candidate.core_hits)/static_cast<long double>(core_count);
            candidate.core_hits=0;
        }
        touched.clear();
    }

    double best_score{};
    for(std::size_t atom=0;atom<candidates.size();++atom)
    {
        const auto & candidate=candidates[atom];
        // Compare scores rounded to binary64 in native atom order to resolve ties.
        const double score=static_cast<double>(candidate.score);
        if(!std::isfinite(score) || score<=best_score || candidate.retained_support<3) continue;
        double minimum=std::numeric_limits<double>::infinity(),maximum{};
        for(const auto & entry:input.support[atom])
        {
            if(blockers[entry.row]!=atom) continue;
            if(!std::isfinite(entry.squared_distance)) { maximum=0; break; }
            minimum=std::min(minimum,entry.squared_distance);
            maximum=std::max(maximum,entry.squared_distance);
        }
        const double tolerance=8*std::numeric_limits<double>::epsilon()*
            std::max(minimum,std::numeric_limits<double>::min());
        if(maximum<=minimum+tolerance) continue;
        result.buffer=atom;
        best_score=score;
    }

    const auto selected=result.buffer.value_or(NoBlocker);
    for(std::size_t row=0;row<blockers.size();++row)
    {
        const auto blocker=blockers[row];
        result.retained_rows[row]=blocker==NoBlocker || blocker==selected;
        if(result.buffer && blocker==selected) ++result.recovered_rows;
    }
    return result;
}
}
