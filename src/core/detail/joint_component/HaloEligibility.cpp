#include "HaloEligibility.hpp"
#include <algorithm>
#include <cmath>
#include <iterator>
#include <stdexcept>

namespace rhbm_gem::core::joint_component {
FullABCHaloEligibility EvaluateFullABCHaloEligibility(
    std::vector<JointSupport> & support,std::size_t minimum_observations)
{
    if(minimum_observations<FullABCHaloMinimumObservations)
        throw std::invalid_argument("FullABC halo minimum observations must be at least three.");
    if(std::any_of(support.begin(),support.end(),[](const auto & entry) {
        return !std::isfinite(entry.squared_distance) || entry.squared_distance<0;
    })) throw std::invalid_argument("Invalid halo support geometry.");

    std::sort(support.begin(),support.end(),[](const auto & left,const auto & right) {
        if(left.row!=right.row) return left.row<right.row;
        return left.squared_distance<right.squared_distance;
    });
    support.erase(std::unique(support.begin(),support.end(),[](const auto & left,const auto & right) {
        return left.row==right.row;
    }),support.end());

    FullABCHaloEligibility out;
    out.support_rows=support.size();
    std::sort(support.begin(),support.end(),[](const auto & left,const auto & right) {
        if(left.squared_distance!=right.squared_distance)
            return left.squared_distance<right.squared_distance;
        return left.row<right.row;
    });
    out.distinct_squared_distances=static_cast<std::size_t>(std::distance(support.begin(),
        std::unique(support.begin(),support.end(),[](const auto & left,const auto & right) {
            return left.squared_distance==right.squared_distance;
        })));
    if(out.support_rows<minimum_observations)
    {
        out.reason="insufficient-observations";
        return out;
    }
    if(out.distinct_squared_distances<FullABCHaloMinimumDistinctRadialInformation)
    {
        out.reason="insufficient-radial-information";
        return out;
    }
    out.eligible=true;
    out.reason="eligible";
    return out;
}
}
