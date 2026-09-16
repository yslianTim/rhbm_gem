#pragma once
#include "support/ObservationMatchedExperiment.hpp"

namespace second_stage_test::matched {

struct FrozenTarget
{
    std::array<Design, 2> design;
    std::array<Eigen::VectorXd, 2> neighbors;
};

// Only checkpoint parameters enter preparation and fitting; there is no truth input.
FrozenTarget PrepareFrozenTarget(const std::vector<Atom> & state, std::size_t target,
    const std::vector<Position> & positions, const std::vector<Stencil> &, double cutoff);
boost::json::array FitFrozenTarget(const FrozenTarget &, const std::vector<double> & distances,
    const Eigen::VectorXd & observations, const Atom & input, double cutoff);
void RunEstimatedNeighborSweep(const std::string & manifest, const std::string & map,
    const std::string & state_index, const std::string & output);

} // namespace second_stage_test::matched
