#pragma once

#include <cstddef>
#include <optional>

namespace rhbm_gem::core::detail {

struct IterationDiagnostics
{
    std::optional<double> accepted_maximum_transformed_change{};
    double proposal_maximum_transformed_change{ 0.0 };
};

struct FinalDependencyPolishDiagnostic
{
    std::size_t component_count{ 0 }, attempted_component_count{ 0 }, accepted_component_count{ 0 };
    std::size_t atom_count{ 0 }, parameter_count{ 0 }, round_count{ 0 }, suspicious_candidate_atom_count{ 0 };
    std::optional<double> objective_before{}, objective_after{};
    double elapsed_milliseconds{ 0.0 };
};

struct SecondStageSeedSummary
{
    std::size_t local_mdpde{ 0 }, global_median{ 0 };
};

} // namespace rhbm_gem::core::detail
