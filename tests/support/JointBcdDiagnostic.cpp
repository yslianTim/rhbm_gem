#include "support/JointBcdDiagnostic.hpp"
#include "core/detail/joint_component/Problem.hpp"

namespace joint_component_test {
namespace n=rhbm_gem::core::joint_component;

n::BcdSearchResult SearchBcdFixture(const rhbm_gem::core::JointProblem & problem,n::VectorRef initial_eta,
    const n::BcdSearchPolicy & policy,std::vector<n::BcdSweepRecord> * sweeps)
{
    const auto & data=rhbm_gem::core::JointProblemAccess::Get(problem);
    auto run_policy=policy;
    run_policy.collect_diagnostics=true;
    if(sweeps)
        run_policy.production_progress=[&](const n::BcdSweepRecord & sweep) {
            sweeps->push_back(sweep);
            if(policy.production_progress) policy.production_progress(sweep);
        };
    return n::SearchBcdComponent(data.domain,data.layout,data.y,initial_eta,data.context,run_policy);
}
}
