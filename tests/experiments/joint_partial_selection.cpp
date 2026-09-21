#include "support/JointPartialSelection.hpp"
#include "data/io/detail/JointResultJson.hpp"
#include <boost/json.hpp>
#include <chrono>
#include <iostream>
#include <sys/resource.h>
int main(int argc,char ** argv)
{
    if(argc!=2) return 2;
    const std::string name=argv[1];
    if(name!="all" && name!="partial" && name!="bridge" && name!="weak") return 2;
    auto fixture=joint_partial_test::Make(name);
    const auto start=std::chrono::steady_clock::now();
    const auto fit=rhbm_gem::core::EstimateJointComponents(*fixture.map,*fixture.model);
    const auto elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    rusage usage{}; getrusage(RUSAGE_SELF,&usage);
#ifdef __APPLE__
    const auto peak=usage.ru_maxrss;
#else
    const auto peak=usage.ru_maxrss*1024;
#endif
    const auto & input=fit.problem->Input(); std::size_t memberships=0,largest=0;
    for(const auto & support:input.support) memberships+=support.size();
    for(const auto & c:fit.components) largest=std::max(largest,c.atoms.size());
    const auto target=input.selection_domain->target_indices.size();
    auto outcome=boost::json::parse(rhbm_gem::joint_result_io::Encode(rhbm_gem::core::CaptureJointAnalysisResult(fit)));
    boost::json::object report{{"case",name},{"targets",target},{"halo",input.atom_ids.size()-target},
        {"rows",input.row_ids.size()},{"memberships",memberships},{"largest_component",largest},
        {"total_seconds",elapsed},{"process_peak_rss_bytes",peak},{"outcome",std::move(outcome)}};
    if(name=="partial")
    {
        const auto matched=rhbm_gem::core::FitJointComponents(*fit.problem,fixture.b);
        auto omitted=input; omitted.atom_ids.resize(1); omitted.support.resize(1);
        const auto missing=rhbm_gem::core::FitJointComponents(rhbm_gem::core::JointProblem(omitted),{fixture.b[0]});
        if(!matched.assembled_state || !missing.assembled_state) return 3;
        report["matched_control"]=boost::json::object{{"complete_objective",*matched.objective},
            {"omitted_objective",*missing.objective},{"complete_target_A_error",matched.assembled_state->ac[0]-fixture.a[0]},
            {"omitted_target_A_error",missing.assembled_state->ac[0]-fixture.a[0]},
            {"complete_target_B_error",matched.assembled_state->b[0]-fixture.b[0]},
            {"omitted_target_B_error",missing.assembled_state->b[0]-fixture.b[0]},
            {"complete_target_C_error",matched.assembled_state->ac[1]-fixture.c[0]},
            {"omitted_target_C_error",missing.assembled_state->ac[1]-fixture.c[0]}};
    }
    std::cout << boost::json::serialize(report) << '\n';
    return 0; // A saved numerical limitation is an outcome, not runner failure.
}
