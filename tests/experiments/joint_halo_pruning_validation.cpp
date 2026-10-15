#include "core/detail/joint_component/HaloEligibility.hpp"
#include "core/detail/joint_component/Numerics.hpp"
#include "core/detail/joint_component/Problem.hpp"
#include <rhbm_gem/core/JointComponentEstimator.hpp>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <initializer_list>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {
namespace core=rhbm_gem::core;
namespace n=core::joint_component;
using Clock=std::chrono::steady_clock;
using Support=std::vector<core::JointSupport>;

struct HaloTruth
{
    std::string id;
    Support support;
    double a{},c{},b{};
};
struct Scenario
{
    std::string name;
    std::vector<HaloTruth> low_support_halos;
};
struct FitRecord
{
    std::string scenario;
    std::size_t threshold{},candidates{},retained{},excluded{},excluded_rows{},excluded_radial{},q_below_three{},rows{},
        components_before{},components_after{},block_search_failures{},failed_components{};
    bool search_completed{},numerically_converged{};
    rhbm_gem::JointCheckStatus target_status{rhbm_gem::JointCheckStatus::Unavailable};
    std::optional<double> objective,target_a_error,target_c_error,target_b_error,residual_rms,global_kkt,global_width_gradient;
    std::optional<std::size_t> construction_peak_bytes;
    double elapsed_seconds{};
};

constexpr std::size_t RowCount=72;
constexpr double TargetA=2.0,TargetC=.18,TargetB=.55;
std::string StatusName(rhbm_gem::JointCheckStatus status)
{
    using S=rhbm_gem::JointCheckStatus;
    switch(status)
    {case S::Passed:return "Passed";case S::Failed:return "Failed";case S::Unavailable:return "Unavailable";case S::NotRun:return "NotRun";}
    return "NotRun";
}
std::optional<double> CheckValue(const core::JointFitResult & fit,const std::string & name)
{
    const auto found=std::find_if(fit.evidence.begin(),fit.evidence.end(),[&](const auto & check) {
        return check.name==name && check.scope==rhbm_gem::JointEvidenceScope::AssembledGlobal;
    });
    if(found==fit.evidence.end()) return {};
    return found->value;
}
Support MakeSupport(std::initializer_list<std::size_t> rows,std::size_t salt)
{
    Support out; out.reserve(rows.size()); std::size_t k{};
    for(const auto row:rows)
    {
        const auto radial=(7*k+5*salt)%30;
        out.push_back({row,.12+.19*static_cast<double>(radial)});
        ++k;
    }
    return out;
}
Support MakeSupport(std::size_t count,std::size_t salt)
{
    Support out; out.reserve(count);
    for(std::size_t k=0;k<count;++k)
    {
        const auto row=(11*salt+7*k)%RowCount;
        const auto radial=(7*k+5*salt)%30;
        out.push_back({row,.12+.19*static_cast<double>(radial)});
    }
    std::sort(out.begin(),out.end(),[](const auto & left,const auto & right){return left.row<right.row;});
    return out;
}
std::vector<HaloTruth> CommonHalos()
{
    return {
        {"m3-weak",MakeSupport(3,1),.015,.003,.62},
        {"m4-overlap",MakeSupport({16,17,18,19},2),.12,.018,.68},
        {"m4-radial-deficient",{{28,.6},{29,.6},{30,1.2},{31,1.2}},.03,.002,.72},
        {"m6-medium",MakeSupport(6,3),.14,.012,.76},
        {"m8-medium",MakeSupport(8,4),.10,-.009,.88},
        {"m12-reference",MakeSupport(12,5),.08,.006,.71}};
}
std::vector<Scenario> Scenarios()
{
    return {
        {"weak-omitted-halo",{{"low-weak",MakeSupport({16,17},6),.03,.004,.65}}},
        {"strong-low-support-halo",{{"low-strong",MakeSupport({16,17},6),8.0,1.2,.75}}},
        {"multiple-overlapping-omitted-halos",{
            {"low-overlap-a",MakeSupport({16,17},6),1.5,.15,.65},
            {"low-overlap-b",MakeSupport({16,18},7),1.2,-.10,.75},
            {"low-overlap-c",MakeSupport({17,18},8),.8,.20,.85}}}};
}
core::JointProblemInput MakeInput(const Scenario & scenario,std::size_t minimum_observations,
    std::vector<double> & initial_b,std::size_t & candidates,std::size_t & retained,
    std::size_t & excluded_rows,std::size_t & excluded_radial)
{
    auto halos=CommonHalos();
    halos.insert(halos.end(),scenario.low_support_halos.begin(),scenario.low_support_halos.end());
    candidates=halos.size(); retained=0; excluded_rows=0; excluded_radial=0;

    core::JointProblemInput full;
    full.atom_ids={"target"}; full.support.resize(1);
    full.row_ids.reserve(RowCount); full.observations.assign(RowCount,0.0);
    full.support.front().reserve(RowCount);
    for(std::size_t row=0;row<RowCount;++row)
    {
        full.row_ids.push_back("r"+std::to_string(row));
        const double square=.08+.078*static_cast<double>(row);
        full.support.front().push_back({row,square});
        const auto basis=n::EvaluateKernel(square,TargetB,2.5);
        full.observations[row]=TargetA*basis.gaussian+TargetC*basis.charge;
    }
    const auto add_truth=[&](const HaloTruth & halo) {
        for(const auto & membership:halo.support)
        {
            const auto basis=n::EvaluateKernel(membership.squared_distance,halo.b,2.5);
            full.observations[membership.row]+=halo.a*basis.gaussian+halo.c*basis.charge;
        }
    };
    for(const auto & halo:halos) add_truth(halo);
    full.selection_domain.emplace(); full.selection_domain->target_indices={0};

    core::JointProblemInput pruned;
    pruned.row_ids=full.row_ids; pruned.observations=full.observations;
    pruned.atom_ids={"target"}; pruned.support={full.support.front()}; initial_b={TargetB};
    pruned.selection_domain.emplace(); pruned.selection_domain->target_indices={0};
    pruned.selection_domain->contributor_policy=std::string(rhbm_gem::FullABCHaloContributorPolicy);
    rhbm_gem::JointHaloEligibilityProvenance provenance;
    provenance.minimum_observations=minimum_observations;
    provenance.candidate_halos=candidates;
    for(const auto & halo:halos)
    {
        auto eligibility_support=halo.support;
        const auto eligibility=n::EvaluateFullABCHaloEligibility(eligibility_support,minimum_observations);
        if(!eligibility.eligible)
        {
            provenance.excluded_halos.push_back({halo.id,eligibility.support_rows,
                eligibility.distinct_squared_distances,eligibility.reason});
            continue;
        }
        ++retained;
        pruned.atom_ids.push_back(halo.id); pruned.support.push_back(halo.support); initial_b.push_back(halo.b);
    }
    provenance.retained_halos=retained;
    for(const auto & excluded:provenance.excluded_halos)
    {
        if(excluded.reason=="insufficient-observations") ++excluded_rows;
        else if(excluded.reason=="insufficient-radial-information") ++excluded_radial;
    }
    std::sort(provenance.excluded_halos.begin(),provenance.excluded_halos.end(),
        [](const auto & left,const auto & right){return left.atom_id<right.atom_id;});
    pruned.selection_domain->halo_eligibility=std::move(provenance);
    return pruned;
}
FitRecord Run(const Scenario & scenario,std::size_t threshold)
{
    std::vector<double> initial_b; std::size_t candidates{},retained{},excluded_rows{},excluded_radial{};
    auto input=MakeInput(scenario,threshold,initial_b,candidates,retained,excluded_rows,excluded_radial);
    std::size_t q_below_three{};
    auto geometry=CommonHalos(); geometry.insert(geometry.end(),scenario.low_support_halos.begin(),scenario.low_support_halos.end());
    for(const auto & halo:geometry)
    {
        auto support=halo.support;
        if(n::EvaluateFullABCHaloEligibility(support,3).distinct_squared_distances<3) ++q_below_three;
    }
    auto full=input;
    full.selection_domain->contributor_policy="all-non-hydrogen";
    full.selection_domain->halo_eligibility.reset();
    // Count topology before selection using every synthetic candidate support.
    auto all=CommonHalos(); all.insert(all.end(),scenario.low_support_halos.begin(),scenario.low_support_halos.end());
    for(const auto & halo:all)
    {
        if(std::find(full.atom_ids.begin(),full.atom_ids.end(),halo.id)!=full.atom_ids.end()) continue;
        full.atom_ids.push_back(halo.id); full.support.push_back(halo.support);
    }
    const core::JointProblem before(std::move(full));
    const core::JointProblem problem(std::move(input));
    const auto & data=core::JointProblemAccess::Get(problem);
    const auto start=Clock::now(); const auto fit=core::FitJointComponents(problem,initial_b);
    FitRecord out; out.scenario=scenario.name; out.threshold=threshold; out.candidates=candidates;
    out.retained=retained; out.excluded=candidates-retained; out.excluded_rows=excluded_rows;
    out.excluded_radial=excluded_radial; out.q_below_three=q_below_three;
    out.rows=problem.Input().observations.size();
    out.components_before=core::JointProblemAccess::Get(before).partition.components.size();
    out.components_after=data.partition.components.size();
    out.search_completed=fit.search_completed;
    out.numerically_converged=fit.RuntimeConvergence()==rhbm_gem::JointCheckStatus::Passed;
    out.target_status=fit.TargetRuntimeConvergence(); out.objective=fit.objective;
    out.global_kkt=CheckValue(fit,"conditional-ac-optimality");
    out.global_width_gradient=CheckValue(fit,"width-stationarity");
    out.construction_peak_bytes=data.construction_memory.estimated_peak_bytes;
    for(const auto & component:fit.components)
    {
        if(component.stop_reason=="block-search-failed") ++out.block_search_failures;
        if(!component.state) ++out.failed_components;
    }
    if(fit.prediction)
    {
        double residual_square{};
        for(std::size_t row=0;row<fit.prediction->size();++row)
        {const double residual=fit.prediction->at(row)-problem.Input().observations[row]; residual_square+=residual*residual;}
        out.residual_rms=std::sqrt(residual_square/static_cast<double>(fit.prediction->size()));
    }
    if(fit.assembled_state)
    {
        out.target_a_error=fit.assembled_state->ac[0]-TargetA;
        out.target_c_error=fit.assembled_state->ac[1]-TargetC;
        out.target_b_error=fit.assembled_state->b[0]-TargetB;
    }
    out.elapsed_seconds=std::chrono::duration<double>(Clock::now()-start).count();
    return out;
}
void PrintOptional(std::optional<double> value)
{if(value) std::cout<<*value; else std::cout<<"NA";}
void PrintOptional(std::optional<std::size_t> value)
{if(value) std::cout<<*value; else std::cout<<"NA";}
void Print(const FitRecord & result)
{
    std::cout<<std::setprecision(10)<<result.scenario<<','<<result.threshold<<','<<result.rows<<",1,"<<result.candidates<<','
        <<result.retained<<','<<result.excluded<<','<<result.excluded_rows<<','
        <<result.excluded_radial<<','<<result.q_below_three<<','
        <<static_cast<double>(result.retained)/static_cast<double>(result.candidates)<<','
        <<result.components_before<<','<<result.components_after<<','
        <<result.block_search_failures<<','<<result.failed_components<<','
        <<(result.search_completed?"true":"false")<<','<<(result.numerically_converged?"true":"false")<<','
        <<StatusName(result.target_status)<<',';
    PrintOptional(result.objective); std::cout<<','; PrintOptional(result.global_kkt); std::cout<<',';
    PrintOptional(result.global_width_gradient); std::cout<<',';
    PrintOptional(result.target_a_error); std::cout<<','; PrintOptional(result.target_c_error); std::cout<<',';
    PrintOptional(result.target_b_error); std::cout<<','; PrintOptional(result.residual_rms); std::cout<<',';
    PrintOptional(result.construction_peak_bytes); std::cout<<','<<std::setprecision(9)<<result.elapsed_seconds<<'\n';
}
}

int main()
{
    try
    {
        std::cout<<"scenario,minimum_observations,observation_rows,target_count,candidate_halos,retained_halos,excluded_halos,excluded_insufficient_rows,excluded_insufficient_radial,q_below_three_candidates,retention_ratio,components_before,components_after,block_search_failures,components_without_state,search_completed,numerical_success,target_status,objective,global_ac_kkt,global_logb_gradient,target_delta_a,target_delta_c,target_delta_b,residual_rms,construction_preflight_peak_bytes,elapsed_seconds\n";
        for(const auto & scenario:Scenarios())
            for(const std::size_t threshold:std::initializer_list<std::size_t>{3,4,6,8,12}) Print(Run(scenario,threshold));
        std::cerr<<"process_peak_rss_bytes=";
        if(const auto rss=n::ProcessPeakRssBytes()) std::cerr<<*rss<<'\n'; else std::cerr<<"Unavailable\n";
        return 0;
    }
    catch(const std::exception & error)
    {
        std::cerr<<"joint halo pruning validation failed: "<<error.what()<<'\n'; return 1;
    }
}
