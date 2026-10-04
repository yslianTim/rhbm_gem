#include "core/detail/joint_component/FixedNeighborBlockCoordinate.hpp"
#include "core/detail/joint_component/Problem.hpp"
#include "support/JointOperatorWorkload.hpp"
#include <boost/json.hpp>
#include <Eigen/Core>
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>
#include <sys/resource.h>

namespace {
namespace n=rhbm_gem::core::joint_component;
namespace j=boost::json;
using Clock=std::chrono::steady_clock;
using Input=rhbm_gem::core::JointProblemInput;
double Seconds(Clock::time_point start) {return std::chrono::duration<double>(Clock::now()-start).count();}
double PeakRssMb()
{
    rusage usage{}; if(getrusage(RUSAGE_SELF,&usage)!=0) throw std::runtime_error("Could not read peak RSS.");
#if defined(__APPLE__)
    return static_cast<double>(usage.ru_maxrss)/(1024.0*1024.0);
#else
    return static_cast<double>(usage.ru_maxrss)/1024.0;
#endif
}
const rhbm_gem::JointCheck * FindCheck(const rhbm_gem::core::JointFitResult & fit,const std::string & name)
{
    const auto found=std::find_if(fit.evidence.begin(),fit.evidence.end(),[&](const auto & check) {
        return check.name==name && check.scope==rhbm_gem::JointEvidenceScope::AssembledGlobal;
    });
    return found==fit.evidence.end() ? nullptr : &*found;
}
std::string CheckName(rhbm_gem::JointCheckStatus status)
{
    using S=rhbm_gem::JointCheckStatus;
    switch(status) {case S::Passed:return "Passed"; case S::Failed:return "Failed";
        case S::Unavailable:return "Unavailable"; case S::NotRun:return "NotRun";}
    return "NotRun";
}
j::value CheckValue(const rhbm_gem::core::JointFitResult & fit,const std::string & name)
{
    const auto * check=FindCheck(fit,name);
    return check ? j::value{{"status",CheckName(check->status)},
        {"value",check->value ? j::value(*check->value) : j::value(nullptr)},
        {"threshold",check->threshold ? j::value(*check->threshold) : j::value(nullptr)}} : j::value(nullptr);
}
double WidthGradient(const rhbm_gem::core::JointFitResult & fit)
{
    if(!fit.assembled_state || fit.assembled_state->width_gradient.empty()) return std::numeric_limits<double>::infinity();
    double out{}; for(double value:fit.assembled_state->width_gradient) out=std::max(out,std::abs(value)); return out;
}
j::object PackFit(const std::string & method,const rhbm_gem::core::JointFitResult & fit,double seconds)
{
    return {{"method",method},{"search_completed",fit.search_completed},
        {"search_seconds",fit.costs.search_seconds},{"total_elapsed_seconds",seconds},
        {"objective",fit.objective ? j::value(*fit.objective) : j::value(nullptr)},
        {"global_kkt",CheckValue(fit,"kkt")},{"width_gradient_inf_norm",WidthGradient(fit)},
        {"assessment_inner",CheckValue(fit,"inner")},{"assessment_gradient",CheckValue(fit,"width-stationarity")},
        {"assessment_local",CheckValue(fit,"local-correction")},{"assessment_identified",CheckValue(fit,"numerical-identifiability")},
        {"runtime_convergence",CheckName(fit.RuntimeConvergence())}};
}
struct MethodState {std::string name; const rhbm_gem::core::JointFitResult * fit{};};
double ScaledAcDifference(const rhbm_gem::core::JointProblem & problem,const MethodState & lhs,const MethodState & rhs)
{
    if(!lhs.fit || !rhs.fit || !lhs.fit->assembled_state || !rhs.fit->assembled_state) return std::numeric_limits<double>::infinity();
    const auto & input=problem.Input(); const auto & layout=problem.ParameterLayout();
    const auto & a=lhs.fit->assembled_state->ac; const auto & b=rhs.fit->assembled_state->ac;
    if(a.size()!=b.size()) return std::numeric_limits<double>::infinity();
    const double scale=problem.ObservationScale(); double difference{};
    for(std::size_t k=0;k<layout.full_atoms.size();++k)
    {
        const auto atom=layout.full_atoms[k]; const double width=rhs.fit->assembled_state->b.at(k);
        double gaussian{},charge{};
        for(const auto & support:input.support.at(atom))
        {
            if(!std::binary_search(layout.informative_rows.begin(),layout.informative_rows.end(),support.row)) continue;
            const auto basis=n::EvaluateKernel(support.squared_distance,width,2.5);
            gaussian+=basis.gaussian*basis.gaussian; charge+=basis.charge*basis.charge;
        }
        const double scales[]{std::sqrt(gaussian)/scale,std::sqrt(charge)/scale};
        for(std::size_t kind=0;kind<2;++kind)
            difference=std::max(difference,std::abs(a[2*k+kind]-b[2*k+kind])*scales[kind]);
    }
    return difference;
}
j::object SweepJson(const n::FixedNeighborBlockSweep & sweep)
{
    return {{"objective_before",sweep.objective_before},{"objective_after",sweep.objective_after},
        {"global_a_feasibility",sweep.global_a_feasibility},{"global_ac_kkt",sweep.global_ac_kkt},
        {"global_width_gradient_inf_norm",sweep.global_width_gradient_inf_norm},
        {"cache_replay_error",sweep.cache_replay_error},{"objective_replay_error",sweep.objective_replay_error},
        {"block_solves",sweep.block_solves},{"profile_evaluations",sweep.profile_evaluations},
        {"maximum_block_rows",sweep.maximum_block_rows},{"maximum_block_columns",sweep.maximum_block_columns},
        {"wall_seconds",sweep.wall_seconds}};
}
void Write(const std::filesystem::path &,const j::value &);
j::object Run(const std::string & topology,int atoms,const std::filesystem::path & output_path,bool compare_global)
{
    auto input=std::make_shared<Input>(second_stage_test::OperatorWorkload(topology,atoms));
    const rhbm_gem::core::JointProblem problem(*input);
    const std::vector<double> initial_b(static_cast<std::size_t>(atoms),.55);
    n::Vector initial_eta= n::Vector::Constant(atoms,std::log(.55)); n::FixedNeighborPolicy neighbor_policy;
    neighbor_policy.core_atoms=128;
    j::array progress_sweeps;
    const auto progress_path=std::filesystem::path(output_path.string()+".progress.json");
    neighbor_policy.sweep_observer=[&](const auto & sweep) {
        progress_sweeps.push_back(SweepJson(sweep));
        Write(progress_path,j::object{{"topology",topology},{"atoms",atoms},
            {"sweeps_completed",progress_sweeps.size()},{"sweep_telemetry",progress_sweeps}});
        std::cerr<<topology<<'-'<<atoms<<" FixedNeighbor sweep "<<progress_sweeps.size()
            <<" KKT="<<sweep.global_ac_kkt<<" width-grad="<<sweep.global_width_gradient_inf_norm
            <<" seconds="<<sweep.wall_seconds<<'\n';
    };
    std::cerr<<topology<<'-'<<atoms<<" FixedNeighbor started\n";
    auto started=Clock::now(); const auto neighbor=n::SearchFixedNeighbor(problem,initial_eta,neighbor_policy);
    const double neighbor_seconds=Seconds(started);
    const double neighbor_search_seconds=std::accumulate(neighbor.sweeps.begin(),neighbor.sweeps.end(),0.0,
        [](double total,const auto & sweep){return total+sweep.wall_seconds;});
    std::cerr<<topology<<'-'<<atoms<<" FixedNeighbor finished in "<<neighbor_seconds<<" s; converged="
        <<neighbor.search_converged<<" sweeps="<<neighbor.sweeps.size()<<" KKT="
        <<(neighbor.sweeps.empty() ? 0.0 : neighbor.sweeps.back().global_ac_kkt)<<" endpoint="
        <<neighbor.endpoint_certified<<" runtime="<<CheckName(neighbor.fit.RuntimeConvergence())<<'\n';
    j::array sweeps;
    for(const auto & sweep:neighbor.sweeps)
        sweeps.push_back(SweepJson(sweep));
    j::array blocks;
    for(const auto & block:neighbor.blocks)
        blocks.push_back({{"sweep",block.sweep},{"block",block.block},{"atoms",block.core_atoms.size()},
            {"affected_rows",block.affected_rows},{"objective_before",block.objective_before},
            {"objective_after",block.objective_after},{"local_objective_before",block.local_objective_before},
            {"local_objective_after",block.local_objective_after},{"global_replay_delta",block.global_replay_delta},
            {"local_global_delta_error",block.local_global_delta_error},
            {"objective_replay_enclosure",block.objective_replay_enclosure},
            {"profile_evaluations",block.profile_evaluations},{"accepted_updates",block.accepted_updates},
            {"search_seconds",block.search_seconds},{"accepted",block.accepted},
            {"status",block.status},{"reason",block.reason},{"local_search_stop_reason",block.local_search_stop_reason}});
    j::object neighbor_json{{"method","FixedNeighbor"},{"search_converged",neighbor.search_converged},
        {"search_reason",neighbor.reason},{"sweeps",neighbor.sweeps.size()},
        {"sweeps_to_stationarity",neighbor.sweeps_to_stationarity},
        {"endpoint_certified",neighbor.endpoint_certified},
        {"endpoint_trust_reason",neighbor.endpoint_trust.reason},
        {"objective",neighbor.fit.objective ? j::value(*neighbor.fit.objective) : j::value(nullptr)},
        {"global_kkt",CheckValue(neighbor.fit,"kkt")},{"width_gradient_inf_norm",WidthGradient(neighbor.fit)},
        {"assessment_inner",CheckValue(neighbor.fit,"inner")},{"assessment_gradient",CheckValue(neighbor.fit,"width-stationarity")},
        {"assessment_local",CheckValue(neighbor.fit,"local-correction")},{"assessment_identified",CheckValue(neighbor.fit,"numerical-identifiability")},
        {"runtime_convergence",CheckName(neighbor.fit.RuntimeConvergence())},
        {"sweep_telemetry",sweeps},{"block_telemetry",blocks},
        {"search_seconds",neighbor_search_seconds},{"total_elapsed_seconds",neighbor_seconds}};
    if(!compare_global)
        return {{"topology",topology},{"atoms",atoms},{"core_atoms",neighbor_policy.core_atoms},
            {"maximum_sweeps",neighbor_policy.maximum_sweeps},{"observation_scale",problem.ObservationScale()},
            {"fixed_neighbor",neighbor_json},{"peak_rss_mb",PeakRssMb()}};
    n::SearchPolicy legacy_policy;
    std::cerr<<topology<<'-'<<atoms<<" global LegacyCompact started\n";
    started=Clock::now(); const auto legacy=n::FitWithSearchPolicy(problem,initial_b,legacy_policy);
    const double legacy_seconds=Seconds(started);
    std::cerr<<topology<<'-'<<atoms<<" global LegacyCompact finished in "<<legacy_seconds<<" s\n";
    n::SearchPolicy operator_policy; operator_policy.method=n::SearchMethod::OperatorPcg;
    std::cerr<<topology<<'-'<<atoms<<" global OperatorPcg started\n";
    started=Clock::now(); const auto operator_fit=n::FitWithSearchPolicy(problem,initial_b,operator_policy);
    const double operator_seconds=Seconds(started);
    std::cerr<<topology<<'-'<<atoms<<" global OperatorPcg finished in "<<operator_seconds<<" s\n";
    const MethodState reference{"LegacyCompact",&legacy};
    j::array comparisons;
    for(const MethodState & candidate:std::array<MethodState,2>{{{"OperatorPcg",&operator_fit},{"FixedNeighbor",&neighbor.fit}}})
    {
        const double objective_difference=legacy.objective && candidate.fit->objective ?
            std::abs(*legacy.objective-*candidate.fit->objective) : std::numeric_limits<double>::infinity();
        const double eta_difference=legacy.assembled_state && candidate.fit->assembled_state ?
            (n::Vector::Map(legacy.assembled_state->log_b.data(),static_cast<Eigen::Index>(legacy.assembled_state->log_b.size()))-
             n::Vector::Map(candidate.fit->assembled_state->log_b.data(),static_cast<Eigen::Index>(candidate.fit->assembled_state->log_b.size()))).lpNorm<Eigen::Infinity>() :
            std::numeric_limits<double>::infinity();
        comparisons.push_back({{"against","LegacyCompact"},{"method",candidate.name},
            {"objective_difference",objective_difference},{"eta_inf_difference",eta_difference},
            {"scaled_ac_inf_difference",ScaledAcDifference(problem,reference,candidate)}});
    }
    return {{"topology",topology},{"atoms",atoms},{"core_atoms",neighbor_policy.core_atoms},
        {"maximum_sweeps",neighbor_policy.maximum_sweeps},{"observation_scale",problem.ObservationScale()},
        {"global_legacy_compact",PackFit("LegacyCompact",legacy,legacy_seconds)},
        {"global_operator_pcg",PackFit("OperatorPcg",operator_fit,operator_seconds)},
        {"fixed_neighbor",neighbor_json},{"comparisons",comparisons},{"peak_rss_mb",PeakRssMb()}};
}
void Write(const std::filesystem::path & path,const j::value & value)
{std::ofstream output(path); if(!output) throw std::runtime_error("Could not open experiment output."); output<<j::serialize(value)<<'\n';}
}
int main(int argc,char ** argv)
{
    try {
        if(argc!=5 || (std::string(argv[1])!="--case" && std::string(argv[1])!="--neighbor-only"))
            throw std::invalid_argument("Usage: joint_fixed_neighbor_experiment --case|--neighbor-only OUTPUT_FILE TOPOLOGY ATOMS");
        Eigen::setNbThreads(1);
        const std::filesystem::path output_path(argv[2]);
        if(output_path.has_parent_path()) std::filesystem::create_directories(output_path.parent_path());
        Write(output_path,Run(argv[3],std::stoi(argv[4]),output_path,std::string(argv[1])=="--case"));
        std::cout<<argv[3]<<'-'<<argv[4]<<" fixed-neighbor endpoint experiment complete\n";
        return 0;
    } catch(const std::exception & error) {std::cerr<<error.what()<<'\n'; return 1;}
}
