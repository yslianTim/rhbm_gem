#include "core/detail/joint_component/Problem.hpp"
#include "core/detail/joint_component/TargetEvidence.hpp"
#include "support/JointSyntheticWorkload.hpp"
#include <boost/json.hpp>
#include <fstream>
#include <chrono>
#include <ctime>
#include <sys/resource.h>
#include <iostream>
#include <numeric>

namespace n=rhbm_gem::core::joint_component;
using namespace rhbm_gem;
namespace j=boost::json;
using Object=j::object;
template<class T> j::value Optional(const std::optional<T> & value)
{return value ? j::value(*value) : j::value(nullptr);}
auto Number=Optional<double>;
auto OptionalSize=Optional<std::size_t>;
auto OptionalBoolean=Optional<bool>;
j::array BlockConditioning(const std::vector<JointBlockConditioning> & values)
{
    j::array out;
    for(const auto & x:values) out.emplace_back(Object{{"block_id",x.block_id},{"status",x.status},{"reason",x.reason},
        {"conditioning_classification",x.conditioning_classification},
        {"core_atoms",j::value_from(x.core_atoms)},{"core_atom_count",x.core_atom_count},
        {"affected_row_count",x.affected_row_count},{"free_ac_columns",x.free_ac_columns},
        {"active_amplitude_count",x.active_amplitude_count},{"free_ac_rank",OptionalSize(x.free_ac_rank)},
        {"projected_width_rank",OptionalSize(x.projected_width_rank)},{"profile_rank",OptionalSize(x.profile_rank)},
        {"active_width_bound_targets",OptionalSize(x.active_width_bound_targets)},
        {"workspace_peak_memory_estimate_bytes",OptionalSize(x.workspace_peak_memory_estimate_bytes)},
        {"free_ac_sigma_min",Number(x.free_ac_sigma_min)},{"free_ac_sigma_max",Number(x.free_ac_sigma_max)},
        {"free_ac_rank_threshold",Number(x.free_ac_rank_threshold)},
        {"free_ac_rank_threshold_margin_relative",Number(x.free_ac_rank_threshold_margin_relative)},
        {"local_ac_projected_kkt",Number(x.local_ac_projected_kkt)},
        {"local_projected_width_kkt",Number(x.local_projected_width_kkt)},
        {"free_ac_threshold_sensitive",OptionalBoolean(x.free_ac_threshold_sensitive)},
        {"free_ac_condition",Number(x.free_ac_condition)},{"projected_width_condition",Number(x.projected_width_condition)},
        {"projected_width_sigma_min",Number(x.projected_width_sigma_min)},
        {"projected_width_sigma_max",Number(x.projected_width_sigma_max)},
        {"projected_width_rank_threshold",Number(x.projected_width_rank_threshold)},
        {"projected_width_rank_threshold_margin_relative",Number(x.projected_width_rank_threshold_margin_relative)},
        {"projected_width_threshold_sensitive",OptionalBoolean(x.projected_width_threshold_sensitive)},
        {"profile_sigma_min",Number(x.profile_sigma_min)},{"profile_sigma_max",Number(x.profile_sigma_max)},
        {"profile_rank_threshold",Number(x.profile_rank_threshold)},
        {"profile_rank_threshold_margin_relative",Number(x.profile_rank_threshold_margin_relative)},
        {"profile_threshold_sensitive",OptionalBoolean(x.profile_threshold_sensitive)},
        {"profile_condition",Number(x.profile_condition)},{"free_ac_status",x.free_ac_status},
        {"free_ac_reason",x.free_ac_reason},{"projected_width_status",x.projected_width_status},
        {"projected_width_reason",x.projected_width_reason},{"free_ac_conditioning",x.free_ac_conditioning},
        {"projected_width_conditioning",x.projected_width_conditioning},{"profile_conditioning",x.profile_conditioning},
        {"profile_status",x.profile_status},
        {"profile_reason",x.profile_reason}});
    return out;
}

j::array Values(const n::Vector & values)
{
    j::array out;
    for(Eigen::Index k=0;k<values.size();++k) out.push_back(values(k));
    return out;
}
int main(int argc,char ** argv)
{
    if(argc!=2) {std::cerr<<"usage: joint cleanup snapshot <output.jsonl>\n"; return 2;}
    std::ofstream output(argv[1]);
    if(!output) return 2;
    using Input=rhbm_gem::core::JointProblemInput;
    struct Fixture {std::string name,digest; Input input; std::size_t core_atoms;};
    std::vector<Fixture> fixtures;
    fixtures.push_back({"high-overlap-cube-8",
        "2a1894dc5dc5f7789fbc802477e77e08f520a6c25f0c7eac6b42a4caf3f717b2",
        second_stage_test::SyntheticJointWorkload("cube",8),4});
    auto rank_deficient=second_stage_test::SyntheticJointWorkload("chain",3);
    for(auto & support:rank_deficient.support)
        for(auto & entry:support) entry.squared_distance=0.0;
    fixtures.push_back({"rank-deficient-constant-support",
        "fd39aa13434e60a2193641f97ff9bdd05669fb6bc1b57a4d52f454bdacaa7151",
        std::move(rank_deficient),2});
    Input near_singular;
    near_singular.atom_ids={"1","2"};
    near_singular.row_ids={"0","1","2","3","4"};
    near_singular.observations.assign(5,0.0);
    near_singular.support={{{0,.1},{1,.5},{2,1.0},{3,1.5},{4,2.0}},
        {{0,.1+1e-8},{1,.5-1e-8},{2,1.0+1e-8},{3,1.5-1e-8},{4,2.0+1e-8}}};
    near_singular.selection_domain.emplace(); near_singular.selection_domain->target_indices={0,1};
    for(std::size_t row=0;row<near_singular.observations.size();++row)
    {
        const auto first=n::EvaluateKernel(near_singular.support[0][row].squared_distance,.7,2.5);
        const auto second=n::EvaluateKernel(near_singular.support[1][row].squared_distance,.7,2.5);
        near_singular.observations[row]=1.0*first.gaussian+.2*first.charge+
            1.4*second.gaussian-.1*second.charge;
    }
    fixtures.push_back({"near-singular-overlapping-support",
        "fc9e6c0ec7f746d9de6221d1b680334bd870d7797012a666c32cdcdb467454ef",
        std::move(near_singular),2});
    auto low_signal=second_stage_test::SyntheticJointWorkload("chain",3);
    for(std::size_t row=0;row<low_signal.observations.size();++row)
        low_signal.observations[row]=1e-5*std::sin(.17*static_cast<double>(row));
    fixtures.push_back({"low-signal-noise-only",
        "0d2ecf601cb21d6ce3d4e23f1aa7f37b78f17fe7ec58a19c49b397374ad99d2f",
        std::move(low_signal),2});
    Input weak_width;
    weak_width.atom_ids={"1","2"};
    for(Eigen::Index row=0;row<12;++row)
    {
        weak_width.row_ids.push_back(std::to_string(row));
        const double first_square=.22*static_cast<double>(row);
        const double second_radius=.12*static_cast<double>(row)+.31;
        weak_width.support.resize(2);
        weak_width.support[0].push_back({static_cast<std::size_t>(row),first_square});
        weak_width.support[1].push_back({static_cast<std::size_t>(row),second_radius*second_radius});
    }
    const std::array<double,2> truth_widths{.72,1.08},truth_coefficients{1.25,.19};
    const std::array<double,2> second_coefficients{.83,.11};
    weak_width.observations.resize(12);
    for(Eigen::Index row=0;row<12;++row)
    {
        const double first_square=weak_width.support[0][static_cast<std::size_t>(row)].squared_distance;
        const double second_square=weak_width.support[1][static_cast<std::size_t>(row)].squared_distance;
        const auto first=n::EvaluateKernel(first_square,truth_widths[0],2.5);
        const auto second=n::EvaluateKernel(second_square,truth_widths[1],2.5);
        weak_width.observations[static_cast<std::size_t>(row)]=truth_coefficients[0]*first.gaussian+
            truth_coefficients[1]*first.charge+second_coefficients[0]*second.gaussian+
            second_coefficients[1]*second.charge;
    }
    fixtures.push_back({"weak-width-direction",
        "aa1d50f8e6dbcdaa22e401be857e849573ca0c6457d714583029003f23f9ef9c",
        std::move(weak_width),2});


    for(const auto count:{3,30,300,1000})
        fixtures.push_back({"chain-"+std::to_string(count),"",
            second_stage_test::SyntheticJointWorkload("chain",count),12});
    auto active=second_stage_test::SyntheticJointWorkload("chain",3);
    fixtures.push_back({"active-width","",std::move(active),2});
    for(auto & fixture:fixtures)
    {
        fixture.input.selection_domain.emplace();
        auto & targets=fixture.input.selection_domain->target_indices;
        targets.resize(fixture.input.atom_ids.size());
        std::iota(targets.begin(),targets.end(),0);
        const auto digest=second_stage_test::SyntheticJointWorkloadHash(fixture.input);
        const rhbm_gem::core::JointProblem problem(std::move(fixture.input));
        const auto & data=rhbm_gem::core::JointProblemAccess::Get(problem);
        n::Vector eta=n::Vector::Constant(static_cast<Eigen::Index>(data.domain.atoms.size()),
            std::log(fixture.name=="active-width" ? .4 : .55));
        if(fixture.name=="weak-width-direction") eta<<std::log(.72),std::log(1.08);
        n::FixedNeighborPolicy policy;
        policy.core_atoms=fixture.core_atoms; policy.maximum_sweeps=100;
        policy.convergence_criteria=n::PracticalJointConvergenceCriteria;
        policy.eta_upper_bound=std::log(fixture.name=="active-width" ? .45 : 2.2);
        std::size_t accepted{},unchanged{},recovery_attempts{},recovery_successes{};
        policy.production_progress=[&](const auto & sweep) {
            accepted+=sweep.accepted_blocks; unchanged+=sweep.unchanged_blocks;
            recovery_attempts+=sweep.recovery_attempts; recovery_successes+=sweep.recovery_successes;
        };
        const auto started=std::chrono::steady_clock::now();
        const auto cpu_started=std::clock();
        const auto result=n::SearchFixedNeighborComponent(data.domain,data.layout,data.y,eta,data.context,policy);
        const double wall=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
        const double cpu=static_cast<double>(std::clock()-cpu_started)/CLOCKS_PER_SEC;
        rusage usage{}; getrusage(RUSAGE_SELF,&usage);
#ifdef __APPLE__
        const auto peak_rss=usage.ru_maxrss;
#else
        const auto peak_rss=usage.ru_maxrss*1024;
#endif
        const auto & state=result.state;
        const bool qualified=result.certification && result.certification->qualified();
        const bool replay=result.certification && result.certification->assessment.independent_replay.passed;
        std::optional<JointState> point;
        if(qualified) point=JointState{{state.beta.data(),state.beta.data()+state.beta.size()},
            {},{state.eta.data(),state.eta.data()+state.eta.size()},{},state.objective};
        if(point) for(Eigen::Index k=0;k<state.eta.size();++k) point->b.push_back(std::exp(state.eta(k)));
        const auto evidence=n::AssessTargetPointState(problem,data.layout,point,
            qualified ? JointCheckStatus::Passed : JointCheckStatus::Unavailable,
            result.local_conditioning,JointEvidenceScope::AssembledGlobal,data.input->observations.size());
        j::array points;
        for(const auto & value:evidence.point_atoms)
            points.emplace_back(Object{{"atom",value.atom},{"quality",value.point_quality},
                {"conditioning",value.local_conditioning},
                {"convergence",static_cast<int>(value.component_convergence)},
                {"validity",static_cast<int>(value.parameter_validity)}});
        j::array widths;
        for(Eigen::Index k=0;k<state.eta.size();++k) widths.push_back(std::exp(state.eta(k)));
        Object record{{"wall_seconds",wall},{"cpu_seconds",cpu},{"peak_rss_bytes",peak_rss},
            {"fixture",fixture.name},{"input_digest",digest},
            {"eta",Values(state.eta)},{"beta",Values(state.beta)},{"b",std::move(widths)},
            {"prediction",Values(state.prediction)},{"residual",Values(state.residual)},
            {"objective",state.objective},{"reason",result.reason},{"sweeps",result.sweep_count},
            {"accepted",accepted},{"unchanged",unchanged},{"recovery_attempts",recovery_attempts},
            {"recovery_successes",recovery_successes},{"integrity",result.accepted_state_integrity},
            {"qualified",qualified},{"replay",replay},{"conditioning",BlockConditioning(result.local_conditioning)},
            {"points",std::move(points)},{"target_convergence",static_cast<int>(n::TargetConvergence(evidence,point.has_value()))}};
        if(result.final_sweep)
        {
            record["ac_kkt"]=result.final_sweep->global_ac_kkt;
            record["width_kkt"]=result.final_sweep->projected_width_gradient_inf_norm;
        }
        output<<j::serialize(record)<<'\n';
        std::cerr<<fixture.name<<": "<<result.reason<<" qualified="<<qualified<<'\n';
    }
}
