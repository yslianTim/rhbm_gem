#include "core/detail/FirstStageInitialization.hpp"
#include "core/detail/joint_component/Numerics.hpp"
#include "core/detail/joint_component/FixedNeighborBlockCoordinate.hpp"
#include "core/detail/joint_component/Problem.hpp"
#include "core/detail/joint_component/StructuralPartition.hpp"
#include "support/DynamicHaloMean.hpp"

#include <rhbm_gem/data/io/ModelMapFileIO.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <boost/json.hpp>
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <sys/resource.h>

namespace {
namespace core=rhbm_gem::core;
namespace joint=core::joint_component;
namespace dyn=joint_test::dynamic_halo;
namespace json=boost::json;
using Clock=std::chrono::steady_clock;
using Matrix=Eigen::MatrixXd;
using Vector=Eigen::VectorXd;

struct Truth {dyn::Parameters value;};
struct Profile
{
    Vector eta,beta,residual;
    Matrix jacobian;
    double objective{},gradient{},kkt{};
    bool valid{};
    std::string reason;
};
struct GroupResult
{
    std::string name;
    dyn::Mode mode{dyn::Mode::ABC};
    dyn::FullState state;
    std::size_t sweeps{},profile_evaluations{},accepted{},rejected{},failed_blocks{};
    double ac_kkt{std::numeric_limits<double>::infinity()},width_gradient{std::numeric_limits<double>::infinity()},
        parameter_change{std::numeric_limits<double>::infinity()},eta_change{std::numeric_limits<double>::infinity()},
        halo_mean_error{std::numeric_limits<double>::infinity()},replay_difference{std::numeric_limits<double>::infinity()},
        replay_objective_difference{std::numeric_limits<double>::infinity()},
        search_seconds{};
    std::string stop_reason{"maximum-sweeps"};
    bool qualified{};
};

json::value ReadJson(const std::filesystem::path & path)
{
    std::ifstream input(path);
    if(!input) throw std::runtime_error("Cannot read "+path.string());
    json::parse_options options; options.numbers=json::number_precision::precise;
    return json::parse(std::string(std::istreambuf_iterator<char>(input),{}),{},options);
}

std::unordered_map<std::string,Truth> ReadTruth(const std::filesystem::path & path)
{
    const auto source=ReadJson(path);
    std::unordered_map<std::string,Truth> result;
    for(const auto & item:source.at("atoms").as_array())
    {
        const auto id=std::to_string(json::value_to<std::uint64_t>(item.at("serial_id")));
        const double b=json::value_to<double>(item.at("effective_gaussian_width"));
        if(b!=json::value_to<double>(item.at("effective_charge_width")))
            throw std::runtime_error("Simulation truth has unequal Gaussian and charge widths at "+id);
        result.emplace(id,Truth{{static_cast<double>(json::value_to<int>(item.at("element"))),b,
            json::value_to<double>(item.at("charge_used"))}});
    }
    return result;
}

std::vector<dyn::Block> BuildBlocks(const core::JointProblem & problem)
{
    const auto & data=core::JointProblemAccess::Get(problem);
    std::vector<dyn::Block> blocks;
    for(std::size_t component=0;component<data.partition.components.size();++component)
    {
        const auto view=joint::MakeComponentProblemView(data.layout,data.partition.components[component]);
        joint::StructuralTopology topology(view.domain,view.atom_ids,view.local_layout);
        const auto partition=topology.BuildCores(12,false);
        for(const auto & core:partition.cores)
        {
            dyn::Block block; block.id=core.id; block.component=component;
            for(const auto local:core.atoms)
            {
                const auto local_atom=static_cast<std::size_t>(view.local_to_parent_atoms.at(static_cast<std::size_t>(local)));
                block.targets.push_back(data.target_to_contributor.at(local_atom));
            }
            blocks.push_back(std::move(block));
        }
    }
    return blocks;
}

dyn::Parameters FitUnits(dyn::Parameters value,double scale)
{value.a*=scale; value.c*=scale; return value;}

double Objective(const std::vector<double> & prediction,const core::JointProblemInput & input,double scale)
{
    long double sum{};
    for(std::size_t row=0;row<prediction.size();++row)
    {
        const long double residual=static_cast<long double>(prediction[row])-input.observations[row];
        sum+=residual*residual;
    }
    return static_cast<double>(0.5L*sum/(static_cast<long double>(scale)*scale));
}

double DenseProjectedKkt(const Matrix & design,joint::VectorRef response,const Vector & beta,double scale)
{
    const Vector residual=design*beta-response;
    const Vector norms=design.colwise().norm();
    if(!residual.allFinite() || !norms.allFinite() || (norms.array()<=0.0).any())
        return std::numeric_limits<double>::infinity();
    const Vector u=norms.array()*beta.array()/scale;
    const Vector gradient=design.transpose()*residual;
    double maximum{};
    for(Eigen::Index column=0;column<beta.size();++column)
    {
        const double projected=u(column)-gradient(column)/(norms(column)*scale);
        const double feasible_projection=column%2==0?std::max(0.0,projected):projected;
        maximum=std::max(maximum,std::abs(u(column)-feasible_projection));
    }
    return maximum;
}

dyn::FullState IndependentState(const core::JointProblemInput & input,const dyn::Ownership & ownership,
    std::span<const dyn::Parameters> targets,dyn::Mode mode,double scale)
{
    dyn::FullState state;
    state.contributors.assign(targets.begin(),targets.end());
    if(mode!=dyn::Mode::StaticABC)
        for(const auto & owner:ownership.halos)
        {
            const auto & members=ownership.blocks[owner.owner_block].targets;
            dyn::Parameters mean{0.0,0.0,0.0};
            for(const auto target:members)
            {
                mean.a+=targets[target].a;
                mean.b+=targets[target].b;
                if(mode==dyn::Mode::ABC) mean.c+=targets[target].c;
            }
            const double count=static_cast<double>(members.size());
            state.contributors[owner.halo]={mean.a/count,mean.b/count,mode==dyn::Mode::ABC?mean.c/count:0.0};
        }
    state.prediction.assign(input.observations.size(),0.0);
    // Independent endpoint replay traverses contributor support directly, without using a block design.
    for(std::size_t atom=0;atom<state.contributors.size();++atom)
    {
        const auto & parameter=state.contributors[atom];
        for(const auto & support:input.support[atom])
        {
            const auto basis=joint::EvaluateKernel(support.squared_distance,parameter.b,2.5);
            state.prediction[support.row]+=parameter.a*basis.gaussian+parameter.c*basis.charge;
        }
    }
    state.objective=Objective(state.prediction,input,scale);
    return state;
}

Eigen::VectorXd Coefficients(std::span<const dyn::Parameters> state,const dyn::Block & block)
{
    Vector beta(static_cast<Eigen::Index>(2*block.targets.size()));
    for(std::size_t slot=0;slot<block.targets.size();++slot)
    {
        const auto & parameter=state[block.targets[slot]];
        beta(static_cast<Eigen::Index>(2*slot))=parameter.a;
        beta(static_cast<Eigen::Index>(2*slot+1))=parameter.c;
    }
    return beta;
}

Profile EvaluateProfile(const core::JointProblemInput & input,const dyn::Ownership & ownership,
    std::size_t block,std::span<const dyn::Parameters> state,joint::VectorRef response,Vector eta,
    dyn::Mode mode,double scale,const joint::EvaluationContext & context,std::size_t & evaluation_count)
{
    Profile result; ++evaluation_count; result.eta=std::move(eta);
    const auto & targets=ownership.blocks.at(block).targets;
    if(result.eta.size()!=static_cast<Eigen::Index>(targets.size()) || response.size()==0)
    {result.reason="invalid-profile-input"; return result;}
    std::vector<dyn::Parameters> trial_state(state.begin(),state.end());
    for(std::size_t slot=0;slot<targets.size();++slot)
    {
        const double width=std::exp(result.eta(static_cast<Eigen::Index>(slot)));
        if(!std::isfinite(width) || !(width>0.0)) {result.reason="invalid-width"; return result;}
        trial_state[targets[slot]].b=width;
    }
    const auto design=dyn::BuildAugmentedDesign(input,ownership,block,trial_state,mode);
    if(static_cast<Eigen::Index>(design.rows.size())!=response.size())
    {result.reason="affected-row-mismatch"; return result;}
    const Vector weights=Vector::Ones(response.size());
    const auto solved=joint::SolveLinear(design.values,response,weights,false,false,nullptr,&context.linear,scale);
    result.beta=solved.beta;
    result.kkt=DenseProjectedKkt(design.values,response,result.beta,scale);
    bool feasible=true;
    for(Eigen::Index column=0;column<result.beta.size();column+=2) feasible&=result.beta(column)>=0.0;
    if(!solved.valid || !std::isfinite(result.kkt) || result.kkt>1e-10 || !feasible)
    {result.reason="local-ac-kkt-failed:"+solved.reason; return result;}
    result.residual=design.values*result.beta-response;
    result.objective=0.5*result.residual.squaredNorm()/(scale*scale);
    const auto width_jacobian=dyn::CoupledWidthJacobian(input,ownership,block,trial_state,result.beta,mode);
    std::vector<Eigen::Index> free_columns;
    for(std::size_t slot=0;slot<targets.size();++slot)
    {
        if(result.beta(static_cast<Eigen::Index>(2*slot))>0.0)
            free_columns.push_back(static_cast<Eigen::Index>(2*slot));
        free_columns.push_back(static_cast<Eigen::Index>(2*slot+1));
    }
    Matrix free_design(design.values.rows(),static_cast<Eigen::Index>(free_columns.size()));
    for(std::size_t column=0;column<free_columns.size();++column)
        free_design.col(static_cast<Eigen::Index>(column))=design.values.col(free_columns[column]);
    Matrix projected=width_jacobian;
    if(free_design.cols())
        projected.noalias()-=free_design*free_design.colPivHouseholderQr().solve(width_jacobian);
    result.jacobian=std::move(projected);
    result.gradient=(result.jacobian.transpose()*result.residual).lpNorm<Eigen::Infinity>()/(scale*scale);
    result.valid=result.residual.allFinite() && result.jacobian.allFinite() && std::isfinite(result.gradient);
    result.reason=result.valid?"ok":"nonfinite-profile";
    return result;
}

Profile OptimizeBlock(const core::JointProblemInput & input,const dyn::Ownership & ownership,
    std::size_t block,std::span<const dyn::Parameters> state,joint::VectorRef response,dyn::Mode mode,double scale,
    const joint::EvaluationContext & context,std::size_t & evaluations)
{
    const auto & targets=ownership.blocks.at(block).targets;
    Vector eta(static_cast<Eigen::Index>(targets.size()));
    for(std::size_t slot=0;slot<targets.size();++slot)
        eta(static_cast<Eigen::Index>(slot))=std::log(state[targets[slot]].b);
    auto current=EvaluateProfile(input,ownership,block,state,response,eta,mode,scale,context,evaluations);
    if(!current.valid) return current;
    double damping=1e-3;
    for(std::size_t iteration=0;iteration<12;++iteration)
    {
        if(current.gradient<=1e-10) break;
        const Matrix hessian=current.jacobian.transpose()*current.jacobian;
        const Vector gradient=current.jacobian.transpose()*current.residual;
        Vector diagonal=hessian.diagonal().cwiseAbs();
        const double largest=diagonal.size()?diagonal.maxCoeff():0.0;
        if(!(largest>0.0) || !std::isfinite(largest)) break;
        diagonal=diagonal.cwiseMax(largest*1e-12);
        bool accepted{};
        for(std::size_t retry=0;retry<8;++retry)
        {
            Matrix damped=hessian;
            damped.diagonal().array()+=damping*diagonal.array();
            Vector step=damped.ldlt().solve(-gradient);
            if(!step.allFinite()) {damping*=10.0; continue;}
            const double infinity=step.lpNorm<Eigen::Infinity>();
            if(infinity>0.25) step*=0.25/infinity;
            auto trial=EvaluateProfile(input,ownership,block,state,response,current.eta+step,mode,scale,context,evaluations);
            if(trial.valid && trial.objective<current.objective)
            {
                current=std::move(trial); damping=std::max(1e-12,damping/3.0); accepted=true; break;
            }
            damping*=10.0;
        }
        if(!accepted) break;
    }
    return current;
}

struct Stationarity {double ac{},width{},amplitude_change{},eta_change{};};
Stationarity GlobalStationarity(const core::JointProblemInput & input,const dyn::Ownership & ownership,
    const dyn::FullState & state,dyn::Mode mode,double scale,std::span<const dyn::Parameters> previous)
{
    Vector residual(static_cast<Eigen::Index>(input.observations.size()));
    for(std::size_t row=0;row<input.observations.size();++row)
        residual(static_cast<Eigen::Index>(row))=state.prediction[row]-input.observations[row];
    Stationarity result;
    for(std::size_t block=0;block<ownership.blocks.size();++block)
    {
        const auto & targets=ownership.blocks[block].targets;
        std::vector<dyn::Parameters> members;
        for(const auto target:targets) members.push_back(state.contributors[target]);
        const auto mean=dyn::Mean(members,mode!=dyn::Mode::AB);
        for(const auto target:targets)
        {
            const auto & parameter=state.contributors[target];
            std::unordered_map<std::size_t,double> column_a,column_c;
            column_a.reserve(input.support[target].size());
            column_c.reserve(input.support[target].size());
            const auto add_amplitude=[&](std::size_t atom,double width,double share,bool halo) {
                for(const auto & support:input.support[atom])
                {
                    const auto basis=joint::EvaluateKernel(support.squared_distance,width,2.5);
                    column_a[support.row]+=share*basis.gaussian;
                    if(!halo || mode!=dyn::Mode::AB) column_c[support.row]+=share*basis.charge;
                }
            };
            add_amplitude(target,parameter.b,1.0,false);
            const double share=1.0/static_cast<double>(targets.size());
            if(mode!=dyn::Mode::StaticABC)
                for(const auto halo:ownership.block_halos[block]) add_amplitude(halo,mean.b,share,true);
            double squared_a{},dot_a{},squared_c{},dot_c{},width_dot{};
            for(const auto & [row,value]:column_a)
            {
                squared_a+=value*value;
                dot_a+=value*residual(static_cast<Eigen::Index>(row));
            }
            for(const auto & [row,value]:column_c)
            {
                squared_c+=value*value;
                dot_c+=value*residual(static_cast<Eigen::Index>(row));
            }
            const double norm_a=std::sqrt(squared_a),norm_c=std::sqrt(squared_c);
            if(norm_a>0.0)
            {
                const double u=norm_a*parameter.a/scale;
                const double grad=dot_a/(norm_a*scale);
                const double projected=std::max(0.0,u-grad);
                result.ac=std::max(result.ac,std::abs(u-projected));
                const auto old=previous[target];
                result.amplitude_change=std::max(result.amplitude_change,std::abs(parameter.a-old.a)*norm_a/scale);
            }
            if(norm_c>0.0)
            {
                const double grad=dot_c/(norm_c*scale);
                result.ac=std::max(result.ac,std::abs(grad));
                const auto old=previous[target];
                result.amplitude_change=std::max(result.amplitude_change,std::abs(parameter.c-old.c)*norm_c/scale);
            }
            for(const auto & support:input.support[target])
            {
                const auto basis=joint::EvaluateKernel(support.squared_distance,parameter.b,2.5);
                width_dot+=(parameter.a*basis.gaussian_log_width+parameter.c*basis.charge_log_width)*
                    residual(static_cast<Eigen::Index>(support.row));
            }
            if(mode!=dyn::Mode::StaticABC)
            {
                const double chain=parameter.b*share/mean.b;
                for(const auto halo:ownership.block_halos[block])
                    for(const auto & support:input.support[halo])
                    {
                        const auto basis=joint::EvaluateKernel(support.squared_distance,mean.b,2.5);
                        const double charge=mode==dyn::Mode::AB ? 0.0 : mean.c*basis.charge_log_width;
                        width_dot+=chain*(mean.a*basis.gaussian_log_width+charge)*
                            residual(static_cast<Eigen::Index>(support.row));
                    }
            }
            result.width=std::max(result.width,std::abs(width_dot)/(scale*scale));
            result.eta_change=std::max(result.eta_change,std::abs(std::log(parameter.b/previous[target].b)));
        }
    }
    return result;
}

double MeanConsistency(const dyn::Ownership & ownership,const dyn::FullState & state,dyn::Mode mode,
    std::span<const dyn::Parameters> fixed_halos)
{
    double maximum{};
    for(const auto & owner:ownership.halos)
    {
        if(mode==dyn::Mode::StaticABC)
        {
            const auto actual=state.contributors[owner.halo],reference=fixed_halos[owner.halo];
            maximum=std::max({maximum,std::abs(actual.a-reference.a),std::abs(actual.b-reference.b),
                std::abs(actual.c-reference.c)});
            continue;
        }
        std::vector<dyn::Parameters> targets;
        for(const auto target:ownership.blocks[owner.owner_block].targets) targets.push_back(state.contributors[target]);
        const auto mean=dyn::Mean(targets,mode!=dyn::Mode::AB);
        const auto actual=state.contributors[owner.halo];
        maximum=std::max({maximum,std::abs(mean.a-actual.a),std::abs(mean.b-actual.b),std::abs(mean.c-actual.c)});
    }
    return maximum;
}

GroupResult RunGroup(const core::JointProblemInput & input,const dyn::Ownership & ownership,
    const std::vector<dyn::Parameters> & initial,dyn::Mode mode,std::string name,double scale,
    const joint::EvaluationContext & context,std::size_t maximum_sweeps)
{
    GroupResult result; result.name=std::move(name); result.mode=mode;
    result.state=IndependentState(input,ownership,initial,mode,scale);
    const auto started=Clock::now();
    std::size_t stable_sweeps{};
    for(std::size_t sweep=0;sweep<maximum_sweeps;++sweep)
    {
        const auto sweep_state=result.state.contributors;
        for(std::size_t block=0;block<ownership.blocks.size();++block)
        {
            const auto & targets=ownership.blocks[block].targets;
            const auto & rows=ownership.block_dynamic_rows[block];
            auto design=dyn::BuildAugmentedDesign(input,ownership,block,result.state.contributors,mode);
            Vector beta=Coefficients(result.state.contributors,ownership.blocks[block]);
            const Vector local_current=design.values*beta;
            Vector response(static_cast<Eigen::Index>(rows.size()));
            for(std::size_t local=0;local<rows.size();++local)
                response(static_cast<Eigen::Index>(local))=input.observations[rows[local]]-
                    result.state.prediction[rows[local]]+local_current(static_cast<Eigen::Index>(local));
            auto fit=OptimizeBlock(input,ownership,block,result.state.contributors,response,mode,scale,context,
                result.profile_evaluations);
            if(!fit.valid) {++result.failed_blocks; continue;}
            const auto profile_residual=design.values*beta-response;
            const double local_before=0.5*profile_residual.squaredNorm()/(scale*scale);
            std::vector<dyn::Parameters> candidate_targets=result.state.contributors;
            for(std::size_t slot=0;slot<targets.size();++slot)
            {
                auto & parameter=candidate_targets[targets[slot]];
                parameter.a=fit.beta(static_cast<Eigen::Index>(2*slot));
                parameter.c=fit.beta(static_cast<Eigen::Index>(2*slot+1));
                parameter.b=std::exp(fit.eta(static_cast<Eigen::Index>(slot)));
            }
            auto candidate=mode==dyn::Mode::StaticABC ?
                IndependentState(input,ownership,candidate_targets,mode,scale) :
                dyn::BuildState(input,ownership,candidate_targets,mode,scale);
            const auto candidate_design=dyn::BuildAugmentedDesign(input,ownership,block,candidate_targets,mode);
            const auto candidate_local=candidate_design.values*fit.beta-response;
            const double local_after=0.5*candidate_local.squaredNorm()/(scale*scale);
            const double replay_error=std::abs((candidate.objective-result.state.objective)-(local_after-local_before));
            const double replay_bound=joint::BlockObjectiveReplayEnclosure(result.state.objective);
            if(replay_error>replay_bound)
            {++result.rejected; continue;}
            if(dyn::CommitMonotone(result.state,std::move(candidate))) ++result.accepted;
            else ++result.rejected;
        }
        result.sweeps=sweep+1;
        const auto independently_replayed=IndependentState(input,ownership,result.state.contributors,mode,scale);
        result.replay_difference=0.0;
        for(std::size_t row=0;row<result.state.prediction.size();++row)
            result.replay_difference=std::max(result.replay_difference,
                std::abs(result.state.prediction[row]-independently_replayed.prediction[row]));
        result.replay_objective_difference=std::abs(result.state.objective-independently_replayed.objective);
        const auto stationarity=GlobalStationarity(input,ownership,result.state,mode,scale,sweep_state);
        result.ac_kkt=stationarity.ac; result.width_gradient=stationarity.width;
        result.parameter_change=std::max(stationarity.amplitude_change,stationarity.eta_change);
        result.eta_change=stationarity.eta_change;
        result.halo_mean_error=MeanConsistency(ownership,result.state,mode,initial);
        const bool passed=result.ac_kkt<=joint::PracticalJointConvergenceCriteria.global_ac_kkt &&
            result.width_gradient<=joint::PracticalJointConvergenceCriteria.global_width_gradient &&
            result.eta_change<=joint::PracticalJointConvergenceCriteria.confirmed_log_width_change &&
            stationarity.amplitude_change<=joint::PracticalJointConvergenceCriteria.confirmed_log_width_change &&
            result.halo_mean_error<=1e-12 && result.replay_difference<=1e-12 &&
            result.replay_objective_difference<=1e-12;
        stable_sweeps=passed?stable_sweeps+1:0;
        if(stable_sweeps>=2)
        {result.qualified=true; result.stop_reason="coupled-practical-v1-confirmed"; break;}
    }
    result.search_seconds=std::chrono::duration<double>(Clock::now()-started).count();
    return result;
}

struct ErrorMetrics {double mae{},rmse{},median{},p90{},p95{},p99{},maximum{};};
ErrorMetrics Errors(std::vector<double> errors)
{
    for(auto & error:errors) error=std::abs(error);
    std::sort(errors.begin(),errors.end());
    long double sum{},squared{};
    for(const auto error:errors) {sum+=error; squared+=static_cast<long double>(error)*error;}
    const auto q=[&](double p){return errors[static_cast<std::size_t>(
        std::ceil(p*static_cast<double>(errors.size()))-1)];};
    return {static_cast<double>(sum/errors.size()),std::sqrt(static_cast<double>(squared/errors.size())),
        errors[(errors.size()-1)/2],q(.9),q(.95),q(.99),errors.back()};
}
json::object ErrorJson(const ErrorMetrics & e)
{return {{"mae",e.mae},{"rmse",e.rmse},{"median_absolute_error",e.median},{"p90_absolute_error",e.p90},
    {"p95_absolute_error",e.p95},{"p99_absolute_error",e.p99},{"maximum_absolute_error",e.maximum}};}

double PeakRssBytes()
{
    struct rusage usage{};
    if(getrusage(RUSAGE_SELF,&usage)!=0) return -1.0;
#if defined(__APPLE__)
    return static_cast<double>(usage.ru_maxrss);
#else
    return static_cast<double>(usage.ru_maxrss)*1024.0;
#endif
}

json::object SummarizeGroup(const core::JointProblemInput & input,const dyn::Ownership & ownership,
    const GroupResult & group,const std::unordered_map<std::string,Truth> & truth,double amplitude_scale,double scale)
{
    std::array<std::vector<double>,3> target_errors,halo_errors;
    std::vector<double> target_prediction(input.observations.size()),halo_prediction(input.observations.size()),
        truth_target(input.observations.size()),truth_halo(input.observations.size());
    const auto target_indices=input.selection_domain->target_indices;
    std::vector<bool> is_target(input.atom_ids.size());
    for(const auto target:target_indices) is_target[target]=true;
    std::size_t target_count{};
    for(std::size_t atom=0;atom<input.atom_ids.size();++atom)
    {
        const auto & p=group.state.contributors[atom];
        const auto it=truth.find(input.atom_ids[atom]);
        if(it==truth.end()) throw std::runtime_error("Missing truth for contributor "+input.atom_ids[atom]);
        const auto t=it->second.value;
        if(is_target[atom])
        {
            ++target_count;
            target_errors[0].push_back(p.a/amplitude_scale-t.a);
            target_errors[1].push_back(p.b-t.b);
            target_errors[2].push_back(p.c/amplitude_scale-t.c);
        }
        else
        {
            halo_errors[0].push_back(p.a/amplitude_scale-t.a);
            halo_errors[1].push_back(p.b-t.b);
            halo_errors[2].push_back(p.c/amplitude_scale-t.c);
        }
        for(const auto & support:input.support[atom])
        {
            const auto basis=joint::EvaluateKernel(support.squared_distance,p.b,2.5);
            const double value=p.a*basis.gaussian+p.c*basis.charge;
            if(is_target[atom]) target_prediction[support.row]+=value; else halo_prediction[support.row]+=value;
            const auto tp=FitUnits(t,amplitude_scale);
            const auto truth_basis=joint::EvaluateKernel(support.squared_distance,tp.b,2.5);
            const double truth_value=tp.a*truth_basis.gaussian+tp.c*truth_basis.charge;
            if(is_target[atom]) truth_target[support.row]+=truth_value; else truth_halo[support.row]+=truth_value;
        }
    }
    const auto field_rmse=[&](const std::vector<double> & p,const std::vector<double> & t) {
        long double sum{};
        for(std::size_t row=0;row<p.size();++row)
        {const long double d=static_cast<long double>(p[row])-t[row]; sum+=d*d;}
        return std::sqrt(static_cast<double>(sum/static_cast<long double>(p.size())));
    };
    std::vector<double> full(group.state.prediction.size()),full_truth(full.size());
    for(std::size_t row=0;row<full.size();++row)
    {full[row]=target_prediction[row]+halo_prediction[row]; full_truth[row]=truth_target[row]+truth_halo[row];}
    return {{"group",group.name},{"status",group.qualified?"qualified":"finite-candidate-unqualified"},
        {"qualification_reason",group.stop_reason},{"target_count",target_count},{"halo_count",ownership.halos.size()},
        {"target_A_raw_units",ErrorJson(Errors(target_errors[0]))},{"target_B_angstrom",ErrorJson(Errors(target_errors[1]))},
        {"target_C_raw_units",ErrorJson(Errors(target_errors[2]))},{"halo_A_raw_units",ErrorJson(Errors(halo_errors[0]))},
        {"halo_B_angstrom",ErrorJson(Errors(halo_errors[1]))},{"halo_C_raw_units",ErrorJson(Errors(halo_errors[2]))},
        {"halo_prediction_rmse_raw_map_units",field_rmse(halo_prediction,truth_halo)*amplitude_scale},
        {"target_prediction_rmse_raw_map_units",field_rmse(target_prediction,truth_target)*amplitude_scale},
        {"full_prediction_rmse_raw_map_units",field_rmse(full,full_truth)*amplitude_scale},
        {"objective",Objective(group.state.prediction,input,scale)},{"full_residual_norm",[&] {
            long double sum{}; for(std::size_t row=0;row<full.size();++row)
            {const long double d=static_cast<long double>(group.state.prediction[row])-input.observations[row];sum+=d*d;}
            return std::sqrt(static_cast<double>(sum));}()},{"sweeps",group.sweeps},
        {"profile_evaluations",group.profile_evaluations},{"accepted_transactions",group.accepted},
        {"rejected_transactions",group.rejected},{"failed_blocks",group.failed_blocks},
        {"global_coupled_ac_kkt",group.ac_kkt},{"global_coupled_width_gradient",group.width_gradient},
        {"complete_sweep_log_width_change",group.eta_change},{"complete_sweep_scaled_amplitude_change",group.parameter_change},
        {"halo_mean_consistency_error",group.halo_mean_error},{"independent_replay_max_prediction_difference",group.replay_difference},
        {"independent_replay_objective_difference",group.replay_objective_difference},
        {"search_seconds",group.search_seconds}};
}

int Run(const std::filesystem::path & model_path,const std::filesystem::path & map_path,
    const std::filesystem::path & baseline_summary,const std::filesystem::path & output_dir,
    const std::string & baseline_commit,std::size_t maximum_sweeps,const std::string & group_filter)
{
    auto model=rhbm_gem::ReadModel(model_path.string()); auto map=rhbm_gem::ReadMap(map_path.string());
    const auto divisor=map->MapValueArrayNormalization();
    if(!divisor) throw std::runtime_error("Map normalization is unavailable.");
    const double amplitude_scale=1.0/(*divisor);
    model->SelectAtoms([](const rhbm_gem::AtomObject & atom) {
        return atom.GetElement()!=Element::HYDROGEN && atom.GetChainID()=="A";
    });
    const auto problem=core::BuildJointProblem(*map,*model);
    const auto initialization=joint::InitializeContributors(*map,*model,problem);
    if(!initialization.valid) throw std::runtime_error("FirstStage initialization is not qualified.");
    const auto & input=problem.Input();
    const auto ownership=dyn::BuildOwnership(input,BuildBlocks(problem));
    if(ownership.cross_component_halos) throw std::runtime_error("unsupported-cross-component-halo-coupling");
    const auto truth=ReadTruth(map_path.string()+".simulation.json");
    const double scale=problem.ObservationScale();
    std::vector<dyn::Parameters> initial(input.atom_ids.size());
    for(std::size_t atom=0;atom<input.atom_ids.size();++atom)
    {
        if(initialization.atoms[atom].id!=input.atom_ids[atom])
            throw std::runtime_error("Contributor initialization identity mismatch.");
        initial[atom]={0.0,initialization.b[atom],0.0};
    }
    const auto context=joint::CreateContext(std::make_shared<core::JointProblemInput>(input));
    const bool run_all=group_filter=="all";
    std::vector<GroupResult> groups;
    if(run_all || group_filter=="E1")
        groups.push_back(RunGroup(input,ownership,initial,dyn::Mode::AB,"E1 Dynamic AB",scale,context,maximum_sweeps));
    if(run_all || group_filter=="E2")
        groups.push_back(RunGroup(input,ownership,initial,dyn::Mode::ABC,"E2 Dynamic ABC",scale,context,maximum_sweeps));
    if(run_all || group_filter=="E3")
    {
        const auto static_initial=dyn::DeriveHalos(input,ownership,initial,dyn::Mode::StaticABC);
        groups.push_back(RunGroup(input,ownership,static_initial,dyn::Mode::StaticABC,
            "E3 Static Mean ABC",scale,context,maximum_sweeps));
    }
    if(groups.empty()) throw std::invalid_argument("Group filter must be all, E1, E2, or E3.");
    const auto baseline=ReadJson(baseline_summary);
    const auto & baseline_p1=baseline.at("p1");
    const auto & e0_prediction=baseline_p1.at("models").as_array().front();
    const auto & ownership_summary=baseline.at("ownership");
    constexpr std::size_t large_observations=5'938'087,large_contributors=37'406,core_targets=12;
    const auto large_global_design=dyn::DenseBytes(large_observations,2*large_contributors).value();
    const auto observed_local_design=dyn::DenseBytes(
        json::value_to<std::size_t>(ownership_summary.at("maximum_dynamic_block_rows")),2*core_targets).value();
    const auto two_global_buffers=dyn::DenseBytes(large_observations,2).value();
    const auto row_block_headers=large_observations*sizeof(std::vector<std::size_t>);
    const auto add_improvements=[&](json::object & group) {
        const auto improvement=[&](const char * key,const char * baseline_key) {
            const double reference=json::value_to<double>(baseline_p1.at(baseline_key));
            const double candidate=json::value_to<double>(group.at(key).at("rmse"));
            return 100.0*(reference-candidate)/reference;
        };
        group["target_A_improvement_percent_vs_E0"]=improvement("target_A_raw_units","target_rmse_A_raw_units");
        group["target_B_improvement_percent_vs_E0"]=improvement("target_B_angstrom","target_rmse_B_angstrom");
        group["target_C_improvement_percent_vs_E0"]=improvement("target_C_raw_units","target_rmse_C_raw_units");
    };
    json::array group_summaries,qualified_groups;
    for(const auto & group:groups)
    {
        auto group_summary=SummarizeGroup(input,ownership,group,truth,amplitude_scale,scale);
        add_improvements(group_summary);
        if(group.qualified) qualified_groups.emplace_back(group.name);
        group_summaries.emplace_back(std::move(group_summary));
    }
    std::filesystem::create_directories(output_dir);
    json::object summary{{"contract","joint-dynamic-halo-mean-coupled-candidate-v1"},
        {"baseline_commit",baseline_commit},{"model_sha256",baseline.at("model_sha256")},
        {"map_sha256",baseline.at("map_sha256")},
        {"simulation_truth_sha256",baseline.at("simulation_truth_sha256")},
        {"normalization_divisor",*divisor},
        {"observation_scale",scale},{"observation_rows",input.observations.size()},
        {"target_atoms",input.selection_domain->target_indices.size()},{"halo_atoms",ownership.halos.size()},
        {"block_count",ownership.blocks.size()},{"maximum_sweeps",maximum_sweeps},
        {"group_filter",group_filter},{"block_order","forward"},
        {"starting_target_state","A=C=0, B=qualified FirstStage width for every target"},
        {"E0_production",json::object{{"status",baseline.at("p0").at("convergence_status")},
            {"objective",e0_prediction.at("objective")},
            {"target_A_rmse_raw_units",baseline_p1.at("target_rmse_A_raw_units")},
            {"target_B_rmse_angstrom",baseline_p1.at("target_rmse_B_angstrom")},
            {"target_C_rmse_raw_units",baseline_p1.at("target_rmse_C_raw_units")},
            {"halo_prediction_rmse_raw_map_units",e0_prediction.at("halo_field_rmse_raw_map_units")},
            {"target_prediction_rmse_raw_map_units",e0_prediction.at("target_prediction_rmse_raw_map_units")},
            {"full_prediction_rmse_raw_map_units",e0_prediction.at("full_prediction_rmse_raw_map_units")},
            {"sweeps",baseline.at("p0").at("sweeps")},
            {"profile_evaluations",baseline.at("p0").at("profile_evaluations")},
            {"accepted_updates",baseline.at("p0").at("accepted_updates")}}},
        {"initialization_truth_used",false},{"ownership_truth_used",false},
        {"halo_mean_scope","whole owner core; physical arithmetic mean of A/B and, for ABC, C"},
        {"ownership",json::object{{"shared_halos",ownership_summary.at("shared_halos")},
            {"cross_component_halos",ownership_summary.at("cross_component_halos")},
            {"row_to_block_memberships",ownership_summary.at("row_to_block_memberships")},
            {"target_block_rows",ownership_summary.at("total_target_block_rows")},
            {"dynamic_affected_rows",ownership_summary.at("total_dynamic_affected_rows")},
            {"maximum_dynamic_block_rows",ownership_summary.at("maximum_dynamic_block_rows")},
            {"maximum_row_expansion",ownership_summary.at("maximum_affected_row_expansion")}}},
        {"groups",std::move(group_summaries)},
        {"qualified_groups",qualified_groups},
        {"memory",json::object{{"peak_rss_bytes",PeakRssBytes()},
            {"local_design_dense","per block only: dynamic affected rows × 2 × core target count"},
            {"global_design_dense",false},
            {"observed_6z6u_local_design_bytes_at_max_rows",observed_local_design},
            {"large_component_arithmetic_preflight",json::object{
                {"observations",large_observations},{"contributors",large_contributors},
                {"forbidden_global_dense_design_bytes",large_global_design},
                {"two_prediction_or_residual_buffers_bytes",two_global_buffers},
                {"row_to_block_vector_headers_bytes",row_block_headers},
                {"status","arithmetic-only; large input unavailable; excludes input/support/allocator storage"}}}}}};
    std::ofstream out(output_dir/"summary.json"); out.exceptions(std::ios::failbit|std::ios::badbit);
    out<<json::serialize(summary)<<'\n';
    std::ofstream params(output_dir/"targets.csv"); params.exceptions(std::ios::failbit|std::ios::badbit);
    params<<"group,contributor_index,atom_id,target_block,A,B,C,truth_A,truth_B,truth_C\n"<<std::setprecision(17);
    for(const auto & group:groups)
    {
        std::vector<std::string> blocks(input.atom_ids.size());
        for(const auto & block:ownership.blocks) for(const auto target:block.targets) blocks[target]=block.id;
        for(const auto target:input.selection_domain->target_indices)
        {
            const auto & value=group.state.contributors[target];
            const auto truth_value=FitUnits(truth.at(input.atom_ids[target]).value,amplitude_scale);
            params<<group.name<<','<<target<<','<<input.atom_ids[target]<<','<<blocks[target]<<','
                <<value.a/amplitude_scale<<','<<value.b<<','<<value.c/amplitude_scale<<','
                <<truth_value.a<<','<<truth_value.b<<','<<truth_value.c<<'\n';
        }
    }
    std::cout<<json::serialize(json::object{{"summary",(output_dir/"summary.json").string()},
        {"observations",input.observations.size()},{"targets",input.selection_domain->target_indices.size()},
        {"halos",ownership.halos.size()},{"blocks",ownership.blocks.size()},
        {"group_filter",group_filter},{"qualified_groups",qualified_groups}})<<'\n';
    return 0;
}
}

int main(int argc,char ** argv)
{
    try
    {
        if(argc<6 || argc>8) throw std::invalid_argument(
            "Usage: joint_dynamic_halo_joint MODEL MAP P1_SUMMARY OUTPUT_DIR BASELINE_COMMIT [MAX_SWEEPS] [all|E1|E2|E3]");
        const auto sweeps=argc>=7?static_cast<std::size_t>(std::stoul(argv[6])):79u;
        const auto group_filter=argc==8?argv[7]:"all";
        return Run(argv[1],argv[2],argv[3],argv[4],argv[5],sweeps,group_filter);
    }
    catch(const std::exception & error) {std::cerr<<error.what()<<'\n'; return 1;}
}
