#include "support/JointHaloObservationExclusion.hpp"
#include "core/detail/joint_component/Problem.hpp"
#include <boost/json.hpp>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <numbers>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>
#ifndef _WIN32
#include <sys/resource.h>
#endif

namespace {
namespace json = boost::json;
namespace core = rhbm_gem::core;
namespace joint = core::joint_component;
namespace audit = joint_test::halo_exclusion;
using Input = core::JointProblemInput;
using Position = audit::Position;

constexpr std::array<double, 3> SimulationRadii{2.5, 3.5, 5.0};
constexpr std::array<double, 4> NoiseSigmas{0.0, 0.005, 0.02, 0.05};
constexpr std::array<std::uint64_t, 10> Seeds{1001,1002,1003,1004,1005,1006,1007,1008,1009,1010};
constexpr std::array<double, 6> Weights{0.0,0.01,0.05,0.1,0.25,1.0};
constexpr std::array<double, 3> TargetTruth{1.4, 0.72, 0.12};
constexpr std::array<double, 2> FirstStageHalo{0.58, 0.82};
constexpr Position HaloPosition{0.0, 0.0, 2.0};

struct Scenario
{
    Input input;
    std::vector<Position> positions;
    std::vector<double> target_truth, halo_truth, clean_observations, halo_correction;
};

struct Parameters { double a{}, b{}, c{}; };
struct Gradient { double a{}, b{}, c{}; };
struct Fit
{
    Parameters parameters;
    double objective{}, max_gradient_error{}, zero_row_delta{}, zero_row_gradient_delta{}, replay_objective_delta{};
};

json::value Number(double value)
{ return std::isfinite(value) ? json::value(value) : json::value(nullptr); }

void Write(const std::filesystem::path & path, const json::value & value)
{
    if (!path.parent_path().empty()) std::filesystem::create_directories(path.parent_path());
    const auto temporary = std::filesystem::path(path.string() + ".tmp");
    std::ofstream output(temporary);
    output.exceptions(std::ios::failbit | std::ios::badbit);
    output << json::serialize(value) << '\n';
    output.close();
    std::filesystem::rename(temporary, path);
}

joint::BasisValues SimulationBasis(double square, double width, double support_radius)
{
    if (square > support_radius * support_radius) return {};
    const double distance = std::sqrt(square);
    auto value = joint::EvaluateKernel(square, width, support_radius);
    value.charge = distance < 1e-5 ? std::sqrt(2.0 / std::numbers::pi) / width :
        std::erf(distance / width / std::sqrt(2.0)) / distance;
    return value;
}

Scenario BuildScenario(double simulation_radius)
{
    Scenario scenario;
    scenario.input.atom_ids = {"target", "halo"};
    scenario.input.support.resize(2);
    scenario.input.selection_domain.emplace();
    scenario.input.selection_domain->target_indices = {0};
    scenario.input.selection_domain->contributor_policy = std::string(rhbm_gem::FixedHaloContributorPolicy);
    scenario.input.selection_domain->fixed_halo_selection = rhbm_gem::JointFixedHaloSelectionProvenance{1,1,1,0};
    for (int z = -5; z <= 5; ++z)
        for (int y = -5; y <= 5; ++y)
            for (int x = -5; x <= 5; ++x)
            {
                const Position position{0.5 * x, 0.5 * y, 0.5 * z};
                const double target_square = position[0]*position[0] + position[1]*position[1] + position[2]*position[2];
                if (target_square > 6.25) continue;
                const auto row = scenario.positions.size();
                scenario.positions.push_back(position);
                scenario.input.row_ids.push_back(std::to_string(x)+","+std::to_string(y)+","+std::to_string(z));
                scenario.input.support[0].push_back({row, target_square});
                const double dx=position[0]-HaloPosition[0], dy=position[1]-HaloPosition[1], dz=position[2]-HaloPosition[2];
                const double halo_square=dx*dx+dy*dy+dz*dz;
                if (halo_square <= 6.25) scenario.input.support[1].push_back({row, halo_square});
                const auto target=SimulationBasis(target_square,TargetTruth[1],2.5);
                const auto halo=SimulationBasis(halo_square,TargetTruth[1],simulation_radius);
                scenario.target_truth.push_back(TargetTruth[0]*target.gaussian+TargetTruth[2]*target.charge);
                scenario.halo_truth.push_back(0.7*halo.gaussian+0.18*halo.charge);
                scenario.input.observations.push_back(scenario.target_truth.back()+scenario.halo_truth.back());
                scenario.halo_correction.push_back(halo_square<=6.25 ? FirstStageHalo[0]*
                    joint::EvaluateKernel(halo_square,FirstStageHalo[1],2.5).gaussian : 0.0);
            }
    scenario.clean_observations=scenario.input.observations;
    return scenario;
}

struct RowBasis { double g{}, q{}, dg{}, dq{}; };
RowBasis Basis(double square, double width)
{
    const auto value=joint::EvaluateKernel(square,width,2.5);
    const double r2=square, r=std::sqrt(square);
    const double dg=value.gaussian*r2/(width*width);
    const double dq=r<1e-5 ? -value.charge :
        -std::sqrt(2.0/std::numbers::pi)*std::exp(-r2/(2.0*width*width))/width;
    return {value.gaussian,value.charge,dg,dq};
}

std::vector<double> TargetSquares(const Scenario & scenario)
{
    std::vector<double> result(scenario.input.observations.size());
    for (const auto & row:scenario.input.support[0]) result.at(row.row)=row.squared_distance;
    return result;
}

double Objective(const Scenario & scenario, const std::vector<double> & observations,
    const std::vector<double> & weights, bool correct_halo, const Parameters & p)
{
    const auto squares=TargetSquares(scenario);
    long double total{};
    for (std::size_t row=0;row<weights.size();++row)
    {
        const auto basis=Basis(squares[row],p.b);
        const double correction=correct_halo ? scenario.halo_correction[row] : 0.0;
        const double residual=p.a*basis.g+p.c*basis.q+correction-observations[row];
        total+=0.5L*weights[row]*residual*residual;
    }
    return static_cast<double>(total);
}

double ReplayObjective(const Scenario & scenario, const std::vector<double> & observations,
    const std::vector<double> & weights, bool correct_halo, const Parameters & p)
{
    long double total{};
    std::vector<double> correction(observations.size());
    if (correct_halo)
        for (const auto & edge:scenario.input.support[1])
            correction.at(edge.row)+=FirstStageHalo[0]*
                joint::EvaluateKernel(edge.squared_distance,FirstStageHalo[1],2.5).gaussian;
    for(const auto & edge:scenario.input.support[0])
    {
        const auto values=joint::EvaluateKernel(edge.squared_distance,p.b,2.5);
        const double residual=p.a*values.gaussian+p.c*values.charge+correction.at(edge.row)-observations.at(edge.row);
        total+=0.5L*weights.at(edge.row)*residual*residual;
    }
    return static_cast<double>(total);
}

Gradient ObjectiveGradient(const Scenario & scenario, const std::vector<double> & observations,
    const std::vector<double> & weights, bool correct_halo, const Parameters & p)
{
    const auto squares=TargetSquares(scenario);
    Gradient out;
    for (std::size_t row=0;row<weights.size();++row)
    {
        const auto basis=Basis(squares[row],p.b);
        const double correction=correct_halo ? scenario.halo_correction[row] : 0.0;
        const double residual=p.a*basis.g+p.c*basis.q+correction-observations[row];
        out.a+=weights[row]*residual*basis.g;
        out.c+=weights[row]*residual*basis.q;
        out.b+=weights[row]*residual*(p.a*basis.dg+p.c*basis.dq);
    }
    return out;
}

Parameters SolveAC(const Scenario & scenario, const std::vector<double> & observations,
    const std::vector<double> & weights, bool correct_halo, double width)
{
    const auto squares=TargetSquares(scenario);
    long double gg{},gq{},qq{},gy{},qy{};
    for (std::size_t row=0;row<weights.size();++row)
    {
        const auto basis=Basis(squares[row],width);
        const double y=observations[row]-(correct_halo ? scenario.halo_correction[row] : 0.0);
        gg+=weights[row]*basis.g*basis.g; gq+=weights[row]*basis.g*basis.q; qq+=weights[row]*basis.q*basis.q;
        gy+=weights[row]*basis.g*y; qy+=weights[row]*basis.q*y;
    }
    const long double determinant=gg*qq-gq*gq;
    std::vector<Parameters> candidates{{0,width,0}};
    if (determinant>1e-24L)
    {
        const double a=static_cast<double>((gy*qq-qy*gq)/determinant);
        const double c=static_cast<double>((qy*gg-gy*gq)/determinant);
        if(a>=0 && c>=0)candidates.push_back({a,width,c});
    }
    if(qq>0)candidates.push_back({0,width,std::max(0.0,static_cast<double>(qy/qq))});
    if(gg>0)candidates.push_back({std::max(0.0,static_cast<double>(gy/gg)),width,0});
    auto best=candidates.front(); double best_value=Objective(scenario,observations,weights,correct_halo,best);
    for(const auto & candidate:candidates)
    {
        const double value=Objective(scenario,observations,weights,correct_halo,candidate);
        if(value<best_value){best=candidate;best_value=value;}
    }
    return best;
}

Fit FitProfile(const Scenario & scenario, const std::vector<double> & observations,
    const std::vector<double> & weights, bool correct_halo)
{
    const double low=std::log(0.25), high=std::log(2.0);
    constexpr int grid=256;
    int best_index{}; double best_value=std::numeric_limits<double>::infinity();
    const auto evaluate=[&](double eta) {
        const auto p=SolveAC(scenario,observations,weights,correct_halo,std::exp(eta));
        return std::pair{Objective(scenario,observations,weights,correct_halo,p),p};
    };
    for(int i=0;i<=grid;++i)
    {
        const double eta=low+(high-low)*i/grid;
        const auto [value,unused]=evaluate(eta); (void)unused;
        if(value<best_value){best_value=value;best_index=i;}
    }
    double left=low+(high-low)*std::max(0,best_index-1)/grid;
    double right=low+(high-low)*std::min(grid,best_index+1)/grid;
    constexpr double ratio=0.6180339887498948482;
    double x1=right-ratio*(right-left),x2=left+ratio*(right-left);
    double f1=evaluate(x1).first,f2=evaluate(x2).first;
    for(int i=0;i<100;++i)
    {
        if(f1<=f2){right=x2;x2=x1;f2=f1;x1=right-ratio*(right-left);f1=evaluate(x1).first;}
        else{left=x1;x1=x2;f1=f2;x2=left+ratio*(right-left);f2=evaluate(x2).first;}
    }
    const auto [objective,parameters]=evaluate((left+right)/2);
    const auto gradient=ObjectiveGradient(scenario,observations,weights,correct_halo,parameters);
    const auto finite=[&](int dimension) {
        auto plus=parameters,minus=parameters;
        const double h=1e-6;
        if(dimension==0){plus.a+=h;minus.a-=h;}
        else if(dimension==1){plus.b*=std::exp(h);minus.b*=std::exp(-h);}
        else {plus.c+=h;minus.c-=h;}
        return (Objective(scenario,observations,weights,correct_halo,plus)-
            Objective(scenario,observations,weights,correct_halo,minus))/(2*h);
    };
    const double max_gradient_error=std::max({std::abs(gradient.a-finite(0)),
        std::abs(gradient.b-finite(1)),std::abs(gradient.c-finite(2))});
    auto changed=observations;
    for(std::size_t row=0;row<weights.size();++row)if(weights[row]==0)changed[row]+=123.0;
    const double zero_row_delta=std::abs(objective-Objective(scenario,changed,weights,correct_halo,parameters));
    const auto changed_gradient=ObjectiveGradient(scenario,changed,weights,correct_halo,parameters);
    const double zero_row_gradient_delta=std::max({std::abs(gradient.a-changed_gradient.a),
        std::abs(gradient.b-changed_gradient.b),std::abs(gradient.c-changed_gradient.c)});
    const double replay_delta=std::abs(objective-ReplayObjective(scenario,observations,weights,correct_halo,parameters));
    if(max_gradient_error>2e-6 || zero_row_delta!=0.0 || zero_row_gradient_delta!=0.0 || replay_delta>2e-12)
        throw std::runtime_error("Weighted objective gradient or zero-row invariant failed.");
    return {parameters,objective,max_gradient_error,zero_row_delta,zero_row_gradient_delta,replay_delta};
}

struct Stats { std::size_t count{}; double bias{},variance{},mse{},rmse{},median_abs{},p95_abs{},max_abs{}; };
Stats Summarize(const std::vector<double> & values)
{
    Stats out; out.count=values.size(); if(values.empty())return out;
    long double sum{},squares{}; std::vector<double> absolute; absolute.reserve(values.size());
    for(double value:values){sum+=value;squares+=static_cast<long double>(value)*value;absolute.push_back(std::abs(value));}
    const auto count=static_cast<long double>(values.size());
    out.bias=static_cast<double>(sum/count);
    out.variance=static_cast<double>(values.size()>1 ? (squares-sum*sum/count)/(count-1) : 0);
    out.mse=static_cast<double>(squares/count); out.rmse=std::sqrt(out.mse);
    std::sort(absolute.begin(),absolute.end());
    out.median_abs=absolute[(absolute.size()-1)/2];
    out.p95_abs=absolute[static_cast<std::size_t>(std::ceil(.95*static_cast<double>(absolute.size())))-1];
    out.max_abs=absolute.back(); return out;
}
json::object StatsJson(const Stats & s)
{ return {{"count",s.count},{"bias",Number(s.bias)},{"sample_variance",Number(s.variance)},
    {"mse",Number(s.mse)},{"rmse",Number(s.rmse)},{"median_absolute_error",Number(s.median_abs)},
    {"p95_absolute_error",Number(s.p95_abs)},{"maximum_absolute_error",Number(s.max_abs)}}; }

json::object FitJson(const Fit & fit)
{ return {{"A",fit.parameters.a},{"B_angstrom",fit.parameters.b},{"C",fit.parameters.c},
    {"objective",fit.objective},{"maximum_gradient_finite_difference_error",fit.max_gradient_error},
    {"zero_weight_row_objective_delta",fit.zero_row_delta},
    {"zero_weight_row_gradient_delta",fit.zero_row_gradient_delta},
    {"independent_replay_objective_delta",fit.replay_objective_delta}}; }

std::size_t PeakRssBytes()
{
#ifdef _WIN32
    return 0;
#else
    rusage usage{}; getrusage(RUSAGE_SELF,&usage); return static_cast<std::size_t>(usage.ru_maxrss)*1024;
#endif
}

json::object RunCondition(const Scenario & scenario,double sigma)
{
    const auto support_mask=audit::HaloSupportExclusionRows(scenario.input);
    std::vector<double> errors[7][3];
    std::vector<double> full_prediction[7],target_prediction[7],retained_halo[7];
    json::array noise_realizations,repetitions[7];
    const std::array<std::string,7> names{"B-pure-hard-mask","W-lambda-0","W-lambda-0.01","W-lambda-0.05","W-lambda-0.1","W-lambda-0.25","W-lambda-1"};
    const auto append=[&](std::size_t group,const Scenario & s,const std::vector<double> & observations,
        const std::vector<double> & weights,bool correct_halo,std::uint64_t seed) {
        const auto fit=FitProfile(s,observations,weights,correct_halo);
        const std::array<double,3> truth{TargetTruth[0],TargetTruth[1],TargetTruth[2]};
        const std::array<double,3> estimate{fit.parameters.a,fit.parameters.b,fit.parameters.c};
        json::array error_json;
        for(std::size_t p=0;p<3;++p){errors[group][p].push_back(estimate[p]-truth[p]);error_json.push_back(estimate[p]-truth[p]);}
        std::vector<double> target_errors,full_errors,halo_values;
        const auto squares=TargetSquares(s);
        for(std::size_t row=0;row<observations.size();++row)
        {
            const auto basis=Basis(squares[row],fit.parameters.b);
            const double target=fit.parameters.a*basis.g+fit.parameters.c*basis.q;
            const double halo=correct_halo ? s.halo_correction[row] : 0.0;
            target_errors.push_back(target-s.target_truth[row]);
            full_errors.push_back(target+halo-s.clean_observations[row]);
            if(support_mask[row])halo_values.push_back(s.halo_truth[row]);
        }
        target_prediction[group].push_back(std::sqrt(std::inner_product(target_errors.begin(),target_errors.end(),target_errors.begin(),0.0)/static_cast<double>(target_errors.size())));
        full_prediction[group].push_back(std::sqrt(std::inner_product(full_errors.begin(),full_errors.end(),full_errors.begin(),0.0)/static_cast<double>(full_errors.size())));
        retained_halo[group].push_back(std::sqrt(std::inner_product(halo_values.begin(),halo_values.end(),halo_values.begin(),0.0)/static_cast<double>(halo_values.size())));
        repetitions[group].emplace_back(json::object{{"seed",seed},{"estimate",json::array{estimate[0],estimate[1],estimate[2]}},
            {"parameter_errors",std::move(error_json)},{"fit",FitJson(fit)}});
    };
    for(auto seed:Seeds)
    {
        std::vector<double> noisy(scenario.input.observations.size());
        const auto noise=audit::FillGaussianNoise(scenario.clean_observations,noisy,sigma,seed);
        noise_realizations.emplace_back(json::object{{"seed",seed},{"realized_mean",noise.mean},
            {"realized_standard_deviation",noise.standard_deviation}});
        std::vector<double> pure_weights(noisy.size()),zero_weights(noisy.size());
        for(std::size_t row=0;row<noisy.size();++row)
        {
            pure_weights[row]=support_mask[row] ? 1.0 : 0.0;
            zero_weights[row]=support_mask[row] ? 1.0 : 0.0;
        }
        append(0,scenario,noisy,pure_weights,false,seed);
        append(1,scenario,noisy,zero_weights,true,seed);
        const auto & pure=repetitions[0].back().at("fit").as_object();
        const auto & zero=repetitions[1].back().at("fit").as_object();
        for(const auto key:{"A","B_angstrom","C","objective"})
            if(pure.at(key).as_double()!=zero.at(key).as_double())
                throw std::runtime_error("lambda=0 differs from exact hard-mask estimate.");
        for(std::size_t i=1;i<Weights.size();++i)
        {
            std::vector<double> weights(noisy.size());
            for(std::size_t row=0;row<noisy.size();++row)weights[row]=support_mask[row] ? 1.0 : Weights[i];
            append(i+1,scenario,noisy,weights,true,seed);
        }
    }
    json::object groups;
    for(std::size_t group=0;group<names.size();++group)
    {
        json::array parameter_stats;
        for(const auto & values:errors[group])parameter_stats.push_back(StatsJson(Summarize(values)));
        groups[names[group]]=json::object{{"attempted_repetitions",Seeds.size()},
            {"parameter_order",json::array{"A","B_angstrom","C"}},
            {"parameter_bias_variance_mse",std::move(parameter_stats)},
            {"target_prediction_rmse",StatsJson(Summarize(target_prediction[group]))},
            {"full_prediction_rmse",StatsJson(Summarize(full_prediction[group]))},
            {"retained_true_halo_rmse",StatsJson(Summarize(retained_halo[group]))},
            {"repetitions",std::move(repetitions[group])}};
    }
    const double row_count=static_cast<double>(support_mask.size());
    const auto exact=std::count(support_mask.begin(),support_mask.end(),std::uint8_t{1});
    return {{"noise_sigma_synthetic_normalized_units",sigma},{"paired_seeds",json::value_from(Seeds)},
        {"noise_realizations",std::move(noise_realizations)},
        {"coverage",json::object{{"original_rows",support_mask.size()},{"exact_retained_rows",exact},
            {"exact_retained_fraction",static_cast<double>(exact)/row_count}}},
        {"groups",std::move(groups)}};
}

void Run(const std::filesystem::path & output)
{
    json::array conditions;
    for(double radius:SimulationRadii)
    {
        const auto scenario=BuildScenario(radius);
        for(double sigma:NoiseSigmas)
            conditions.emplace_back(json::object{{"simulation_halo_support_radius_angstrom",radius},
                {"joint_support_radius_angstrom",2.5},{"condition",RunCondition(scenario,sigma)}});
    }
    Write(output,json::object{{"contract","joint-halo-support-weighting-reference-v1"},
        {"noise_layer","paired iid Gaussian noise added to normalized synthetic Joint observation rows"},
        {"objective","0.5 * sum_v w_v * (Target + fixed FirstStage Gaussian Halo - noisy observation)^2; s=1"},
        {"weight_policy","weights fixed before optimization; w=lambda on exact Halo-support rows, 1 elsewhere"},
        {"first_stage_halo_prior",json::object{{"A",FirstStageHalo[0]},{"B_angstrom",FirstStageHalo[1]},{"C",0.0},
            {"source","fixed synthetic FirstStage-like prior, same for all seeds"},
            {"qualification_status","assumed-qualified synthetic control; FirstStage qualification implementation not run"}}},
        {"lambda_values",json::value_from(Weights)},
        {"profile_solver","reference constrained A/C active-set with log-width grid plus golden-section profile search; not Production LM"},
        {"gradient_validation","analytic A/C/logB gradients compared with central finite differences; zero-weight-row perturbation asserted invariant"},
        {"support_mismatch","Halo simulation Gaussian and charge terms extend to R_sim; Joint support remains 2.5 A"},
        {"target_count",1},{"halo_count",1},{"target_truth",json::value_from(TargetTruth)},
        {"simulation_support_regimes_angstrom",json::value_from(SimulationRadii)},
        {"noise_sigmas_synthetic_normalized_units",json::value_from(NoiseSigmas)},
        {"peak_rss_bytes",PeakRssBytes()},{"threads",json::object{{"Eigen",1},
            {"OpenMP",std::getenv("OMP_NUM_THREADS")?std::getenv("OMP_NUM_THREADS"):"unset"},
            {"OpenBLAS",std::getenv("OPENBLAS_NUM_THREADS")?std::getenv("OPENBLAS_NUM_THREADS"):"unset"},
            {"Accelerate",std::getenv("VECLIB_MAXIMUM_THREADS")?std::getenv("VECLIB_MAXIMUM_THREADS"):"unset"}}},
        {"conditions",std::move(conditions)}});
}
}

int main(int argc,char ** argv)
{
    try
    {
        if(argc!=2)throw std::invalid_argument("Usage: joint_halo_support_weighting_reference OUTPUT.json");
        Run(argv[1]); return 0;
    }
    catch(const std::exception & error){std::cerr<<error.what()<<'\n';return 1;}
}
