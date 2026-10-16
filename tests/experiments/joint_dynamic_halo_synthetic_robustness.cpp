#include "core/detail/joint_component/Numerics.hpp"
#include "support/DynamicHaloMean.hpp"

#include <boost/json.hpp>
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
namespace core=rhbm_gem::core;
namespace joint=rhbm_gem::core::joint_component;
namespace dyn=joint_test::dynamic_halo;
namespace json=boost::json;

struct Case
{
    std::string name;
    std::array<dyn::Parameters,2> targets;
    dyn::Parameters halo,first_stage;
    std::array<std::string,2> elements;
    std::size_t halo_support_rows{},seed{};
    double noise_sd{};
};

std::vector<double> Prediction(const rhbm_gem::core::JointProblemInput & input,
    std::span<const dyn::Parameters> parameters)
{
    return dyn::ReplayPrediction(input,parameters);
}

double Rmse(const std::vector<double> & left,const std::vector<double> & right)
{
    if(left.size()!=right.size() || left.empty()) throw std::invalid_argument("Synthetic RMSE size mismatch.");
    long double sum{};
    for(std::size_t i=0;i<left.size();++i)
    {
        const long double difference=static_cast<long double>(left[i])-right[i];
        sum+=difference*difference;
    }
    return std::sqrt(static_cast<double>(sum/static_cast<long double>(left.size())));
}

double Objective(const std::vector<double> & prediction,const std::vector<double> & observations,double scale)
{
    long double sum{};
    for(std::size_t row=0;row<prediction.size();++row)
    {
        const long double residual=static_cast<long double>(prediction[row])-observations[row];
        sum+=residual*residual;
    }
    return static_cast<double>(0.5L*sum/(static_cast<long double>(scale)*scale));
}

json::object Evaluate(const Case & test_case)
{
    using core::JointProblemInput;
    JointProblemInput input;
    constexpr std::size_t row_count=41;
    input.observations.assign(row_count,0.0);
    input.support.resize(3);
    input.atom_ids={"target-0","target-1","halo"};
    input.selection_domain.emplace();
    input.selection_domain->target_indices={0,1};
    for(std::size_t row=0;row<row_count;++row) input.row_ids.push_back("row-"+std::to_string(row));

    constexpr std::array centers{-0.7,0.7,0.0};
    std::vector<std::pair<double,std::size_t>> halo_candidates;
    for(std::size_t row=0;row<row_count;++row)
    {
        const double coordinate=-2.4+0.12*static_cast<double>(row);
        for(std::size_t atom=0;atom<3;++atom)
        {
            const double distance=coordinate-centers[atom];
            const double square=distance*distance;
            if(square<=6.25)
            {
                if(atom<2) input.support[atom].push_back({row,square});
                else halo_candidates.emplace_back(square,row);
            }
        }
    }
    std::sort(halo_candidates.begin(),halo_candidates.end());
    for(std::size_t i=0;i<std::min(test_case.halo_support_rows,halo_candidates.size());++i)
        input.support[2].push_back({halo_candidates[i].second,halo_candidates[i].first});
    std::sort(input.support[2].begin(),input.support[2].end(),[](const auto & a,const auto & b) {
        return a.row<b.row;
    });
    const auto ownership=dyn::BuildOwnership(input,{{"synthetic-core",0,{0,1}}});

    std::vector<dyn::Parameters> truth{
        test_case.targets[0],test_case.targets[1],test_case.halo};
    const auto truth_full=Prediction(input,truth);
    std::mt19937_64 generator(test_case.seed);
    std::normal_distribution<double> noise(0.0,test_case.noise_sd);
    for(std::size_t row=0;row<row_count;++row) input.observations[row]=truth_full[row]+noise(generator);
    double observation_norm{};
    for(const auto value:input.observations) observation_norm+=value*value;
    const double scale=std::max(1.0,std::sqrt(observation_norm));
    const auto context=joint::CreateContext(std::make_shared<JointProblemInput>(input));
    const Eigen::VectorXd weights=Eigen::VectorXd::Ones(static_cast<Eigen::Index>(row_count));

    std::vector<dyn::Parameters> target_seed=truth;
    for(std::size_t target=0;target<2;++target)
    {target_seed[target].a=0.0; target_seed[target].c=0.0;}
    target_seed[2]={test_case.first_stage.a,test_case.first_stage.b,0.0};
    auto fixed_halo_prediction=Prediction(input,target_seed);
    const auto target_only_design=dyn::BuildAugmentedDesign(input,ownership,0,target_seed,dyn::Mode::StaticABC);
    Eigen::VectorXd target_response(static_cast<Eigen::Index>(target_only_design.rows.size()));
    for(std::size_t i=0;i<target_only_design.rows.size();++i)
        target_response(static_cast<Eigen::Index>(i))=input.observations[target_only_design.rows[i]]-
            fixed_halo_prediction[target_only_design.rows[i]];
    const auto fixed_fit=joint::SolveLinear(target_only_design.values,target_response,weights,false,false,
        nullptr,&context.linear,scale);
    if(!fixed_fit.valid) throw std::runtime_error("Synthetic fixed-Halo solve failed: "+fixed_fit.reason);
    auto fixed_state=target_seed;
    for(std::size_t target=0;target<2;++target)
    {
        fixed_state[target].a=fixed_fit.beta(static_cast<Eigen::Index>(2*target));
        fixed_state[target].c=fixed_fit.beta(static_cast<Eigen::Index>(2*target+1));
    }
    const auto fixed_full=Prediction(input,fixed_state);
    auto fixed_targets_only=fixed_state; fixed_targets_only[2]={0.0,test_case.halo.b,0.0};
    auto fixed_halo_only=fixed_state; fixed_halo_only[0]={0.0,test_case.targets[0].b,0.0};
    fixed_halo_only[1]={0.0,test_case.targets[1].b,0.0};

    const auto mean_design=dyn::BuildAugmentedDesign(input,ownership,0,target_seed,dyn::Mode::ABC);
    Eigen::VectorXd mean_response(static_cast<Eigen::Index>(mean_design.rows.size()));
    for(std::size_t i=0;i<mean_design.rows.size();++i)
        mean_response(static_cast<Eigen::Index>(i))=input.observations[mean_design.rows[i]];
    const auto mean_fit=joint::SolveLinear(mean_design.values,mean_response,weights,false,false,
        nullptr,&context.linear,scale);
    if(!mean_fit.valid) throw std::runtime_error("Synthetic dynamic-mean solve failed: "+mean_fit.reason);
    auto mean_state=target_seed;
    for(std::size_t target=0;target<2;++target)
    {
        mean_state[target].a=mean_fit.beta(static_cast<Eigen::Index>(2*target));
        mean_state[target].c=mean_fit.beta(static_cast<Eigen::Index>(2*target+1));
    }
    mean_state=dyn::DeriveHalos(input,ownership,mean_state,dyn::Mode::ABC);
    const auto mean_full=Prediction(input,mean_state);
    auto mean_targets_only=mean_state; mean_targets_only[2]={0.0,mean_state[2].b,0.0};
    auto mean_halo_only=mean_state; mean_halo_only[0]={0.0,mean_state[0].b,0.0};
    mean_halo_only[1]={0.0,mean_state[1].b,0.0};

    std::vector<dyn::Parameters> truth_targets=truth;
    truth_targets[2]={0.0,truth[2].b,0.0};
    std::vector<dyn::Parameters> truth_halo=truth;
    truth_halo[0]={0.0,truth[0].b,0.0}; truth_halo[1]={0.0,truth[1].b,0.0};
    const auto truth_target_prediction=Prediction(input,truth_targets);
    const auto truth_halo_prediction=Prediction(input,truth_halo);
    const auto fixed_target_prediction=Prediction(input,fixed_targets_only);
    const auto fixed_halo_field=Prediction(input,fixed_halo_only);
    const auto mean_target_prediction=Prediction(input,mean_targets_only);
    const auto mean_halo_field=Prediction(input,mean_halo_only);
    double fixed_a_error{},fixed_c_error{},mean_a_error{},mean_c_error{};
    for(std::size_t target=0;target<2;++target)
    {
        fixed_a_error+=std::pow(fixed_state[target].a-truth[target].a,2);
        fixed_c_error+=std::pow(fixed_state[target].c-truth[target].c,2);
        mean_a_error+=std::pow(mean_state[target].a-truth[target].a,2);
        mean_c_error+=std::pow(mean_state[target].c-truth[target].c,2);
    }
    const double target_prediction_truth_rmse_fixed=Rmse(fixed_target_prediction,truth_target_prediction);
    const double target_prediction_truth_rmse_mean=Rmse(mean_target_prediction,truth_target_prediction);
    json::object metrics{
        {"target_A_rmse_fixed",std::sqrt(fixed_a_error/2.0)},
        {"target_A_rmse_dynamic",std::sqrt(mean_a_error/2.0)},
        {"target_C_rmse_fixed",std::sqrt(fixed_c_error/2.0)},
        {"target_C_rmse_dynamic",std::sqrt(mean_c_error/2.0)},
        {"halo_prediction_rmse_fixed",Rmse(fixed_halo_field,truth_halo_prediction)},
        {"halo_prediction_rmse_dynamic",Rmse(mean_halo_field,truth_halo_prediction)},
        {"target_prediction_rmse_fixed",target_prediction_truth_rmse_fixed},
        {"target_prediction_rmse_dynamic",target_prediction_truth_rmse_mean},
        {"full_prediction_rmse_fixed",Rmse(fixed_full,truth_full)},
        {"full_prediction_rmse_dynamic",Rmse(mean_full,truth_full)},
        {"objective_fixed",Objective(fixed_full,input.observations,scale)},
        {"objective_dynamic",Objective(mean_full,input.observations,scale)},
        {"dynamic_halo_mean",json::object{{"A",mean_state[2].a},{"B",mean_state[2].b},{"C",mean_state[2].c}}}};
    return json::object{
        {"name",test_case.name},{"seed",test_case.seed},{"noise_sd",test_case.noise_sd},
        {"target_elements",json::array{test_case.elements[0],test_case.elements[1]}},
        {"halo_support_rows",input.support[2].size()},
        {"target_truth",json::array{
            json::object{{"A",truth[0].a},{"B",truth[0].b},{"C",truth[0].c}},
            json::object{{"A",truth[1].a},{"B",truth[1].b},{"C",truth[1].c}}}},
        {"halo_truth",json::object{{"A",truth[2].a},{"B",truth[2].b},{"C",truth[2].c}}},
        {"first_stage_halo",json::object{{"A",test_case.first_stage.a},
            {"B",test_case.first_stage.b},{"C",0.0}}},
        {"metrics",std::move(metrics)}};
}
}

int main(int argc,char ** argv)
{
    try
    {
        if(argc!=2) throw std::invalid_argument("Usage: joint_dynamic_halo_synthetic_robustness OUTPUT_JSON");
        const std::vector<Case> cases{
            {"homogeneous-matched-low-noise",{{{1.0,0.45,-0.2},{1.0,0.45,-0.2}}},
                {1.0,0.45,-0.2},{0.85,0.45,0},{"C","C"},8,11,0.005},
            {"heterogeneous-matched-charge",{{{0.4,0.4,-0.8},{1.6,0.8,0.4}}},
                {1.0,0.6,-0.2},{0.75,0.6,0},{"C","O"},16,17,0.005},
            {"heterogeneous-halo-mismatch",{{{0.3,0.42,-0.6},{1.9,0.7,0.8}}},
                {0.2,0.38,-0.6},{0.45,0.38,0},{"N","O"},16,19,0.02},
            {"mixed-elements-short-support-charge-match",{{{0.6,0.45,0.7},{1.3,0.56,-0.1}}},
                {0.95,0.505,0.3},{0.8,0.505,0},{"C","O"},2,23,0.02},
            {"mixed-sign-charge-short-support-mismatch",{{{0.9,0.5,-0.9},{1.2,0.55,0.3}}},
                {1.05,0.525,0.6},{0.9,0.525,0},{"N","O"},1,29,0.05},
            {"higher-noise-matched",{{{0.8,0.48,-0.3},{1.2,0.62,0.5}}},
                {1.0,0.55,0.1},{0.9,0.55,0},{"C","N"},8,31,0.1}};
        json::array results;
        for(const auto & item:cases) results.emplace_back(Evaluate(item));
        std::ofstream output(argv[1]);
        output.exceptions(std::ios::failbit|std::ios::badbit);
        output<<json::serialize(json::object{{"contract","dynamic-halo-synthetic-fixed-width-v1"},
            {"method","one-block constrained linear A/C fit; Target widths held at truth; E0 fixed Halo and E2 whole-core mean"},
            {"noise","std::mt19937_64 with fixed per-case seeds; Gaussian observation noise"},
            {"kernel","production EvaluateKernel, cutoff 2.5"},{"cases",std::move(results)}})<<'\n';
        std::cout<<"Wrote "<<argv[1]<<" with 6 deterministic synthetic cases.\n";
        return 0;
    }
    catch(const std::exception & error) {std::cerr<<error.what()<<'\n'; return 1;}
}
