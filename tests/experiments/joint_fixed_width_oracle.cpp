#include "core/detail/joint_component/Problem.hpp"
#include <rhbm_gem/data/io/ModelMapFileIO.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <boost/json.hpp>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <iostream>
#ifndef _WIN32
#include <sys/resource.h>
#endif
#include <unordered_map>
#include <vector>

namespace {
namespace core=rhbm_gem::core;
namespace joint=core::joint_component;
namespace json=boost::json;
using Clock=std::chrono::steady_clock;
struct Truth {double amplitude{},charge{},width{};};
struct Measurement
{
    double a_squared{},c_squared{},b_squared{},a_truth_squared{},c_truth_squared{};
    std::size_t a_count{},c_count{},b_count{};
    void Add(double a,double c,double b,const Truth & truth)
    {
        a_squared+=std::pow(a-truth.amplitude,2); c_squared+=std::pow(c-truth.charge,2);
        b_squared+=std::pow(b-truth.width,2); a_truth_squared+=truth.amplitude*truth.amplitude;
        c_truth_squared+=truth.charge*truth.charge;
        ++a_count; ++c_count; ++b_count;
    }
};
json::value Number(double value) {return std::isfinite(value) ? json::value(value) : json::value(nullptr);}
json::value Read(const std::filesystem::path & path)
{
    std::ifstream input(path);
    if(!input) throw std::runtime_error("Cannot read simulation metadata: "+path.string());
    json::parse_options options; options.numbers=json::number_precision::precise;
    return json::parse(std::string(std::istreambuf_iterator<char>(input),{}),{},options);
}
std::unordered_map<std::string,Truth> Truths(const std::filesystem::path & path)
{
    const auto metadata=Read(path); std::unordered_map<std::string,Truth> out;
    for(const auto & atom:metadata.at("atoms").as_array())
    {
        const auto id=std::to_string(json::value_to<std::uint64_t>(atom.at("serial_id")));
        const double width=json::value_to<double>(atom.at("effective_gaussian_width"));
        const double charge_width=json::value_to<double>(atom.at("effective_charge_width"));
        if(width!=charge_width) throw std::runtime_error("Map truth uses different Gaussian and charge widths at atom "+id);
        out.emplace(id,Truth{static_cast<double>(json::value_to<int>(atom.at("element"))),
            json::value_to<double>(atom.at("charge_used")),width});
    }
    return out;
}
double Rmse(const Measurement & values,bool amplitude)
{
    const auto count=amplitude ? values.a_count : values.c_count;
    return count ? std::sqrt((amplitude ? values.a_squared : values.c_squared)/count) : joint::unavailable;
}
double WidthRmse(const Measurement & values)
{return values.b_count ? std::sqrt(values.b_squared/values.b_count) : joint::unavailable;}
double RelativeRmse(const Measurement & values,bool amplitude)
{
    const auto count=amplitude ? values.a_count : values.c_count;
    const double truth=std::sqrt((amplitude ? values.a_truth_squared : values.c_truth_squared)/std::max<std::size_t>(1,count));
    return truth>0 ? Rmse(values,amplitude)/truth : joint::unavailable;
}
long PeakRssBytes()
{
#ifdef _WIN32
    return -1;
#else
    rusage usage{};
    if(getrusage(RUSAGE_SELF,&usage)!=0) return -1;
#ifdef __APPLE__
    return usage.ru_maxrss;
#else
    return usage.ru_maxrss*1024;
#endif
#endif
}
json::value PeakRss()
{const auto bytes=PeakRssBytes(); return bytes>=0 ? json::value(bytes) : json::value(nullptr);}
json::object RunCase(const std::string & name,const joint::ProblemData & data,
    const std::unordered_map<std::string,Truth> & truths,
    const std::vector<bool> & targets,const std::vector<double> & initial_widths,
    std::vector<std::string> fixed_width_ids,std::size_t maximum_sweeps,const std::filesystem::path & output)
{
    if(data.partition.components.size()!=1)
        throw std::runtime_error("The 6Z6U oracle control expects one structural component.");
    joint::Vector initial(static_cast<Eigen::Index>(initial_widths.size()));
    for(std::size_t atom=0;atom<initial_widths.size();++atom) initial(static_cast<Eigen::Index>(atom))=initial_widths[atom];
    const auto & component=data.partition.components.front();
    const auto prepared=joint::MakeComponentProblemView(data.layout,component);
    joint::FixedNeighborSearchPolicy policy;
    policy.core_atoms=12; policy.maximum_sweeps=maximum_sweeps;
    policy.order=joint::FixedNeighborBlockOrder::Forward;
    policy.configure_search_for_testing=[&](joint::FixedNeighborPolicy & local) {
        local.certify_endpoint=false; local.collect_diagnostics=true; local.collect_telemetry=true;
        local.assess_local_conditioning=false; local.assess_parameter_stability=false;
        local.fixed_width_atom_ids=fixed_width_ids;
    };
    const auto started=Clock::now(); const auto fit=joint::SolveFixedNeighborComponent(
        prepared,initial,data.context,policy);
    const double seconds=std::chrono::duration<double>(Clock::now()-started).count();
    const auto & state=fit.search.state;
    joint::Vector eta=joint::Vector::Zero(static_cast<Eigen::Index>(data.input->support.size()));
    joint::Vector beta=joint::Vector::Zero(2*eta.size());
    std::vector<Eigen::Index> local_by_parent(data.input->support.size(),-1);
    for(std::size_t k=0;k<prepared.local_to_parent_atoms.size();++k)
    {
        const auto parent=prepared.local_to_parent_atoms[k]; local_by_parent.at(static_cast<std::size_t>(parent))=static_cast<Eigen::Index>(k);
        eta(parent)=state.eta(static_cast<Eigen::Index>(k));
        beta.segment<2>(2*parent)=state.beta.segment<2>(2*static_cast<Eigen::Index>(k));
    }
    const auto primary=joint::EvaluateStateStreaming(data.domain,data.y,eta,beta,data.context);
    const auto replay=joint::CheckIndependentReplay(data.domain,data.y,primary,data.context);
    const auto objective=joint::NormalizedProfileObjective(primary,data.context.scale);
    double target_gradient=joint::unavailable,free_target_gradient=joint::unavailable,
        all_width_gradient=joint::unavailable;
    if(primary.gradient.size()==static_cast<Eigen::Index>(data.input->support.size()))
    {
        target_gradient=0.0; free_target_gradient=0.0; all_width_gradient=0.0;
        for(std::size_t atom=0;atom<data.input->support.size();++atom)
        {
            const double gradient=std::abs(primary.gradient(static_cast<Eigen::Index>(atom)));
            all_width_gradient=std::max(all_width_gradient,gradient);
            if(targets[atom])
            {
                target_gradient=std::max(target_gradient,gradient);
                if(std::find(fixed_width_ids.begin(),fixed_width_ids.end(),data.input->atom_ids[atom])==
                    fixed_width_ids.end()) free_target_gradient=std::max(free_target_gradient,gradient);
            }
        }
    }
    Measurement target,halo;
    joint::Vector target_prediction=joint::Vector::Zero(static_cast<Eigen::Index>(data.input->observations.size()));
    joint::Vector halo_prediction=joint::Vector::Zero(target_prediction.size());
    std::ofstream csv(output/(name+"-estimates.csv"));
    csv.exceptions(std::ios::failbit|std::ios::badbit);
    csv<<"atom_id,role,A_truth,A_estimate,C_truth,C_estimate,B_truth,B_estimate\n";
    for(std::size_t atom=0;atom<data.input->support.size();++atom)
    {
        const auto found=truths.find(data.input->atom_ids[atom]);
        if(found==truths.end()) throw std::runtime_error("No simulation truth for contributor "+data.input->atom_ids[atom]);
        const auto & truth=found->second; const auto local=local_by_parent.at(atom);
        if(local<0) throw std::runtime_error("A contributor is absent from the component state.");
        const double a=beta(2*static_cast<Eigen::Index>(atom)),c=beta(2*static_cast<Eigen::Index>(atom)+1),b=std::exp(eta(static_cast<Eigen::Index>(atom)));
        (targets[atom] ? target : halo).Add(a,c,b,truth);
        csv<<data.input->atom_ids[atom]<<','<<(targets[atom] ? "target" : "halo")<<','
            <<truth.amplitude<<','<<a<<','<<truth.charge<<','<<c<<','<<truth.width<<','<<b<<'\n';
        auto & contribution=targets[atom] ? target_prediction : halo_prediction;
        for(const auto & support:data.input->support[atom])
        {
            const auto basis=joint::EvaluateKernel(support.squared_distance,b,2.5);
            contribution(static_cast<Eigen::Index>(support.row))+=a*basis.gaussian+c*basis.charge;
        }
    }
    const joint::Vector prediction=target_prediction+halo_prediction;
    const joint::Vector residual=prediction-data.y;
    std::size_t block_successes{},block_failures{};
    for(const auto & block:fit.search.blocks)
        if(block.status=="accepted" || block.status=="unchanged") ++block_successes;
        else ++block_failures;
    const auto block_visits=block_successes+block_failures;
    const bool ac_kkt_passed=primary.certificate.available && primary.certificate.projected_kkt<=1e-10;
    const bool free_width_gradient_passed=free_target_gradient<=1e-12;
    const bool confirmed=fit.search.search_converged && ac_kkt_passed && free_width_gradient_passed && replay.passed;
    return { {"name",name},{"fixed_width_atoms",fixed_width_ids.size()},
        {"restricted_search_converged",fit.search.search_converged},{"formal_joint_converged",false},
        {"restricted_fit_qualified",confirmed},{"stop_reason",fit.stop_reason},
        {"sweeps",fit.search.sweep_count},{"block_visits",fit.search.total_block_solves},
        {"block_successes",block_successes},{"block_failures",block_failures},
        {"profile_evaluations",fit.search.total_profile_evaluations},
        {"accepted_width_updates",fit.search.total_accepted_local_updates},
        {"block_success_rate",block_visits ? Number(static_cast<double>(block_successes)/block_visits) : json::value(nullptr)},
        {"objective",objective ? Number(*objective) : json::value(nullptr)},
        {"global_ac_kkt",Number(primary.certificate.projected_kkt)},
        {"all_width_gradient_inf_norm",Number(all_width_gradient)},
        {"target_width_gradient_inf_norm",Number(target_gradient)},
        {"free_target_width_gradient_inf_norm",Number(free_target_gradient)},
        {"independent_replay_passed",replay.passed},{"independent_replay_reason",replay.reason},
        {"target_count",target.a_count},{"halo_count",halo.a_count},
        {"target_rmse_A",Number(Rmse(target,true))},{"target_relative_rmse_A",Number(RelativeRmse(target,true))},
        {"target_rmse_C",Number(Rmse(target,false))},{"target_relative_rmse_C",Number(RelativeRmse(target,false))},
        {"target_rmse_B",Number(WidthRmse(target))},{"halo_rmse_A",Number(Rmse(halo,true))},
        {"halo_relative_rmse_A",Number(RelativeRmse(halo,true))},
        {"halo_rmse_C",Number(Rmse(halo,false))},{"halo_relative_rmse_C",Number(RelativeRmse(halo,false))},
        {"halo_rmse_B",Number(WidthRmse(halo))},
        {"target_prediction_norm",Number(target_prediction.norm())},
        {"halo_prediction_norm",Number(halo_prediction.norm())},
        {"prediction_rmse",Number(residual.norm()/std::sqrt(static_cast<double>(residual.size())))},
        {"estimated_peak_bytes",fit.search.memory_estimate.estimated_peak_bytes},
        {"wall_seconds",seconds},{"process_peak_rss_bytes",PeakRss()},
        {"qr_factorization_count",nullptr}};
}
void Run(const std::filesystem::path & model_path,const std::filesystem::path & map_path,
    const std::filesystem::path & output,std::size_t maximum_sweeps)
{
    if(std::filesystem::exists(output)) throw std::runtime_error("Use a fresh experiment output directory.");
    std::filesystem::create_directories(output);
    auto model=rhbm_gem::ReadModel(model_path.string()); auto map=rhbm_gem::ReadMap(map_path.string());
    model->SelectAtoms([](const rhbm_gem::AtomObject & atom) {
        return atom.GetElement()!=Element::HYDROGEN && atom.GetChainID()=="A";
    });
    const auto problem=core::BuildJointProblem(*map,*model);
    const auto initialization=joint::InitializeContributors(*map,*model,problem);
    if(!initialization.valid) throw std::runtime_error("First-stage initialization failed: "+initialization.reason);
    const auto & data=core::JointProblemAccess::Get(problem);
    if(!data.input->selection_domain) throw std::runtime_error("Selected target identities are unavailable.");
    const auto truths=Truths(map_path.string()+".simulation.json");
    std::vector<bool> targets(data.input->support.size());
    for(const auto atom:data.input->selection_domain->target_indices) targets.at(atom)=true;
    if(initialization.b.size()!=static_cast<Eigen::Index>(targets.size()))
        throw std::runtime_error("First-stage width count does not match FullABC contributors.");
    std::vector<double> true_widths(targets.size());
    std::vector<std::string> all_ids,halo_ids;
    for(std::size_t atom=0;atom<targets.size();++atom)
    {
        const auto id=data.input->atom_ids[atom]; const auto found=truths.find(id);
        if(found==truths.end()) throw std::runtime_error("No simulation truth for contributor "+id);
        true_widths[atom]=found->second.width; all_ids.push_back(id);
        if(!targets[atom]) halo_ids.push_back(id);
    }
    std::vector<double> oracle_halo_widths=initialization.b;
    for(std::size_t atom=0;atom<targets.size();++atom)
        if(!targets[atom]) oracle_halo_widths[atom]=true_widths[atom];
    const std::vector<double> seed_halo_widths=initialization.b;
    json::array experiments;
    experiments.push_back(RunCase("all-widths-oracle",data,truths,targets,true_widths,
        std::move(all_ids),maximum_sweeps,output));
    experiments.push_back(RunCase("oracle-halo-widths",data,truths,targets,oracle_halo_widths,
        std::move(halo_ids),maximum_sweeps,output));
    std::vector<std::string> seed_halo_ids;
    for(std::size_t atom=0;atom<targets.size();++atom)
        if(!targets[atom]) seed_halo_ids.push_back(data.input->atom_ids[atom]);
    experiments.push_back(RunCase("seed-halo-widths",data,truths,targets,seed_halo_widths,
        std::move(seed_halo_ids),maximum_sweeps,output));
    std::ofstream summary(output/"summary.json"); summary.exceptions(std::ios::failbit|std::ios::badbit);
    summary<<json::serialize(json::object{{"model",model_path.string()},{"map",map_path.string()},
        {"simulation_metadata",map_path.string()+".simulation.json"},{"selected_targets",data.input->selection_domain->target_indices.size()},
        {"contributors",data.input->support.size()},{"observations",data.input->observations.size()},
        {"maximum_sweeps",maximum_sweeps},{"experiments",std::move(experiments)},
        {"process_peak_rss_bytes",PeakRss()}})<<'\n';
}
}
int main(int argc,char ** argv)
{
    try
    {
        if(argc<4 || argc>5) throw std::invalid_argument(
            "Usage: joint_fixed_width_oracle MODEL MAP OUTPUT_DIR [MAX_SWEEPS]");
        const auto sweeps=argc==5 ? static_cast<std::size_t>(std::stoull(argv[4])) : 200;
        if(sweeps==0 || sweeps>200) throw std::invalid_argument("MAX_SWEEPS must be between 1 and 200.");
        Run(argv[1],argv[2],argv[3],sweeps); return 0;
    }
    catch(const std::exception & error) {std::cerr<<error.what()<<'\n'; return 1;}
}
