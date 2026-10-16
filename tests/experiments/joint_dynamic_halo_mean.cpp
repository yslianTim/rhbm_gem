#include "core/detail/FirstStageInitialization.hpp"
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
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace {
namespace core=rhbm_gem::core;
namespace joint=core::joint_component;
namespace json=boost::json;
using Parameters=joint_test::dynamic_halo::Parameters;

struct Truth {Parameters value;};
struct Metrics
{
    double mae{},rmse{},median{},p90{},p95{},p99{},maximum{};
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

std::vector<double> NumberArray(const json::value & value)
{
    std::vector<double> result;
    for(const auto & item:value.as_array()) result.push_back(json::value_to<double>(item));
    return result;
}

std::vector<std::string> StringArray(const json::value & value)
{
    std::vector<std::string> result;
    for(const auto & item:value.as_array()) result.push_back(json::value_to<std::string>(item));
    return result;
}

void RequireEqual(const std::vector<std::string> & a,const std::vector<std::string> & b,const char * name)
{
    if(a!=b) throw std::runtime_error(std::string("Production result ")+name+" identity mismatch.");
}

std::string Csv(std::string_view value)
{
    std::string result="\"";
    for(const auto c:value) {if(c=='\"') result+='\"'; result+=c;}
    result+='\"'; return result;
}

std::vector<joint_test::dynamic_halo::Block> BuildBlocks(const core::JointProblem & problem)
{
    const auto & data=core::JointProblemAccess::Get(problem);
    std::vector<joint_test::dynamic_halo::Block> blocks;
    for(std::size_t component_index=0;component_index<data.partition.components.size();++component_index)
    {
        const auto view=joint::MakeComponentProblemView(data.layout,data.partition.components[component_index]);
        joint::StructuralTopology topology(view.domain,view.atom_ids,view.local_layout);
        const auto partition=topology.BuildCores(12,false);
        for(const auto & core:partition.cores)
        {
            joint_test::dynamic_halo::Block block;
            block.id=core.id; block.component=component_index;
            block.targets.reserve(core.atoms.size());
            for(const auto local_atom:core.atoms)
            {
                const auto target_local=static_cast<std::size_t>(view.local_to_parent_atoms.at(
                    static_cast<std::size_t>(local_atom)));
                block.targets.push_back(data.target_to_contributor.at(target_local));
            }
            blocks.push_back(std::move(block));
        }
    }
    return blocks;
}

Metrics ErrorMetrics(std::vector<double> errors)
{
    if(errors.empty()) throw std::runtime_error("Cannot summarize an empty parameter error set.");
    for(auto & value:errors) value=std::abs(value);
    std::sort(errors.begin(),errors.end());
    long double sum{},sum_squared{};
    for(const auto value:errors) {sum+=value; sum_squared+=static_cast<long double>(value)*value;}
    const auto quantile=[&](double p) {return errors[static_cast<std::size_t>(
        std::ceil(p*static_cast<double>(errors.size()))-1)];};
    return {static_cast<double>(sum/errors.size()),
        std::sqrt(static_cast<double>(sum_squared/errors.size())),
        errors[(errors.size()-1)/2],quantile(.90),quantile(.95),quantile(.99),errors.back()};
}

json::object ToJson(const Metrics & value)
{
    return {{"mae",value.mae},{"rmse",value.rmse},{"median_absolute_error",value.median},
        {"p90_absolute_error",value.p90},{"p95_absolute_error",value.p95},
        {"p99_absolute_error",value.p99},{"maximum_absolute_error",value.maximum}};
}

void Add(std::vector<double> & field,const core::JointProblemInput & input,std::size_t atom,
    const Parameters & parameter)
{
    for(const auto & support:input.support.at(atom))
    {
        const auto basis=joint::EvaluateKernel(support.squared_distance,parameter.b,2.5);
        field.at(support.row)+=parameter.a*basis.gaussian+parameter.c*basis.charge;
    }
}

Parameters FitUnits(Parameters parameter,double amplitude_scale)
{
    parameter.a*=amplitude_scale; parameter.c*=amplitude_scale;
    return parameter;
}

double Rmse(const std::vector<double> & a,const std::vector<double> & b)
{
    if(a.size()!=b.size() || a.empty()) throw std::invalid_argument("Invalid RMSE vectors.");
    long double sum{};
    for(std::size_t row=0;row<a.size();++row)
    {
        const long double difference=static_cast<long double>(a[row])-b[row];
        sum+=difference*difference;
    }
    return std::sqrt(static_cast<double>(sum/static_cast<long double>(a.size())));
}

double Objective(const std::vector<double> & prediction,const core::JointProblemInput & input,double scale)
{
    long double sum{};
    for(std::size_t row=0;row<prediction.size();++row)
    {
        const long double difference=static_cast<long double>(prediction[row])-input.observations[row];
        sum+=difference*difference;
    }
    return static_cast<double>(.5L*sum/(static_cast<long double>(scale)*scale));
}

void WriteText(const std::filesystem::path & path,const std::string & text)
{
    std::ofstream output(path); output.exceptions(std::ios::failbit|std::ios::badbit); output<<text;
}

struct InputState
{
    std::vector<Parameters> initial,production,first_stage,truth,static_initial_mean;
    std::vector<std::string> target_block;
    std::vector<std::string> halo_owner;
    std::vector<std::size_t> halo_overlap,halo_intersections;
};

InputState ReadStates(const core::JointProblem & problem,const core::JointInitialization & initialization,
    const joint_test::dynamic_halo::Ownership & ownership,const json::value & production,
    const std::unordered_map<std::string,Truth> & truth)
{
    const auto & input=problem.Input();
    const auto target_indices=input.selection_domain->target_indices;
    const auto initial_width=initialization.b;
    if(initial_width.size()!=input.atom_ids.size())
        throw std::runtime_error("FirstStage widths do not cover the contributor identities.");
    InputState result;
    result.initial.resize(input.atom_ids.size()); result.production.resize(input.atom_ids.size());
    result.first_stage.resize(input.atom_ids.size()); result.truth.resize(input.atom_ids.size());
    result.static_initial_mean.resize(input.atom_ids.size());
    result.target_block.resize(input.atom_ids.size()); result.halo_owner.resize(input.atom_ids.size());
    result.halo_overlap.resize(input.atom_ids.size()); result.halo_intersections.resize(input.atom_ids.size());

    const auto & target_state=production.at("assembled_state");
    const auto target_ac=NumberArray(target_state.at("ac"));
    const auto target_b=NumberArray(target_state.at("b"));
    if(target_ac.size()!=target_indices.size()*2 || target_b.size()!=target_indices.size())
        throw std::runtime_error("Production target state has inconsistent dimensions.");
    const auto & fixed_halos=production.at("fixed_halos").as_array();
    std::vector<bool> fixed_seen(input.atom_ids.size());
    for(const auto & halo:fixed_halos)
    {
        const auto contributor=json::value_to<std::size_t>(halo.at("contributor_index"));
        if(contributor>=input.atom_ids.size() || fixed_seen[contributor])
            throw std::runtime_error("Production fixed-halo mapping is duplicated or out of range.");
        fixed_seen[contributor]=true;
        if(json::value_to<std::string>(halo.at("atom_id"))!=input.atom_ids[contributor] ||
            json::value_to<double>(halo.at("charge"))!=0.0 ||
            !json::value_to<bool>(halo.at("first_stage_fit_qualified")))
            throw std::runtime_error("Production fixed-halo contract does not match this research input.");
        result.production[contributor]={json::value_to<double>(halo.at("amplitude")),
            json::value_to<double>(halo.at("width")),json::value_to<double>(halo.at("charge"))};
    }

    std::vector<std::size_t> target_local(input.atom_ids.size(),std::numeric_limits<std::size_t>::max());
    for(std::size_t local=0;local<target_indices.size();++local) target_local[target_indices[local]]=local;
    for(std::size_t atom=0;atom<input.atom_ids.size();++atom)
    {
        const auto found=truth.find(input.atom_ids[atom]);
        if(found==truth.end()) throw std::runtime_error("Missing truth for contributor "+input.atom_ids[atom]);
        const auto & init=initialization.atoms.at(atom);
        if(init.id!=input.atom_ids[atom]) throw std::runtime_error("FirstStage identity mapping mismatch.");
        result.initial[atom]={0.0,initial_width[atom],0.0};
        result.first_stage[atom]=result.production[atom];
        result.truth[atom]=found->second.value;
        if(target_local[atom]!=std::numeric_limits<std::size_t>::max())
        {
            if(fixed_seen[atom]) throw std::runtime_error("Production classified a target as a fixed halo.");
            const auto local=target_local[atom];
            result.production[atom]={target_ac[2*local],target_b[local],target_ac[2*local+1]};
            result.first_stage[atom]={init.mdpde[0],initial_width[atom],0.0};
        }
        else if(!fixed_seen[atom]) throw std::runtime_error("Production result is missing a supported halo.");
    }
    for(const auto & block:ownership.blocks)
        for(const auto target:block.targets) result.target_block[target]=block.id;
    for(const auto & owner:ownership.halos)
    {
        result.halo_owner[owner.halo]=ownership.blocks[owner.owner_block].id;
        result.halo_overlap[owner.halo]=owner.owner_overlap;
        result.halo_intersections[owner.halo]=owner.intersecting_blocks.size();
        std::vector<Parameters> initial_targets;
        for(const auto target:ownership.blocks[owner.owner_block].targets)
            initial_targets.push_back(result.initial[target]);
        result.static_initial_mean[owner.halo]=joint_test::dynamic_halo::Mean(initial_targets,true);
    }
    return result;
}

json::object FieldMetrics(const std::string & name,const std::vector<double> & halo,
    const std::vector<double> & target,const std::vector<double> & halo_truth,
    const std::vector<double> & target_truth,const core::JointProblemInput & input,double scale,
    double normalization_divisor)
{
    std::vector<double> full(halo.size()),full_truth(halo.size());
    for(std::size_t row=0;row<full.size();++row)
    {full[row]=halo[row]+target[row]; full_truth[row]=halo_truth[row]+target_truth[row];}
    return {{"name",name},{"halo_field_rmse_normalized_map_units",Rmse(halo,halo_truth)},
        {"halo_field_rmse_raw_map_units",Rmse(halo,halo_truth)*normalization_divisor},
        {"target_prediction_rmse_normalized_map_units",Rmse(target,target_truth)},
        {"target_prediction_rmse_raw_map_units",Rmse(target,target_truth)*normalization_divisor},
        {"full_prediction_rmse_normalized_map_units",Rmse(full,full_truth)},
        {"full_prediction_rmse_raw_map_units",Rmse(full,full_truth)*normalization_divisor},
        {"objective",Objective(full,input,scale)},
        {"full_residual_norm",[&] {
            long double sum{};
            for(std::size_t row=0;row<full.size();++row)
            {const long double d=static_cast<long double>(full[row])-input.observations[row]; sum+=d*d;}
            return std::sqrt(static_cast<double>(sum));
        }()}};
}

Metrics BlockSpread(const std::vector<Parameters> & values,const std::vector<joint_test::dynamic_halo::Block> & blocks,
    int parameter,double amplitude_scale)
{
    std::vector<double> spreads;
    for(const auto & block:blocks)
    {
        std::vector<double> sample;
        for(const auto atom:block.targets)
            sample.push_back(parameter==0 ? values[atom].a/amplitude_scale :
                parameter==1 ? values[atom].b : values[atom].c/amplitude_scale);
        const double mean=std::accumulate(sample.begin(),sample.end(),0.0)/static_cast<double>(sample.size());
        long double squared{};
        for(const auto value:sample) squared+=(value-mean)*(value-mean);
        spreads.push_back(std::sqrt(static_cast<double>(squared/static_cast<long double>(sample.size()))));
    }
    return ErrorMetrics(std::move(spreads));
}

int Run(const std::filesystem::path & model_path,const std::filesystem::path & map_path,
    const std::filesystem::path & production_path,const std::filesystem::path & output_dir,
    const std::string & baseline_commit)
{
    auto model=rhbm_gem::ReadModel(model_path.string());
    auto map=rhbm_gem::ReadMap(map_path.string());
    const auto divisor=map->MapValueArrayNormalization();
    const double amplitude_scale=1.0/divisor.value_or(1.0);
    if(!divisor) throw std::runtime_error("The production CLI normalization could not be reproduced.");
    model->SelectAtoms([](const rhbm_gem::AtomObject & atom) {
        return atom.GetElement()!=Element::HYDROGEN && atom.GetChainID()=="A";
    });
    const auto problem=core::BuildJointProblem(*map,*model);
    const auto initialization=joint::InitializeContributors(*map,*model,problem);
    if(!initialization.valid) throw std::runtime_error("Production FirstStage initialization is not qualified.");
    const auto production=ReadJson(production_path);
    const auto & input=problem.Input();
    RequireEqual(input.atom_ids,StringArray(production.at("atom_ids")),"contributor");
    RequireEqual(input.row_ids,StringArray(production.at("row_ids")),"observation");
    if(production.at("runtime_convergence")!="passed" || production.at("target_runtime_convergence")!="passed")
        throw std::runtime_error("The supplied Production endpoint is not qualified.");
    if(std::abs(json::value_to<double>(production.at("metadata").at("map_normalization").at("divisor"))-
        *divisor)>1e-14)
        throw std::runtime_error("Production and research map normalization do not match.");
    if(std::abs(json::value_to<double>(production.at("observation_scale"))-problem.ObservationScale())>
        1e-10*problem.ObservationScale())
        throw std::runtime_error("Production and research observation scales do not match.");

    // Ownership is derived before simulation truth is opened or inspected.
    const auto ownership=joint_test::dynamic_halo::BuildOwnership(input,BuildBlocks(problem));
    if(ownership.cross_component_halos)
        throw std::runtime_error("unsupported-cross-component-halo-coupling");
    const auto truth=ReadTruth(map_path.string()+".simulation.json");
    const auto state=ReadStates(problem,initialization,ownership,production,truth);
    const auto rows=input.observations.size();
    std::vector<double> target_prediction(rows),target_truth(rows),halo_first(rows),halo_ab(rows),
        halo_abc(rows),halo_truth(rows);
    std::vector<std::vector<Parameters>> block_values(ownership.blocks.size());
    for(std::size_t block=0;block<ownership.blocks.size();++block)
        for(const auto target:ownership.blocks[block].targets)
        {
            block_values[block].push_back(state.production[target]);
        }
    std::vector<Parameters> halo_ab_parameters(input.atom_ids.size()),halo_abc_parameters(input.atom_ids.size());
    for(std::size_t atom=0;atom<input.atom_ids.size();++atom)
    {
        if(state.target_block[atom].size())
        {
            Add(target_prediction,input,atom,state.production[atom]);
            Add(target_truth,input,atom,FitUnits(state.truth[atom],amplitude_scale));
        }
    }
    std::array<std::array<std::vector<double>,3>,3> support_pred;
    std::array<std::vector<double>,3> support_truth;
    std::array<std::vector<bool>,3> support_rows;
    for(std::size_t bucket=0;bucket<3;++bucket)
    {
        for(auto & field:support_pred[bucket]) field.assign(rows,0);
        support_truth[bucket].assign(rows,0); support_rows[bucket].assign(rows,false);
    }
    std::array<std::size_t,3> support_counts{};
    std::array<std::array<std::vector<double>,3>,3> halo_parameter_errors;
    for(const auto & owner:ownership.halos)
    {
        const auto halo=owner.halo;
        const auto bucket=std::min<std::size_t>(2,input.support[halo].size()-1);
        ++support_counts[bucket];
        const auto & members=block_values.at(owner.owner_block);
        halo_ab_parameters[halo]=joint_test::dynamic_halo::Mean(members,false);
        halo_abc_parameters[halo]=joint_test::dynamic_halo::Mean(members,true);
        Add(halo_first,input,halo,state.first_stage[halo]);
        Add(halo_ab,input,halo,halo_ab_parameters[halo]);
        Add(halo_abc,input,halo,halo_abc_parameters[halo]);
        Add(halo_truth,input,halo,FitUnits(state.truth[halo],amplitude_scale));
        Add(support_pred[bucket][0],input,halo,state.first_stage[halo]);
        Add(support_pred[bucket][1],input,halo,halo_ab_parameters[halo]);
        Add(support_pred[bucket][2],input,halo,halo_abc_parameters[halo]);
        Add(support_truth[bucket],input,halo,FitUnits(state.truth[halo],amplitude_scale));
        for(const auto & support:input.support[halo]) support_rows[bucket][support.row]=true;
        halo_parameter_errors[0][0].push_back(state.first_stage[halo].a/amplitude_scale-state.truth[halo].a);
        halo_parameter_errors[0][1].push_back(state.first_stage[halo].b-state.truth[halo].b);
        halo_parameter_errors[0][2].push_back(state.first_stage[halo].c/amplitude_scale-state.truth[halo].c);
        halo_parameter_errors[1][0].push_back(halo_ab_parameters[halo].a/amplitude_scale-state.truth[halo].a);
        halo_parameter_errors[1][1].push_back(halo_ab_parameters[halo].b-state.truth[halo].b);
        halo_parameter_errors[1][2].push_back(halo_ab_parameters[halo].c/amplitude_scale-state.truth[halo].c);
        halo_parameter_errors[2][0].push_back(halo_abc_parameters[halo].a/amplitude_scale-state.truth[halo].a);
        halo_parameter_errors[2][1].push_back(halo_abc_parameters[halo].b-state.truth[halo].b);
        halo_parameter_errors[2][2].push_back(halo_abc_parameters[halo].c/amplitude_scale-state.truth[halo].c);
    }

    const std::array<std::string,3> bucket_names{"one","two","three-plus"};
    json::array support_metrics;
    for(std::size_t bucket=0;bucket<3;++bucket)
    {
        std::array<long double,3> squared{};
        std::size_t count{};
        for(std::size_t row=0;row<rows;++row) if(support_rows[bucket][row])
        {
            for(std::size_t group=0;group<3;++group)
            {
                const long double d=static_cast<long double>(support_pred[bucket][group][row])-support_truth[bucket][row];
                squared[group]+=d*d;
            }
            ++count;
        }
        json::array model_errors;
        const std::array<std::string,3> model_names{"production-first-stage", "postfit-AB", "postfit-ABC"};
        for(std::size_t group=0;group<3;++group)
            model_errors.emplace_back(json::object{{"model",model_names[group]},
                {"halo_field_rmse_normalized_map_units",count>0 ?
                    std::sqrt(static_cast<double>(squared[group]/static_cast<long double>(count))) : 0.0},
                {"halo_field_rmse_raw_map_units",count>0 ?
                    std::sqrt(static_cast<double>(squared[group]/static_cast<long double>(count)))*(*divisor) : 0.0}});
        support_metrics.emplace_back(json::object{{"support_rows",bucket_names[bucket]},
            {"halo_count",support_counts[bucket]},{"active_observation_rows",count},
            {"models",std::move(model_errors)}});
    }

    const auto production_fields=FieldMetrics("E0 Production FirstStage A/B, C=0",halo_first,
        target_prediction,halo_truth,target_truth,input,problem.ObservationScale(),*divisor);
    const auto ab_fields=FieldMetrics("P1 post-fit block mean AB",halo_ab,target_prediction,
        halo_truth,target_truth,input,problem.ObservationScale(),*divisor);
    const auto abc_fields=FieldMetrics("P1 post-fit block mean ABC",halo_abc,target_prediction,
        halo_truth,target_truth,input,problem.ObservationScale(),*divisor);
    json::array block_means;
    for(std::size_t block=0;block<ownership.blocks.size();++block)
    {
        const auto mean_ab=joint_test::dynamic_halo::Mean(block_values[block],false);
        const auto mean_abc=joint_test::dynamic_halo::Mean(block_values[block],true);
        std::vector<Parameters> initial_members;
        for(const auto target:ownership.blocks[block].targets) initial_members.push_back(state.initial[target]);
        const auto mean_initial=joint_test::dynamic_halo::Mean(initial_members,true);
        block_means.emplace_back(json::object{{"block_id",ownership.blocks[block].id},
            {"component",ownership.blocks[block].component},{"target_count",ownership.blocks[block].targets.size()},
            {"halo_count",ownership.block_halos[block].size()},
            {"initial_mean_A",mean_initial.a/amplitude_scale},{"initial_mean_B",mean_initial.b},
            {"initial_mean_C",mean_initial.c/amplitude_scale},
            {"mean_A",mean_abc.a/amplitude_scale},{"mean_B",mean_abc.b},{"mean_C",mean_abc.c/amplitude_scale},
            {"mean_A_AB",mean_ab.a/amplitude_scale},{"mean_B_AB",mean_ab.b}});
    }
    json::array ownership_halos;
    std::size_t shared_halos{};
    std::vector<std::string> target_block(input.atom_ids.size());
    for(std::size_t block=0;block<ownership.blocks.size();++block)
        for(const auto target:ownership.blocks[block].targets) target_block[target]=ownership.blocks[block].id;
    for(const auto & owner:ownership.halos)
    {
        if(owner.intersecting_blocks.size()>1) ++shared_halos;
        json::array intersecting;
        for(const auto block:owner.intersecting_blocks) intersecting.emplace_back(ownership.blocks[block].id);
        ownership_halos.emplace_back(json::object{{"atom_id",input.atom_ids[owner.halo]},
            {"owner_block",ownership.blocks[owner.owner_block].id},{"owner_overlap_rows",owner.owner_overlap},
            {"observation_support_rows",input.support[owner.halo].size()},
            {"intersecting_blocks",std::move(intersecting)},
            {"overlapping_target_count",owner.overlapping_targets.size()}});
    }
    std::size_t total_dynamic_rows{},total_target_rows{},maximum_dynamic_rows{},maximum_expansion{};
    for(std::size_t block=0;block<ownership.blocks.size();++block)
    {
        total_dynamic_rows+=ownership.block_dynamic_rows[block].size();
        maximum_dynamic_rows=std::max(maximum_dynamic_rows,ownership.block_dynamic_rows[block].size());
        std::set<std::size_t> target_rows;
        for(const auto target:ownership.blocks[block].targets)
            for(const auto & support:input.support[target]) target_rows.insert(support.row);
        total_target_rows+=target_rows.size();
        maximum_expansion=std::max(maximum_expansion,
            ownership.block_dynamic_rows[block].size()-target_rows.size());
    }

    std::filesystem::create_directories(output_dir);
    {
        std::ofstream output(output_dir/"rows.csv"); output.exceptions(std::ios::failbit|std::ios::badbit);
        output<<"row_index,row_id,observation\n"<<std::setprecision(17);
        for(std::size_t row=0;row<rows;++row)
            output<<row<<','<<Csv(input.row_ids[row])<<','<<input.observations[row]<<'\n';
    }
    {
        std::ofstream output(output_dir/"supports.csv"); output.exceptions(std::ios::failbit|std::ios::badbit);
        output<<"contributor_index,atom_id,role,row_index,row_id,squared_distance\n"<<std::setprecision(17);
        const auto targets=input.selection_domain->target_indices;
        for(std::size_t atom=0;atom<input.atom_ids.size();++atom)
            for(const auto & support:input.support[atom])
            {
                const bool target=std::binary_search(targets.begin(),targets.end(),atom);
                output<<atom<<','<<Csv(input.atom_ids[atom])<<','<<(target?"target":"halo")<<','
                    <<support.row<<','<<Csv(input.row_ids[support.row])<<','<<support.squared_distance<<'\n';
            }
    }
    {
        std::ofstream output(output_dir/"contributors.csv"); output.exceptions(std::ios::failbit|std::ios::badbit);
        output<<"contributor_index,atom_id,role,target_block,halo_owner,owner_overlap_rows,intersecting_block_count,support_rows,"
            "initial_A,initial_B,initial_C,production_A,production_B,production_C,first_stage_qualified,truth_A,truth_B,truth_C,"
            "static_initial_mean_A,static_initial_mean_B,static_initial_mean_C,"
            "postfit_AB_A,postfit_AB_B,postfit_AB_C,postfit_ABC_A,postfit_ABC_B,postfit_ABC_C\n"<<std::setprecision(17);
        for(std::size_t atom=0;atom<input.atom_ids.size();++atom)
        {
            const bool target=!target_block[atom].empty();
            const auto & init=initialization.atoms[atom];
            output<<atom<<','<<Csv(input.atom_ids[atom])<<','<<(target?"target":"halo")<<','
                <<Csv(target?target_block[atom]:"")<<','<<Csv(target?"":state.halo_owner[atom])<<','
                <<state.halo_overlap[atom]<<','<<state.halo_intersections[atom]<<','<<input.support[atom].size()<<','
                <<state.initial[atom].a/amplitude_scale<<','<<state.initial[atom].b<<','<<state.initial[atom].c/amplitude_scale<<','
                <<state.production[atom].a/amplitude_scale<<','<<state.production[atom].b<<','
                <<state.production[atom].c/amplitude_scale<<','<<(target?"":init.first_stage_fit_qualified?"true":"false")<<','
                <<state.truth[atom].a<<','<<state.truth[atom].b<<','<<state.truth[atom].c;
            if(target) output<<",,,,,,,,,\n";
            else output<<','<<state.static_initial_mean[atom].a/amplitude_scale<<','
                <<state.static_initial_mean[atom].b<<','<<state.static_initial_mean[atom].c/amplitude_scale<<','
                <<halo_ab_parameters[atom].a/amplitude_scale<<','<<halo_ab_parameters[atom].b<<",0,"
                <<halo_abc_parameters[atom].a/amplitude_scale<<','<<halo_abc_parameters[atom].b<<','
                <<halo_abc_parameters[atom].c/amplitude_scale<<'\n';
        }
    }

    const auto & production_component=production.at("components").as_array().front();
    const auto component_evidence=[&](const std::string & name) {
        for(const auto & evidence:production_component.at("evidence").as_array())
            if(json::value_to<std::string>(evidence.at("name"))==name)
                return evidence.at("value").is_null() ? json::value(nullptr) : evidence.at("value");
        return json::value(nullptr);
    };
    const auto param_error_set=[&](std::size_t group) {
        return json::object{{"A_raw_units",ToJson(ErrorMetrics(halo_parameter_errors[group][0]))},
            {"B_angstrom",ToJson(ErrorMetrics(halo_parameter_errors[group][1]))},
            {"C_raw_units",ToJson(ErrorMetrics(halo_parameter_errors[group][2]))}};
    };
    const auto baseline_halo_rmse=json::value_to<double>(production_fields.at("halo_field_rmse_normalized_map_units"));
    const auto ab_halo_rmse=json::value_to<double>(ab_fields.at("halo_field_rmse_normalized_map_units"));
    const auto abc_halo_rmse=json::value_to<double>(abc_fields.at("halo_field_rmse_normalized_map_units"));
    json::array metrics{production_fields,ab_fields,abc_fields};
    const auto & fixed_halos=production.at("fixed_halos").as_array();
    const auto qualified_halos=std::count_if(fixed_halos.begin(),fixed_halos.end(),[](const auto & halo) {
        return json::value_to<bool>(halo.at("first_stage_fit_qualified"));
    });
    json::object summary{
        {"contract","joint-dynamic-halo-mean-postfit-v1"},
        {"baseline_commit",baseline_commit},
        {"production_result",production_path.string()},
        {"model_sha256",production.at("metadata").at("model_sha256")},
        {"map_sha256",production.at("metadata").at("map_sha256")},
        {"simulation_truth_sidecar",map_path.string()+".simulation.json"},
        {"normalization_divisor",*divisor},{"observation_scale",problem.ObservationScale()},
        {"observation_rows",rows},{"target_atoms",input.selection_domain->target_indices.size()},
        {"halo_atoms",ownership.halos.size()},{"block_count",ownership.blocks.size()},
        {"p0",json::object{{"branch","develop"},{"commit",baseline_commit},{"cli_exit_status",0},
            {"targets",input.selection_domain->target_indices.size()},{"qualified_targets",production_component.at("atoms").as_array().size()},
            {"fixed_halos",ownership.halos.size()},{"qualified_fixed_halos",qualified_halos},
            {"observation_rows",rows},{"objective",production.at("objective")},
            {"sweeps",production_component.at("search_sweeps")},
            {"profile_evaluations",production_component.at("profile_evaluations")},
            {"accepted_updates",production_component.at("accepted_updates")},
            {"global_ac_kkt",component_evidence("conditional-ac-optimality")},
            {"global_width_gradient",component_evidence("width-stationarity")},
            {"confirmed_log_width_change",production_component.at("candidate_diagnostic").at("candidate_log_width_change_inf_norm")},
            {"criteria_id",production_component.at("candidate_diagnostic").at("convergence_criteria_id")},
            {"convergence_status",production.at("runtime_convergence")},
            {"independent_replay",[&] {
            for(const auto & check:production.at("evidence").as_array())
                if(json::value_to<std::string>(check.at("name"))=="independent-replay") return check.at("status");
            return json::value("unavailable");
            }()},
            {"full_observation_replay",[&] {
                for(const auto & check:production.at("evidence").as_array())
                    if(json::value_to<std::string>(check.at("name"))=="full-observation-replay") return check.at("status");
                return json::value("unavailable");
            }()},
            {"search_seconds",production.at("costs").at("search_seconds")}}},
        {"ownership",json::object{{"supported_halos_with_one_owner",ownership.halos.size()},
            {"cross_component_halos",ownership.cross_component_halos},{"shared_halos",shared_halos},
            {"row_to_block_memberships",ownership.row_block_entries},{"total_target_block_rows",total_target_rows},
            {"total_dynamic_affected_rows",total_dynamic_rows},{"maximum_dynamic_block_rows",maximum_dynamic_rows},
            {"maximum_affected_row_expansion",maximum_expansion},
            {"halos",std::move(ownership_halos)},{"blocks",std::move(block_means)}}},
        {"p1",json::object{{"interpretation","Post-fit Halo Mean keeps every Production Target A/B/C fixed; it evaluates halo-field reconstruction only."},
            {"common_initial_target_state","The Production FixedNeighbor state starts with A=C=0 and FirstStage B; the saved static mean uses exactly that state."},
            {"models",std::move(metrics)},
            {"halo_field_rmse_improvement_percent_vs_production",json::object{
                {"AB",100.0*(baseline_halo_rmse-ab_halo_rmse)/baseline_halo_rmse},
                {"ABC",100.0*(baseline_halo_rmse-abc_halo_rmse)/baseline_halo_rmse}}},
            {"target_rmse_A_raw_units",ErrorMetrics([&] {std::vector<double> e; for(const auto atom:input.selection_domain->target_indices) e.push_back(state.production[atom].a/amplitude_scale-state.truth[atom].a); return e;}()).rmse},
            {"target_rmse_B_angstrom",ErrorMetrics([&] {std::vector<double> e; for(const auto atom:input.selection_domain->target_indices) e.push_back(state.production[atom].b-state.truth[atom].b); return e;}()).rmse},
            {"target_rmse_C_raw_units",ErrorMetrics([&] {std::vector<double> e; for(const auto atom:input.selection_domain->target_indices) e.push_back(state.production[atom].c/amplitude_scale-state.truth[atom].c); return e;}()).rmse},
            {"halo_parameter_errors",json::object{{"production_first_stage",param_error_set(0)},
                {"postfit_AB_mean",param_error_set(1)},{"postfit_ABC_mean",param_error_set(2)}}},
            {"support_strata",std::move(support_metrics)},
            {"target_block_parameter_spread",json::object{
                {"A_raw_units",ToJson(BlockSpread(state.production,ownership.blocks,0,amplitude_scale))},
                {"B_angstrom",ToJson(BlockSpread(state.production,ownership.blocks,1,amplitude_scale))},
                {"C_raw_units",ToJson(BlockSpread(state.production,ownership.blocks,2,amplitude_scale))}}},
            {"target_rmse_is_unchanged_across_postfit_models",true},
            {"target_jacobian_projection","not-computed"}}}
    };
    WriteText(output_dir/"summary.json",json::serialize(summary)+"\n");
    json::object checkpoint{{"contract","joint-dynamic-halo-mean-checkpoint-v1"},
        {"branch","develop"},{"baseline_commit",baseline_commit},
        {"model_sha256",production.at("metadata").at("model_sha256")},
        {"map_sha256",production.at("metadata").at("map_sha256")},
        {"observation_scale",problem.ObservationScale()},{"normalization_divisor",*divisor},
        {"row_identity_file","rows.csv"},{"support_mapping_file","supports.csv"},
        {"parameter_and_owner_file","contributors.csv"},{"observation_rows",rows},
        {"contributors",input.atom_ids.size()},{"targets",input.selection_domain->target_indices.size()},
        {"halos",ownership.halos.size()},{"blocks",ownership.blocks.size()},
        {"support_entries",std::accumulate(input.support.begin(),input.support.end(),std::size_t{},
            [](std::size_t sum,const auto & support){return sum+support.size();})}};
    WriteText(output_dir/"checkpoint.json",json::serialize(checkpoint)+"\n");
    std::cout<<json::serialize(json::object{{"summary",(output_dir/"summary.json").string()},
        {"checkpoint",(output_dir/"checkpoint.json").string()},{"observations",rows},
        {"targets",input.selection_domain->target_indices.size()},{"halos",ownership.halos.size()},
        {"blocks",ownership.blocks.size()},{"cross_component_halos",ownership.cross_component_halos}})<<'\n';
    return 0;
}
}

int main(int argc,char ** argv)
{
    try
    {
        if(argc!=6) throw std::invalid_argument(
            "Usage: joint_dynamic_halo_mean MODEL MAP QUALIFIED_PRODUCTION_RESULT OUTPUT_DIR BASELINE_COMMIT");
        return Run(argv[1],argv[2],argv[3],argv[4],argv[5]);
    }
    catch(const std::exception & error)
    {std::cerr<<error.what()<<'\n'; return 1;}
}
