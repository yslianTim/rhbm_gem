#include "core/detail/HaloMedianReconstruction.hpp"
#include "core/detail/PostFitPeeling.hpp"
#include "core/detail/MapInterpolation.hpp"
#include "core/detail/PotentialFittingWorkflow.hpp"
#include "core/detail/GroupPotentialFitting.hpp"
#include "core/detail/joint_component/Problem.hpp"
#include "core/command/detail/SimulationGeometry.hpp"
#include "data/detail/AtomClassifier.hpp"
#include "data/detail/JointStageAdapter.hpp"
#include <rhbm_gem/core/JointComponentEstimator.hpp>
#include <rhbm_gem/data/io/ModelMapFileIO.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <rhbm_gem/data/object/ModelAnalysisView.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/math/GaussianModel3D.hpp>
#include <boost/json.hpp>
#include <Eigen/Core>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <set>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace {
namespace rg=rhbm_gem;
namespace core=rg::core;
namespace detail=core::detail;
namespace joint=core::joint_component;
namespace json=boost::json;

struct Truth
{
    double a{};
    double b{};
    double c{};
};

struct Summary
{
    std::size_t count{};
    double bias{},rmse{},median_absolute{},p95_absolute{},maximum_absolute{};
};

json::value Number(double value)
{
    return std::isfinite(value) ? json::value(value) : json::value(nullptr);
}

json::value ReadJson(const std::filesystem::path & path)
{
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Cannot read " + path.string());
    return json::parse(std::string(std::istreambuf_iterator<char>(input),{}));
}

std::unordered_map<int,Truth> ReadTruth(const std::filesystem::path & path)
{
    const auto metadata=ReadJson(path);
    std::unordered_map<int,Truth> truth;
    for (const auto & atom:metadata.at("atoms").as_array())
    {
        const int id=static_cast<int>(json::value_to<std::uint64_t>(atom.at("serial_id")));
        const double gaussian_width=json::value_to<double>(atom.at("effective_gaussian_width"));
        const double charge_width=json::value_to<double>(atom.at("effective_charge_width"));
        if (gaussian_width!=charge_width)
            throw std::runtime_error("Simulation truth has unequal Gaussian/charge widths for atom " + std::to_string(id));
        truth.emplace(id,Truth{
            static_cast<double>(json::value_to<int>(atom.at("element"))),
            gaussian_width,
            json::value_to<double>(atom.at("charge_used"))});
    }
    return truth;
}

double Percentile(std::vector<double> values,double probability)
{
    if (values.empty()) return std::numeric_limits<double>::quiet_NaN();
    std::sort(values.begin(),values.end());
    const auto index=static_cast<std::size_t>(std::ceil(probability * static_cast<double>(values.size()-1)));
    return values.at(index);
}

Summary Summarize(const std::vector<double> & errors)
{
    Summary result; result.count=errors.size();
    if (errors.empty())
    {
        result.bias=result.rmse=result.median_absolute=result.p95_absolute=result.maximum_absolute=
            std::numeric_limits<double>::quiet_NaN();
        return result;
    }
    long double sum{},square{}; std::vector<double> absolute; absolute.reserve(errors.size());
    for (const double error:errors)
    {
        sum+=error; square+=static_cast<long double>(error)*error; absolute.push_back(std::abs(error));
    }
    result.bias=static_cast<double>(sum/static_cast<long double>(errors.size()));
    result.rmse=std::sqrt(static_cast<double>(square/static_cast<long double>(errors.size())));
    result.median_absolute=Percentile(absolute,.5);
    result.p95_absolute=Percentile(absolute,.95);
    result.maximum_absolute=*std::max_element(absolute.begin(),absolute.end());
    return result;
}

json::object SummaryJson(const Summary & value)
{
    return {{"count",value.count},{"bias",Number(value.bias)},{"rmse",Number(value.rmse)},
        {"median_absolute_error",Number(value.median_absolute)},
        {"p95_absolute_error",Number(value.p95_absolute)},
        {"maximum_absolute_error",Number(value.maximum_absolute)}};
}

template<std::size_t N>
json::object StatsJson(const std::array<std::vector<double>,N> & errors,
    const std::array<std::string,N> & names)
{
    json::object output;
    for (std::size_t i=0;i<errors.size();++i) output[names[i]]=SummaryJson(Summarize(errors[i]));
    return output;
}

std::string EvidenceStatusName(rg::EvidenceStatus status)
{
    switch (status)
    {
    case rg::EvidenceStatus::NotRun: return "NotRun";
    case rg::EvidenceStatus::Available: return "Available";
    case rg::EvidenceStatus::Unavailable: return "Unavailable";
    case rg::EvidenceStatus::Ineligible: return "Ineligible";
    }
    return "Unknown";
}

bool IsTarget(const core::JointProblemInput & input,std::size_t atom)
{
    return input.selection_domain && std::binary_search(
        input.selection_domain->target_indices.begin(),input.selection_domain->target_indices.end(),atom);
}

std::vector<std::size_t> TargetIndices(const core::JointProblemInput & input)
{
    if (input.selection_domain) return input.selection_domain->target_indices;
    std::vector<std::size_t> all(input.atom_ids.size());
    std::iota(all.begin(),all.end(),0);
    return all;
}

std::map<std::size_t,rg::GaussianModel3D> SnapshotPoints(const rg::JointAnalysisResult & result)
{
    std::map<std::size_t,rg::GaussianModel3D> points;
    for (const auto & component:result.components)
    {
        if (!component.state) continue;
        const auto & full=component.layout ? component.layout->full_atoms : component.atoms;
        if (component.state->ac.size()!=2*full.size() || component.state->b.size()!=full.size())
            throw std::runtime_error("Invalid Joint state while building truth comparison.");
        for (std::size_t k=0;k<full.size();++k)
            points[full[k]]=rg::GaussianModel3D{
                component.state->ac[2*k],component.state->b[k],component.state->ac[2*k+1]};
    }
    return points;
}

json::object TargetAccuracy(const core::JointProblemInput & input,
    const rg::JointAnalysisResult & result,const std::unordered_map<int,Truth> & truth,double scale)
{
    const auto points=SnapshotPoints(result);
    std::array<std::vector<double>,3> errors;
    for (const auto atom:TargetIndices(input))
    {
        const auto found=points.find(atom);
        const auto t=truth.find(std::stoi(input.atom_ids[atom]));
        if (found==points.end() || t==truth.end()) continue;
        const auto & point=found->second;
        errors[0].push_back(point.GetAmplitude()-scale*t->second.a);
        errors[1].push_back(point.GetWidth()-t->second.b);
        errors[2].push_back(point.GetOffset()-scale*t->second.c);
    }
    return {{"status",errors[0].empty() ? "NotAssessed" : "Measured"},
        {"parameter_order",json::array{"A","B_angstrom","C_charge_coefficient"}},
        {"errors",StatsJson(errors,{"A","B_angstrom","C"})}};
}

std::map<std::size_t,rg::GaussianModel3D> HaloPoints(
    const core::JointProblemInput & input,const rg::JointAnalysisResult & result)
{
    std::map<std::size_t,rg::GaussianModel3D> points;
    for (const auto & fixed:result.fixed_halos) points[fixed.contributor_index]=
        rg::GaussianModel3D{fixed.amplitude,fixed.width,fixed.charge};
    if (result.halo_median_reconstruction)
        for (const auto & halo:result.halo_median_reconstruction->halos)
            points[halo.halo_contributor_index]=rg::GaussianModel3D{halo.amplitude,halo.width,halo.charge};
    (void)input;
    return points;
}

json::object HaloFieldAccuracy(const core::JointProblemInput & input,
    const rg::JointAnalysisResult & result,const std::unordered_map<int,Truth> & truth,
    double scale,std::string_view strategy)
{
    const auto points=HaloPoints(input,result);
    std::vector<double> all_errors;
    std::map<std::string,std::vector<double>> by_source;
    json::array per_halo;
    std::size_t missing{};
    for (std::size_t atom=0;atom<input.atom_ids.size();++atom)
    {
        if (IsTarget(input,atom)) continue;
        const auto point=points.find(atom);
        const auto t=truth.find(std::stoi(input.atom_ids[atom]));
        if (point==points.end() || t==truth.end()) {++missing; continue;}
        std::vector<double> errors;
        for (const auto & support:input.support[atom])
        {
            const auto truth_basis=joint::EvaluateKernel(support.squared_distance,t->second.b,2.5);
            const auto fit_basis=joint::EvaluateKernel(support.squared_distance,point->second.GetWidth(),2.5);
            const double truth_value=scale*t->second.a*truth_basis.gaussian+scale*t->second.c*truth_basis.charge;
            const double fitted_value=point->second.GetAmplitude()*fit_basis.gaussian+
                point->second.GetOffset()*fit_basis.charge;
            errors.push_back(fitted_value-truth_value);
        }
        const auto summary=Summarize(errors);
        const auto source=(result.halo_median_reconstruction &&
            std::find_if(result.halo_median_reconstruction->halos.begin(),
                result.halo_median_reconstruction->halos.end(),[atom](const auto & halo) {
                    return halo.halo_contributor_index==atom;
                }) != result.halo_median_reconstruction->halos.end()) ?
            [&] {
                const auto & halos=*result.halo_median_reconstruction;
                const auto found=std::find_if(halos.halos.begin(),halos.halos.end(),[atom](const auto & halo) {
                    return halo.halo_contributor_index==atom;
                });
                if (found->source==rg::HaloMedianSource::GroupKeyMedian) return std::string("same-group");
                return std::string("global-fallback");
            }() : std::string(strategy=="A" ? "fixed-first-stage" : "oracle");
        all_errors.insert(all_errors.end(),errors.begin(),errors.end());
        by_source[source].insert(by_source[source].end(),errors.begin(),errors.end());
        per_halo.push_back(json::object{{"atom_id",input.atom_ids[atom]},{"source",source},
            {"support_count",errors.size()},{"error",SummaryJson(summary)}});
    }
    json::object sources;
    for (const auto & [name,errors]:by_source) sources[name]=SummaryJson(Summarize(errors));
    return {{"status",all_errors.empty() ? "NotAssessed" : "Measured"},
        {"strategy",strategy},{"missing_halo_reconstructions",missing},
        {"aggregate",SummaryJson(Summarize(all_errors))},{"by_source",std::move(sources)},
        {"per_halo",std::move(per_halo)}};
}

double TargetTruthAtGrid(const rg::MapObject & map,const rg::ModelObject & model,
    const core::JointProblemInput & input,std::size_t target,
    const std::array<int,3> & node,const std::unordered_map<int,Truth> & truth,double scale)
{
    const auto origin=map.GetOrigin(),spacing=map.GetGridSpacing();
    const std::array<double,3> position{
        origin[0]+spacing[0]*node[0],origin[1]+spacing[1]*node[1],origin[2]+spacing[2]*node[2]};
    const auto * atom=model.FindAtomPtr(std::stoi(input.atom_ids[target]));
    const auto t=truth.at(atom->GetSerialID());
    const auto atom_position=atom->GetPosition();
    const std::array<double,3> delta{
        position[0]-atom_position[0],position[1]-atom_position[1],position[2]-atom_position[2]};
    const double squared=delta[0]*delta[0]+delta[1]*delta[1]+delta[2]*delta[2];
    if (squared>6.25) return 0.0;
    const auto basis=joint::EvaluateKernel(squared,t.b,2.5);
    return scale*t.a*basis.gaussian+scale*t.c*basis.charge;
}

struct PeelingData
{
    std::vector<double> errors;
    std::map<std::pair<int,std::size_t>,double> by_sample;
    json::object report;
};

PeelingData PeelingAccuracy(const rg::MapObject & map,const rg::ModelObject & model,
    const core::JointProblemInput & input,const std::map<int,rg::PostFitPeelingResult> & peelings,
    const std::unordered_map<int,Truth> & truth,double scale)
{
    PeelingData output;
    std::size_t raw_count{},available_count{},unavailable_count{},low_coverage{};
    json::array by_target;
    for (const auto target:TargetIndices(input))
    {
        const int id=std::stoi(input.atom_ids[target]);
        const auto * atom=model.FindAtomPtr(id);
        const auto local=rg::AtomLocalPotentialView::For(*atom);
        const auto peeling=peelings.find(id);
        if (peeling==peelings.end()) continue;
        const auto & raw=local.GetRawSamplingEntriesView();
        raw_count+=raw.size();
        std::vector<double> target_errors;
        for (std::size_t i=0;i<raw.size() && i<peeling->second.samples.size();++i)
        {
            const auto & sample=peeling->second.samples[i];
            if (!sample.response) {++unavailable_count; continue;}
            const auto stencil=detail::MakeTricubicStencil(map,raw[i].point.position);
            const double reference=detail::InterpolateTricubic(stencil,[&](const auto & node) {
                return TargetTruthAtGrid(map,model,input,target,node,truth,scale);
            });
            const double error=*sample.response-reference;
            output.errors.push_back(error); target_errors.push_back(error);
            output.by_sample.emplace(std::make_pair(id,i),error); ++available_count;
        }
        if (raw.size()>0 && target_errors.size()*2<raw.size()) ++low_coverage;
        by_target.push_back(json::object{{"atom_id",id},{"raw_sample_count",raw.size()},
            {"available_sample_count",target_errors.size()},{"unavailable_sample_count",raw.size()-target_errors.size()},
            {"neighbor_count",peeling->second.neighbor_count},{"error",SummaryJson(Summarize(target_errors))}});
    }
    output.report={{"status",output.errors.empty() ? "NotAssessed" : "Measured"},
        {"raw_sample_count",raw_count},{"available_sample_count",available_count},
        {"unavailable_sample_count",unavailable_count},{"low_coverage_target_count",low_coverage},
        {"error",SummaryJson(Summarize(output.errors))},{"by_target",std::move(by_target)}};
    return output;
}

json::object LocalAccuracy(const rg::ModelObject & model,const core::JointProblemInput & input,
    const std::unordered_map<int,Truth> & truth,double scale)
{
    std::array<std::vector<double>,2> errors;
    std::map<std::string,std::size_t> statuses;
    std::size_t effective{};
    for (const auto target:TargetIndices(input))
    {
        const int id=std::stoi(input.atom_ids[target]);
        const auto local=rg::AtomLocalPotentialView::For(*model.FindAtomPtr(id));
        const auto result=local.GetPostJointLocalResult();
        if (!result) {++statuses["missing"]; continue;}
        ++statuses[EvidenceStatusName(result->status)]; effective+=result->effective_sample_count;
        if (result->status!=rg::EvidenceStatus::Available) continue;
        const auto t=truth.at(id); const auto point=result->fit.mdpde.GetModel();
        errors[0].push_back(point.GetAmplitude()-scale*t.a);
        errors[1].push_back(point.GetWidth()-t.b);
    }
    json::object status;
    for (const auto & [name,count]:statuses) status[name]=count;
    return {{"status",errors[0].empty() ? "NotAssessed" : "Measured"},
        {"available_count",errors[0].size()},{"effective_sample_count",effective},
        {"status_counts",std::move(status)},
        {"errors",StatsJson(std::array<std::vector<double>,2>{errors[0],errors[1]},
            std::array<std::string,2>{"A","B_angstrom"})}};
}

json::object GroupAccuracy(const rg::ModelObject & model,const core::JointProblemInput & input,
    const std::unordered_map<int,Truth> & truth,double scale)
{
    (void)input;
    const auto view=model.GetAnalysisView();
    std::size_t available{},unavailable{},single_member{},singular{};
    json::array groups; std::array<std::vector<double>,2> member_errors,group_errors;
    for (const auto key:view.CollectAtomGroupKeys())
    {
        const auto members=view.GetAtomObjectList(key);
        const auto summary=view.GetGroupParameterSummary(key);
        const std::string status=summary ? EvidenceStatusName(summary->status) : "missing";
        if (status=="Available") ++available; else ++unavailable;
        if (members.size()==1) ++single_member;
        if (summary && summary->reason.find("singular")!=std::string::npos) ++singular;
        json::object group{{"group_key",static_cast<std::uint64_t>(key)},
            {"member_count",members.size()},{"status",status},
            {"reason",summary ? summary->reason : "missing"}};
        json::array member_list;
        for (const auto * atom:members)
        {
            const auto member=rg::AtomLocalPotentialView::For(*atom).GetGroupMemberResult();
            json::object entry{{"atom_id",atom->GetSerialID()},{"available",member.has_value()}};
            const auto t=truth.find(atom->GetSerialID());
            if (member && t!=truth.end())
            {
                const auto point=member->posterior.GetModel();
                const double a_error=point.GetAmplitude()-scale*t->second.a;
                const double b_error=point.GetWidth()-t->second.b;
                member_errors[0].push_back(a_error); member_errors[1].push_back(b_error);
                entry["A_error"]=a_error; entry["B_error_angstrom"]=b_error;
            }
            member_list.push_back(std::move(entry));
        }
        group["members"]=std::move(member_list);
        if (summary && summary->inference)
        {
            const auto point=summary->inference->mdpde;
            group["group_A"]=point.GetAmplitude(); group["group_B_angstrom"]=point.GetWidth();
            std::vector<double> truth_a,truth_b;
            for (const auto * atom:members)
            {
                const auto t=truth.find(atom->GetSerialID());
                if (t!=truth.end()) {truth_a.push_back(scale*t->second.a); truth_b.push_back(t->second.b);}
            }
            if (!truth_a.empty())
            {
                const auto mean=[](const std::vector<double> & values) {
                    double sum{}; for (const double value:values) sum+=value;
                    return sum/static_cast<double>(values.size());
                };
                group_errors[0].push_back(point.GetAmplitude()-mean(truth_a));
                group_errors[1].push_back(point.GetWidth()-mean(truth_b));
            }
        }
        groups.push_back(std::move(group));
    }
    return {{"status",groups.empty() ? "NotAssessed" : "Measured"},
        {"group_count",groups.size()},{"available_group_count",available},
        {"unavailable_group_count",unavailable},{"single_member_group_count",single_member},
        {"singular_covariance_group_count",singular},
        {"member_parameter_error",StatsJson(member_errors,{"A","B_angstrom"})},
        {"group_parameter_error_against_member_truth_mean",StatsJson(group_errors,{"A","B_angstrom"})},
        {"groups",std::move(groups)}};
}

struct VariantRun
{
    std::string name;
    bool completed{};
    std::string failure;
    std::unique_ptr<rg::ModelObject> model;
    std::map<std::pair<int,std::size_t>,double> peeling_errors;
    std::map<int,std::array<double,2>> local_points;
    json::object result;
};

VariantRun RunVariant(std::string name,const rg::MapObject & map,const rg::ModelObject & base_model,
    const core::JointProblem & problem,const rg::JointAnalysisResult & base_snapshot,
    rg::JointHaloMedianSnapshot halo_snapshot,const std::unordered_map<int,Truth> & truth,
    double scale,const core::FitOptions & options)
{
    VariantRun output; output.name=std::move(name);
    try
    {
        output.model=std::make_unique<rg::ModelObject>(base_model);
        auto snapshot=base_snapshot;
        snapshot.halo_median_reconstruction=std::move(halo_snapshot);
        rg::data_internal::ApplyJointStageEstimates(*output.model,snapshot,"p0-"+output.name);
        output.model->EditAnalysis().SetJointResult(snapshot);
        const auto & stored=*output.model->GetAnalysisView().GetJointResult();
        detail::JointPostprocessingMetrics metrics;
        const detail::JointPostprocessingContext context(map,*output.model,problem,&stored,&metrics);
        const auto targets=TargetIndices(problem.OriginalInput());
        const auto peelings=context.BuildPostFitPeelingSamples(targets);
        const auto peeling=PeelingAccuracy(map,*output.model,problem.OriginalInput(),peelings,truth,scale);
        output.peeling_errors=peeling.by_sample;
        std::vector<std::pair<int,rg::PostJointLocalResult>> local_updates;
        local_updates.reserve(targets.size());
        for (const auto target:targets)
        {
            const int id=std::stoi(problem.OriginalInput().atom_ids[target]);
            const auto found=peelings.find(id);
            if (found==peelings.end()) continue;
            local_updates.emplace_back(id,detail::FitPostJointLocalTargetForResearch(
                context,*output.model,id,found->second,options,metrics));
        }
        std::vector<std::pair<int,rg::PostFitPeelingResult>> peeling_updates;
        peeling_updates.reserve(peelings.size());
        for (const auto & [id,value]:peelings) peeling_updates.emplace_back(id,value);
        output.model->EditAnalysis().ApplyPostFitPeelingsBatch(std::move(peeling_updates));
        output.model->EditAnalysis().ApplyPostJointLocalResultsBatch(std::move(local_updates));
        detail::RunJointGroupPotentialFitting(*output.model,options,&metrics);
        for (const auto target:targets)
        {
            const int id=std::stoi(problem.OriginalInput().atom_ids[target]);
            const auto local=rg::AtomLocalPotentialView::For(*output.model->FindAtomPtr(id)).GetPostJointLocalResult();
            if (local && local->status==rg::EvidenceStatus::Available)
            {
                const auto point=local->fit.mdpde.GetModel();
                output.local_points[id]={point.GetAmplitude(),point.GetWidth()};
            }
        }
        output.result={
            {"status","Measured"},{"strategy",output.name},
            {"joint_target_accuracy",TargetAccuracy(problem.OriginalInput(),stored,truth,scale)},
            {"halo_field_accuracy",HaloFieldAccuracy(problem.OriginalInput(),stored,truth,scale,output.name)},
            {"peeling_accuracy",peeling.report},
            {"local_mdPDE_accuracy",LocalAccuracy(*output.model,problem.OriginalInput(),truth,scale)},
            {"group_fitting_accuracy",GroupAccuracy(*output.model,problem.OriginalInput(),truth,scale)},
            {"diagnostics",json::object{
                {"raw_peeling_sample_count",metrics.raw_peeling_sample_count},
                {"nonzero_stencil_node_count",metrics.nonzero_stencil_node_count},
                {"available_peeling_sample_count",metrics.available_peeling_sample_count},
                {"unavailable_peeling_sample_count",metrics.unavailable_peeling_sample_count},
                {"map_geometry_boundary_samples",metrics.map_geometry_boundary_samples},
                {"original_support_boundary_samples",metrics.original_support_boundary_samples},
                {"fitting_mask_only_intersection_samples",metrics.fitting_mask_only_intersection_samples},
                {"mapping_mismatch_samples",metrics.mapping_mismatch_samples},
                {"missing_contributor_samples",metrics.missing_contributor_samples},
                {"other_peeling_failure_samples",metrics.other_peeling_failure_samples}}}};
        output.completed=true;
    }
    catch (const std::exception & error)
    {
        output.failure=error.what();
        output.result={{"status","NotAssessed"},{"strategy",output.name},{"reason",output.failure}};
    }
    return output;
}

rg::JointHaloMedianSnapshot GlobalOnly(const rg::JointHaloMedianSnapshot & source)
{
    auto result=source;
    for (auto & halo:result.halos)
    {
        halo.amplitude=result.global_median[0]; halo.width=result.global_median[1]; halo.charge=result.global_median[2];
        halo.source=rg::HaloMedianSource::GlobalTargetMedian;
        halo.donor_count=result.valid_target_count;
    }
    return result;
}

rg::JointHaloMedianSnapshot Oracle(const rg::JointHaloMedianSnapshot & source,
    const core::JointProblemInput & input,const std::unordered_map<int,Truth> & truth,double scale)
{
    auto result=source;
    for (auto & halo:result.halos)
    {
        const auto t=truth.find(std::stoi(input.atom_ids[halo.halo_contributor_index]));
        if (t==truth.end()) throw std::runtime_error("Missing halo truth for " + halo.atom_id);
        halo.amplitude=scale*t->second.a; halo.width=t->second.b; halo.charge=scale*t->second.c;
        // The wire contract has only production median source values. Oracle is
        // research-only and retains a valid GroupKey source tag for the parser.
        halo.source=rg::HaloMedianSource::GroupKeyMedian;
        halo.donor_count=1;
    }
    return result;
}

json::object FixedBaseline(const rg::MapObject & map,const rg::ModelObject & source,
    const std::unordered_map<int,Truth> & truth,double scale)
{
    auto fixed_map=std::make_unique<rg::MapObject>(map);
    rg::ModelObject fixed_model(source);
    core::FitOptions options; options.estimator=core::PotentialEstimator::JOINT_COMPONENTS;
    options.quiet_mode=true; options.thread_size=1;
    detail::JointPostprocessingMetrics metrics;
    detail::RunJointComponentWorkflow(*fixed_map,fixed_model,options,&metrics,false);
    const auto result=fixed_model.GetAnalysisView().GetJointResult();
    if (!result) throw std::runtime_error("Fixed Halo baseline did not produce a Joint result.");
    const auto problem=core::BuildJointProblem(*fixed_map,fixed_model);
    std::map<int,rg::PostFitPeelingResult> peelings;
    for (const auto target:problem.Input().selection_domain->target_indices)
    {
        const int id=std::stoi(problem.Input().atom_ids[target]);
        const auto peeling=rg::AtomLocalPotentialView::For(*fixed_model.FindAtomPtr(id)).GetPostFitPeeling();
        if (peeling) peelings.emplace(id,*peeling);
    }
    const auto peeling=PeelingAccuracy(*fixed_map,fixed_model,problem.Input(),peelings,truth,scale);
    return {{"status","Measured"},{"strategy","A"},
        {"joint_target_accuracy",TargetAccuracy(problem.Input(),*result,truth,scale)},
        {"halo_field_accuracy",HaloFieldAccuracy(problem.Input(),*result,truth,scale,"A")},
        {"peeling_accuracy",peeling.report},
        {"local_mdPDE_accuracy",LocalAccuracy(fixed_model,problem.Input(),truth,scale)},
        {"group_fitting_accuracy",GroupAccuracy(fixed_model,problem.Input(),truth,scale)},
        {"diagnostics",json::object{
            {"raw_peeling_sample_count",metrics.raw_peeling_sample_count},
            {"nonzero_stencil_node_count",metrics.nonzero_stencil_node_count},
            {"available_peeling_sample_count",metrics.available_peeling_sample_count},
            {"unavailable_peeling_sample_count",metrics.unavailable_peeling_sample_count},
            {"map_geometry_boundary_samples",metrics.map_geometry_boundary_samples},
            {"original_support_boundary_samples",metrics.original_support_boundary_samples},
            {"fitting_mask_only_intersection_samples",metrics.fitting_mask_only_intersection_samples},
            {"mapping_mismatch_samples",metrics.mapping_mismatch_samples},
            {"missing_contributor_samples",metrics.missing_contributor_samples},
            {"other_peeling_failure_samples",metrics.other_peeling_failure_samples}}}};
}

json::object Threads()
{
    const auto env=[](const char * name) {
        const auto * value=std::getenv(name); return value ? std::string(value) : std::string("unset");
    };
    return {{"Eigen",Eigen::nbThreads()},{"OpenMP",env("OMP_NUM_THREADS")},
        {"OpenBLAS",env("OPENBLAS_NUM_THREADS")},{"Accelerate",env("VECLIB_MAXIMUM_THREADS")},
        {"workflow_thread_size",1}};
}

void Write(const std::filesystem::path & path,const json::value & value)
{
    std::ofstream output(path); output.exceptions(std::ios::failbit|std::ios::badbit);
    output<<json::serialize(value)<<'\n';
}
}

int main(int argc,char ** argv)
{
    try
    {
        if (argc<4 || argc>5) throw std::invalid_argument(
            "Usage: joint_halo_support_exclusion_validation MODEL MAP OUTPUT.json [--cli-normalization]");
        Eigen::setNbThreads(1);
        const std::filesystem::path model_path=argv[1],map_path=argv[2],output_path=argv[3];
        bool normalize=false;
        if (argc==5)
        {
            if (std::string(argv[4])!="--cli-normalization") throw std::invalid_argument("Unknown research option.");
            normalize=true;
        }
        auto model=rg::ReadModel(model_path.string());
        auto map=rg::ReadMap(map_path.string());
        double normalization_divisor=1.0;
        if (normalize)
            if (const auto divisor=map->MapValueArrayNormalization()) normalization_divisor=*divisor;
        const double truth_scale=1.0/normalization_divisor;
        model->SelectAtoms([](const rg::AtomObject & atom) {
            return atom.GetElement()!=Element::HYDROGEN && atom.GetChainID()=="A";
        });
        model->EditAnalysis().InitializeFromSelection();
        const auto truth=ReadTruth(map_path.string()+".simulation.json");
        const auto problem=core::BuildPureHaloSupportExclusionProblem(*map,*model);
        const auto initialization=joint::InitializeContributors(*map,*model,problem);
        if (!initialization.valid) throw std::runtime_error("Pure Exclusion initialization failed: "+initialization.reason);
        const auto fit=core::FitJointComponentsInitialized(problem,initialization);
        auto base_snapshot=core::CaptureJointAnalysisResult(fit);
        const auto median=detail::BuildHaloMedianReconstruction(*model,problem,base_snapshot,"p0-joint-endpoint");
        rg::ModelObject base_model(*model);
        core::FitOptions options; options.estimator=core::PotentialEstimator::JOINT_COMPONENTS;
        options.quiet_mode=true; options.thread_size=1;

        std::array<VariantRun,3> variants{
            RunVariant("B-group-key-global-fallback",*map,base_model,problem,base_snapshot,median.snapshot,truth,truth_scale,options),
            RunVariant("C-global-median-only",*map,base_model,problem,base_snapshot,GlobalOnly(median.snapshot),truth,truth_scale,options),
            RunVariant("D-oracle-halo-reconstruction",*map,base_model,problem,base_snapshot,Oracle(median.snapshot,problem.OriginalInput(),truth,truth_scale),truth,truth_scale,options)};

        for (auto & variant:variants)
        {
            if (!variant.completed) continue;
            std::vector<double> paired;
            const auto & reference=variants[2].peeling_errors;
            for (const auto & [key,error]:variant.peeling_errors)
                if (reference.contains(key)) paired.push_back(error);
            variant.result["peeling_accuracy"].as_object()["paired_common_with_oracle"]=
                SummaryJson(Summarize(paired));
            std::array<std::vector<double>,2> local_deviation;
            for (const auto & [id,point]:variant.local_points)
            {
                const auto oracle=variants[2].local_points.find(id);
                if (oracle==variants[2].local_points.end()) continue;
                local_deviation[0].push_back(point[0]-oracle->second[0]);
                local_deviation[1].push_back(point[1]-oracle->second[1]);
            }
            variant.result["local_mdPDE_accuracy"].as_object()["deviation_from_oracle_pipeline"]=
                StatsJson(local_deviation,std::array<std::string,2>{"A","B_angstrom"});
        }

        json::array variant_json;
        for (const auto & variant:variants) variant_json.push_back(variant.result);
        const auto fixed=FixedBaseline(*map,*model,truth,truth_scale);
        const auto & original=problem.OriginalInput();
        Write(output_path,json::object{
            {"contract","pure-halo-support-exclusion-production-validation-v1"},
            {"status","Measured"},{"model",model_path.string()},{"map",map_path.string()},
            {"simulation_truth",map_path.string()+".simulation.json"},
            {"map_normalization",json::object{{"requested",normalize},{"divisor",normalization_divisor},
                {"truth_amplitude_and_charge_scale",truth_scale}}},
            {"kernel","single-gaus-v1 plus charge kernel; C is charge coefficient"},
            {"support","sphere-fma-v1; squared distance <= 6.25 A^2"},
            {"joint_policy","Pure Halo-Support Exclusion; R=2.5 A; one shared Joint endpoint for B/C/D"},
            {"fixed_seed_configuration","deterministic sampling; no random noise; Eigen=1; workflow thread_size=1"},
            {"threads",Threads()},{"observations",original.observations.size()},
            {"contributors",original.atom_ids.size()},
            {"targets",original.selection_domain ? original.selection_domain->target_indices.size() : 0},
            {"retained_halos",original.selection_domain && original.selection_domain->pure_halo_support_exclusion ?
                original.selection_domain->pure_halo_support_exclusion->retained_halos : 0},
            {"pure_halo_median_audit",json::object{
                {"valid_target_count",median.audit.valid_target_count},
                {"group_count",median.audit.group_count},
                {"same_group_halos",median.audit.same_group_halos},
                {"fallback_halos",median.audit.fallback_halos},
                {"fallback_support_rows",median.audit.fallback_support_rows},
                {"global_median",json::array{median.audit.global_median[0],median.audit.global_median[1],median.audit.global_median[2]}}}},
            {"joint_endpoint",json::object{{"runtime_convergence",static_cast<int>(fit.RuntimeConvergence())},
                {"target_runtime_convergence",static_cast<int>(fit.TargetRuntimeConvergence())},
                {"independent_replay",std::any_of(fit.evidence.begin(),fit.evidence.end(),[](const auto & check) {
                    return check.name=="independent-replay" && check.status==rg::JointCheckStatus::Passed;
                })},{"target_accuracy",TargetAccuracy(original,base_snapshot,truth,truth_scale)}}},
            {"A_fixed_halo_baseline",fixed},{"B_C_D_controlled_variants",std::move(variant_json)},
            {"interpretation","B/C/D share the same Pure Joint target point estimates, fitting rows, map and raw sample geometry; only halo peeling reconstruction differs. Local/group truth summaries are A/B only because the Group estimator has no C parameter."}});
        std::cout<<json::serialize(json::object{{"output",output_path.string()},
            {"observations",original.observations.size()},{"targets",original.selection_domain->target_indices.size()},
            {"halo_variants",3}})<<'\n';
        return 0;
    }
    catch (const std::exception & error)
    {
        std::cerr<<error.what()<<'\n';
        return 1;
    }
}
