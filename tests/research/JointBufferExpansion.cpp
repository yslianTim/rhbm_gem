#include "core/estimator/InitialStage.hpp"
#include "core/estimator/PotentialFittingWorkflow.hpp"
#include "core/estimator/gaussian_fit/PreparedLocalGaussianFit.hpp"
#include "core/estimator/joint_component/HaloMedianReconstruction.hpp"
#include "core/estimator/joint_component/JointCliProgressReporter.hpp"
#include "core/estimator/joint_component/Numerics.hpp"
#include "core/estimator/joint_component/PostJointQScore.hpp"
#include "data/detail/JointStageAdapter.hpp"
#include "data/io/detail/JointResultJson.hpp"
#include <rhbm_gem/core/QScoreHelper.hpp>
#include <rhbm_gem/data/io/ModelMapFileIO.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>

#include <boost/json.hpp>
#if defined(__APPLE__) || defined(__linux__)
#include <sys/resource.h>
#endif

#include <algorithm>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>

namespace {
using namespace rhbm_gem;
namespace fs = std::filesystem;
namespace core = rhbm_gem::core;
namespace detail = core::detail;
namespace n = core::joint_component;
namespace j = boost::json;
using Clock = std::chrono::steady_clock;

std::vector<std::string> Split(const std::string & value)
{
    std::vector<std::string> tokens;
    if (value == "-") return tokens;
    std::istringstream input(value);
    for (std::string x; std::getline(input, x, ',');) tokens.push_back(std::move(x));
    return tokens;
}

j::value Number(double value)
{
    if (!std::isfinite(value)) return nullptr;
    return value;
}

j::array Point(const GaussianModel3D & value)
{
    return {Number(value.GetAmplitude()), Number(value.GetWidth()), Number(value.GetOffset())};
}

j::array Checks(const std::vector<JointCheck> & evidence)
{
    j::array checks;
    for (const auto & e : evidence)
        checks.push_back(j::object{{"name",e.name},{"status",std::string(joint_result_io::StatusText(e.status))},
            {"value",e.value ? Number(*e.value) : j::value{}},
            {"threshold",e.threshold ? Number(*e.threshold) : j::value{}},{"reason",e.reason}});
    return checks;
}

void Write(const fs::path & path, const j::value & data)
{
    std::ofstream file(path);
    file.exceptions(std::ios::failbit | std::ios::badbit);
    file << j::serialize(data) << '\n';
}

j::object DiagnoseLocal(const LocalPotentialSampleList & samples)
{
    std::vector<double> y;
    for (const auto & x : samples) y.push_back(x.response);
    const detail::PreparedLocalGaussianDesign design(samples, 0.0, PostJointLocalFittingRadius);
    const auto data = design.BuildSelectedDataset(y, {0,1,0});
    j::object out{{"effective_samples",data.y.size()}};
    if (data.y.size() <= 2) return out;
    const auto alpha = detail::TrainLocalGaussianAlpha(data, 1);
    const auto result = detail::EstimateLocalGaussianDataset(data, alpha, 1, 0.0);
    if (result.fit_result)
    {
        const auto & fit = *result.fit_result;
        out["status"] = static_cast<int>(fit.status);
        out["alpha"] = alpha;
        out["beta_ols"] = j::array{Number(fit.beta_ols(0)), Number(fit.beta_ols(1))};
        out["beta_mdpde"] = j::array{Number(fit.beta_mdpde(0)), Number(fit.beta_mdpde(1))};
        out["ols"] = Point(result.ols.GetModel());
        out["mdpde"] = Point(result.mdpde.GetModel());
    }
    return out;
}

void ReuseOutputSeedFallback(JointInitialization & initialization, std::span<const std::size_t> output)
{
    JointInitialization seed;
    for (const auto i : output)
    {
        seed.atoms.push_back(initialization.atoms.at(i));
        seed.tau.push_back(seed.atoms.back().seed_source == "median-fallback" ?
            std::numeric_limits<double>::quiet_NaN() : initialization.tau.at(i));
    }
    detail::ApplyJointSeedFallback(seed);
    if (!seed.valid) throw std::invalid_argument("Invalid output-only width initialization.");
    for (std::size_t i=0; i<output.size(); ++i)
    {
        initialization.tau[output[i]]=seed.tau[i];
        initialization.atoms[output[i]]=std::move(seed.atoms[i]);
    }
}

void ExportFocalGrid(const fs::path & path, const core::JointProblem & problem,
    const JointAnalysisResult & snapshot, const JointHaloMedianSnapshot * halos, int focal)
{
    const auto & original = problem.OriginalInput();
    const auto id = std::find(original.atom_ids.begin(), original.atom_ids.end(), std::to_string(focal));
    if (id == original.atom_ids.end()) return;
    const auto & support = original.support.at(static_cast<std::size_t>(id-original.atom_ids.begin()));
    std::vector<std::size_t> slots(original.observations.size(), support.size());
    for (std::size_t i=0; i<support.size(); ++i) slots[support[i].row]=i;
    std::vector<std::optional<GaussianModel3D>> points(original.atom_ids.size());
    const auto & state=*snapshot.assembled_state;
    for (std::size_t i=0; i<snapshot.fitting_to_original_contributors.size(); ++i)
        points[snapshot.fitting_to_original_contributors[i]]=GaussianModel3D{state.ac[2*i],state.tau[i],state.ac[2*i+1]};
    if (halos)
        for (const auto & halo : halos->halos)
            points[halo.halo_contributor_index]=GaussianModel3D{halo.amplitude,halo.width,halo.charge};
    std::vector<double> total(support.size()), fitting(support.size());
    for (std::size_t i=0; i<points.size(); ++i)
        if (const auto & point=points[i])
            for (const auto & incidence : original.support[i])
            {
                const auto slot=slots[incidence.row];
                if (slot==support.size()) continue;
                const auto basis=n::EvaluateKernel(incidence.squared_distance,point->GetWidth(),2.5);
                const double value=point->GetAmplitude()*basis.gaussian+point->GetOffset()*basis.charge;
                total[slot]+=value;
                if (std::binary_search(snapshot.fitting_to_original_contributors.begin(),
                    snapshot.fitting_to_original_contributors.end(),i)) fitting[slot]+=value;
            }
    std::ofstream file(path);
    file.exceptions(std::ios::failbit | std::ios::badbit);
    file << std::setprecision(17) << "grid\tsquared_distance\tfitting\tobserved\tprediction\tfitted_contributors\n";
    for (std::size_t i=0; i<support.size(); ++i)
    {
        const auto row=support[i].row;
        const bool retained=std::binary_search(snapshot.fitting_to_original_rows.begin(),snapshot.fitting_to_original_rows.end(),row);
        file << original.row_ids[row] << '\t' << support[i].squared_distance << '\t' << retained << '\t'
            << original.observations[row] << '\t' << total[i] << '\t' << fitting[i] << '\n';
    }
}

void Run(const char * model_path, const char * map_path, const std::string & chain_argument,
    const std::string & buffer_argument, const std::string & domain_argument,
    const std::string & reference_argument, int focal, const fs::path & out)
{
    const auto start=Clock::now();
    const auto elapsed=[](auto since) {return std::chrono::duration<double>(Clock::now()-since).count();};
    const auto chains=Split(chain_argument);
    std::vector<int> buffer;
    for (const auto & value:Split(buffer_argument)) buffer.push_back(std::stoi(value));
    if (domain_argument!="fixed" && domain_argument!="expanded") throw std::invalid_argument("Unknown domain mode.");
    const auto domain=domain_argument=="fixed" ? detail::BufferObservationDomain::Fixed : detail::BufferObservationDomain::Expanded;
    std::optional<double> reference_norm;
    if (reference_argument!="-") reference_norm=std::stod(reference_argument);

    auto model=ReadModel(model_path);
    model->SelectAtoms([&](const AtomObject & atom) {
        return atom.GetElement()!=Element::HYDROGEN &&
            (chains.empty() || std::find(chains.begin(),chains.end(),atom.GetChainID())!=chains.end());
    });
    std::set<int> output_ids;
    for (const auto * atom:model->GetSelectedAtoms()) output_ids.insert(atom->GetSerialID());
    auto map=ReadMap(map_path);
    const auto divisor=map->MapValueArrayNormalization().value_or(1);
    const auto [height,offset]=core::GetReferenceGaussianParameters(*map);
    model->SetReferenceHeight(height);
    model->SetReferenceOffset(offset);
    fs::create_directories(out);
    const auto construction_start=Clock::now();
    auto experiment=detail::BuildJointBufferProblem(*map,*model,buffer,domain,reference_norm);
    const auto & problem=experiment.problem;
    const auto & fitting=problem.Input();
    const double construction_seconds=elapsed(construction_start);
    if (experiment.output_fitting_indices.size()!=output_ids.size() ||
        fitting.atom_ids.size()!=output_ids.size()+buffer.size())
        throw std::logic_error("Joint buffer population size mismatch.");
    for (auto i:experiment.output_fitting_indices)
        if (!output_ids.contains(std::stoi(fitting.atom_ids.at(i))))
            throw std::logic_error("Joint buffer output identity mismatch.");

    std::cerr << "Original rows " << problem.OriginalInput().row_ids.size() << ", fitting rows " << fitting.row_ids.size()
        << ", output targets " << output_ids.size() << ", buffer atoms " << buffer.size()
        << ", observation scale " << std::setprecision(17) << problem.ObservationScale() << '\n';
    core::FitOptions options;
    options.joint_initialization_thread_size=4;
    const auto init_start=Clock::now();
    auto workset=detail::BuildInitialStageWorkset(*model,problem);
    model->EditAnalysis().InitializeFromSelection();
    auto init=detail::RunInitialStage(*map,*model,workset,options);
    if (reference_norm) ReuseOutputSeedFallback(init,experiment.output_fitting_indices);
    const double init_seconds=elapsed(init_start);
    detail::JointCliProgressReporter progress;
    const auto calibrated=reference_norm ? std::span<const std::size_t>(experiment.output_fitting_indices) :
        std::span<const std::size_t>{};
    auto fit=detail::FitJointComponentsInitializedWithProgress(problem,init,progress.Observer(),
        n::AssessWidthGridSpacing(*map),calibrated);
    fit.costs.initialization_seconds=init_seconds;
    fit.costs.construction_seconds=construction_seconds;
    JointAnalysisMetadata metadata;
    metadata.grid_size=map->GetGridSize();
    metadata.grid_spacing=map->GetGridSpacing();
    metadata.origin=map->GetOrigin();
    auto snapshot=core::CaptureJointAnalysisResult(std::move(fit),std::move(metadata));
    j::array components;
    for (const auto & c:snapshot.components)
        components.push_back(j::object{{"id",c.id},{"sweeps",c.search_sweeps},{"evaluations",c.profile_evaluations},
            {"stop",c.stop_reason},{"state",c.state.has_value()},{"evidence",Checks(c.evidence)}});
    j::object summary{{"fitting_targets",fitting.atom_ids.size()},{"output_targets",output_ids.size()},
        {"original_rows",problem.OriginalInput().row_ids.size()},{"fitting_rows",fitting.row_ids.size()},
        {"buffer",j::value_from(buffer)},{"domain",domain_argument},
        {"calibration",reference_norm ? "fixed-output" : "fitting-population"},
        {"map_divisor",divisor},{"observation_scale",snapshot.observation_scale},
        {"runtime",std::string(joint_result_io::StatusText(snapshot.runtime_convergence))},
        {"target_runtime",std::string(joint_result_io::StatusText(snapshot.target_runtime_convergence))},
        {"objective",snapshot.objective ? Number(*snapshot.objective) : j::value{}},
        {"components",components},{"evidence",Checks(snapshot.evidence)},
        {"initialization_seconds",init_seconds},{"construction_seconds",construction_seconds},
        {"search_seconds",snapshot.costs.search_seconds},{"seconds_to_solve",elapsed(start)}};
    Write(out/"summary.json",summary);
    if (!IsQualifiedJointAnalysisResult(snapshot)) throw std::runtime_error("Joint experiment endpoint did not qualify.");
    if (model->GetSelectedAtomCount()!=output_ids.size()) throw std::logic_error("Model output selection changed.");

    const auto & state=*snapshot.assembled_state;
    std::ofstream parameters(out/"parameters.tsv");
    parameters.exceptions(std::ios::failbit | std::ios::badbit);
    parameters << std::setprecision(17) << "serial\toutput\tA\ttau\tC\tinit_tau\teta\teta_gradient\tlocal_conditioning\tpoint_quality\n";
    std::vector<const JointAtomPointValidation *> point_validation(problem.OriginalInput().atom_ids.size());
    for (const auto & p:snapshot.target_evidence->point_atoms) point_validation.at(p.atom)=&p;
    for (std::size_t i=0; i<fitting.atom_ids.size(); ++i)
    {
        const auto * p=point_validation.at(snapshot.fitting_to_original_contributors.at(i));
        if (!p) throw std::logic_error("Missing fitted atom point validation.");
        parameters << fitting.atom_ids[i] << '\t' << output_ids.contains(std::stoi(fitting.atom_ids[i])) << '\t'
            << state.ac[2*i] << '\t' << state.tau[i] << '\t' << state.ac[2*i+1] << '\t' << init.tau[i]
            << '\t' << state.log_tau[i] << '\t' << state.eta_gradient[i] << '\t' << p->local_conditioning << '\t' << p->point_quality << '\n';
    }

    // Keep research outputs on the original population. The private in-memory
    // model still retains every fitted contribution for the Peeling replay.
    const std::string run_id="buffer-research";
    auto halo_build=detail::BuildHaloMedianReconstruction(*model,problem,snapshot,run_id);
    ExportFocalGrid(out/"focal-grid.tsv",problem,snapshot,&*halo_build.snapshot,focal);
    model->EditAnalysis().ApplyJointStageEstimates(data_internal::BuildJointStageEstimates(snapshot,run_id));
    if (problem.OriginalInput().atom_ids.size()>fitting.atom_ids.size())
        snapshot.halo_median_reconstruction=halo_build.TakeSnapshot();
    detail::JointPostprocessingMetrics metrics;
    detail::JointPostprocessingContext post(*map,*model,problem,snapshot,&metrics,&halo_build,
        detail::BuildJointPeelingContributorUniverse(*model,problem,snapshot));
    std::vector<std::size_t> focal_indices;
    for (const auto i:experiment.output_fitting_indices)
        if (const auto id=std::stoi(fitting.atom_ids[i]); id==focal || id==8607 || id==2768)
            focal_indices.push_back(snapshot.fitting_to_original_contributors[i]);
    auto qscore=detail::BuildPostJointQScores(*map,*model,problem,snapshot,&halo_build,metrics);
    j::object points;
    for (const auto & [id,peeling]:post.BuildPostFitPeelingBatch(focal_indices))
    {
        auto local=detail::FitPostJointLocalTarget(post,*model,id,peeling,options,metrics);
        auto diagnostics=DiagnoseLocal(post.BuildPostJointLocalSamples(id,peeling).BuildGaussianSamples());
        diagnostics["reason"]=local.reason;
        diagnostics["peeling_samples"]=peeling.samples.size();
        diagnostics["peeling_available"]=std::count_if(peeling.samples.begin(),peeling.samples.end(),
            [](const auto & sample) {return sample.response.has_value();});
        const auto q=std::find_if(qscore.begin(),qscore.end(),[id](const auto & p){return p.first==id;});
        if (q!=qscore.end())
        {
            if (q->second.scores) diagnostics["qscore"]=j::array{q->second.scores->standard,q->second.scores->joint};
            else diagnostics["qscore_reason"]=q->second.reason;
        }
        points[std::to_string(id)]=std::move(diagnostics);
    }
    summary["post_joint"]=std::move(points);
    summary["selected_after_run"]=model->GetSelectedAtomCount();
    summary["elapsed_seconds"]=elapsed(start);
#if defined(__APPLE__) || defined(__linux__)
    rusage usage{};
    if (getrusage(RUSAGE_SELF,&usage)==0)
    {
        summary["peak_rss_platform_units"]=usage.ru_maxrss;
    }
#endif
    Write(out/"summary.json",summary);
}

} // namespace

int main(int argc, char ** argv)
{
    try
    {
        if (argc==3 && std::string_view(argv[1])=="local")
        {
            std::ifstream samples(argv[2]);
            if (!samples) throw std::invalid_argument("Unable to read corrected samples.");
            LocalPotentialSampleList input;
            for (double radius{},value{}; samples >> radius >> value;)
                input.push_back({value,{radius,{0,0,0},true}});
            std::cout << boost::json::serialize(DiagnoseLocal(input)) << '\n';
        }
        else if (argc==9) Run(argv[1],argv[2],argv[3],argv[4],argv[5],argv[6],std::stoi(argv[7]),argv[8]);
        else throw std::invalid_argument("Usage: joint_buffer_expansion MODEL MAP CHAINS BUFFER_IDS fixed|expanded REFERENCE_NORM|- FOCAL OUTPUT_DIR");
        return 0;
    }
    catch (const std::exception & error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
