#include "core/estimator/InitialStage.hpp"
#include "core/estimator/GroupingStage.hpp"
#include "core/estimator/PotentialFittingWorkflow.hpp"
#include "core/classifier/LocalFittingFeatures.hpp"
#include "core/classifier/UmapEmbeddingAlgorithm.hpp"
#include "core/classifier/UmapEmbeddingFeatures.hpp"
#include "core/estimator/gaussian_fit/PreparedLocalGaussianFit.hpp"
#include "core/estimator/joint_component/HaloMedianReconstruction.hpp"
#include "core/estimator/joint_component/JointCliProgressReporter.hpp"
#include "core/estimator/joint_component/Numerics.hpp"
#include "core/estimator/joint_component/PostJointQScore.hpp"
#include "core/estimator/joint_component/Problem.hpp"
#include "data/detail/JointStageAdapter.hpp"
#include "data/io/detail/JointResultJson.hpp"
#include <rhbm_gem/core/QScoreHelper.hpp>
#include <rhbm_gem/data/io/ModelMapFileIO.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <rhbm_gem/data/object/ModelAnalysisView.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>

#include <boost/json.hpp>
#if defined(__APPLE__) || defined(__linux__)
#include <sys/resource.h>
#endif

#include <algorithm>
#include <chrono>
#include <cstdint>
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

#ifdef RHBM_GEM_ENABLE_UMAP
void WriteUmapMatrix(const fs::path & path,
    const std::vector<core::classifier::LocalFittingFeatureRow> & rows,
    const core::classifier::UmapFeatureMatrix & matrix)
{
    std::ofstream file(path);
    file.exceptions(std::ios::failbit | std::ios::badbit);
    file << std::setprecision(17);
    if (matrix.values.empty()) return;
    for (std::size_t i=0;i<rows.size();++i)
    {
        file << rows[i].serial_id;
        for (std::size_t column=0;column<matrix.column_count;++column)
            file << '\t' << matrix.values[i*matrix.column_count+column];
        file << '\n';
    }
}
#endif

void SelectOutput(ModelObject & model, const std::string & argument)
{
    const auto chains=Split(argument);
    model.SelectAtoms([&](const AtomObject & atom) {
        return atom.GetElement()!=Element::HYDROGEN &&
            (chains.empty() || std::find(chains.begin(),chains.end(),atom.GetChainID())!=chains.end());
    });
}

void ExportGeometry(const MapObject & map, const ModelObject & model,
    const core::JointProblem & problem, const fs::path & out)
{
    fs::create_directories(out);
    const auto & original=problem.OriginalInput();
    j::array atoms;
    std::ofstream incidences(out/"support.bin",std::ios::binary), rows(out/"rows.tsv");
    incidences.exceptions(std::ios::failbit | std::ios::badbit);
    rows.exceptions(std::ios::failbit | std::ios::badbit);
    if (original.row_ids.size()>std::numeric_limits<std::uint32_t>::max())
        throw std::length_error("Research incidence row index exceeds uint32.");
    for (std::size_t i=0; i<original.atom_ids.size(); ++i)
    {
        const auto id=std::stoi(original.atom_ids[i]);
        const auto & atom=*model.FindAtomPtr(id);
        const auto & support=original.support[i];
        atoms.push_back(j::object{{"id",id},{"chain",atom.GetChainID()},{"residue",atom.GetComponentID()},
            {"sequence",atom.GetSequenceID()},{"atom",atom.GetAtomID()},{"element",static_cast<int>(atom.GetElement())},
            {"position",j::value_from(atom.GetPositionRef())},{"support",support.size()},
            {"output",std::binary_search(original.selection_domain->target_indices.begin(),
                original.selection_domain->target_indices.end(),i)}});
        for (const auto & s:support)
        {
            const auto row=static_cast<std::uint32_t>(s.row);
            incidences.write(reinterpret_cast<const char *>(&row),sizeof(row));
            incidences.write(reinterpret_cast<const char *>(&s.squared_distance),sizeof(s.squared_distance));
        }
    }
    rows << "voxel\tfitting\n";
    auto fitting=problem.FittingToOriginalRows().begin();
    for (std::size_t i=0; i<original.row_ids.size(); ++i)
    {
        const bool retained=fitting!=problem.FittingToOriginalRows().end() && *fitting==i;
        rows << original.row_ids[i] << '\t' << retained << '\n';
        if (retained) ++fitting;
    }
    const auto spacing=n::AssessWidthGridSpacing(map);
    const auto & memory=core::JointProblemAccess::Get(problem).construction_memory;
    Write(out/"geometry.json",j::object{{"grid",j::value_from(map.GetGridSize())},
        {"spacing",j::value_from(map.GetGridSpacing())},{"origin",j::value_from(map.GetOrigin())},
        {"isotropic_spacing",spacing.spacing ? Number(*spacing.spacing) : j::value{}},
        {"spacing_reason",std::string(spacing.disabled_reason)},{"atoms",std::move(atoms)},
        {"fitting_rows",problem.Input().row_ids.size()},{"original_rows",original.row_ids.size()},
        {"construction_peak_estimate_bytes",memory.estimated_peak_bytes}});
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
    if (std::none_of(initialization.atoms.begin(),initialization.atoms.end(),
        [](const auto & atom){return atom.seed_source=="median-fallback";})) return;
    JointInitialization seed;
    for (const auto i : output)
    {
        seed.atoms.push_back(initialization.atoms.at(i));
        seed.tau.push_back(seed.atoms.back().seed_source == "median-fallback" ?
            std::numeric_limits<double>::quiet_NaN() : initialization.tau.at(i));
    }
    // Use the existing fallback calculation for an additional fitting atom
    // with no valid width, without admitting that atom to the donor pool.
    seed.atoms.emplace_back().target=true;
    seed.tau.push_back(std::numeric_limits<double>::quiet_NaN());
    detail::ApplyJointSeedFallback(seed);
    if (!seed.valid) throw std::invalid_argument("Invalid output-only width initialization.");
    for (std::size_t i=0; i<output.size(); ++i)
    {
        initialization.tau[output[i]]=seed.tau[i];
        initialization.atoms[output[i]]=std::move(seed.atoms[i]);
    }
    for (std::size_t i=0;i<initialization.atoms.size();++i)
        if (auto & atom=initialization.atoms[i]; atom.seed_source=="median-fallback")
        {
            initialization.tau[i]=seed.tau.back();
            atom.used_tau=seed.tau.back();
            atom.donor_count=seed.atoms.back().donor_count;
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

void ExportPrediction(const fs::path & out, const core::JointProblem & problem,
    const JointAnalysisResult & snapshot, const ModelObject & model,
    const detail::HaloMedianBuildResult & halos, const std::set<int> & output_ids)
{
    const auto & original=problem.OriginalInput();
    const auto & state=*snapshot.assembled_state;
    std::vector<std::array<double,3>> predicted(original.row_ids.size());
    auto add=[&](std::size_t atom, const std::array<double,3> & point, std::size_t col) {
        for (const auto & s:original.support[atom])
        {
            const auto b=n::EvaluateKernel(s.squared_distance,point[1],2.5);
            predicted[s.row][col]+=point[0]*b.gaussian+point[2]*b.charge;
        }
    };
    for (std::size_t i=0;i<state.tau.size();++i)
        add(snapshot.fitting_to_original_contributors[i],{state.ac[2*i],state.tau[i],state.ac[2*i+1]},0);
    for (auto & row:predicted) row[1]=row[2]=row[0];
    std::vector<detail::HaloMedianTargetDonor> donors;
    for (const auto & donor:halos.target_donors)
        if (output_ids.contains(std::stoi(donor.atom_id))) donors.push_back(donor);
    const detail::GeometryLocalTargetDonorIndex output_only(donors);
    const auto & medians=*halos.snapshot;
    j::array halo_diagnostics;
    for (const auto & halo:medians.halos)
    {
        const auto * atom=model.FindAtomPtr(std::stoi(halo.atom_id));
        const auto group=output_only.Query(atom->GetPositionRef());
        const auto alternative=detail::GeometryLocalTargetMedian(group);
        const std::array point{halo.amplitude,halo.width,halo.charge};
        add(halo.halo_contributor_index,point,1);
        add(halo.halo_contributor_index,alternative,2);
        std::vector<std::string> alternative_ids;
        for (const auto * donor:group) alternative_ids.push_back(donor->atom_id);
        halo_diagnostics.push_back(j::object{{"id",halo.atom_id},{"fitting_donors",j::value_from(halo.donor_atom_ids)},
            {"output_donors",j::value_from(alternative_ids)},{"fitting_point",j::value_from(point)},
            {"output_point",j::value_from(alternative)}});
    }
    Write(out/"donors.json",halo_diagnostics);
    std::ofstream file(out/"prediction.tsv");
    file.exceptions(std::ios::failbit | std::ios::badbit);
    file << std::setprecision(17) << "voxel\tfitting\tfit_only\tfull\tfull_output_donors\n";
    for (std::size_t i=0;i<original.row_ids.size();++i)
        file << original.row_ids[i] << '\t' << snapshot.available_row_mask[i] << '\t'
            << predicted[i][0] << '\t' << predicted[i][1] << '\t' << predicted[i][2] << '\n';
}

#ifdef RHBM_GEM_ENABLE_UMAP
void ExportResearchEmbedding(const fs::path & out, const core::classifier::PreparedUmapInput & prepared)
{
    namespace c=core::classifier;
    try
    {
        const auto first=c::RunUmap(prepared.features,prepared.effective_neighbors,.1,0,42,1);
        const auto shifted=c::RunUmapMeanShift(first,prepared.feature_rows,1);
        auto binning=c::BinMeanShiftEmbedding(prepared.feature_rows,shifted.coordinates,shifted.scale);
        c::PopulateUmapBinFeatureMedians(prepared.feature_rows,binning);
        const auto second_features=c::BuildSecondUmapFeatures(prepared.feature_rows,binning);
        const auto second=c::RunUmap(second_features,prepared.effective_neighbors,.1,0,42,1);
        if (!c::HasFiniteCoordinates(first) || !c::HasFiniteCoordinates(second))
            throw std::runtime_error("Nonfinite UMAP coordinates.");
        std::ofstream embedding(out/"umap.tsv");
        embedding.exceptions(std::ios::failbit | std::ios::badbit);
        embedding << std::setprecision(17);
        for (std::size_t i=0;i<prepared.feature_rows.size();++i)
            embedding << prepared.feature_rows[i].serial_id << '\t' << first[2*i] << '\t' << first[2*i+1]
                << '\t' << second[2*i] << '\t' << second[2*i+1] << '\n';
        Write(out/"umap.json",j::object{{"status","passed"},{"rows",prepared.feature_rows.size()},
            {"excluded",prepared.excluded_rows},{"exclusions",j::value_from(prepared.exclusions)}});
    }
    catch (const std::exception & e)
    {
        Write(out/"umap.json",j::object{{"status","failed"},{"reason",e.what()}});
    }
}
#endif

void ExportOutputAnalysis(const fs::path & out, ModelObject & model, const core::FitOptions & options)
{
    namespace c=core::classifier;
    detail::RunGroupingStage(model,options);
    const auto view=model.GetAnalysisView();
    j::array groups;
    for (const auto key:view.CollectAtomGroupKeys())
    {
        const auto s=view.GetGroupParameterSummary(key);
        if (s) groups.push_back(j::object{{"key",key},{"status",static_cast<int>(s->status)},
            {"reason",s->reason},{"eligible",s->eligible_count},{"excluded",s->excluded_count},
            {"mdpde",s->inference ? j::value(Point(s->inference->mdpde)) : j::value{}}});
    }
    Write(out/"groups.json",groups);
#ifdef RHBM_GEM_ENABLE_UMAP
    const auto rows=c::BuildLocalFittingFeatureRows(model);
    std::ofstream features(out/"features.tsv");
    features.exceptions(std::ios::failbit | std::ios::badbit);
    features << std::setprecision(17);
    for (const auto & r:rows)
    {
        features << r.serial_id;
        for (const double v:r.features) features << '\t' << v;
        features << '\n';
    }
    std::string reason;
    auto prepared=c::BuildAndStandardizeInput(model,rows,"buffer-research",15,reason);
    if (!prepared)
    {
        Write(out/"umap.json",j::object{{"status","unavailable"},{"reason",reason}});
        return;
    }
    WriteUmapMatrix(out/"standardized.tsv",prepared->feature_rows,prepared->features);
    ExportResearchEmbedding(out,*prepared);
#endif
}

#ifdef RHBM_GEM_ENABLE_UMAP
void ReplayUmapFeatures(const ModelObject & model, const fs::path & input, const fs::path & out)
{
    namespace c = core::classifier;
    std::ifstream stream(input);
    if (!stream) throw std::runtime_error("Unable to read recorded UMAP inputs.");
    std::vector<c::LocalFittingFeatureRow> rows;
    std::vector<double> coordinates;
    c::UmapBinning bins;
    std::size_t recorded_bins{};
    for (std::string line; std::getline(stream,line);)
    {
        std::istringstream fields(line);
        c::LocalFittingFeatureRow row;
        if (!(fields >> row.serial_id)) throw std::runtime_error("Missing recorded UMAP atom ID.");
        for (double & feature:row.features)
        {
            std::string token;
            if (!(fields >> token)) throw std::runtime_error("Incomplete recorded UMAP features.");
            feature=std::stod(token);
        }
        double x{},y{};
        std::int64_t bin{};
        if (!(fields >> x >> y >> bin) || !model.FindAtomPtr(row.serial_id))
            throw std::runtime_error("Invalid recorded UMAP row.");
        if (bin>0)
        {
            bins.bins[{bin,0}].points.push_back(rows.size());
            ++recorded_bins;
        }
        rows.push_back(row);
        coordinates.insert(coordinates.end(),{x,y});
    }
    const auto row_count=rows.size();
    if (recorded_bins!=0 && recorded_bins!=row_count)
        throw std::runtime_error("Recorded UMAP bins must cover all rows.");
    std::string reason;
    const auto prepared=c::BuildAndStandardizeInput(model,std::move(rows),"buffer-audit",15,reason);
    if (!prepared || prepared->feature_rows.size()!=row_count)
        throw std::runtime_error("Recorded UMAP comparison rows were excluded during feature replay. " + reason);
    if (bins.bins.empty())
    {
        const auto shifted=c::RunUmapMeanShift(coordinates,prepared->feature_rows,1);
        bins=c::BinMeanShiftEmbedding(prepared->feature_rows,shifted.coordinates,shifted.scale);
    }
    c::PopulateUmapBinFeatureMedians(prepared->feature_rows,bins);
    const auto second=c::BuildSecondUmapFeatures(prepared->feature_rows,bins);
    WriteUmapMatrix(out/"audited-first-features.tsv",prepared->feature_rows,prepared->features);
    WriteUmapMatrix(out/"audited-second-features.tsv",prepared->feature_rows,second);
    Write(out/"audited-features.json",j::object{{"first",j::value_from(c::GetUmapFeatureNames())},
        {"second",j::value_from(c::GetSecondUmapFeatureNames())},{"bins",bins.bins.size()}});
    ExportResearchEmbedding(out,*prepared);
}
#endif

j::object RunPostprocessing(const fs::path & out, const MapObject & map, ModelObject & model,
    const core::JointProblem & problem, JointAnalysisResult & snapshot, const detail::HaloMedianBuildResult & halo,
    const std::set<int> & output_ids, const std::vector<std::size_t> & focal_indices,
    const core::FitOptions & options, bool full, bool export_samples)
{
    const auto start=Clock::now();
    const auto elapsed=[](auto since) {return std::chrono::duration<double>(Clock::now()-since).count();};
    if (problem.OriginalInput().atom_ids.size()>problem.Input().atom_ids.size())
        snapshot.halo_median_reconstruction=*halo.snapshot;
    model.EditAnalysis().ApplyJointResult(snapshot,halo.snapshot->joint_run_id);
    model.EditAnalysis().BindJointMapObservations(map);
    detail::JointPostprocessingMetrics metrics;
    detail::JointPostprocessingContext post(map,model,problem,snapshot,&metrics,&halo,
        detail::BuildJointPeelingContributorUniverse(model,problem,snapshot));
    auto qscore=detail::BuildPostJointQScores(map,model,problem,snapshot,&halo,metrics);
    const auto peeling_start=Clock::now();
    auto peelings=post.BuildPostFitPeelingBatch(focal_indices);
    const double peeling_seconds=elapsed(peeling_start);
    const auto rss_peeling=n::CurrentRssBytes();
    std::vector<std::pair<int,PostJointLocalResult>> local_results;
    std::erase_if(qscore,[&](const auto & p){return !output_ids.contains(p.first);});
    const std::map<int,PostJointQScoreResult> scores(qscore.begin(),qscore.end());
    std::ofstream samples;
    if (export_samples)
    {
        samples.exceptions(std::ios::failbit | std::ios::badbit);
        samples.open(out/"peeling-responses.tsv");
        samples << std::setprecision(17) << "atom\tvoxel\tresponse\treason\n";
    }
    j::object points;
    for (const auto & [id,peeling]:peelings)
    {
        if (export_samples)
            for (const auto & sample:peeling.samples)
            {
                samples << id << '\t' << sample.observation.global_grid_index << '\t';
                if (sample.response) samples << *sample.response;
                samples << '\t' << sample.reason << '\n';
            }
        auto local=detail::FitPostJointLocalTarget(post,model,id,peeling,options,metrics);
        j::object diagnostics{{"effective_samples",local.effective_sample_count},
            {"mdpde",Point(local.fit.mdpde.GetModel())},{"available",local.information.has_value()}};
        diagnostics["reason"]=local.reason;
        diagnostics["peeling_samples"]=peeling.samples.size();
        diagnostics["peeling_available"]=std::count_if(peeling.samples.begin(),peeling.samples.end(),
            [](const auto & sample) {return sample.response.has_value();});
        const auto q=scores.find(id);
        if (q!=scores.end())
        {
            if (q->second.scores) diagnostics["qscore"]=j::array{q->second.scores->standard,q->second.scores->joint};
            else diagnostics["qscore_reason"]=q->second.reason;
        }
        points[std::to_string(id)]=std::move(diagnostics);
        if (full) local_results.emplace_back(id,std::move(local));
    }
    const auto output_start=Clock::now();
    if (full)
    {
        auto editor=model.EditAnalysis();
        editor.ApplyPostFitPeelingsBatch(std::move(peelings));
        editor.ApplyPostJointLocalResultsBatch(std::move(local_results));
        ExportOutputAnalysis(out,model,options);
    }
    return {{"points",std::move(points)}, {"donors",halo.target_donors.size()},
        {"runtime_seconds",elapsed(start)}, {"grouping_umap_seconds",elapsed(output_start)},
        {"qscore_seconds",metrics.qscore_seconds},{"peeling_seconds",peeling_seconds},
        {"qscore_estimated_bytes",metrics.estimated_qscore_peak_bytes},
        {"peeling_estimated_bytes",metrics.estimated_peeling_workspace_bytes},
        {"peeling_result_estimated_bytes",metrics.estimated_peeling_result_bytes},
        {"rss_after_peeling_bytes",j::value_from(rss_peeling)}};
}

void Run(const char * model_path, const char * map_path, const std::string & chain_argument,
    const std::string & buffer_argument, const std::string & domain_argument,
    const std::string & reference_argument, int focal, const fs::path & out, const std::string & stage)
{
    const auto start=Clock::now();
    const auto elapsed=[](auto since) {return std::chrono::duration<double>(Clock::now()-since).count();};
    std::vector<int> buffer;
    for (const auto & value:Split(buffer_argument)) buffer.push_back(std::stoi(value));
    if (domain_argument!="fixed" && domain_argument!="expanded") throw std::invalid_argument("Unknown domain mode.");
    const auto domain=domain_argument=="fixed" ? detail::BufferObservationDomain::Fixed : detail::BufferObservationDomain::Expanded;
    std::optional<double> reference_norm;
    const bool output_calibration=reference_argument!="-";
    if (output_calibration && reference_argument!="output") reference_norm=std::stod(reference_argument);
    const bool simulation=stage=="full-simulation";
    const bool donor_audit=stage=="full-donor-audit";
    const bool full=stage=="full" || simulation || donor_audit;
    if (!full && stage!="focal") throw std::invalid_argument("Unknown research stage.");

    j::object resources;
    const auto sample_rss=[&](std::string_view key){resources[key]=j::value_from(n::CurrentRssBytes());};
    sample_rss("rss_before_inputs");
    auto model=ReadModel(model_path);
    sample_rss("rss_after_model");
    resources["model_atoms"]=model->GetNumberOfAtom();
    resources["model_bonds"]=model->GetNumberOfBond();
    SelectOutput(*model,chain_argument);
    std::set<int> output_ids;
    for (const auto * atom:model->GetSelectedAtoms()) output_ids.insert(atom->GetSerialID());
    auto map=ReadMap(map_path);
    sample_rss("rss_after_map");
    resources["map_payload_bytes"]=map->GetMapValueArraySize()*sizeof(double);
    if (simulation) model->ApplySimulationMetadata(.5);
    const auto divisor=simulation ? 1.0 : map->MapValueArrayNormalization().value_or(1);
    const auto [height,offset]=core::GetReferenceGaussianParameters(*map);
    model->SetReferenceHeight(height);
    model->SetReferenceOffset(offset);
    fs::create_directories(out);
    const auto construction_start=Clock::now();
    auto experiment=detail::BuildJointBufferProblem(*map,*model,buffer,domain,reference_norm);
    sample_rss("rss_after_construction");
    const auto & problem=experiment.problem;
    const auto & fitting=problem.Input();
    resources["original_snapshot_bytes"]=j::value_from(n::EstimateSnapshotBytes(problem.OriginalInput()));
    resources["fitting_snapshot_bytes"]=j::value_from(n::EstimateSnapshotBytes(fitting));
    resources["original_support_incidence"]=j::value_from(n::CountSupportEntries(problem.OriginalInput()));
    resources["fitting_support_storage_bytes"]=n::CountSupportEntries(fitting).value()*sizeof(core::JointSupport);
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
    if (output_calibration) ReuseOutputSeedFallback(init,experiment.output_fitting_indices);
    const double init_seconds=elapsed(init_start);
    j::array seeds;
    for (const auto & atom:init.atoms)
        seeds.push_back(j::object{{"id",atom.id},{"initial_mdpde",Point(GaussianModel3D{
            atom.mdpde[0],atom.mdpde[1],atom.mdpde[2]})},{"seed_source",atom.seed_source},
            {"used_tau",atom.used_tau ? Number(*atom.used_tau) : j::value{}},{"donor_count",atom.donor_count},
            {"qualified",atom.initial_stage_fit_qualified},{"reason",atom.reason}});
    Write(out/"initialization.json",seeds);
    detail::JointCliProgressReporter progress;
    const auto calibrated=output_calibration ? std::span<const std::size_t>(experiment.output_fitting_indices) :
        std::span<const std::size_t>{};
    auto fit=detail::FitJointComponentsInitializedWithProgress(problem,init,progress.Observer(),
        n::AssessWidthGridSpacing(*map),calibrated);
    fit.costs.initialization_seconds=init_seconds;
    fit.costs.construction_seconds=construction_seconds;
    sample_rss("rss_after_fit");
    JointAnalysisMetadata metadata;
    metadata.grid_size=map->GetGridSize();
    metadata.grid_spacing=map->GetGridSpacing();
    metadata.origin=map->GetOrigin();
    metadata.map_normalization=JointMapNormalization{!simulation,!simulation,divisor};
    auto snapshot=core::CaptureJointAnalysisResult(std::move(fit),std::move(metadata));
    j::array components;
    std::size_t conditioning_peak{};
    for (const auto & c:snapshot.components)
    {
        components.push_back(j::object{{"id",c.id},{"sweeps",c.search_sweeps},{"evaluations",c.profile_evaluations},
            {"atoms",c.atoms.size()},{"rows",c.rows.size()},{"conditioning_blocks",c.local_conditioning.size()},
            {"stop",c.stop_reason},{"state",c.state.has_value()},{"evidence",Checks(c.evidence)}});
        conditioning_peak=std::max(conditioning_peak,c.local_conditioning_peak_memory_estimate_bytes.value_or(0));
    }
    j::object summary{{"fitting_targets",fitting.atom_ids.size()},{"output_targets",output_ids.size()},
        {"original_rows",problem.OriginalInput().row_ids.size()},{"fitting_rows",fitting.row_ids.size()},
        {"buffer",j::value_from(buffer)},{"domain",domain_argument},
        {"calibration",output_calibration ? "fixed-output" : "fitting-population"},
        {"fitting_incidence",n::CountSupportEntries(fitting).value()},
        {"construction_estimated_peak_bytes",core::JointProblemAccess::Get(problem).construction_memory.estimated_peak_bytes},
        {"conditioning_estimated_peak_bytes",conditioning_peak},
        {"map_divisor",divisor},{"observation_scale",snapshot.observation_scale},
        {"runtime",std::string(joint_result_io::StatusText(snapshot.runtime_convergence))},
        {"target_runtime",std::string(joint_result_io::StatusText(snapshot.target_runtime_convergence))},
        {"objective",snapshot.objective ? Number(*snapshot.objective) : j::value{}},
        {"components",components},{"evidence",Checks(snapshot.evidence)},
        {"initialization_seconds",init_seconds},{"construction_seconds",construction_seconds},
        {"assessment_seconds",snapshot.costs.assessment_seconds},{"assembly_seconds",snapshot.costs.assembly_seconds},
        {"search_seconds",snapshot.costs.search_seconds},{"seconds_to_solve",elapsed(start)}};
    const auto write_summary=[&]
    {
        summary["resource_samples"]=resources;
        summary["selected_after_run"]=model->GetSelectedAtomCount();
        summary["elapsed_seconds"]=elapsed(start);
#if defined(__APPLE__) || defined(__linux__)
        rusage usage{};
        if (getrusage(RUSAGE_SELF,&usage)==0) summary["peak_rss_platform_units"]=usage.ru_maxrss;
#endif
        Write(out/"summary.json",summary);
    };
    write_summary();
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
    const auto halo_start=Clock::now();
    auto halo_build=detail::BuildHaloMedianReconstruction(*model,problem,snapshot,run_id);
    summary["halo_build_seconds"]=elapsed(halo_start);
    sample_rss("rss_after_halo");
    const auto export_start=Clock::now();
    ExportFocalGrid(out/"focal-grid.tsv",problem,snapshot,&*halo_build.snapshot,focal);
    ExportPrediction(out,problem,snapshot,*model,halo_build,output_ids);
    summary["research_prediction_export_seconds"]=elapsed(export_start);
    sample_rss("rss_after_prediction_exports");
    std::vector<std::size_t> focal_indices;
    for (const auto i:experiment.output_fitting_indices)
        if (const auto id=std::stoi(fitting.atom_ids[i]); full || id==focal || id==8607 || id==2768)
            focal_indices.push_back(snapshot.fitting_to_original_contributors[i]);
    auto post=RunPostprocessing(out,*map,*model,problem,snapshot,halo_build,
        output_ids,focal_indices,options,full,donor_audit);
    summary["post_joint"]=std::move(post.at("points"));
    post.erase("points");
    summary["postprocessing_metrics"]=std::move(post);
    summary["wall_through_fitting_donors_seconds"]=elapsed(start);
    sample_rss("rss_after_grouping_umap");
    if (donor_audit && !buffer.empty())
    {
        // This adapter masks only global donor adequacy in the local research
        // copy. Component certification, fitted points and contributor roles
        // remain the inputs for prediction and Peeling.
        for (auto & point:snapshot.target_evidence->point_atoms)
            if (!output_ids.contains(std::stoi(snapshot.atom_ids.at(point.atom))))
                point.local_conditioning="Warning";
        const auto alternate_run="buffer-research-output-donors";
        auto alternate=detail::BuildHaloMedianReconstruction(*model,problem,snapshot,alternate_run);
        const auto alternate_out=out/"output-only-donors";
        fs::create_directories(alternate_out);
        auto audit=RunPostprocessing(alternate_out,*map,*model,problem,snapshot,alternate,
            output_ids,focal_indices,options,true,true);
        Write(alternate_out/"postprocessing.json",audit);
        sample_rss("rss_after_output_only_donors");
    }
    write_summary();
}

} // namespace

void RunBufferTruth(const std::filesystem::path &, const std::filesystem::path &);

int main(int argc, char ** argv)
{
    try
    {
        if (argc==4 && std::string_view(argv[1])=="truth") RunBufferTruth(argv[2],argv[3]);
#ifdef RHBM_GEM_ENABLE_UMAP
        else if (argc==5 && std::string_view(argv[1])=="umap-audit")
            ReplayUmapFeatures(*ReadModel(argv[2]),argv[3],argv[4]);
#endif
        else if (argc==3 && std::string_view(argv[1])=="local")
        {
            std::ifstream samples(argv[2]);
            if (!samples) throw std::invalid_argument("Unable to read corrected samples.");
            LocalPotentialSampleList input;
            for (double radius{},value{}; samples >> radius >> value;)
                input.push_back({value,{radius,{0,0,0},true}});
            std::cout << boost::json::serialize(DiagnoseLocal(input)) << '\n';
        }
        else if (argc==6 && std::string_view(argv[1])=="geometry")
        {
            const auto model=ReadModel(argv[2]);
            SelectOutput(*model,argv[4]);
            const auto map=ReadMap(argv[3]);
            const auto problem=core::BuildPureHaloSupportExclusionProblem(*map,*model);
            ExportGeometry(*map,*model,problem,argv[5]);
        }
        else if (argc==9 || argc==10) Run(argv[1],argv[2],argv[3],argv[4],argv[5],argv[6],
            std::stoi(argv[7]),argv[8],argc==10 ? argv[9] : "focal");
        else throw std::invalid_argument("Usage: joint_buffer_expansion MODEL MAP CHAINS BUFFER_IDS fixed|expanded REFERENCE_NORM|- FOCAL OUTPUT_DIR");
        return 0;
    }
    catch (const std::exception & error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
