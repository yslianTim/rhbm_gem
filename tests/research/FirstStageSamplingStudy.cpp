// Reproducible, single-process FirstStage/JCE study. No production instrumentation.
#include "core/detail/FirstStageInitialization.hpp"
#include "support/FirstStageSamplingFixture.hpp"
#include "core/detail/FittingWorkset.hpp"
#include "core/detail/PotentialFittingWorkflow.hpp"
#include "core/detail/gaussian_fit/PreparedLocalGaussianFit.hpp"
#include "core/detail/joint_component/Problem.hpp"
#include "core/detail/joint_component/WidthCalibration.hpp"
#include <rhbm_gem/core/QScoreHelper.hpp>
#include <rhbm_gem/data/io/ModelMapFileIO.hpp>
#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <rhbm_gem/data/object/ModelAnalysisView.hpp>
#include <rhbm_gem/utils/math/SphereSampler.hpp>
#include <rhbm_gem/utils/domain/SampleFilter.hpp>
#include <boost/json.hpp>
#include <Eigen/QR>
#include <Eigen/Cholesky>
#include <Eigen/Eigenvalues>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include <memory>
#include <sys/resource.h>

namespace {
using namespace rhbm_gem;
using namespace rhbm_gem::core;
namespace d = rhbm_gem::core::detail;
namespace n = rhbm_gem::core::joint_component;
namespace j = boost::json;
using Clock = std::chrono::steady_clock;
double Seconds(Clock::time_point start)
{ return std::chrono::duration<double>(Clock::now()-start).count(); }
j::value Number(double x) { return std::isfinite(x) ? j::value(x) : j::value(nullptr); }
j::array Values(const std::vector<double> & values)
{ j::array out; for (double v : values) out.push_back(Number(v)); return out; }
j::array Checks(const std::vector<JointCheck> & checks)
{
    j::array out;
    for (const auto & c : checks)
    {
        j::object check{{"name",c.name},{"status",static_cast<int>(c.status)},{"reason",c.reason}};
        if (c.value) check["value"]=Number(*c.value);
        if (c.threshold) check["threshold"]=Number(*c.threshold);
        out.push_back(check);
    }
    return out;
}
j::object DescribeFit(const JointFitResult & fit)
{
    j::array components;
    for (const auto & c : fit.components)
        components.push_back(j::object{{"id",c.id},{"sweeps",c.search_sweeps},
            {"stop",c.stop_reason},{"checks",Checks(c.evidence)}});
    j::object out{{"runtime_qualification",static_cast<int>(fit.RuntimeConvergence())},
        {"target_qualification",static_cast<int>(fit.TargetRuntimeConvergence())},
        {"checks",Checks(fit.evidence)},{"components",components}};
    if (fit.objective) out["objective"]=Number(*fit.objective);
    if (fit.assembled_state)
    {
        out["tau"]=Values(fit.assembled_state->tau);
        out["ac"]=Values(fit.assembled_state->ac);
        out["eta_gradient"]=Values(fit.assembled_state->eta_gradient);
    }
    return out;
}

void GeometryControl(MapObject & map,ModelObject & model,const d::FittingWorkset & workset,
    d::MapInterpolationMethod interpolation)
{
    auto editor=model.EditAnalysis();
    for (const auto * atom:workset.contributors)
    {
        // Match F0 shell radii/count/center/neighbor policy; change sphere geometry only.
        const auto fibonacci=sphere_sampler::GenerateFibonacciDeterministic(atom->GetPosition());
        SamplingPointList points;
        for (std::size_t i=0;i<fibonacci.size();i+=10)
            for (const auto & p:d::GenerateSpiralSpherePoints(atom->GetPositionRef(),fibonacci[i].distance,10))
                points.push_back({fibonacci[i].distance,p,true});
        std::vector<std::array<double,3>> neighbors;
        for (const auto * neighbor:atom->FindNeighborAtoms())
            if (neighbor->GetElement()!=Element::HYDROGEN) neighbors.push_back(neighbor->GetPosition());
        sample_filter::FilterSamplingPointList(points,atom->GetPosition(),neighbors);
        editor.SetAtomLocalRawSamplingEntries(*atom,d::BuildLocalPotentialSampleList(map,points,interpolation));
    }
}

j::object WorkingInformation(const RHBMMemberDataset & dataset,const LocalGaussianResult & local)
{
    j::object out;
    if (!local.fit_result || local.fit_result->data_covariance.rows()!=dataset.y.size()) return out;
    const auto diagonal=local.fit_result->data_covariance.diagonal();
    if (!diagonal.allFinite() || !(diagonal.array()>0.).all()) return out;
    const Eigen::Matrix2d information=dataset.X.transpose()*diagonal.cwiseInverse().asDiagonal()*dataset.X;
    Eigen::LLT<Eigen::Matrix2d> factor(information);
    if (factor.info()!=Eigen::Success) return out;
    const Eigen::Matrix2d covariance=factor.solve(Eigen::Matrix2d::Identity());
    const auto model=local.mdpde.GetModel(); const double a=model.GetAmplitude(),t=model.GetWidth();
    Eigen::Matrix2d jacobian; jacobian << a,-1.5*a*t*t,0.,-.5*t*t*t;
    const Eigen::Matrix2d physical=jacobian*covariance*jacobian.transpose();
    const Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> spectrum(information);
    out["eigenvalues"]=Values({spectrum.eigenvalues()[0],spectrum.eigenvalues()[1]});
    out["beta_covariance"]=j::array{covariance(0,0),covariance(0,1),covariance(1,0),covariance(1,1)};
    out["amplitude_width_covariance"]=j::array{physical(0,0),physical(0,1),physical(1,0),physical(1,1)};
    return out;
}

void SelectRealTargets(ModelObject & model)
{
    model.SelectAllAtoms(false); std::size_t count=0;
    for (const auto & atom:model.GetAtomList())
        if (atom->GetElement()!=Element::HYDROGEN && count++<120) model.SetAtomSelected(atom->GetSerialID(),true);
}

j::object AuditEndpoints(const JointProblem & problem,const std::filesystem::path & directory,
    const std::string & prefix)
{
    const auto read=[](const std::filesystem::path & path) {
        std::ifstream stream(path);
        return j::parse(std::string(std::istreambuf_iterator<char>(stream),{})).as_object();
    };
    const auto reference=read(directory/(prefix+"-F0-production.json"));
    const auto & reference_fit=reference.at("joint").as_object();
    const auto predict=[&](const j::object & fit) {
        std::vector<double> prediction(problem.Input().observations.size(),0.);
        const auto & tau=fit.at("tau").as_array(), & ac=fit.at("ac").as_array();
        for (std::size_t atom=0;atom<tau.size();++atom)
            for (const auto & support:problem.Input().support[atom])
            {
                const auto kernel=n::EvaluateKernel(support.squared_distance,j::value_to<double>(tau[atom]),2.5);
                prediction[support.row]+=j::value_to<double>(ac[2*atom])*kernel.gaussian+
                    j::value_to<double>(ac[2*atom+1])*kernel.charge;
            }
        return prediction;
    };
    const auto reference_prediction=predict(reference_fit);
    j::array states;
    for (const auto & entry:std::filesystem::directory_iterator(directory))
    {
        const auto name=entry.path().filename().string();
        if (!name.starts_with(prefix+"-") || entry.path().extension()!=".json") continue;
        const auto data=read(entry.path()); const auto & fit=data.at("joint").as_object();
        if (!fit.contains("tau")) continue;
        const auto prediction=predict(fit); double prediction_difference=0.,objective=0.;
        for (std::size_t row=0;row<prediction.size();++row)
        {
            prediction_difference=std::max(prediction_difference,std::abs(prediction[row]-reference_prediction[row]));
            const double residual=(prediction[row]-problem.Input().observations[row])/problem.ObservationScale();
            objective+=residual*residual/2.;
        }
        const auto & tau=fit.at("tau").as_array(), & ac=fit.at("ac").as_array();
        j::array changed; double maximum_tau_difference=0.; std::size_t active_widths=0;
        const double bound=data.contains("common_bound") ? j::value_to<double>(data.at("common_bound")) :
            j::value_to<double>(data.at("calibration").as_object().at("bound"));
        for (std::size_t atom=0;atom<tau.size();++atom)
        {
            const double t=j::value_to<double>(tau[atom]),a=j::value_to<double>(ac[2*atom]),c=j::value_to<double>(ac[2*atom+1]);
            if (t>=bound*(1.-1e-10)) ++active_widths;
            const double difference=std::abs(t-j::value_to<double>(reference_fit.at("tau").as_array()[atom]));
            maximum_tau_difference=std::max(maximum_tau_difference,difference);
            if (difference<=1e-6) continue;
            double norm_square=0.,minimum_square=std::numeric_limits<double>::infinity();
            for (const auto & support:problem.Input().support[atom])
            {
                const auto kernel=n::EvaluateKernel(support.squared_distance,t,2.5);
                const double derivative=(a*kernel.d_gaussian_d_log_tau+c*kernel.d_charge_d_log_tau)/problem.ObservationScale();
                norm_square+=derivative*derivative; minimum_square=std::min(minimum_square,support.squared_distance);
            }
            changed.push_back(j::object{{"id",problem.Input().atom_ids[atom]},{"tau",t},{"difference",difference},
                {"amplitude",a},{"charge",c},{"minimum_support_radius",std::sqrt(minimum_square)},
                {"normalized_width_column_norm",std::sqrt(norm_square)}});
        }
        states.push_back(j::object{{"file",name},{"objective_recomputed",objective},
            {"objective_recompute_error",objective-j::value_to<double>(fit.at("objective"))},
            {"prediction_max_difference",prediction_difference},
            {"normalized_prediction_max_difference",prediction_difference/problem.ObservationScale()},
            {"maximum_tau_difference",maximum_tau_difference},{"active_widths",active_widths},{"changed_widths",changed}});
    }
    return {{"reference",prefix+"-F0-production.json"},{"states",states}};
}
}

int main(int argc,char ** argv)
{
    try
    {
        if (argc!=5 && argc!=7) throw std::invalid_argument(
            "Usage: first_stage_sampling_study output.json scenario[:replicate] F0|F1|Q0|Q1|Q8|S0|S1 production|common_bound [model.cif map.map]");
        if (std::string(argv[2])=="audit" && argc==7)
        {
            auto model=ReadModel(argv[5]); auto map=ReadMap(argv[6]); SelectRealTargets(*model);
            const auto problem=BuildPureHaloSupportExclusionProblem(*map,*model);
            std::ofstream stream(argv[1]); stream << j::serialize(AuditEndpoints(problem,argv[4],argv[3])) << '\n';
            return 0;
        }
        const auto start=Clock::now();
        const std::string scenario_arg=argv[2],id=argv[3],bound_mode=argv[4];
        const auto colon=scenario_arg.find(':'); const auto scenario=scenario_arg.substr(0,colon);
        const int replicate=colon==std::string::npos ? 0 : std::stoi(scenario_arg.substr(colon+1));
        if (replicate<0 || replicate>2) throw std::invalid_argument("Replicate must be 0..2.");
        if (scenario!="fine" && scenario!="coarse" && scenario!="noisy" && scenario!="crowded" && scenario!="real")
            throw std::invalid_argument("Unknown study scenario.");
        d::JointSeedSamplingPolicy policy;
        if (id=="Q0" || id=="Q1" || id=="Q8") policy.sampling=d::JointSeedSamplingMethod::QScoreSpiralShells;
        else if (id!="F0" && id!="F1" && id!="S0" && id!="S1") throw std::invalid_argument("Unknown study policy.");
        if (id=="F1" || id=="Q1" || id=="S1") policy.interpolation=d::MapInterpolationMethod::Trilinear;
        const bool real=argc==7,production=bound_mode=="production";
        if (real!=(scenario=="real")) throw std::invalid_argument("Real scenario requires model/map paths.");
        auto model=real ? ReadModel(argv[5]) : second_stage_test::FirstStageStudyModel(scenario=="crowded" ? 1.4 : 5.5,replicate);
        auto map=real ? ReadMap(argv[6]) : second_stage_test::FirstStageStudyMap(*model,scenario=="fine" ? .2 : .4,scenario=="noisy" ? .05 : 1e-4,replicate);
        if (real)
            SelectRealTargets(*model);
        const auto problem=BuildPureHaloSupportExclusionProblem(*map,*model);
        const auto workset=d::MakeJointFittingWorkset(*model,problem);
        model->EditAnalysis().InitializeFromSelection();
        FitOptions options; options.quiet_mode=true;
        const auto init_start=Clock::now();
        d::JointSeedSamplingReport report;
        const auto init=[&] {
            if (id!="Q8" && id!="S0" && id!="S1")
                return d::PrepareJointWidthInitialization(*map,*model,workset,options,policy,&report);
            if (id=="Q8")
            {
                d::RunJointSeedSamplingWorkflow(*map,*model,workset.contributors,1,{},true,policy,&report);
                for (const auto * atom:workset.contributors)
                {
                    auto raw=AtomLocalPotentialView::For(*atom).GetRawSamplingEntries(false);
                    const auto center=raw.front(); raw.insert(raw.begin(),7,center);
                    model->EditAnalysis().SetAtomLocalRawSamplingEntries(*atom,std::move(raw));
                }
            }
            else GeometryControl(*map,*model,workset,policy.interpolation);
            return d::RunJointFirstStageInitializationFromPreparedSamples(*model,workset,options);
        }();
        j::object output{{"policy",id},{"scenario",scenario_arg},{"bound_mode",production ? "production" : "common-research"},
            {"sampling",id=="S0" || id=="S1" ? "matched-spiral-control" : std::string(d::JointSeedSamplingName(policy.sampling))},
            {"interpolation",d::MapInterpolationName(policy.interpolation)},
            {"initialization_seconds",Seconds(init_start)},
            {"grid_spacing",Values({map->GetGridSpacing()[0],map->GetGridSpacing()[1],map->GetGridSpacing()[2]})}};
        j::array records;
        for (std::size_t i=0;i<init.atoms.size();++i)
        {
            const auto & record=init.atoms[i];
            const auto view=AtomLocalPotentialView::For(*workset.contributors[i]);
            const auto & raw=view.GetRawSamplingEntriesView();
            const auto selected=view.GetSamplingEntries(FittingStage::First);
            std::vector<double> responses; std::map<double,std::size_t> radii;
            std::size_t positive=0,fit_count=0,positive_selected=0;
            std::map<double,std::size_t> fitting_radii;
            for (const auto & s:raw) { ++radii[s.point.distance]; if (s.response>0) ++positive; }
            for (const auto & s:selected)
            {
                responses.push_back(s.response); if (s.response>0.) ++positive_selected;
                if (s.point.distance<=1.) { ++fit_count; if (s.response>0.) ++fitting_radii[s.point.distance]; }
            }
            const d::PreparedLocalGaussianDesign design(selected,0.,1.);
            const auto dataset=design.BuildSelectedDataset(responses,GaussianModel3D{0.,1.,0.});
            j::array radial;
            for (const auto & [radius,count]:radii) radial.push_back(j::array{radius,count});
            j::object saved{{"id",record.id},{"raw",raw.size()},{"selected",selected.size()},
                {"positive_raw",positive},{"fitting_candidates",fit_count},{"retained",dataset.y.size()},
                {"positive_selected",positive_selected},{"distinct_fitting_radii",fitting_radii.size()},
                {"rank",dataset.y.size() ? Eigen::ColPivHouseholderQR<RHBMDesignMatrix>(dataset.X).rank() : 0},
                {"radial_counts",radial},{"alpha",Number(record.alpha)},
                {"ols",Values({record.ols[0],record.ols[1],record.ols[2]})},
                {"mdpde",Values({record.mdpde[0],record.mdpde[1],record.mdpde[2]})},
                {"qualified",record.first_stage_fit_qualified},{"reason",record.reason},{"seed_source",record.seed_source},
                {"working_information",WorkingInformation(dataset,view.GetGaussianResult(FittingStage::First))}};
            if (!real) saved["truth"]=Values({1.5+.02*static_cast<double>(i),.35+.005*static_cast<double>(i)});
            if (report.atoms.size()==init.atoms.size())
            {
                j::array shells;
                for (const auto & shell:report.atoms[i].shells)
                    shells.push_back(j::object{{"radius",shell.radius},{"attempts",shell.attempt_count},
                        {"requested",shell.requested_count},{"candidates",shell.candidate_count},
                        {"accepted",shell.accepted_count},{"successful",shell.successful}});
                saved["shells"]=shells;
            }
            records.push_back(saved);
        }
        output["atoms"]=records;
        const auto calibration=n::CalibrateGlobalWidth(init,"study");
        j::object calibrated{{"status",calibration.status},{"reason",calibration.unavailable_reason},
            {"donors",calibration.qualified_donor_count},{"coverage",calibration.qualified_donor_coverage},
            {"fallback",calibration.fallback_count},{"native",calibration.native_success_count},
            {"refined",calibration.refined_success_count},{"unqualified",calibration.unqualified_count},
            {"initial_violations",calibration.initial_violations4}};
        if (calibration.tau_max4) calibrated["bound"]=*calibration.tau_max4;
        if (calibration.p95) { calibrated["p50"]=*calibration.median_width; calibrated["p90"]=*calibration.p90;
            calibrated["p95"]=*calibration.p95; calibrated["p99"]=*calibration.p99; }
        output["calibration"]=calibrated;
        double initial_objective=0.;
        for (double y:problem.Input().observations) initial_objective+=(y/problem.ObservationScale())*(y/problem.ObservationScale())/2.;
        output["initial_objective"]=initial_objective; // Production BCD starts with A=C=0.
        const auto fit_start=Clock::now();
        const auto fit=[&] {
            if (production) return FitJointComponentsInitialized(problem,init);
            const double bound=std::stod(bound_mode);
            if (!(std::isfinite(bound) && bound>0.)) throw std::invalid_argument("Invalid common bound.");
            auto common=calibration; common.status="Available"; common.tau_max4=bound;
            const auto bounded=n::MakeBoundedInitialization(init,common);
            n::BcdExecutionPolicy solver_policy; solver_policy.eta_upper_bound=std::log(bound);
            output["common_bound"]=bound; output["common_projections"]=bounded.projection_count;
            // Research execution, deliberately never labeled Production calibration qualification.
            if (!init.valid) return JointFitResult{};
            return n::FitBcdComponents(problem,bounded.tau_values,solver_policy);
        }();
        output["joint_seconds"]=Seconds(fit_start); output["joint"]=DescribeFit(fit);
        const auto & atom=*model->GetAtomList().front();
        output["qscore"]=CalculateQScoreForAtom(atom,*map,*model,1.,0.);
        output["total_seconds"]=Seconds(start);
        // Validate the unchanged downstream workflow independently of timing above.
        if (production && (id=="F0" || id=="F1" || id=="Q0" || id=="Q1"))
        {
            const auto downstream_start=Clock::now(); d::JointPostprocessingMetrics metrics;
            d::RunJointComponentWorkflow(*map,*model,options,&metrics,policy);
            const auto & snapshot=*model->GetAnalysisView().GetJointResult();
            std::size_t local_count=0,group_count=0;
            for (const auto * target:workset.contributors)
            {
                const auto view=AtomLocalPotentialView::For(*target);
                if (view.GetPostJointLocalResult() && view.GetPostJointLocalResult()->status==EvidenceStatus::Available) ++local_count;
                if (view.GetGroupMemberResult()) ++group_count;
            }
            output["downstream"]=j::object{{"local_available",local_count},{"group_available",group_count},
                {"same_objective",snapshot.objective==fit.objective},{"seconds",Seconds(downstream_start)},
                {"peeling_samples",metrics.grid_peeling_sample_count},{"peeling_available",metrics.available_peeling_sample_count}};
        }
        rusage usage{}; getrusage(RUSAGE_SELF,&usage);
#ifdef __APPLE__
        output["peak_rss_bytes"]=usage.ru_maxrss;
#else
        output["peak_rss_bytes"]=1024*usage.ru_maxrss;
#endif
        std::ofstream stream(argv[1]); stream << j::serialize(output) << '\n';
        return 0;
    }
    catch (const std::exception & error) { std::cerr << error.what() << '\n'; return 1; }
}
