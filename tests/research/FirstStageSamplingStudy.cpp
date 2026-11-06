// Reproducible, single-process FirstStage/JCE study. No production instrumentation.
#include "core/detail/FirstStageInitialization.hpp"
#include "core/detail/FittingWorkset.hpp"
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
#include <boost/json.hpp>
#include <Eigen/QR>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
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
std::unique_ptr<ModelObject> SyntheticModel(double separation)
{
    std::vector<std::unique_ptr<AtomObject>> atoms;
    for (int i=0;i<24;++i)
    {
        auto atom=std::make_unique<AtomObject>();
        atom->SetSerialID(i+1); atom->SetElement(Element::CARBON);
        atom->SetPosition(.13+separation*(i%4),.17+separation*((i/4)%3),.19+separation*(i/12));
        atoms.push_back(std::move(atom));
    }
    auto model=std::make_unique<ModelObject>(std::move(atoms));
    model->SelectAllAtoms(); return model;
}
std::unique_ptr<MapObject> SyntheticMap(const ModelObject & model,double spacing,double noise)
{
    const auto last=model.GetAtomList().back()->GetPosition();
    const std::array<double,3> origin{-3.,-3.,-3.};
    std::array<int,3> dims{};
    for (std::size_t k=0;k<3;++k) dims[k]=static_cast<int>(std::ceil((last[k]+6.)/spacing));
    const auto count=static_cast<std::size_t>(dims[0])*static_cast<std::size_t>(dims[1])*static_cast<std::size_t>(dims[2]);
    auto values=std::make_unique<double[]>(count);
    for (int z=0;z<dims[2];++z) for (int y=0;y<dims[1];++y) for (int x=0;x<dims[0];++x)
    {
        const std::array<double,3> p{origin[0]+x*spacing,origin[1]+y*spacing,origin[2]+z*spacing};
        double value=noise*std::sin(2.13*x+1.57*y+.73*z);
        for (std::size_t i=0;i<model.GetAtomList().size();++i)
        {
            const auto q=model.GetAtomList()[i]->GetPosition(); double square=0;
            for (std::size_t k=0;k<3;++k) square+=(p[k]-q[k])*(p[k]-q[k]);
            value+=(1.5+.02*static_cast<double>(i))*n::EvaluateKernel(square,.35+.005*static_cast<double>(i),2.5).gaussian;
        }
        const auto index=static_cast<std::size_t>(x)+static_cast<std::size_t>(dims[0])*(static_cast<std::size_t>(y)+static_cast<std::size_t>(dims[1])*static_cast<std::size_t>(z));
        values[index]=value;
    }
    return std::make_unique<MapObject>(dims,std::array<double,3>{spacing,spacing,spacing},origin,std::move(values));
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
}

int main(int argc,char ** argv)
{
    try
    {
        if (argc!=3 && argc!=5) throw std::invalid_argument("Usage: first_stage_sampling_study output.json spacing [model.cif map.map]");
        const auto start=Clock::now();
        auto model=argc==5 ? ReadModel(argv[3]) : SyntheticModel(5.5);
        auto map=argc==5 ? ReadMap(argv[4]) : SyntheticMap(*model,std::stod(argv[2]),1e-4);
        if (argc==5)
        {
            model->SelectAllAtoms(false); std::size_t count=0;
            for (const auto & atom:model->GetAtomList())
                if (atom->GetElement()!=Element::HYDROGEN && count++<120) model->SetAtomSelected(atom->GetSerialID(),true);
        }
        const auto problem=BuildPureHaloSupportExclusionProblem(*map,*model);
        const auto workset=d::MakeJointFittingWorkset(*model,problem);
        model->EditAnalysis().InitializeFromSelection();
        FitOptions options; options.quiet_mode=true;
        const auto init_start=Clock::now();
        const auto init=d::PrepareJointWidthInitialization(*map,*model,workset,options);
        j::object output{{"policy","F0"},{"initialization_seconds",Seconds(init_start)},
            {"grid_spacing",Values({map->GetGridSpacing()[0],map->GetGridSpacing()[1],map->GetGridSpacing()[2]})}};
        j::array records;
        for (std::size_t i=0;i<init.atoms.size();++i)
        {
            const auto & record=init.atoms[i];
            const auto view=AtomLocalPotentialView::For(*workset.contributors[i]);
            const auto & raw=view.GetRawSamplingEntriesView();
            const auto selected=view.GetSamplingEntries(FittingStage::First);
            std::vector<double> responses; std::map<double,std::size_t> radii;
            std::size_t positive=0,fit_count=0;
            for (const auto & s:raw) { ++radii[s.point.distance]; if (s.response>0) ++positive; }
            for (const auto & s:selected) { responses.push_back(s.response); if (s.point.distance<=1.) ++fit_count; }
            const d::PreparedLocalGaussianDesign design(selected,0.,1.);
            const auto dataset=design.BuildSelectedDataset(responses,GaussianModel3D{0.,1.,0.});
            j::array radial;
            for (const auto & [radius,count]:radii) radial.push_back(j::array{radius,count});
            records.push_back(j::object{{"id",record.id},{"raw",raw.size()},{"selected",selected.size()},
                {"positive_raw",positive},{"fitting_candidates",fit_count},{"retained",dataset.y.size()},
                {"rank",dataset.y.size() ? Eigen::ColPivHouseholderQR<RHBMDesignMatrix>(dataset.X).rank() : 0},
                {"radial_counts",radial},{"alpha",Number(record.alpha)},
                {"ols",Values({record.ols[0],record.ols[1],record.ols[2]})},
                {"mdpde",Values({record.mdpde[0],record.mdpde[1],record.mdpde[2]})},
                {"qualified",record.first_stage_fit_qualified},{"reason",record.reason},{"seed_source",record.seed_source}});
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
        const auto fit_start=Clock::now();
        const auto fit=FitJointComponentsInitialized(problem,init);
        output["joint_seconds"]=Seconds(fit_start); output["joint"]=DescribeFit(fit);
        const auto & atom=*model->GetAtomList().front();
        output["qscore"]=CalculateQScoreForAtom(atom,*map,*model,1.,0.);
        output["total_seconds"]=Seconds(start);
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
