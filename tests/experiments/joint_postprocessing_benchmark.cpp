#include "support/JointPartialSelection.hpp"
#include "core/detail/PostFitPeeling.hpp"
#include "core/detail/LocalConditionalCovariance.hpp"
#include "core/detail/joint_component/ResourceWork.hpp"
#include "data/detail/JointStageAdapter.hpp"
#include <rhbm_gem/core/GaussianEstimator.hpp>
#include <rhbm_gem/core/MapSampler.hpp>
#include <rhbm_gem/data/io/DataRepository.hpp>
#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <rhbm_gem/data/object/ModelAnalysisView.hpp>
#include <boost/json.hpp>
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <iostream>
#include <sys/resource.h>
#include <thread>
namespace {
namespace rg = rhbm_gem;
namespace c = rg::core;
namespace n = c::joint_component;
namespace j = boost::json;
using Clock = std::chrono::steady_clock;
double Seconds(Clock::time_point t) { return std::chrono::duration<double>(Clock::now()-t).count(); }
class RssMonitor
{
    std::atomic<bool> m_running{true};
    std::atomic<std::size_t> m_peak{};
    std::thread m_thread;
    void Observe()
    {
        const auto rss=n::CurrentRssBytes();
        if(!rss) return;
        auto peak=m_peak.load();
        while(peak<*rss && !m_peak.compare_exchange_weak(peak,*rss)) {}
    }
public:
    RssMonitor():m_thread([this] {
        while(m_running.load())
        {
            Observe();
            std::this_thread::sleep_for(std::chrono::milliseconds(5));
        }
        Observe();
    }) { Observe(); }
    ~RssMonitor() { m_running=false; m_thread.join(); }
    std::size_t Peak() { Observe(); return m_peak.load(); }
};
int RunLargeLocalCovariance()
{
    constexpr std::size_t observations=5'938'087, contributors=37'406, targets=37'406;
    constexpr double amplitude=1.8, width=.65, charge=-.3;
    rg::LocalStageEstimate stage;
    stage.point=rg::GaussianModel3D(amplitude,width,charge);
    stage.source={rg::EstimateMethod::JointComponents,"1","large-local","large-local-run",rg::FittingRole::Target};
    stage.convergence=rg::JointCheckStatus::Passed;
    LocalPotentialSampleList raw(observations);
    rg::PostFitPeelingResult peeling;
    peeling.source=stage.source;
    peeling.samples.resize(observations);
    for(std::size_t i=0;i<observations;++i)
    {
        const double distance=.05+2.3*static_cast<double>(i%1024)/1023.0;
        raw[i].point.distance=distance;
        const auto basis=n::EvaluateKernel(distance*distance,width,2.5);
        peeling.samples[i].response=amplitude*basis.gaussian+charge*basis.charge+
            .001*std::sin(.13*static_cast<double>(i%997));
    }

    c::detail::LocalConditionalCovariance covariance;
    const auto covariance_start=Clock::now();
    std::size_t covariance_peak_rss_bytes{};
    {
        RssMonitor rss;
        covariance=c::detail::ComputeLocalConditionalCovariance(stage,raw,peeling,true);
        covariance_peak_rss_bytes=rss.Peak();
    }
    rusage usage{};
    getrusage(RUSAGE_SELF,&usage);
#ifdef __APPLE__
    const auto process_peak_rss_bytes=static_cast<std::size_t>(usage.ru_maxrss);
#else
    const auto process_peak_rss_bytes=static_cast<std::size_t>(usage.ru_maxrss)*1024;
#endif
    const auto local_workspace=c::detail::EstimateLocalConditionalCovarianceWorkspaceBytes(observations);
    const auto global_columns=3*contributors;
    const auto global_covariance_bytes=global_columns*global_columns*sizeof(double);
    const auto all_target_covariance_bytes=targets*sizeof(Eigen::Matrix2d);
    std::cout<<j::serialize(j::object{
        {"case","large-local"},{"mode","local-covariance"},{"observations",observations},
        {"global_shape_contributors",contributors},{"global_shape_targets",targets},{"local_targets_processed",1},
        {"local_observation_count",covariance.observation_count},
        {"local_covariance_status",static_cast<int>(covariance.status)},{"local_covariance_reason",covariance.reason},
        {"local_covariance_seconds",Seconds(covariance_start)},
        {"local_covariance_peak_rss_bytes",covariance_peak_rss_bytes},
        {"local_covariance_workspace_estimate_bytes",local_workspace ? *local_workspace : 0},
        {"raw_sample_storage_bytes",raw.size()*sizeof(LocalPotentialSample)},
        {"peeling_response_storage_bytes",peeling.samples.size()*sizeof(rg::PeelingSampleEstimate)},
        {"per_target_covariance_storage_bytes",sizeof(Eigen::Matrix2d)},
        {"all_target_covariance_storage_bytes",all_target_covariance_bytes},
        {"hypothetical_global_covariance_bytes",global_covariance_bytes},
        {"process_peak_rss_bytes",process_peak_rss_bytes}})<<'\n';
    return covariance.status==rg::EvidenceStatus::Available ? 0 : 1;
}
joint_partial_test::Fixture Make(const std::string & name)
{
    joint_partial_test::Fixture f;
    const int components = name == "multi" ? 3 : 1, count = name == "halo" ? 5 : 3;
    std::vector<std::unique_ptr<rg::AtomObject>> atoms;
    for (int group=0; group<components; ++group) for (int k=0; k<count; ++k)
    {
        auto a=std::make_unique<rg::AtomObject>();
        a->SetSerialID(group*count+k+1); a->SetElement(Element::CARBON);
        a->SetPosition(12.*group+.7*k,0,0); a->SetChainID("A"); a->SetComponentID("ALA"); a->SetAtomID("CA");
        atoms.push_back(std::move(a));
        f.a.push_back(2+.2*k); f.b.push_back(.5+.03*k); f.c.push_back(.1-.02*k);
    }
    f.model=std::make_unique<rg::ModelObject>(std::move(atoms));
    f.model->SelectAtoms([&](const auto & a) { return name == "full" || (a.GetSerialID()-1)%count == 0; });
    const std::array<int,3> dims{static_cast<int>((12.*(components-1)+.7*(count-1)+6.4)/.16)+1,41,41};
    const std::array<double,3> spacing{.16,.16,.16}, origin{-3.2,-3.2,-3.2};
    auto values=std::make_unique<double[]>(static_cast<std::size_t>(dims[0]*dims[1]*dims[2]));
    for(int z=0;z<dims[2];++z) for(int y=0;y<dims[1];++y) for(int x=0;x<dims[0];++x)
    {
        const auto p=c::simulation::GridPosition({x,y,z},spacing,origin);
        double value=.0001*std::sin(.7*x+.3*y+.2*z);
        for(std::size_t a=0;a<f.a.size();++a)
        {
            const auto kernel=c::joint_component::EvaluateKernel(c::simulation::SupportSquare(p,f.model->GetAtomList()[a]->GetPosition()),f.b[a],2.5);
            value+=f.a[a]*kernel.gaussian+f.c[a]*kernel.charge;
        }
        values[static_cast<std::size_t>(x+dims[0]*(y+dims[1]*z))]=value;
    }
    f.map=std::make_unique<rg::MapObject>(dims,spacing,origin,std::move(values));
    return f;
}
}
int main(int argc,char ** argv)
{
    try
    {
        if(argc==2 && std::string(argv[1])=="local-large") return RunLargeLocalCovariance();
        if(argc!=4) throw std::invalid_argument("Usage: joint_postprocessing_benchmark full|halo|multi workflow|post DATABASE, or local-large");
        const std::string name=argv[1], mode=argv[2];
        if(name!="full" && name!="halo" && name!="multi") return 2;
        auto f=Make(name); auto editor=f.model->EditAnalysis();
        c::FitOptions options; options.estimator=c::PotentialEstimator::JOINT_COMPONENTS; options.quiet_mode=true; options.thread_size=1;
        const auto start=Clock::now(); j::object phases;
        if(mode=="workflow")
        {
            c::RunPotentialFittingWorkflow(*f.map,*f.model,options);
            phases["workflow_seconds"]=Seconds(start);
        }
        else if(mode=="post")
        {
            const auto problem=c::BuildJointProblem(*f.map,*f.model);
            auto snapshot=c::CaptureJointAnalysisResult(c::FitJointComponents(problem,f.b));
            editor.InitializeFromSelection();
            for(const auto & atom:f.model->GetAtomList())
                editor.SetAtomLocalRawSamplingEntries(*atom,c::SampleAtomMapValues(*f.map,*atom,SphereSamplingMethod::FibonacciDeterministic));
            rg::data_internal::ApplyJointStageEstimates(*f.model,snapshot,"benchmark");
            const auto peeling_start=Clock::now();
            RssMonitor peeling_rss;
            for(auto & [id, p]:c::detail::BuildPostFitPeelingSamples(*f.map,*f.model,problem)) editor.SetAtomPostFitPeeling(*f.model->FindAtomPtr(id),std::move(p));
            phases["peeling_seconds"]=Seconds(peeling_start);
            phases["post_fit_peeling_peak_rss_bytes"]=peeling_rss.Peak();
            const auto covariance_start=Clock::now();
            RssMonitor covariance_rss;
            std::size_t local_workspace_estimate_bytes=0;
            for(const auto * selected:f.model->GetSelectedAtoms())
            {
                const auto local=rg::AtomLocalPotentialView::For(*selected);
                const auto & stage=local.GetStageEstimate(rg::FittingStage::Second);
                const auto & peeling=local.GetPostFitPeeling();
                const auto raw=local.GetRawSamplingEntries(false);
                c::detail::LocalConditionalCovariance covariance;
                if(peeling)
                {
                    covariance=c::detail::ComputeLocalConditionalCovariance(stage,raw,*peeling,local.HasSampleGeometry());
                    const auto workspace=c::detail::EstimateLocalConditionalCovarianceWorkspaceBytes(
                        covariance.observation_count);
                    local_workspace_estimate_bytes=std::max(local_workspace_estimate_bytes,
                        raw.size()*sizeof(LocalPotentialSample)+(workspace ? *workspace : 0));
                }
                covariance.source=stage.source;
                editor.SetAtomGroupEvidence(*selected,c::detail::BuildLocalGroupEvidence(
                    stage,peeling ? &*peeling : nullptr,covariance));
            }
            phases["local_covariance_seconds"]=Seconds(covariance_start);
            phases["local_covariance_peak_rss_bytes"]=covariance_rss.Peak();
            phases["local_covariance_workspace_estimate_bytes"]=local_workspace_estimate_bytes;
            phases["per_target_covariance_storage_bytes"]=f.model->GetSelectedAtoms().size()*sizeof(Eigen::Matrix2d);
            editor.SetJointResult(std::move(snapshot));
            const auto group_start=Clock::now(); RssMonitor group_rss; c::RunGroupPotentialFitting(*f.model,options);
            phases["group_seconds"]=Seconds(group_start);
            phases["group_fitting_peak_rss_bytes"]=group_rss.Peak();
        }
        else return 2;
        const auto save_start=Clock::now();
        {rg::DataRepository repository(argv[3]); repository.SaveModel(*f.model,"model");}
        phases["save_seconds"]=Seconds(save_start);
        phases["total_seconds"]=Seconds(start);
        rusage usage{}; getrusage(RUSAGE_SELF,&usage);
        j::array endpoints, targets;
        const auto & snapshot=*f.model->GetAnalysisView().GetJointResult();
        for(const auto & component:snapshot.components)
            endpoints.push_back(j::object{{"id",component.id},{"convergence",static_cast<int>(component.runtime_convergence)},
                {"ac",component.state ? j::value_from(component.state->ac):j::value{}}, {"b",component.state ? j::value_from(component.state->b):j::value{}}});
        for(const auto * atom:f.model->GetSelectedAtoms())
        {
            const auto view=rg::AtomLocalPotentialView::For(*atom); const auto & stage=view.GetStageEstimate(rg::FittingStage::Second);
            j::array covariance, estimate, peeled, point;
            if(stage.point)
            {
                point={stage.point->GetAmplitude(),stage.point->GetWidth(),stage.point->GetOffset()};
            }
            if(view.GetGroupEvidence())
            {
                const auto & evidence=*view.GetGroupEvidence();
                if(evidence.covariance) for(int r=0;r<2;++r) for(int col=0;col<2;++col) covariance.push_back((*evidence.covariance)(r,col));
                if(evidence.estimate) for(int k=0;k<2;++k) estimate.push_back((*evidence.estimate)(k));
            }
            if(view.GetPostFitPeeling()) for(const auto & s:view.GetPostFitPeeling()->samples) peeled.push_back(s.response ? j::value(*s.response):j::value(s.reason));
            const auto & evidence=view.GetGroupEvidence();
            targets.push_back(j::object{{"id",atom->GetSerialID()},{"point",point},
                {"estimate",estimate},{"covariance",covariance},{"evidence_status",evidence ? static_cast<int>(evidence->status):-1},
                {"uncertainty_method",evidence ? evidence->uncertainty_method:""},
                {"uncertainty_reason",evidence ? evidence->reason:"missing-evidence"},{"peeled",peeled}});
        }
#ifdef __APPLE__
        const auto peak=usage.ru_maxrss;
#else
        const auto peak=usage.ru_maxrss*1024;
#endif
        std::cout<<j::serialize(j::object{{"case",name},{"mode",mode},{"phases",phases},{"process_peak_rss_bytes",peak},
            {"objective",snapshot.objective ? j::value(*snapshot.objective):j::value{}},{"endpoints",endpoints},{"targets",targets}})<<'\n';
    }
    catch(const std::exception & e) {std::cerr<<e.what()<<'\n'; return 1;}
}
