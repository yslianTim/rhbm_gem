#include "support/JointPartialSelection.hpp"
#include "core/detail/PostFitPeeling.hpp"
#include "core/detail/JointUncertainty.hpp"
#include "data/detail/JointStageAdapter.hpp"
#include <rhbm_gem/core/GaussianEstimator.hpp>
#include <rhbm_gem/core/MapSampler.hpp>
#include <rhbm_gem/data/io/DataRepository.hpp>
#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <rhbm_gem/data/object/ModelAnalysisView.hpp>
#include <boost/json.hpp>
#include <chrono>
#include <iostream>
#include <sys/resource.h>
namespace {
namespace rg = rhbm_gem;
namespace c = rg::core;
namespace j = boost::json;
using Clock = std::chrono::steady_clock;
double Seconds(Clock::time_point t) { return std::chrono::duration<double>(Clock::now()-t).count(); }
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
        if(argc!=4) throw std::invalid_argument("Usage: joint_postprocessing_benchmark full|halo|multi workflow|post DATABASE");
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
            for(auto & [id, p]:c::detail::BuildPostFitPeelingSamples(*f.map,*f.model,problem)) editor.SetAtomPostFitPeeling(*f.model->FindAtomPtr(id),std::move(p));
            phases["peeling_seconds"]=Seconds(peeling_start);
            const auto uncertainty_start=Clock::now();
            for(auto & [id,u]:c::detail::ComputeJointUncertainty(problem,snapshot))
            {
                auto & atom=*f.model->FindAtomPtr(id); auto stage=rg::AtomLocalPotentialView::For(atom).GetStageEstimate(rg::FittingStage::Second);
                stage.uncertainty=std::move(u); editor.SetAtomStageEstimate(rg::FittingStage::Second,atom,stage);
                editor.SetAtomGroupEvidence(atom,c::detail::BuildJointParameterEvidence(stage));
            }
            phases["uncertainty_seconds"]=Seconds(uncertainty_start);
            editor.SetJointResult(std::move(snapshot));
            const auto group_start=Clock::now(); c::RunGroupPotentialFitting(*f.model,options);
            phases["group_seconds"]=Seconds(group_start);
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
            j::array covariance, peeled;
            if(stage.uncertainty.covariance) for(int r=0;r<3;++r) for(int col=0;col<3;++col) covariance.push_back((*stage.uncertainty.covariance)(r,col));
            if(view.GetPostFitPeeling()) for(const auto & s:view.GetPostFitPeeling()->samples) peeled.push_back(s.response ? j::value(*s.response):j::value(s.reason));
            targets.push_back(j::object{{"id",atom->GetSerialID()},{"covariance",covariance},{"uncertainty_reason",stage.uncertainty.reason},{"peeled",peeled}});
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
