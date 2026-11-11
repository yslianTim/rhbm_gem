#include <gtest/gtest.h>
#include "core/estimator/InitialStage.hpp"
#include "core/sampling/MapSampler.hpp"
#include <rhbm_gem/core/QScoreHelper.hpp>
#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <algorithm>
#include <memory>

namespace {
using namespace rhbm_gem;
using namespace rhbm_gem::core;
namespace d = rhbm_gem::core::detail;
std::unique_ptr<ModelObject> Model()
{
    std::vector<std::unique_ptr<AtomObject>> atoms;
    for (int i=0;i<2;++i)
    {
        auto atom=std::make_unique<AtomObject>();
        atom->SetSerialID(i+1); atom->SetElement(Element::CARBON);
        atom->SetPosition(.13+i*1.3,.17,.19); atoms.push_back(std::move(atom));
    }
    auto model=std::make_unique<ModelObject>(std::move(atoms));
    model->SelectAllAtoms(); model->EditAnalysis().InitializeFromSelection(); return model;
}
MapObject Map()
{
    auto values=std::make_unique<double[]>(32*32*32);
    for (int z=0;z<32;++z) for (int y=0;y<32;++y) for (int x=0;x<32;++x)
    {
        const double a=-3.+x*.25,b=-3.+y*.25,c=-3.+z*.25;
        values[static_cast<std::size_t>(x+32*(y+32*z))]=.01+std::exp(-.5*(a*a+b*b+c*c));
    }
    return MapObject{{32,32,32},{.25,.25,.25},{-3.,-3.,-3.},std::move(values)};
}
void SameSamples(const LocalPotentialSampleList & a,const LocalPotentialSampleList & b)
{
    ASSERT_EQ(a.size(),b.size());
    for (std::size_t i=0;i<a.size();++i)
    {
        EXPECT_EQ(a[i].point.position,b[i].point.position);
        EXPECT_EQ(a[i].point.distance,b[i].point.distance);
        EXPECT_EQ(a[i].point.is_selected,b[i].point.is_selected);
        EXPECT_EQ(a[i].response,b[i].response);
    }
}
void SameInitialization(const JointInitialization & a,const JointInitialization & b)
{
    EXPECT_EQ(a.tau,b.tau); EXPECT_EQ(a.valid,b.valid); ASSERT_EQ(a.atoms.size(),b.atoms.size());
    for (std::size_t i=0;i<a.atoms.size();++i)
    {
        EXPECT_EQ(a.atoms[i].ols,b.atoms[i].ols); EXPECT_EQ(a.atoms[i].mdpde,b.atoms[i].mdpde);
        EXPECT_EQ(a.atoms[i].alpha,b.atoms[i].alpha);
        EXPECT_EQ(a.atoms[i].sample_count,b.atoms[i].sample_count);
        EXPECT_EQ(a.atoms[i].initial_stage_qualification,b.atoms[i].initial_stage_qualification);
    }
}
}

TEST(InitialStageSamplingTest, DefaultExactlyMatchesExplicitInitialStageSampling)
{
    auto map=Map(); auto explicit_sampling=Model(); auto current=Model();
    FitOptions options; options.quiet_mode=true;
    const auto explicit_workset=d::BuildInitialStageWorkset(*explicit_sampling,BuildPureHaloSupportExclusionProblem(map,*explicit_sampling));
    RunPotentialSamplingWorkflow(map,*explicit_sampling,explicit_workset.contributors,1,{},true);
    const auto explicit_init=d::FitInitialStageFromPreparedSamples(*explicit_sampling,explicit_workset,options);
    const auto new_workset=d::BuildInitialStageWorkset(*current,BuildPureHaloSupportExclusionProblem(map,*current));
    const auto new_init=d::RunInitialStage(map,*current,new_workset,options);
    SameInitialization(explicit_init,new_init);
    for (std::size_t i=0;i<2;++i)
        SameSamples(AtomLocalPotentialView::For(*explicit_workset.contributors[i]).GetRawSamplingEntriesView(),
            AtomLocalPotentialView::For(*new_workset.contributors[i]).GetRawSamplingEntriesView());
}

TEST(InitialStageSamplingTest, DeterministicAcrossWorkersAndRepeatedRuns)
{
    auto map=Map();
    for (const auto method : { SphereSamplingMethod::FibonacciDeterministic,
             SphereSamplingMethod::SpiralDeterministic })
    {
        FitOptions options; options.quiet_mode=true; options.initial_stage_sampling_method=method;
        auto serial=Model(); auto parallel=Model();
        const auto problem=BuildPureHaloSupportExclusionProblem(map,*serial);
        const auto sw=d::BuildInitialStageWorkset(*serial,problem);
        const auto pw=d::BuildInitialStageWorkset(*parallel,problem);
        const auto baseline=d::RunInitialStage(map,*serial,sw,options);
        ASSERT_TRUE(baseline.valid);
        for (const int workers : { 2, 4, 4 })
        {
            SCOPED_TRACE(workers);
            options.joint_initialization_thread_size=workers;
            const auto result=d::RunInitialStage(map,*parallel,pw,options);
            SameInitialization(baseline,result);
            for (std::size_t i=0;i<sw.contributors.size();++i)
                SameSamples(AtomLocalPotentialView::For(*sw.contributors[i]).GetRawSamplingEntriesView(),
                    AtomLocalPotentialView::For(*pw.contributors[i]).GetRawSamplingEntriesView());
        }
    }
}

TEST(InitialStageSamplingTest, ParallelSamplingFailurePreservesEveryPreviousRawEntry)
{
    auto model=Model(); auto map=Map(); const auto atoms=model->GetSelectedAtoms();
    RunPotentialSamplingWorkflow(map,*model,atoms,1,{},true);
    std::vector<LocalPotentialSampleList> before;
    for(const auto * atom:atoms) before.push_back(AtomLocalPotentialView::For(*atom).GetRawSamplingEntriesView());
    AtomObject unattached;
    auto batch=atoms; batch.insert(batch.begin()+1,&unattached);
    auto values=std::make_unique<double[]>(32*32*32);
    std::fill_n(values.get(),32*32*32,2.);
    MapObject replacement{{32,32,32},{.25,.25,.25},{-3.,-3.,-3.},std::move(values)};
    for (const auto method : { SphereSamplingMethod::FibonacciDeterministic,
             SphereSamplingMethod::SpiralDeterministic })
    {
        EXPECT_THROW(RunPotentialSamplingWorkflow(replacement,*model,batch,4,{},true,method),std::runtime_error);
        for(std::size_t i=0;i<atoms.size();++i)
            SameSamples(before[i],AtomLocalPotentialView::For(*atoms[i]).GetRawSamplingEntriesView());
    }
}

TEST(InitialStageSamplingTest, SpiralInitializesWidthsAndLeavesQScoreUnchanged)
{
    auto map=Map(); auto model=Model(); const auto atoms=model->GetSelectedAtoms();
    const auto score=CalculateQScoreForAtom(*atoms[0],map,*model,1.,0.);
    const auto default_samples=SampleAtomMapValues(map,*atoms[0]);
    FitOptions options; options.quiet_mode=true;
    options.initial_stage_sampling_method=SphereSamplingMethod::SpiralDeterministic;
    const auto workset=d::BuildInitialStageWorkset(*model,BuildPureHaloSupportExclusionProblem(map,*model));
    const auto initialization=d::RunInitialStage(map,*model,workset,options);
    ASSERT_TRUE(initialization.valid);
    ASSERT_EQ(initialization.atoms.size(),atoms.size());
    EXPECT_EQ(initialization.atoms.front().seed_source,"fitted");
    for (const auto * atom : workset.contributors)
        SameSamples(SampleAtomMapValues(map,*atom,options.initial_stage_sampling_method),
            AtomLocalPotentialView::For(*atom).GetRawSamplingEntriesView());
    const auto & spiral=AtomLocalPotentialView::For(*atoms[0]).GetRawSamplingEntriesView();
    ASSERT_FALSE(spiral.empty());
    ASSERT_FALSE(default_samples.empty());
    EXPECT_NE(spiral.front().point.position,default_samples.front().point.position);
    EXPECT_EQ(CalculateQScoreForAtom(*atoms[0],map,*model,1.,0.),score);
}
