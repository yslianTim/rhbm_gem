#include <gtest/gtest.h>
#include "core/detail/FirstStageInitialization.hpp"
#include "core/detail/FittingWorkset.hpp"
#include <rhbm_gem/core/MapSampler.hpp>
#include <rhbm_gem/core/QScoreHelper.hpp>
#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
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
        EXPECT_DOUBLE_EQ(a[i].response,b[i].response);
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
        EXPECT_EQ(a.atoms[i].first_stage_qualification,b.atoms[i].first_stage_qualification);
    }
}
}

TEST(JointFirstStageSamplingTest, DefaultExactlyMatchesLegacySamplingAndEstimator)
{
    auto map=Map(); auto legacy=Model(); auto current=Model();
    FitOptions options; options.quiet_mode=true;
    const auto old_workset=d::MakeTwoStageFittingWorkset(*legacy);
    RunPotentialSamplingWorkflow(map,*legacy,old_workset.contributors,1,{},true);
    const auto old_init=d::RunJointFirstStageInitializationFromPreparedSamples(*legacy,old_workset,options);
    const auto new_workset=d::MakeTwoStageFittingWorkset(*current);
    const auto new_init=d::PrepareJointWidthInitialization(map,*current,new_workset,options);
    SameInitialization(old_init,new_init);
    for (std::size_t i=0;i<2;++i)
        SameSamples(AtomLocalPotentialView::For(*old_workset.contributors[i]).GetRawSamplingEntriesView(),
            AtomLocalPotentialView::For(*new_workset.contributors[i]).GetRawSamplingEntriesView());
}

TEST(JointFirstStageSamplingTest, DeterministicAcrossWorkersAndRepeatedRuns)
{
    auto map=Map(); FitOptions options; options.quiet_mode=true;
    auto serial=Model(); auto parallel=Model();
    auto sw=d::MakeTwoStageFittingWorkset(*serial),pw=d::MakeTwoStageFittingWorkset(*parallel);
    const auto a=d::PrepareJointWidthInitialization(map,*serial,sw,options);
    options.joint_initialization_thread_size=4;
    const auto b=d::PrepareJointWidthInitialization(map,*parallel,pw,options);
    SameInitialization(a,b);
    for (std::size_t i=0;i<2;++i)
        SameSamples(AtomLocalPotentialView::For(*sw.contributors[i]).GetRawSamplingEntriesView(),
            AtomLocalPotentialView::For(*pw.contributors[i]).GetRawSamplingEntriesView());
    const auto c=d::PrepareJointWidthInitialization(map,*parallel,pw,options);
    SameInitialization(b,c);
}

TEST(JointFirstStageSamplingTest, ParallelSamplingFailurePreservesEveryPreviousRawEntry)
{
    auto model=Model(); auto map=Map(); const auto atoms=model->GetSelectedAtoms();
    RunPotentialSamplingWorkflow(map,*model,atoms,1,{},true);
    std::vector<LocalPotentialSampleList> before;
    for(const auto * atom:atoms) before.push_back(AtomLocalPotentialView::For(*atom).GetRawSamplingEntries(false));
    AtomObject unattached;
    auto batch=atoms; batch.insert(batch.begin()+1,&unattached);
    EXPECT_THROW(RunPotentialSamplingWorkflow(map,*model,batch,4,{},true),std::runtime_error);
    for(std::size_t i=0;i<atoms.size();++i)
        SameSamples(before[i],AtomLocalPotentialView::For(*atoms[i]).GetRawSamplingEntriesView());
}

TEST(JointFirstStageSamplingTest, FirstStageLeavesIndependentQScoreUnchanged)
{
    auto map=Map(); auto model=Model(); const auto atoms=model->GetSelectedAtoms();
    const auto score=CalculateQScoreForAtom(*atoms[0],map,*model,1.,0.);
    FitOptions options; options.quiet_mode=true;
    d::PrepareJointWidthInitialization(map,*model,d::MakeTwoStageFittingWorkset(*model),options);
    EXPECT_EQ(CalculateQScoreForAtom(*atoms[0],map,*model,1.,0.),score);
}
