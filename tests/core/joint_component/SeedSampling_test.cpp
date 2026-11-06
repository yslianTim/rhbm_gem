#include <gtest/gtest.h>
#include "core/detail/FirstStageInitialization.hpp"
#include "core/detail/FittingWorkset.hpp"
#include <rhbm_gem/core/MapSampler.hpp>
#include <rhbm_gem/core/QScoreHelper.hpp>
#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <memory>
#include "support/FibonacciNeighborSelection.hpp"
#include <rhbm_gem/utils/math/ArrayHelper.hpp>

namespace {
using namespace rhbm_gem;
using namespace rhbm_gem::core;
namespace d = rhbm_gem::core::detail;
constexpr d::JointSeedSamplingPolicy f0{},
    f1{d::JointSeedSamplingMethod::FibonacciShells,d::MapInterpolationMethod::Trilinear},
    q0{d::JointSeedSamplingMethod::QScoreSpiralShells,d::MapInterpolationMethod::Tricubic},
    q1{d::JointSeedSamplingMethod::QScoreSpiralShells,d::MapInterpolationMethod::Trilinear};
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
void SameSamples(const LocalPotentialSampleList & a,const LocalPotentialSampleList & b,bool responses=true)
{
    ASSERT_EQ(a.size(),b.size());
    for (std::size_t i=0;i<a.size();++i)
    {
        EXPECT_EQ(a[i].point.position,b[i].point.position);
        EXPECT_EQ(a[i].point.distance,b[i].point.distance);
        EXPECT_EQ(a[i].point.is_selected,b[i].point.is_selected);
        if (responses) EXPECT_DOUBLE_EQ(a[i].response,b[i].response);
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

TEST(JointSeedSamplingTest, DefaultExactlyMatchesLegacySamplingAndEstimator)
{
    auto map=Map(); auto legacy=Model(); auto current=Model();
    FitOptions options; options.quiet_mode=true;
    const auto old_workset=d::MakeTwoStageFittingWorkset(*legacy);
    RunPotentialSamplingWorkflow(map,*legacy,old_workset.contributors,1,{},true);
    const auto old_init=d::RunJointFirstStageInitializationFromPreparedSamples(*legacy,old_workset,options);
    const auto new_workset=d::MakeTwoStageFittingWorkset(*current);
    d::JointSeedSamplingReport report;
    const auto new_init=d::PrepareJointWidthInitialization(map,*current,new_workset,options,{},&report);
    SameInitialization(old_init,new_init);
    ASSERT_EQ(report.atoms.size(),2u);
    EXPECT_EQ(report.policy.sampling,f0.sampling); EXPECT_EQ(report.policy.interpolation,f0.interpolation);
    for (std::size_t i=0;i<2;++i)
        SameSamples(AtomLocalPotentialView::For(*old_workset.contributors[i]).GetRawSamplingEntriesView(),
            AtomLocalPotentialView::For(*new_workset.contributors[i]).GetRawSamplingEntriesView());
}

TEST(JointSeedSamplingTest, FactorialSwitchesAreIndependentAndQScoreIsUnaffected)
{
    auto map=Map(); auto model=Model(); const auto atoms=model->GetSelectedAtoms();
    const auto score=CalculateQScoreForAtom(*atoms[0],map,*model,1.,0.);
    std::vector<LocalPotentialSampleList> samples;
    for (auto policy:{f0,f1,q0,q1})
    {
        d::JointSeedSamplingReport report;
        d::RunJointSeedSamplingWorkflow(map,*model,atoms,1,{},true,policy,&report);
        samples.push_back(AtomLocalPotentialView::For(*atoms[0]).GetRawSamplingEntries(false));
        EXPECT_EQ(report.policy.sampling,policy.sampling);
        EXPECT_EQ(report.policy.interpolation,policy.interpolation);
        for (const auto & sample:samples.back())
            EXPECT_DOUBLE_EQ(sample.response,d::InterpolateMapValue(map,sample.point.position,policy.interpolation));
        if (policy.sampling==d::JointSeedSamplingMethod::QScoreSpiralShells)
        {
            ASSERT_EQ(report.atoms[0].shells.size(),20u);
            EXPECT_EQ(report.atoms[0].positive_fitting_count,report.atoms[0].fitting_count);
        }
        EXPECT_DOUBLE_EQ(CalculateQScoreForAtom(*atoms[0],map,*model,1.,0.),score);
    }
    SameSamples(samples[0],samples[1],false); SameSamples(samples[2],samples[3],false);
    EXPECT_NE(samples[0].size(),samples[2].size());
    EXPECT_NE(samples[0][0].response,samples[1][0].response);
}

TEST(JointSeedSamplingTest, PoliciesAreDeterministicAcrossWorkersAndRepeatedRuns)
{
    auto map=Map(); FitOptions options; options.quiet_mode=true;
    for (auto policy:{f0,f1,q0,q1})
    {
        auto serial=Model(); auto parallel=Model();
        auto sw=d::MakeTwoStageFittingWorkset(*serial),pw=d::MakeTwoStageFittingWorkset(*parallel);
        const auto a=d::PrepareJointWidthInitialization(map,*serial,sw,options,policy);
        options.joint_initialization_thread_size=4;
        const auto b=d::PrepareJointWidthInitialization(map,*parallel,pw,options,policy);
        SameInitialization(a,b);
        for (std::size_t i=0;i<2;++i)
            SameSamples(AtomLocalPotentialView::For(*sw.contributors[i]).GetRawSamplingEntriesView(),
                AtomLocalPotentialView::For(*pw.contributors[i]).GetRawSamplingEntriesView());
        d::PrepareJointWidthInitialization(map,*parallel,pw,options,q1);
        const auto c=d::PrepareJointWidthInitialization(map,*parallel,pw,options,policy);
        SameInitialization(b,c);
        options.joint_initialization_thread_size=1;
    }
}

TEST(JointSeedSamplingTest, ParallelInterpolationFailurePreservesPreviousRawEntries)
{
    auto model=Model(); auto map=Map(); const auto atoms=model->GetSelectedAtoms();
    d::RunJointSeedSamplingWorkflow(map,*model,atoms,1,{},true,f0);
    const auto before=AtomLocalPotentialView::For(*atoms[0]).GetRawSamplingEntries(false);
    auto values=std::make_unique<double[]>(1); values[0]=1.;
    MapObject small{{1,1,1},{.1,.1,.1},{0.,0.,0.},std::move(values)};
    EXPECT_THROW(d::RunJointSeedSamplingWorkflow(small,*model,atoms,4,{},true,q1),std::out_of_range);
    SameSamples(before,AtomLocalPotentialView::For(*atoms[0]).GetRawSamplingEntriesView());
}


TEST(JointSeedSamplingTest, NeighborControlsPreserveRawGeometryResponsesAndParallelEstimates)
{
    auto map=Map(); FitOptions options; options.quiet_mode=true;
    for (bool qscore:{false,true})
    {
        auto serial=Model(),parallel=Model();
        std::vector<JointInitialization> initializations;
        for (auto * model:{serial.get(),parallel.get()})
        {
            const int workers=model==serial.get() ? 1 : 4;
            const auto workset=d::MakeTwoStageFittingWorkset(*model);
            d::RunJointSeedSamplingWorkflow(map,*model,workset.contributors,workers,{},true,f0);
            for (const auto * atom:workset.contributors)
            {
                const auto baseline=AtomLocalPotentialView::For(*atom).GetRawSamplingEntries(false);
                const auto masks=second_stage_test::FibonacciNeighborMasks(*atom,*model,baseline);
                auto controlled=baseline;
                second_stage_test::ApplyFibonacciNeighborSelection(controlled,masks,qscore);
                ASSERT_EQ(controlled.size(),200u);
                for (std::size_t k=0;k<controlled.size();++k)
                {
                    const auto & m=masks.rows[k];
                    EXPECT_EQ(m.f0,baseline[k].point.is_selected);
                    EXPECT_EQ(m.f0,!(m.ownership || m.cone));
                    EXPECT_EQ(controlled[k].point.position,baseline[k].point.position);
                    EXPECT_EQ(controlled[k].point.distance,baseline[k].point.distance);
                    EXPECT_EQ(controlled[k].response,baseline[k].response);
                    EXPECT_EQ(controlled[k].point.is_selected,!(qscore ? m.qscore : m.ownership));
                    if (m.shared_qscore) EXPECT_TRUE(m.ownership);
                    if (baseline[k].point.distance<=1. && m.qscore) EXPECT_TRUE(m.ownership);
                    bool independently_rejected=false;
                    const double r=baseline[k].point.distance;
                    for (const auto * neighbor:model->FindNeighborAtoms(*atom,2.*r))
                    {
                        const double distance=array_helper::ComputeNorm(baseline[k].point.position,neighbor->GetPositionRef());
                        independently_rejected |= neighbor->GetElement()!=Element::HYDROGEN && distance<.9*r;
                    }
                    EXPECT_EQ(m.qscore,independently_rejected);
                }
                auto repeated=baseline;
                second_stage_test::ApplyFibonacciNeighborSelection(repeated,masks,qscore);
                SameSamples(controlled,repeated);
                model->EditAnalysis().SetAtomLocalRawSamplingEntries(*atom,std::move(controlled));
            }
            options.joint_initialization_thread_size=workers;
            initializations.push_back(d::RunJointFirstStageInitializationFromPreparedSamples(*model,workset,options));
        }
        SameInitialization(initializations[0],initializations[1]);
    }
}

namespace {
std::unique_ptr<ModelObject> NeighborControlModel(double neighbor_position,Element element=Element::CARBON)
{
    std::vector<std::unique_ptr<AtomObject>> atoms;
    for (int i=0;i<2;++i)
    {
        auto atom=std::make_unique<AtomObject>(); atom->SetSerialID(i+1);
        atom->SetElement(i==0 ? Element::CARBON : element);
        atom->SetPosition(i==0 ? -1. : neighbor_position,0.,0.);
        atoms.push_back(std::move(atom));
    }
    return std::make_unique<ModelObject>(std::move(atoms));
}
}

TEST(JointSeedSamplingTest, NeighborControlStrictBoundaries)
{
    const LocalPotentialSampleList raw{{1.,SamplingPoint{1.,{0.,0.,0.},true}}};
    auto ownership=NeighborControlModel(1.);
    EXPECT_FALSE(second_stage_test::FibonacciNeighborMasks(*ownership->FindAtomPtr(1),*ownership,raw).rows[0].ownership);
    for (double neighbor:{std::nextafter(.9,0.),.9,std::nextafter(.9,1.)})
    {
        auto model=NeighborControlModel(neighbor);
        const auto masks=second_stage_test::FibonacciNeighborMasks(*model->FindAtomPtr(1),*model,raw);
        EXPECT_EQ(masks.rows[0].qscore,neighbor<.9);
    }
    auto model=NeighborControlModel(1.); // center-neighbor distance equals 2r
    const auto masks=second_stage_test::FibonacciNeighborMasks(*model->FindAtomPtr(1),*model,raw);
    EXPECT_EQ(masks.shells.at(1.).qscore_ids,(std::vector<int>{2}));
}

TEST(JointSeedSamplingTest, NeighborControlsIgnoreHydrogenAndEmptyNeighborSets)
{
    const LocalPotentialSampleList raw{{1.,SamplingPoint{1.,{0.,0.,0.},true}}};
    for (double position:{.1,10.})
    {
        auto model=NeighborControlModel(position,position==.1 ? Element::HYDROGEN : Element::CARBON);
        const auto masks=second_stage_test::FibonacciNeighborMasks(*model->FindAtomPtr(1),*model,raw);
        EXPECT_TRUE(masks.rows[0].f0); EXPECT_FALSE(masks.rows[0].ownership);
        EXPECT_FALSE(masks.rows[0].qscore); EXPECT_FALSE(masks.rows[0].cone);
        EXPECT_TRUE(masks.shells.at(1.).fixed_ids.empty());
    }
}

TEST(JointSeedSamplingTest, NeighborUniverseDifferenceCanReverseRawDomainInclusion)
{
    auto model=NeighborControlModel(2.); // distance 3 A: outside fixed 2.5 A query
    const LocalPotentialSampleList raw{{1.,SamplingPoint{1.,{0.,0.,0.},true}},
        {1.,SamplingPoint{1.95,{.95,0.,0.},true}}};
    const auto masks=second_stage_test::FibonacciNeighborMasks(*model->FindAtomPtr(1),*model,raw);
    EXPECT_FALSE(masks.rows[0].qscore);
    EXPECT_TRUE(masks.rows[1].f0); EXPECT_FALSE(masks.rows[1].ownership);
    EXPECT_TRUE(masks.rows[1].qscore); EXPECT_FALSE(masks.rows[1].shared_qscore);
    EXPECT_TRUE(masks.shells.at(1.95).fixed_ids.empty());
    EXPECT_EQ(masks.shells.at(1.95).qscore_ids,(std::vector<int>{2}));
}
