#include <gtest/gtest.h>
#include "support/JointABCCoverage.hpp"
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <cmath>

namespace {
namespace c=second_stage_test::matched::coverage;
namespace j=boost::json;
}

TEST(JointABCCoverageTest, SyntheticControlsKeepTwelveOverlappingAtomsAndExplicitTruth)
{
    const auto base=c::SyntheticAtoms("baseline"); ASSERT_EQ(base.size(),12);
    EXPECT_EQ(base[0].width,.5); EXPECT_EQ(base[1].width,.45); EXPECT_EQ(base[2].width,.4);
    for (const auto k:{1u,5u,9u})
    {
        EXPECT_DOUBLE_EQ(c::SyntheticAtoms("weak-1e-4")[k].amplitude,base[k].amplitude*1e-4);
        EXPECT_EQ(c::SyntheticAtoms("active-a")[k].amplitude,0);
        EXPECT_EQ(c::SyntheticAtoms("active-a")[k].charge,.2);
    }
    const auto zero=c::SyntheticAtoms("zero-signal"); EXPECT_EQ(zero[5].amplitude+zero[5].charge,0);
    const auto duplicate=c::SyntheticAtoms("duplicate");
    EXPECT_EQ(duplicate[0].position,duplicate[3].position); EXPECT_EQ(duplicate[0].width,duplicate[3].width);
    EXPECT_NEAR(c::SyntheticAtoms("near-0.02")[3].position[0],.02,1e-15);
    EXPECT_THROW(c::SyntheticAtoms("unknown"),std::invalid_argument);
}

TEST(JointABCCoverageTest, StartsDependOnlyOnFirstStageWidthsAndSerialParity)
{
    const Eigen::Vector3d b(.6,.7,.8); const j::array ids{j::object{{"serial_id",2}},j::object{{"serial_id",5}},j::object{{"serial_id",8}}};
    EXPECT_TRUE(c::InitialWidths(b,ids,"first-stage").isApprox(b,0));
    EXPECT_TRUE(c::InitialWidths(b,ids,"mixed").isApprox(Eigen::Vector3d(.72,.56,.96),1e-15));
    EXPECT_THROW(c::InitialWidths(Eigen::Vector3d(.6,0,.8),ids,"first-stage"),std::invalid_argument);
    EXPECT_THROW(c::InitialWidths(b,ids,"truth"),std::invalid_argument);
}

TEST(JointABCCoverageTest, FirstStageRunsOnMapWithoutPeelingAndRejectsWrongIdentity)
{
    auto atom=std::make_unique<rhbm_gem::AtomObject>(); atom->SetSerialID(17);
    atom->SetElement(Element::CARBON); atom->SetPosition(0,0,0);
    const j::array identities{c::Identity(*atom)};
    std::vector<std::unique_ptr<rhbm_gem::AtomObject>> atoms; atoms.push_back(std::move(atom));
    rhbm_gem::ModelObject model(std::move(atoms));
    rhbm_gem::MapObject map({25,25,25},{.3,.3,.3},{-3.6,-3.6,-3.6});
    auto values=std::make_unique<double[]>(map.GetMapValueArraySize());
    for (std::size_t k=0;k<map.GetMapValueArraySize();++k)
    {
        const auto p=map.GetGridPosition(k); const double square=p[0]*p[0]+p[1]*p[1]+p[2]*p[2];
        values[k]=static_cast<float>(6*std::pow(2*M_PI*.25,-1.5)*std::exp(-square/.5));
    }
    map.SetMapValueArray(std::move(values));
    const auto first=c::Initialize(model,map,identities);
    ASSERT_TRUE(first.valid)<<j::serialize(first.evidence);
    EXPECT_NEAR(first.b(0),.5,.01); EXPECT_EQ(first.evidence.at("uses_peeling"),false);
    EXPECT_EQ(j::value_to<double>(first.evidence.at("atoms").at(0).at("mdpde").at(2)),0);
    const auto second=c::Initialize(model,map,identities);
    EXPECT_EQ(first.b(0),second.b(0));
    auto wrong=identities; wrong[0].at("position").as_array()[0]=.1;
    EXPECT_THROW(c::Initialize(model,map,wrong),std::runtime_error);
}
