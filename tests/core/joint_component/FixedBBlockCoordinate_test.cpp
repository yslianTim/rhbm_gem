#include <gtest/gtest.h>
#include "core/detail/joint_component/FixedBBlockCoordinate.hpp"
#include "core/detail/joint_component/Problem.hpp"
#include "support/JointOperatorWorkload.hpp"
#include <algorithm>

namespace {
namespace n=rhbm_gem::core::joint_component;
using rhbm_gem::core::JointProblemInput;
n::Vector Select(n::VectorRef values,const std::vector<std::size_t> & indices)
{
    n::Vector out(static_cast<Eigen::Index>(indices.size()));
    for(std::size_t k=0;k<indices.size();++k) out(static_cast<Eigen::Index>(k))=values(static_cast<Eigen::Index>(indices[k]));
    return out;
}
void CheckFixedBCase(const std::string & topology,int atoms)
{
    auto input=std::make_shared<JointProblemInput>(second_stage_test::OperatorWorkload(topology,atoms));
    const auto layout=n::BuildParameterLayout(*input);
    const auto parent_context=n::CreateContext(input);
    const n::Domain parent_domain(input); const auto profile_domain=n::ProfileDomain(parent_domain,layout);
    const auto context=n::ProfileContext(parent_context,layout,parent_domain.rows);
    n::Vector observations=Eigen::Map<const n::Vector>(input->observations.data(),static_cast<Eigen::Index>(input->observations.size()));
    const n::Vector eta=n::Vector::Constant(static_cast<Eigen::Index>(input->atom_ids.size()),std::log(.6));
    const auto global=n::EvaluateProfile(profile_domain,Select(observations,layout.informative_rows),
        Select(eta,layout.full_atoms),false,&context);
    ASSERT_TRUE(global.valid)<<global.reason;
    n::FixedBBlockPolicy policy; policy.core_atoms=128; policy.maximum_sweeps=200;
    const auto block=n::SearchFixedBBlocks(*input,layout,observations,eta,context,policy);
    ASSERT_TRUE(block.success)<<block.reason<<" blocks="<<block.blocks.size()
        <<" before="<<(block.blocks.empty()?0:block.blocks.back().objective_before)
        <<" after="<<(block.blocks.empty()?0:block.blocks.back().objective_after)
        <<" local-reduction="<<(block.blocks.empty()?0:block.blocks.back().objective_reduction)
        <<" kkt="<<(block.sweeps.empty()?0:block.sweeps.back().global_ac_kkt);
    EXPECT_LE(global.certificate.projected_kkt,1e-10);
    EXPECT_LT(block.sweeps_to_global_kkt,policy.maximum_sweeps+1);
    for(auto atom:layout.full_atoms) EXPECT_GE(block.state.beta(2*static_cast<Eigen::Index>(atom)),0.0);
    const double global_objective=global.certificate.objective/(context.scale*context.scale);
    EXPECT_NEAR(block.state.objective,global_objective,1e-11*(1+std::abs(global_objective)));
    EXPECT_LE(block.sweeps.back().global_ac_kkt,1e-10);
    const auto global_prediction=(global.x*global.beta).eval();
    for(const auto & sweep:block.sweeps)
    {
        EXPECT_LE(sweep.objective_after,sweep.objective_before);
        EXPECT_LE(sweep.cache_replay_error,2e-12);
        EXPECT_LE(sweep.objective_replay_error,1e-12);
    }
    for(std::size_t k=0;k<layout.informative_rows.size();++k)
    {
        const auto row=static_cast<Eigen::Index>(layout.informative_rows[k]);
        EXPECT_NEAR(block.state.prediction(row),global_prediction(static_cast<Eigen::Index>(k)),
            2e-10*(1+std::abs(block.state.prediction(row))));
    }
    for(std::size_t k=0;k<layout.full_atoms.size();++k)
    {
        const auto atom=static_cast<Eigen::Index>(layout.full_atoms[k]);
        EXPECT_NEAR(block.state.beta(2*atom),global.beta(2*static_cast<Eigen::Index>(k)),
            2e-8*(1+std::abs(global.beta(2*static_cast<Eigen::Index>(k)))));
        EXPECT_NEAR(block.state.beta(2*atom+1),global.beta(2*static_cast<Eigen::Index>(k)+1),
            2e-8*(1+std::abs(global.beta(2*static_cast<Eigen::Index>(k)+1))));
        EXPECT_DOUBLE_EQ(block.state.eta(atom),eta(atom));
    }
}
}
TEST(JointFixedBBlockCoordinateTest, FixedBBlockSweepsRecoverSmallGlobalConstrainedProfiles)
{
    CheckFixedBCase("chain",8);
    CheckFixedBCase("cube",8);
}
