#include <gtest/gtest.h>
#include "core/detail/joint_component/Problem.hpp"
#include "core/detail/joint_component/StructuralPartition.hpp"
#include "support/JointBlockCoordinate.hpp"
#include "support/JointSyntheticWorkload.hpp"
#include <algorithm>
#include <numeric>
#include <random>
#include <tuple>

namespace {
namespace n=rhbm_gem::core::joint_component;
namespace c=rhbm_gem::core;
namespace b=second_stage_test::joint_block;
n::Vector ObservationVector(const c::JointProblemInput & input)
{return Eigen::Map<const n::Vector>(input.observations.data(),static_cast<Eigen::Index>(input.observations.size()));}
n::Vector InitialEta(const c::JointProblemInput & input)
{return n::Vector::Constant(static_cast<Eigen::Index>(input.atom_ids.size()),std::log(.6));}
n::Vector InitialBeta(const c::JointProblemInput & input)
{
    n::Vector beta= n::Vector::Zero(2*static_cast<Eigen::Index>(input.atom_ids.size()));
    for(std::size_t a=0;a<input.atom_ids.size();++a)
    {
        const double index=static_cast<double>(a);
        beta(2*static_cast<Eigen::Index>(a))=.15+.01*index; beta(2*static_cast<Eigen::Index>(a)+1)=.04-.003*index;
    }
    return beta;
}
void AssertPartition(const n::StructuralBlockPartition & partition,const rhbm_gem::JointParameterLayout & layout,std::size_t atoms)
{
    std::vector<unsigned> membership(atoms);
    for(const auto & core:partition.cores) for(auto atom:core.atoms) ++membership.at(static_cast<std::size_t>(atom));
    for(auto atom:layout.full_atoms) EXPECT_EQ(membership.at(atom),1U);
    for(std::size_t atom=0;atom<atoms;++atom)
        EXPECT_EQ(membership[atom],std::find(layout.full_atoms.begin(),layout.full_atoms.end(),atom)!=layout.full_atoms.end()?1U:0U);
}
void CoreParameters(const b::BlockCoordinateState & state,const n::StructuralCore & core,n::Vector & eta,n::Vector & beta)
{
    eta.resize(static_cast<Eigen::Index>(core.atoms.size())); beta.resize(2*eta.size());
    for(std::size_t k=0;k<core.atoms.size();++k)
    {
        const auto atom=static_cast<Eigen::Index>(core.atoms[k]); eta(static_cast<Eigen::Index>(k))=state.eta(atom);
        beta.segment<2>(2*static_cast<Eigen::Index>(k))=state.beta.segment<2>(2*atom);
    }
}
void ExerciseConditionalReplacements(const c::JointProblemInput & input,const rhbm_gem::JointParameterLayout & layout,
    const n::StructuralBlockPartition & partition,unsigned seed)
{
    AssertPartition(partition,layout,input.atom_ids.size());
    const auto observations=ObservationVector(input);
    const double parent_scale=b::ParentScale(observations);
    auto state=b::Initialize(input,layout,observations,InitialEta(input),InitialBeta(input),parent_scale);
    std::mt19937 generator(seed); std::uniform_real_distribution<double> width_step(-.15,.15),ac_step(-.08,.08);
    for(const auto & core:partition.cores)
    {
        const auto conditional=b::BuildConditionalProblem(input,core,state);
        n::Vector eta,beta; CoreParameters(state,core,eta,beta);
        ASSERT_EQ(conditional.scale,parent_scale);
        ASSERT_EQ(conditional.rows.size(),core.affected_rows.size());
        for(std::size_t k=0;k<conditional.rows.size();++k)
        {
            const auto row=conditional.rows[k];
            EXPECT_DOUBLE_EQ(conditional.effective_response(static_cast<Eigen::Index>(k)),
                conditional.old_core_prediction(static_cast<Eigen::Index>(k))-state.residual(row));
        }
        const double local_before=b::ConditionalObjective(input,core,conditional,eta,beta);
        const double slice_before=.5*std::accumulate(core.affected_rows.begin(),core.affected_rows.end(),0.0,
            [&](double sum,Eigen::Index row){return sum+state.residual(row)*state.residual(row);})/(parent_scale*parent_scale);
        EXPECT_NEAR(local_before,slice_before,2e-12*(1+std::abs(local_before)));

        n::Vector next_eta=eta,next_beta=beta;
        for(Eigen::Index k=0;k<next_eta.size();++k)
        {
            next_eta(k)+=width_step(generator);
            next_beta(2*k)=k==0 ? 0.0 : std::max(0.0,next_beta(2*k)+ac_step(generator));
            next_beta(2*k+1)+=ac_step(generator);
        }
        const double local_after=b::ConditionalObjective(input,core,conditional,next_eta,next_beta);
        const double objective_before=state.objective;
        const auto prior=state;
        b::ReplaceBlock(input,core,next_eta,next_beta,state);
        EXPECT_NEAR((state.objective-objective_before)-(local_after-local_before),
            0.0,4e-12*(1+std::abs(objective_before)));
        for(std::size_t atom=0;atom<input.atom_ids.size();++atom)
            if(std::find(core.atoms.begin(),core.atoms.end(),static_cast<Eigen::Index>(atom))==core.atoms.end())
            {
                EXPECT_DOUBLE_EQ(state.eta(static_cast<Eigen::Index>(atom)),prior.eta(static_cast<Eigen::Index>(atom)));
                EXPECT_TRUE((state.beta.segment<2>(2*static_cast<Eigen::Index>(atom)).array()==
                    prior.beta.segment<2>(2*static_cast<Eigen::Index>(atom)).array()).all());
            }
        for(Eigen::Index row=0;row<state.residual.size();++row)
            if(std::find(core.affected_rows.begin(),core.affected_rows.end(),row)==core.affected_rows.end())
            {
                EXPECT_DOUBLE_EQ(state.prediction(row),prior.prediction(row));
                EXPECT_DOUBLE_EQ(state.residual(row),prior.residual(row));
            }
        const auto replay=b::Replay(input,layout,observations,state.eta,state.beta,parent_scale);
        EXPECT_LT((state.prediction-replay.prediction).lpNorm<Eigen::Infinity>(),2e-12);
        EXPECT_LT((state.residual-replay.residual).lpNorm<Eigen::Infinity>(),2e-12);
        EXPECT_NEAR(state.objective,replay.objective,2e-12*(1+std::abs(replay.objective)));
    }
}
c::JointProblemInput BridgeInput()
{
    c::JointProblemInput input; input.atom_ids={"A","B-bridge","C","D"};
    input.row_ids={"r0","r1","r2","r3","r4"}; input.observations={2.0,-1.0,3.0,.5,-.25};
    input.support={{{0,0.0},{1,.25}},{{1,.25},{2,.5}},{{2,.5},{3,.25}},{{3,.25},{4,0.0}}};
    return input;
}
c::JointProblemInput HaloInput()
{
    c::JointProblemInput input; input.atom_ids={"target-A","halo-multi","target-B","halo-single"};
    input.row_ids={"r0","r1","r2","lambda","constant"}; input.observations={2.0,-3.0,4.0,5.0,-.5};
    input.support={{{0,0.0},{1,.25}},{{1,.25},{2,.5}},{{2,.5}},{{3,0.0}}};
    rhbm_gem::JointSelectionDomain selection; selection.target_indices={0,2}; input.selection_domain=selection;
    return input;
}
}

TEST(JointBlockCoordinateTest, ConditionalResidualAndObjectiveMatchOnChainAndCubeFixtures)
{
    for(const auto & [topology,atoms,seed]:{std::tuple<std::string,int,unsigned>{"chain",8,8},
        {"chain",32,32},{"cube",8,108},{"cube",32,132}})
    {
        auto input=second_stage_test::SyntheticJointWorkload(topology,atoms); const auto layout=n::BuildParameterLayout(input);
        const auto partition=n::BuildStructuralBlockPartition(input,layout,4);
        ExerciseConditionalReplacements(input,layout,partition,seed);
    }
}
TEST(JointBlockCoordinateTest, BridgeAtomRetainsEveryBoundaryCoupling)
{
    const auto input=BridgeInput(); const auto layout=n::BuildParameterLayout(input);
    const auto partition=n::BuildStructuralBlockPartition(input,layout,1);
    const auto bridge=std::find_if(partition.cores.begin(),partition.cores.end(),[](const auto & core){return core.atoms==n::Indices{1};});
    ASSERT_NE(bridge,partition.cores.end());
    EXPECT_EQ(bridge->affected_rows,(n::Indices{1,2}));
    EXPECT_EQ(bridge->neighbor_atoms,(n::Indices{0,2}));
    ExerciseConditionalReplacements(input,layout,partition,501);
}
TEST(JointBlockCoordinateTest, HaloAndContributionOnlyLambdaRowsStayOutOfBlockRows)
{
    const auto input=HaloInput(); const auto layout=n::BuildParameterLayout(input);
    ASSERT_EQ(layout.full_atoms,(std::vector<std::size_t>{0,1,2}));
    ASSERT_EQ(layout.groups.size(),1U); EXPECT_EQ(layout.groups.front().row,3U);
    const auto partition=n::BuildStructuralBlockPartition(input,layout,1);
    const auto first=std::find_if(partition.cores.begin(),partition.cores.end(),[](const auto & core){return core.atoms==n::Indices{0};});
    ASSERT_NE(first,partition.cores.end());
    EXPECT_EQ(first->neighbor_atoms,(n::Indices{1}));
    EXPECT_TRUE(std::find(first->affected_rows.begin(),first->affected_rows.end(),3)==first->affected_rows.end());
    for(const auto & core:partition.cores)
        EXPECT_TRUE(std::find(core.affected_rows.begin(),core.affected_rows.end(),3)==core.affected_rows.end());
    ExerciseConditionalReplacements(input,layout,partition,777);
}
