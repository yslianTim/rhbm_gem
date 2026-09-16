#include <gtest/gtest.h>
#include "support/ObservationMatchedExperiment.hpp"
#include "support/EstimatedNeighborSweep.hpp"
#include "support/ForwardModelExperiment.hpp"
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/utils/math/ElectricPotential.hpp>
#include <cmath>
#include <set>

namespace {
namespace m = second_stage_test::matched;
m::Design ExampleDesign()
{
    m::Design design;
    for (int i=0;i<50;++i) design.push_back({{std::pow(0.04*i,2),1.0}});
    return design;
}
}

TEST(ObservationMatchedTest, CompactStencilMatchesSamplerWithAnisotropicGridAndClampedEdges)
{
    rhbm_gem::MapObject grid({8,9,10},{0.1,0.13,0.17},{-0.713,-0.827,-0.919});
    auto values{std::make_unique<double[]>(grid.GetMapValueArraySize())};
    for (std::size_t i=0;i<grid.GetMapValueArraySize();++i) values[i]=std::sin(static_cast<double>(i));
    grid.SetMapValueArray(std::move(values));
    SamplingPointList points{{0.0,grid.GetOrigin(),true},{0.0,grid.GetGridPosition(377),true},
        {0.0,{-0.702,-0.80,-0.90},true},{0.0,{-0.122,0.15,0.55},true},
        {0.0,{-0.71300000001,-0.827,-0.919},true}};
    const auto samples{second_stage_test::SampleExperimentPoints(grid,points)};
    for (std::size_t i=0;i<points.size();++i)
    {
        const auto stencil{m::MakeStencil(grid,grid,points[i].position)};
        double prediction{}, total{};
        for (const auto & s : stencil.slots) {prediction += s.coefficient*grid.GetMapValue(s.index); total += s.coefficient;}
        EXPECT_NEAR(prediction,samples[i].response,2e-14); EXPECT_NEAR(total,1.0,2e-14);
    }
    const auto edge{m::MakeStencil(grid,grid,points[2].position)};
    std::set<std::size_t> unique; bool negative{};
    for (const auto & s : edge.slots) {unique.insert(s.index); negative |= s.coefficient<0.0;}
    EXPECT_TRUE(edge.boundary); EXPECT_LT(unique.size(),64u); EXPECT_TRUE(negative);
}

TEST(ObservationMatchedTest, KernelMatchesGeneratorIncludingBothCutoffsAndChargeSigns)
{
    ElectricPotential generator; generator.SetBlurringWidth(0.5);
    for (double cutoff : {2.5,3.0}) for (double r : {0.0,0.000001,0.00001,0.3,2.499999,2.5,2.500001,3.01})
        for (double c : {-0.3,0.0,0.3})
        {
            const auto b{m::EvaluateBasis(r*r,0.5,cutoff)};
            const double expected{r*r<=cutoff*cutoff ? generator.GetPotentialValue(Element::CARBON,r,c) : 0.0};
            EXPECT_NEAR(6*b.gaussian+c*b.charge,expected,2e-15);
        }
}

TEST(ObservationMatchedTest, StencilIncludesContributorsBeyondSampleCenterCutoff)
{
    rhbm_gem::MapObject grid({61,61,61},{0.1,0.1,0.1},{-3.0,-3.0,-3.0});
    const auto stencil{m::MakeStencil(grid,grid,{2.51,0.03,0.02})};
    const std::vector<m::Atom> atoms{{{0,0,0},6,0.5,0.3}};
    EXPECT_DOUBLE_EQ(m::EvaluateBasis(m::SquareDistance({2.51,0.03,0.02},{0,0,0}),0.5,2.5).charge,0.0);
    EXPECT_GT(std::abs(m::Predict(stencil,atoms,2.5)),1e-4);
}

TEST(ObservationMatchedTest, QuantizationHappensAfterContributorSum)
{
    rhbm_gem::MapObject grid({4,4,4},{0.1,0.1,0.1},{0,0,0});
    const auto stencil{m::MakeStencil(grid,grid,{0,0,0})};
    const std::vector<m::Atom> atoms{{{0,0,0},6,0.5,0.3},{{0,0,0},1e-7,0.5,-0.3}};
    double magnitudes{}, bound{};
    const double value{m::Predict(stencil,atoms,2.5,false,&magnitudes,&bound)};
    const double rounded{m::Predict(stencil,atoms,2.5,true)};
    EXPECT_DOUBLE_EQ(rounded,static_cast<double>(static_cast<float>(value)));
    EXPECT_LE(std::abs(rounded-value),bound+1e-15); EXPECT_GE(magnitudes,std::abs(value));
}

TEST(ObservationMatchedTest, WidthDerivativesIncludeChargeCoupling)
{
    for (double square : {0.0,1e-12,0.25,6.25,6.26})
    {
        const auto b{m::EvaluateBasis(square,0.5,3.0)};
        const auto plus{m::EvaluateBasis(square,0.5*std::exp(1e-5),3.0)};
        const auto minus{m::EvaluateBasis(square,0.5*std::exp(-1e-5),3.0)};
        EXPECT_NEAR((plus.gaussian-minus.gaussian)/2e-5,b.gaussian_log_width,1e-8);
        EXPECT_NEAR((plus.charge-minus.charge)/2e-5,b.charge_log_width,1e-8);
    }
}

TEST(ObservationMatchedTest, IdentifiableExactFitsRecoverPositiveNegativeAndZeroCharge)
{
    const auto design{ExampleDesign()}; const auto basis{m::EvaluateDesign(design,0.47,2.5)};
    for (double c : {-0.3,0.0,0.3}) for (bool free_charge : {false,true})
    {
        const Eigen::VectorXd y{6.3*basis.col(0)+c*basis.col(1)};
        const auto fit{m::Fit(design,y,c,free_charge,2.5)};
        ASSERT_TRUE(fit.at("qualified").as_bool()) << fit.at("reason");
        EXPECT_NEAR(fit.at("A").as_double(),6.3,6.3e-6);
        EXPECT_NEAR(fit.at("B").as_double(),0.47,1e-6);
        EXPECT_NEAR(fit.at("C").as_double(),c,1e-6);
    }
}

TEST(ObservationMatchedTest, DegeneracyBoundsAndBudgetAreNotQualified)
{
    const m::Design repeated(10,{{0.0,1.0}});
    EXPECT_EQ(m::Fit(repeated,Eigen::VectorXd::Ones(10),0.0,true,2.5).at("reason"),"rank-deficient");
    const auto design{ExampleDesign()};
    const Eigen::VectorXd wide{m::EvaluateDesign(design,3.0,2.5).col(0)};
    EXPECT_EQ(m::Fit(design,wide,0.0,false,2.5).at("reason"),"parameter-boundary");
    const Eigen::VectorXd y{m::EvaluateDesign(design,0.47,2.5).col(0)};
    EXPECT_EQ(m::Fit(design,y,0.0,false,2.5,0).at("reason"),"budget-exhausted");
    EXPECT_FALSE(m::Fit(design,-y,0.0,false,2.5).at("qualified").as_bool());
}

TEST(ObservationMatchedTest, EstimatedNeighborsAreFrozenAndTargetOrderDoesNotChangeFits)
{
    rhbm_gem::MapObject grid({61,61,61},{0.1,0.1,0.1},{-3,-3,-3});
    const std::vector<m::Atom> state{{{0,0,0},6.2,0.47,-0.2},{{1.1,0.2,0.1},7.1,0.52,0.3}};
    std::vector<m::Position> positions;
    std::vector<m::Stencil> stencils;
    std::vector<double> distances;
    Eigen::VectorXd observations(80);
    for (int i=0;i<80;++i)
    {
        const double r{0.025*i};
        positions.push_back({r,0,0}); distances.push_back(r);
        stencils.push_back(m::MakeStencil(grid,grid,positions.back()));
        observations(i)=m::Predict(stencils.back(),state,2.5);
    }
    const auto prepared{m::PrepareFrozenTarget(state,0,positions,stencils,2.5)};
    auto first = m::FitFrozenTarget(prepared,distances,observations,state[0],2.5);
    const auto other{m::PrepareFrozenTarget(state,1,positions,stencils,2.5)};
    static_cast<void>(m::FitFrozenTarget(other,distances,observations,state[1],2.5));
    auto repeated = m::FitFrozenTarget(m::PrepareFrozenTarget(state,0,positions,stencils,2.5),distances,observations,state[0],2.5);
    for (auto & row : first) row.as_object().erase("seconds");
    for (auto & row : repeated) row.as_object().erase("seconds");
    EXPECT_EQ(first,repeated);
    for (const auto & row : first)
    {
        if (row.at("parameters")=="AB") EXPECT_DOUBLE_EQ(row.at("C").as_double(),state[0].charge);
        if (row.at("prediction")!="matched") continue;
        ASSERT_TRUE(row.at("qualified").as_bool()) << row.at("reason");
        EXPECT_NEAR(row.at("A").as_double(),state[0].amplitude,1e-5);
        EXPECT_NEAR(row.at("B").as_double(),state[0].width,1e-6);
        EXPECT_NEAR(row.at("C").as_double(),state[0].charge,1e-6);
    }
    auto changed{state}; changed[1].charge+=0.1;
    const auto perturbed{m::PrepareFrozenTarget(changed,0,positions,stencils,2.5)};
    EXPECT_GT((prepared.neighbors[0]-perturbed.neighbors[0]).norm(),0.1);
    EXPECT_GT((prepared.neighbors[1]-perturbed.neighbors[1]).norm(),0.1);
    changed=state; changed[0].amplitude+=1; changed[0].width+=0.1; changed[0].charge+=0.2;
    const auto target_changed{m::PrepareFrozenTarget(changed,0,positions,stencils,2.5)};
    EXPECT_EQ(prepared.neighbors[0],target_changed.neighbors[0]);
    EXPECT_EQ(prepared.neighbors[1],target_changed.neighbors[1]);
    const auto abc = m::FitFrozenTarget(target_changed,distances,observations,changed[0],2.5).back();
    for (const char * p : {"A","B","C"}) EXPECT_EQ(abc.at(p),first.back().at(p));
}

TEST(ObservationMatchedTest, FrozenTargetRejectsInvalidCheckpointAndObservations)
{
    rhbm_gem::MapObject grid({4,4,4},{0.1,0.1,0.1},{0,0,0});
    const std::vector<m::Position> positions{{0.1,0.1,0.1}};
    const std::vector<m::Stencil> stencils{m::MakeStencil(grid,grid,positions[0])};
    std::vector<m::Atom> state{{{0,0,0},6,0.5,0}};
    EXPECT_THROW(m::PrepareFrozenTarget(state,1,positions,stencils,2.5),std::invalid_argument);
    const auto prepared{m::PrepareFrozenTarget(state,0,positions,stencils,2.5)};
    EXPECT_THROW(m::FitFrozenTarget(prepared,{0.1},Eigen::VectorXd::Constant(1,NAN),state[0],2.5),std::invalid_argument);
    state[0].charge=NAN;
    EXPECT_THROW(m::PrepareFrozenTarget(state,0,positions,stencils,2.5),std::invalid_argument);
}
