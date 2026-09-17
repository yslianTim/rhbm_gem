#include <gtest/gtest.h>
#include "support/UniqueStencilGrid.hpp"
#include "support/ForwardModelExperiment.hpp"
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/utils/math/ElectricPotential.hpp>
#include <cmath>
#include <set>
#include <random>

namespace {
namespace m=second_stage_test::matched;
namespace g=m::unique_grid;
namespace ac=m::joint_ac;
Eigen::MatrixXd Design()
{
    Eigen::MatrixXd x(120,4);
    for (int i=0;i<120;++i) {const double t{static_cast<double>(i)/120}; x.row(i)<<1+t,std::sin(7*t),std::cos(4*t),t*t;}
    return x;
}
Eigen::VectorXd Response(const Eigen::MatrixXd & x)
{
    Eigen::VectorXd y{x*Eigen::Vector4d(2,-.4,1,.2)};
    for (Eigen::Index i=0;i<y.size();++i) y(i)+=.05*std::sin(3.1*static_cast<double>(i));
    return y;
}
}

TEST(UniqueStencilGridTest, DeduplicatesAllSlotsIncludingZerosClampsAndRepeatedSamples)
{
    rhbm_gem::MapObject map({8,9,10},{.1,.13,.17},{-.713,-.827,-.919});
    auto values{std::make_unique<double[]>(map.GetMapValueArraySize())};
    for (std::size_t i=0;i<map.GetMapValueArraySize();++i) values[i]=i%3==0 ? 0 : std::sin(static_cast<double>(i));
    map.SetMapValueArray(std::move(values));
    SamplingPointList points{{0,map.GetOrigin(),true},{0,{-.702,-.80,-.90},false},{0,map.GetOrigin(),false}};
    std::vector<m::Stencil> stencils;
    for (const auto & p:points) stencils.push_back(m::MakeStencil(map,map,p.position));
    const auto grid{g::BuildGrid(stencils,map)};
    std::set<std::size_t> expected; bool zero{},negative{},negative_response{},zero_response{}; std::size_t slots{};
    for (const auto & s:stencils) for (const auto & v:s.slots)
    {expected.insert(v.index); zero|=v.coefficient==0; negative|=v.coefficient<0;}
    ASSERT_EQ(grid.voxels.size(),expected.size()); EXPECT_LT(grid.voxels.size(),64u);
    Eigen::VectorXd y(static_cast<Eigen::Index>(grid.voxels.size()));
    for (std::size_t p=0;p<grid.voxels.size();++p)
    {
        const auto & v{grid.voxels[p]}; slots+=v.multiplicity;
        y(static_cast<Eigen::Index>(p))=v.observed; negative_response|=v.observed<0; zero_response|=v.observed==0;
        if (p) EXPECT_LT(grid.voxels[p-1].index,v.index);
    }
    EXPECT_EQ(slots,192u); EXPECT_TRUE(zero); EXPECT_TRUE(negative); EXPECT_TRUE(negative_response); EXPECT_TRUE(zero_response);
    const auto replay{g::Project(grid,stencils,y)};
    const auto native{second_stage_test::SampleExperimentPoints(map,points)};
    for (std::size_t p=0;p<native.size();++p) EXPECT_NEAR(replay(static_cast<Eigen::Index>(p)),native[p].response,2e-14);
    EXPECT_EQ(grid.sample_rows[0],grid.sample_rows[2]);
}

TEST(UniqueStencilGridTest, KeepsGenerationCoordinatesSeparateFromSamplingHeader)
{
    rhbm_gem::MapObject generation({8,9,10},{.1,.13,.17},{-.713,-.827,-.919});
    rhbm_gem::MapObject sampling({8,9,10},{.10000001,.12999999,.17000001},{-.71300001,-.82699999,-.91900001});
    sampling.ClearMapValueArray();
    const std::vector<m::Stencil> stencils{m::MakeStencil(generation,sampling,{-.69,-.79,-.85}),
        m::MakeStencil(generation,sampling,{-2,-2,-2})};
    const auto grid{g::BuildGrid(stencils,sampling)};
    EXPECT_TRUE(stencils[1].boundary);
    for (const auto & v:grid.voxels) {EXPECT_EQ(v.position,generation.GetGridPosition(v.index)); EXPECT_NE(v.position,sampling.GetGridPosition(v.index));}
    auto changed=stencils; changed[0].slots[1]=changed[0].slots[0]; changed[0].slots[1].position[0]+=1;
    EXPECT_THROW(g::BuildGrid(changed,sampling),std::invalid_argument);
    changed=stencils; changed[0].slots[0].index=sampling.GetMapValueArraySize();
    EXPECT_THROW(g::BuildGrid(changed,sampling),std::invalid_argument);
}

TEST(UniqueStencilGridTest, DirectDesignMatchesIndependentGeneratorAndUsesCheckpointWidths)
{
    ElectricPotential generator;
    for (double width:{.3,.5,.8}) for (double charge:{-.3,0.,.3}) for (double cutoff:{2.5,3.0})
    {
        generator.SetBlurringWidth(width); g::Grid grid;
        for (double r:{0.,1e-6,1e-5,.3,2.499999,2.5,2.500001,3.01}) grid.voxels.push_back({0,1,{r,0,0},0});
        const std::vector<m::Atom> atoms{{{0,0,0},6,width,charge}};
        const auto x{g::BuildDesign(grid,atoms,cutoff)}; const Eigen::VectorXd prediction{x*Eigen::Vector2d(6,charge)};
        for (std::size_t p=0;p<grid.voxels.size();++p)
        {
            const double r{grid.voxels[p].position[0]};
            const double expected{r*r<=cutoff*cutoff ? generator.GetPotentialValue(Element::CARBON,r,charge) : 0};
            EXPECT_NEAR(prediction(static_cast<Eigen::Index>(p)),expected,1e-12);
            EXPECT_NEAR(g::Direct(grid.voxels[p].position,atoms,cutoff),expected,1e-12);
        }
        auto changed=atoms; changed[0].amplitude=0; changed[0].charge=0;
        EXPECT_EQ(g::BuildDesign(grid,changed,cutoff),x);
    }
}

TEST(UniqueStencilGridTest, CommonAlphaFitsPreserveQualificationAndBothHundredStepBudgets)
{
    const auto x{Design()}; const auto y{Response(x)};
    const auto reference{ac::WeightedSolve(x,y,Eigen::VectorXd::Ones(y.size()),true)}; ASSERT_TRUE(reference.valid);
    for (double alpha:g::alphas)
    {
        const auto fit{g::Fit(x,y,Eigen::Vector4d(2,-.4,1,.2),alpha)};
        EXPECT_EQ(fit.at("iteration_budget"),100); EXPECT_EQ(fit.at("refinement_budget"),100);
        ASSERT_EQ(fit.at("blocks").as_array().size(),1u); EXPECT_EQ(fit.at("blocks").at(0).at("scope"),"global");
        EXPECT_FALSE(fit.at("blocks").at(0).as_object().contains("owner"));
        if (alpha<1) EXPECT_TRUE(fit.at("qualified").as_bool())<<fit.at("reason");
        else {EXPECT_FALSE(fit.at("qualified").as_bool()); EXPECT_EQ(fit.at("reason"),"budget-exhausted");}
        for (const auto & b:fit.at("branches").as_array())
        {
            EXPECT_LE(b.at("primary").at("iterations").as_int64(),100);
            EXPECT_LE(b.at("reference").at("iterations").as_int64(),100);
            EXPECT_LE(b.at("trace").as_array().size(),101u); EXPECT_LE(b.at("reference_trace").as_array().size(),101u);
            EXPECT_EQ(b.at("endpoint_diagnostics").at("rank"),4);
            if (alpha==1) {EXPECT_EQ(b.at("primary").at("iterations"),100); EXPECT_TRUE(b.at("reference_trace").as_array().empty());}
        }
        if (alpha==0)
        {
            for (std::size_t k=0;k<4;++k) EXPECT_NEAR(boost::json::value_to<double>(fit.at("beta").at(k)),reference.beta(static_cast<Eigen::Index>(k)),1e-12);
            EXPECT_NEAR(boost::json::value_to<double>(fit.at("variances").at(0)),(y-x*reference.beta).squaredNorm()/120,1e-14);
        }
    }
}

TEST(UniqueStencilGridTest, BoundariesFailuresAndStagnationRemainUnqualified)
{
    const auto x{Design()}; const Eigen::Vector4d beta(2,-.4,1,.2);
    auto fit{g::Fit(x,x*beta,beta,.1)};
    EXPECT_FALSE(fit.at("qualified").as_bool()); EXPECT_EQ(fit.at("reason"),"exact-fit-boundary");
    EXPECT_EQ(fit.at("failure_block"),"global");
    auto duplicate=x; duplicate.col(2)=duplicate.col(0);
    EXPECT_EQ(g::Fit(duplicate,Response(x),beta,.1).at("reason"),"rank-deficient");
    EXPECT_EQ(ac::Evaluate(x,Response(x),beta,Eigen::VectorXd::Constant(1,1e-100),g::GlobalBlock(120,1)).reason,"invalid-denominator");
    EXPECT_EQ(ac::Fit(x,Response(x),Eigen::Vector4d(1,0,2,0),g::GlobalBlock(120,.1),0,0).at("reason"),"budget-exhausted");
    EXPECT_EQ(ac::Fit(x,Response(x),beta,g::GlobalBlock(120,.1),100,0).at("reason"),"reference-unverified");
    Eigen::MatrixXd z(24,2); Eigen::VectorXd y(24);
    for (int i=0;i<24;++i) {z(i,0)=1; z(i,1)=std::sin(.71*i); y(i)=1e9+.1*std::sin(2.17*i);}
    fit=g::Fit(z,y,Eigen::Vector2d(1e9,0),0);
    EXPECT_FALSE(fit.at("qualified").as_bool()); EXPECT_EQ(fit.at("reason"),"stalled");
}

TEST(UniqueStencilGridTest, ZeroAmplitudeAndSignedChargesAreAllowed)
{
    Eigen::MatrixXd x(31,2); Eigen::VectorXd error(31);
    for (int i=0;i<31;++i) {const double t{static_cast<double>(i)-15}; x(i,0)=1; x(i,1)=t; error(i)=.01*(t*t-80);}
    for (double c:{-.4,0.,.4})
    {
        const Eigen::VectorXd y{x*Eigen::Vector2d(-1,c)+error};
        const auto fit{g::Fit(x,y,Eigen::Vector2d(1,0),0)};
        ASSERT_TRUE(fit.at("qualified").as_bool()); EXPECT_DOUBLE_EQ(boost::json::value_to<double>(fit.at("beta").at(0)),0);
        EXPECT_NEAR(boost::json::value_to<double>(fit.at("beta").at(1)),c,1e-12);
    }
}

TEST(UniqueStencilGridTest, BlockedSvdAgreesWithPivotedSvdIncludingConstraintsAndRankFailure)
{
    auto x{Design()}; const auto y{Response(x)};
    Eigen::VectorXd weights(120);
    for (int i=0;i<120;++i) weights(i)=i%7==0 ? 0 : .2+std::abs(std::sin(.1*i));
    for (double sign:{-1.,1.})
    {
        const auto old{ac::WeightedSolve(x,sign*y,weights,true)};
        const auto blocked{ac::WeightedSolve(x,sign*y,weights,true,true)};
        ASSERT_TRUE(old.valid); ASSERT_TRUE(blocked.valid); EXPECT_EQ(old.rank,blocked.rank);
        EXPECT_LT((old.beta-blocked.beta).norm(),1e-11);
    }
    const auto old{ac::Fit(x,y,Eigen::Vector4d(2,-.4,1,.2),g::GlobalBlock(120,.1),100,100)};
    const auto blocked{g::Fit(x,y,Eigen::Vector4d(2,-.4,1,.2),.1)};
    EXPECT_EQ(old.at("qualified"),blocked.at("qualified")); EXPECT_EQ(old.at("beta"),blocked.at("beta"));
    EXPECT_EQ(old.at("design_spectrum").at("rank"),blocked.at("design_spectrum").at("rank"));
    EXPECT_NEAR(boost::json::value_to<double>(old.at("design_spectrum").at("condition")),
        boost::json::value_to<double>(blocked.at("design_spectrum").at("condition")),1e-11);
    x.col(2)=x.col(0);
    EXPECT_EQ(ac::WeightedSolve(x,y,weights,true,true).reason,"rank-deficient");
    Eigen::MatrixXd tall(2000,12); Eigen::VectorXd response(2000),w(2000);
    for (int i=0;i<2000;++i)
    {
        for (int k=0;k<12;++k) tall(i,k)=std::sin(.013*(i+1)*(k+1));
        response(i)=2*tall(i,0)-tall(i,2)+.3*tall(i,5)+.01*std::cos(.31*i);
        w(i)=.1+std::abs(std::cos(.017*i));
    }
    const auto direct{ac::WeightedSolve(tall,response,w,true)},reduced{ac::WeightedSolve(tall,response,w,true,true)};
    ASSERT_TRUE(direct.valid); ASSERT_TRUE(reduced.valid); EXPECT_EQ(direct.rank,reduced.rank);
    EXPECT_LT((direct.beta-reduced.beta).norm(),1e-11);
}

TEST(UniqueStencilGridTest, SparseQrPreservesEveryNonzeroConstraintsAndRankFailures)
{
    std::mt19937 generator(17); std::normal_distribution<double> normal; bool released{};
    for (int seed=0;seed<200;++seed)
    {
        Eigen::MatrixXd x(60,6); Eigen::VectorXd y(60),weights(60);
        for (int i=0;i<60;++i)
        {
            y(i)=normal(generator); weights(i)=i%7==0 ? 0 : .1+std::abs(normal(generator));
            for (int k=0;k<6;++k) x(i,k)=i%3==k%3 ? 0 : normal(generator)+.6*std::cos(.2*i);
        }
        x(0,0)=1e-200;
        Eigen::SparseMatrix<double> sparse{x.sparseView(0.0,0.0)};
        EXPECT_EQ(sparse.nonZeros(),(x.array()!=0).count()); EXPECT_EQ(sparse.coeff(0,0),1e-200);
        const auto dense{ac::WeightedSolve(x,y,weights)};
        const auto fast{ac::WeightedSolve(x,y,weights,false,false,&sparse)};
        const auto check{ac::WeightedSolve(x,y,weights,true,true)};
        ASSERT_TRUE(dense.valid); ASSERT_TRUE(fast.valid); ASSERT_TRUE(check.valid);
        EXPECT_LT((dense.beta-fast.beta).norm(),1e-11); EXPECT_LT((check.beta-fast.beta).norm(),1e-11);
        EXPECT_EQ(dense.rank,fast.rank); released|=fast.releases>0;
        x.col(2)=x.col(0); sparse=x.sparseView(0.0,0.0);
        EXPECT_EQ(ac::WeightedSolve(x,y,weights,false,false,&sparse).reason,"rank-deficient");
        EXPECT_EQ(ac::WeightedSolve(x,y,Eigen::VectorXd::Zero(60),false,false,&sparse).reason,"rank-deficient");
    }
    EXPECT_TRUE(released);
}

TEST(UniqueStencilGridTest, SparseQrCommonAlphaAgreesWithDenseAndIndependentSvd)
{
    const auto x{Design()}; const auto y{Response(x)};
    const Eigen::SparseMatrix<double> sparse{x.sparseView(0.0,0.0)};
    for (double alpha:g::alphas)
    {
        const auto dense{g::Fit(x,y,Eigen::Vector4d(2,-.4,1,.2),alpha)};
        const auto fast{g::Fit(x,y,Eigen::Vector4d(2,-.4,1,.2),alpha,&sparse)};
        EXPECT_EQ(dense.at("qualified"),fast.at("qualified")); EXPECT_EQ(dense.at("reason"),fast.at("reason"));
        EXPECT_EQ(fast.at("linear_solver"),"sparse-qr");
        for (std::size_t k=0;k<4;++k) EXPECT_NEAR(boost::json::value_to<double>(dense.at("beta").at(k)),
            boost::json::value_to<double>(fast.at("beta").at(k)),1e-8);
        EXPECT_NEAR(boost::json::value_to<double>(dense.at("variances").at(0)),
            boost::json::value_to<double>(fast.at("variances").at(0)),1e-10);
    }
}

TEST(UniqueStencilGridTest, SparseRowReductionAcrossTilesPreservesWeightedConstrainedSolution)
{
    Eigen::MatrixXd x{Eigen::MatrixXd::Zero(5003,24)}; Eigen::VectorXd y(5003),weights(5003);
    for (int i=0;i<5003;++i)
    {
        for (int k=0;k<24;++k) if ((i/800+k)%4==0 || k<4) x(i,k)=std::sin(.013*(i+1)*(k+1));
        y(i)=2*x(i,0)-x(i,2)+.3*x(i,5)+.01*std::cos(.31*i);
        weights(i)=i%17==0 ? 0 : .1+std::abs(std::cos(.017*i));
    }
    const Eigen::SparseMatrix<double> sparse{x.sparseView(0.0,0.0)};
    const auto dense{ac::WeightedSolve(x,y,weights,true,true)};
    const auto reduced{ac::WeightedSolve(x,y,weights,false,false,&sparse)};
    ASSERT_TRUE(dense.valid); ASSERT_TRUE(reduced.valid); EXPECT_EQ(dense.rank,reduced.rank);
    EXPECT_LT((dense.beta-reduced.beta).norm(),1e-11);
    EXPECT_NEAR((weights.array()*(y-x*dense.beta).array().square()).sum(),
        (weights.array()*(y-x*reduced.beta).array().square()).sum(),1e-10);
}
