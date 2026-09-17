#include <gtest/gtest.h>
#include "support/MatchedJointAC.hpp"
#include <rhbm_gem/data/object/MapObject.hpp>
#include <cmath>
#include <random>

namespace {
namespace m=second_stage_test::matched;
namespace ac=m::joint_ac;
ac::Blocks Blocks(Eigen::Index n, std::vector<double> alphas)
{
    ac::Blocks blocks;
    for (std::size_t i=0;i<alphas.size();++i) blocks.push_back({i,alphas[i],{}});
    for (Eigen::Index p=0;p<n;++p) blocks[static_cast<std::size_t>(p)%blocks.size()].rows.push_back(p);
    return blocks;
}
Eigen::MatrixXd Design()
{
    Eigen::MatrixXd x(120,4);
    for (int i=0;i<120;++i) {const double t=static_cast<double>(i)/120; x.row(i)<<1+t,std::sin(7*t),std::cos(4*t),t*t;}
    return x;
}
Eigen::VectorXd Response(const Eigen::MatrixXd & x)
{
    Eigen::Vector4d b; b<<2,-.4,1,-.2;
    Eigen::VectorXd y{x*b}; for (Eigen::Index i=0;i<y.size();++i) y(i)+=.05*std::sin(3.1*static_cast<double>(i));
    return y;
}
}

TEST(MatchedJointACTest, ZeroAlphaEqualsConstrainedLeastSquaresAndKeepsChargeSigns)
{
    const auto x{Design()}; const auto y{Response(x)};
    const auto ls{ac::WeightedSolve(x,y,Eigen::VectorXd::Ones(y.size()))}; ASSERT_TRUE(ls.valid);
    const auto fit{ac::Fit(x,y,Eigen::Vector4d(2,-.4,1,-.2),Blocks(120,{0}))};
    ASSERT_TRUE(fit.at("qualified").as_bool())<<fit;
    for (std::size_t i=0;i<4;++i) EXPECT_NEAR(boost::json::value_to<double>(fit.at("beta").at(i)),ls.beta(static_cast<Eigen::Index>(i)),1e-12);
    EXPECT_NEAR(boost::json::value_to<double>(fit.at("variances").at(0)),(y-x*ls.beta).squaredNorm()/120,1e-14);
}

TEST(MatchedJointACTest, MixedActiveSetMatchesExhaustiveSolutionsAndReleasesConstraints)
{
    bool released{};
    std::mt19937 generator(17); std::normal_distribution<double> normal;
    for (int seed=1;seed<=200;++seed)
    {
        Eigen::MatrixXd x(30,6); Eigen::VectorXd y(30);
        for (int i=0;i<30;++i)
        {
            y(i)=normal(generator);
            for (int k=0;k<6;++k) x(i,k)=normal(generator)+.6*std::cos(.2*i);
        }
        const auto solution{ac::WeightedSolve(x,y,Eigen::VectorXd::Ones(30))}; ASSERT_TRUE(solution.valid);
        released |= solution.releases>0;
        double best{std::numeric_limits<double>::infinity()};
        for (int mask=0;mask<8;++mask)
        {
            std::vector<int> cols;
            for (int k=0;k<6;++k) if (k%2 || (mask&(1<<(k/2)))) cols.push_back(k);
            Eigen::MatrixXd subset(30,static_cast<Eigen::Index>(cols.size()));
            for (std::size_t k=0;k<cols.size();++k) subset.col(static_cast<Eigen::Index>(k))=x.col(cols[k]);
            const Eigen::VectorXd b{subset.colPivHouseholderQr().solve(y)};
            bool valid=true;
            for (std::size_t k=0;k<cols.size();++k) if (cols[k]%2==0 && b(static_cast<Eigen::Index>(k))<0) valid=false;
            if (valid) best=std::min(best,(y-subset*b).squaredNorm());
        }
        EXPECT_NEAR((y-x*solution.beta).squaredNorm(),best,1e-10);
    }
    EXPECT_TRUE(released);
}

TEST(MatchedJointACTest, MixedDpdGradientAgreesWithCoefficientAndEveryScaleEquation)
{
    const auto x{Design()}; const auto y{Response(x)}; const auto blocks{Blocks(120,{0,.2,1})};
    const Eigen::Vector4d b(2.01,-.39,1.01,-.19); const Eigen::Vector3d v(.03,.02,.01); const double h{1e-6};
    const auto e{ac::Evaluate(x,y,b,v,blocks)}; ASSERT_TRUE(e.valid);
    const Eigen::VectorXd r{y-x*b}; Eigen::VectorXd weights(120);
    for (std::size_t i=0;i<blocks.size();++i)
    {
        const double a{blocks[i].alpha},vi{v(static_cast<Eigen::Index>(i))};
        const double q{(1+a)*std::pow(2*std::acos(-1.0)*vi,-a/2)/(3*40*vi)};
        for (auto p:blocks[i].rows) weights(p)=q*e.weights(p);
    }
    for (int k=0;k<4;++k)
    {
        Eigen::Vector4d plus=b,minus=b; plus(k)+=h; minus(k)-=h;
        const double numeric{(ac::Objective(y-x*plus,v,blocks)-ac::Objective(y-x*minus,v,blocks))/(2*h)};
        EXPECT_NEAR(numeric,-x.col(k).dot((weights.array()*r.array()).matrix()),1e-6);
    }
    for (int i=0;i<3;++i)
    {
        Eigen::Vector3d plus=v,minus=v; plus(i)*=std::exp(h); minus(i)*=std::exp(-h);
        const double a{blocks[static_cast<std::size_t>(i)].alpha};
        const double numeric{(ac::Objective(r,plus,blocks)-ac::Objective(r,minus,blocks))/(2*h)};
        EXPECT_NEAR(numeric,-(1+a)*std::pow(2*std::acos(-1.0)*v(i),-a/2)/6*e.scaled(4+i),1e-6);
    }
    EXPECT_GT((weights/weights.maxCoeff()-e.weights/e.weights.maxCoeff()).norm(),.1);
    EXPECT_NEAR((weights/weights.maxCoeff()-e.linear_weights/e.linear_weights.maxCoeff()).norm(),0,1e-14);
}

TEST(MatchedJointACTest, EqualBlocksReduceToCommonObjectiveAndZeroAlphaStillUsesScale)
{
    const auto x{Design()}; const auto y{Response(x)}; const Eigen::Vector4d b(2,-.4,1,-.2);
    const Eigen::VectorXd r{y-x*b}; const double v{.02},a{.2};
    const double old{std::pow(2*std::acos(-1.0)*v,-a/2)*(1/std::sqrt(1+a)-
        (1+a)/a*(-.5*a*r.array().square()/v).exp().mean())};
    EXPECT_NEAR(ac::Objective(r,Eigen::Vector3d::Constant(v),Blocks(120,{a,a,a})),old+1/a,1e-13);
    const auto e{ac::Evaluate(x,y,b,Eigen::Vector3d(.01,.02,.04),Blocks(120,{0,0,0}))}; ASSERT_TRUE(e.valid);
    EXPECT_TRUE((e.weights.array()==1).all());
    EXPECT_NEAR(e.linear_weights(0)/e.linear_weights(1),2,1e-14);
    EXPECT_NEAR(e.linear_weights(0)/e.linear_weights(2),4,1e-14);
}

TEST(MatchedJointACTest, PositiveAlphaHasFreshQualifiedEndpointAndIndependentReference)
{
    const auto x{Design()}; const auto y{Response(x)};
    const auto fit{ac::Fit(x,y,Eigen::Vector4d(2,-.4,1,-.2),Blocks(120,{0,.1,.2}))};
    ASSERT_TRUE(fit.at("qualified").as_bool())<<fit;
    EXPECT_LE(boost::json::value_to<double>(fit.at("stationarity")),1e-8);
    EXPECT_EQ(fit.at("branches").as_array().size(),2u);
    for (const auto & branch:fit.at("branches").as_array())
    {
        EXPECT_TRUE(branch.at("qualified").as_bool());
        EXPECT_LE(boost::json::value_to<double>(branch.at("reference").at("stationarity")),1e-10);
    }
}

TEST(MatchedJointACTest, BoundaryAmplitudeAndNegativeResponsesAreLegitimate)
{
    Eigen::MatrixXd x(30,2); x.col(0).setOnes();
    for (int i=0;i<30;++i) x(i,1)=static_cast<double>(i)-15;
    const Eigen::VectorXd y{-Eigen::VectorXd::Ones(30)-.4*x.col(1)};
    const auto fit{ac::Fit(x,y,Eigen::Vector2d(1,-.4),Blocks(30,{0}))};
    ASSERT_TRUE(fit.at("qualified").as_bool())<<fit;
    EXPECT_DOUBLE_EQ(boost::json::value_to<double>(fit.at("beta").at(0)),0);
    EXPECT_LT(boost::json::value_to<double>(fit.at("beta").at(1)),0);
}

TEST(MatchedJointACTest, ExactFitRankFailureDenominatorAndBudgetCannotPass)
{
    const auto x{Design()}; const Eigen::Vector4d b(2,-.4,1,-.2);
    const auto exact{ac::Fit(x,x*b,b,Blocks(120,{.1,.2}))};
    EXPECT_FALSE(exact.at("qualified").as_bool()); EXPECT_EQ(exact.at("reason"),"exact-fit-boundary");
    auto duplicate=x; duplicate.col(2)=duplicate.col(0);
    EXPECT_EQ(ac::Fit(duplicate,Response(x),b,Blocks(120,{.1})).at("reason"),"rank-deficient");
    EXPECT_EQ(ac::Evaluate(x,Response(x),b,Eigen::VectorXd::Constant(1,1e-100),Blocks(120,{.5})).reason,"invalid-denominator");
    const auto exhausted{ac::Fit(x,Response(x),Eigen::Vector4d(1,0,2,0),Blocks(120,{.1}),0,0)};
    EXPECT_FALSE(exhausted.at("qualified").as_bool());
}

TEST(MatchedJointACTest, StructuralComponentsAndDesignPreserveCompleteStencilContributions)
{
    rhbm_gem::MapObject grid({101,61,61},{.1,.1,.1},{-3,-3,-3});
    const std::vector<m::Atom> atoms{{{0,0,0},6,.5,-.2},{{4,0,0},7,.6,.1}};
    const std::vector<m::Stencil> stencils{m::MakeStencil(grid,grid,{2.51,.03,.02}),m::MakeStencil(grid,grid,{0,0,0})};
    const auto c{ac::BuildComponents(stencils,atoms,2.5,{0,0})};
    ASSERT_EQ(c.atoms.size(),1u); EXPECT_EQ(c.atoms[0].size(),2u);
    const auto x{ac::BuildDesign(stencils,atoms,c.rows[0],c.atoms[0],2.5)};
    const Eigen::Vector4d b(6,-.2,7,.1);
    for (int p=0;p<2;++p) EXPECT_NEAR((x*b)(p),m::Predict(stencils[static_cast<std::size_t>(p)],atoms,2.5),1e-14);
    auto changed=atoms; changed[0].amplitude=0; changed[0].charge=0;
    EXPECT_EQ(ac::BuildComponents(stencils,changed,2.5,{0,0}).contributors,c.contributors);
    EXPECT_EQ(ac::BuildDesign(stencils,changed,c.rows[0],c.atoms[0],2.5),x);
    EXPECT_DOUBLE_EQ(atoms[0].width,.5); EXPECT_DOUBLE_EQ(atoms[1].width,.6);
}

TEST(MatchedJointACTest, FrozenFitsUseAllRowsAndDoNotDependOnTargetOrder)
{
    const auto x{Design()}; const auto y{Response(x)}; const Eigen::Vector4d initial(2.1,-.3,.9,-.1);
    const Eigen::VectorXd predicted{x*initial};
    std::array<boost::json::object,2> first;
    for (int k=0;k<2;++k)
    {
        const Eigen::MatrixXd local{x.middleCols(2*k,2)};
        first[static_cast<std::size_t>(k)]=ac::Fit(local,y-predicted+local*initial.segment(2*k,2),initial.segment(2*k,2),Blocks(120,{0,.1,.2}));
    }
    for (int k=1;k>=0;--k)
    {
        const Eigen::MatrixXd local{x.middleCols(2*k,2)};
        auto repeated=ac::Fit(local,y-predicted+local*initial.segment(2*k,2),initial.segment(2*k,2),Blocks(120,{0,.1,.2}));
        repeated.erase("seconds"); first[static_cast<std::size_t>(k)].erase("seconds");
        EXPECT_EQ(repeated,first[static_cast<std::size_t>(k)]);
        EXPECT_EQ(repeated.at("rows"),120);
    }
}

TEST(MatchedJointACTest, SharedOwnerMergesComponentsWithoutChangingStructuralContributors)
{
    rhbm_gem::MapObject grid({161,61,61},{.1,.1,.1},{-3,-3,-3});
    const std::vector<m::Atom> atoms{{{0,0,0},6,.5,-.2},{{8,0,0},7,.6,.1}};
    const std::vector<m::Stencil> stencils{m::MakeStencil(grid,grid,{0,0,0}),m::MakeStencil(grid,grid,{8,0,0})};
    EXPECT_EQ(ac::BuildComponents(stencils,atoms,2.5,{0,1}).atoms.size(),2u);
    const auto linked{ac::BuildComponents(stencils,atoms,2.5,{0,0})};
    EXPECT_EQ(linked.atoms.size(),1u); EXPECT_EQ(linked.contributors[0],std::vector<std::size_t>{0});
    EXPECT_EQ(linked.contributors[1],std::vector<std::size_t>{1});
}

TEST(MatchedJointACTest, InvalidBlockOwnershipAndPartialExactFitCannotPass)
{
    const auto x{Design()}; const Eigen::Vector4d b(2,-.4,1,-.2); auto y{Response(x)};
    auto blocks{Blocks(120,{0,.2})};
    blocks[0].rows.push_back(blocks[1].rows[0]);
    EXPECT_THROW(ac::Fit(x,y,b,blocks),std::invalid_argument);
    blocks=Blocks(120,{0,-.2}); EXPECT_THROW(ac::Fit(x,y,b,blocks),std::invalid_argument);
    blocks=Blocks(120,{0,.2}); for (auto p:blocks[1].rows) y(p)=x.row(p).dot(b);
    const auto e{ac::Evaluate(x,y,b,Eigen::Vector2d(.1,.1),blocks)};
    EXPECT_FALSE(e.valid); EXPECT_EQ(e.reason,"exact-fit-boundary"); EXPECT_EQ(e.failure_owner,1);
}

TEST(MatchedJointACTest, SmallAlphaLimitNonfiniteAndUnderflowAreExplicit)
{
    const auto x{Design()}; auto y{Response(x)}; const Eigen::Vector4d b(2,-.4,1,-.2);
    const auto zero{Blocks(120,{0})}; const Eigen::VectorXd v{Eigen::VectorXd::Constant(1,.02)};
    EXPECT_NEAR(ac::Objective(y-x*b,v,Blocks(120,{1e-10})),ac::Objective(y-x*b,v,zero),1e-8);
    auto e{ac::Evaluate(x,y,b,Eigen::VectorXd::Constant(1,1e-100),Blocks(120,{1}))};
    EXPECT_FALSE(e.valid); EXPECT_EQ(e.reason,"invalid-denominator");
    EXPECT_GT((e.weights.array()==0).count(),0); EXPECT_EQ(e.failure_owner,0);
    y(0)=std::numeric_limits<double>::quiet_NaN();
    EXPECT_EQ(ac::Evaluate(x,y,b,v,zero).reason,"nonfinite");
}

TEST(MatchedJointACTest, ZeroTargetColumnsStillContributeToBlockScale)
{
    Eigen::MatrixXd x(40,2); x.setZero();
    for (int i=0;i<20;++i) {x(i,0)=1; x(i,1)=static_cast<double>(i)-10;}
    Eigen::VectorXd y=x*Eigen::Vector2d(2,-.2);
    for (int i=0;i<40;++i) y(i)+=.03*std::sin(1.7*i);
    const auto fit{ac::Fit(x,y,Eigen::Vector2d(2,-.2),Blocks(40,{0,.1}))};
    ASSERT_TRUE(fit.at("qualified").as_bool())<<fit;
    EXPECT_EQ(fit.at("rows"),40); EXPECT_EQ(fit.at("blocks").at(0).at("rows"),20);
    EXPECT_GT(boost::json::value_to<double>(fit.at("variances").at(1)),0);
}

TEST(MatchedJointACTest, FloatingPointStagnationCannotQualifyAsStationarity)
{
    Eigen::MatrixXd x(24,2); Eigen::VectorXd y(24);
    for (int i=0;i<24;++i)
    {
        x(i,0)=1; x(i,1)=std::sin(.71*i);
        y(i)=1e9+.1*std::sin(2.17*i);
    }
    const auto fit{ac::Fit(x,y,Eigen::Vector2d(1e9,0),Blocks(24,{0}),20,20)};
    EXPECT_FALSE(fit.at("qualified").as_bool());
    EXPECT_EQ(fit.at("reason"),"stalled");
    EXPECT_GT(boost::json::value_to<double>(fit.at("stationarity")),1e-8);
    for (const auto & branch:fit.at("branches").as_array())
    {
        EXPECT_FALSE(branch.at("qualified").as_bool());
        EXPECT_EQ(branch.at("primary").at("stop"),"stalled");
    }
}

TEST(MatchedJointACTest, PositiveNegativeAndZeroChargeRemainUnconstrained)
{
    Eigen::MatrixXd x(31,2); Eigen::VectorXd error(31);
    for (int i=0;i<31;++i)
    {
        const double t{static_cast<double>(i)-15};
        x(i,0)=1; x(i,1)=t; error(i)=.01*(t*t-80);
    }
    for (double charge:{-.4,0.,.4})
    {
        const Eigen::VectorXd y{x*Eigen::Vector2d(2,charge)+error};
        const auto fit{ac::Fit(x,y,Eigen::Vector2d(1,.2),Blocks(31,{0}))};
        ASSERT_TRUE(fit.at("qualified").as_bool())<<fit;
        EXPECT_NEAR(boost::json::value_to<double>(fit.at("beta").at(0)),2,1e-12);
        EXPECT_NEAR(boost::json::value_to<double>(fit.at("beta").at(1)),charge,1e-12);
    }
}
