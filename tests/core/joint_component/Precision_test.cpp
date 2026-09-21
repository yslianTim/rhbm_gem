#include <gtest/gtest.h>
#include "support/JointOfflineAudit.hpp"
#include "support/JointPrecisionAudit.hpp"
#include "support/InstrumentedLM.hpp"
#include "support/CommandTestHelpers.hpp"
#include "core/command/detail/SimulationGeometry.hpp"
#include <rhbm_gem/utils/math/ElectricPotential.hpp>

namespace {
namespace m=second_stage_test::matched;
namespace p=m::joint_abc;
namespace j=boost::json;
using Vector=Eigen::VectorXd;
struct Sample
{
    m::unique_grid::Grid grid;
    std::vector<m::Atom> atoms{{{-.6,0,0},2,.42,-.3},{{.7,.2,0},1.5,.67,.4}};
    Vector y;
    Sample()
    {
        for(int x=-5;x<=5;++x) for(int z=-2;z<=2;++z) for(int y0=-3;y0<=3;++y0)
        {m::unique_grid::Voxel v; v.position={x*.4,y0*.4,z*.4}; v.index=grid.voxels.size(); grid.voxels.push_back(v);}
        y.resize(static_cast<Eigen::Index>(grid.voxels.size()));
        for(Eigen::Index k=0;k<y.size();++k) y(k)=m::unique_grid::Direct(grid.voxels[static_cast<std::size_t>(k)].position,atoms,2.5);
    }
};
struct RejectOnce
{
    int calls{},rejected{}; bool always{},reject_all{},tiny_step{}; std::string failure;
    bool guarded() const {return true;}
    bool retry() const {return calls<40;}
    int values() const {return 2;}
    int operator()(const Vector & x,Vector & r) {++calls; r=Vector::Constant(2,tiny_step ? -1 : x(0)-1); return 0;}
    int df(const Vector &,Eigen::MatrixXd & d) {d=Eigen::MatrixXd::Constant(2,1,tiny_step ? 1e100 : 1.); return 0;}
    bool Trial(const Vector &,const Vector &,const Vector &,double,double,double,double,double,bool proposed)
    {if(reject_all || (proposed && (always || rejected==0))) {++rejected; return false;} return true;}
};
}

TEST(JointOfflineAuditTest, HighPrecisionReferenceIsIndependentAndStable)
{
    Sample s; const p::Domain domain(s.grid,s.atoms); const Vector eta=Eigen::Vector2d(.48,.59).array().log();
    const auto e=p::Evaluate(domain,s.y,eta); ASSERT_TRUE(e.valid);
    const auto result=m::certification::PrecisionAudit(domain,s.y,e,Eigen::MatrixXd::Identity(2,2));
    EXPECT_TRUE(result.at("agreement_passed").as_bool())<<j::serialize(result);
    EXPECT_TRUE(result.at("derivative_passed").as_bool())<<j::serialize(result);
}

TEST(JointOfflineAuditTest, PrecisionProfileChangeUsesTheSameObservations)
{
    Sample s; const p::Domain domain(s.grid,s.atoms);
    const Vector eta=Eigen::Vector2d(.42,.67).array().log(),original=s.y;
    const auto zero=m::certification::PrecisionProfileChange(domain,s.y,eta,eta);
    const auto shifted=m::certification::PrecisionProfileChange(domain,s.y,eta,eta+Vector::Constant(2,.01));
    for(const char * precision:{"precision50","precision100"})
    {
        ASSERT_TRUE(zero.at(precision).at("valid").as_bool());
        EXPECT_EQ(std::stod(j::value_to<std::string>(zero.at(precision).at("objective_change"))),0.);
        ASSERT_TRUE(shifted.at(precision).at("valid").as_bool());
        EXPECT_GT(std::stod(j::value_to<std::string>(shifted.at(precision).at("objective_change"))),0.);
    }
    EXPECT_EQ((original-s.y).norm(),0.);
}

TEST(JointOfflineAuditTest, BoundaryReductionPreservesFeasibleCompensationScans)
{
    m::unique_grid::Grid grid; std::vector<m::Atom> atoms;
    Vector eta=Vector::Constant(12,std::log(.5)),beta=Vector::Zero(24);
    for(int a=0;a<12;++a)
    {
        atoms.push_back({{5.0*a,0,0},0,.5,.2}); beta(2*a+1)=.2;
        for(int k=0;k<8;++k)
        {m::unique_grid::Voxel v; v.position={5.0*a,.2*k,0}; v.index=grid.voxels.size(); grid.voxels.push_back(v);}
    }
    const p::Domain domain(grid,atoms);
    Vector y(static_cast<Eigen::Index>(grid.voxels.size()));
    for(Eigen::Index k=0;k<y.size();++k) y(k)=m::unique_grid::Direct(grid.voxels[static_cast<std::size_t>(k)].position,atoms,2.5);
    const auto audit=m::certification::BoundaryAudit(domain,y,eta,beta);
    ASSERT_TRUE(audit.at("agreement_passed").as_bool());
    for(const auto & row:audit.at("precision100").at("rows").as_array())
    {
        EXPECT_TRUE(row.at("path_feasible").as_bool());
        EXPECT_GT(j::value_to<int>(row.at("near_zero_rows")),0);
        EXPECT_LT(std::stod(j::value_to<std::string>(row.at("first_order_prediction_norm"))),1e-20);
        ASSERT_TRUE(row.at("profile_valid").as_bool());
        EXPECT_LE(std::stod(j::value_to<std::string>(row.at("profile_objective"))),
            std::stod(j::value_to<std::string>(row.at("path_objective")))+1e-40);
    }
}
