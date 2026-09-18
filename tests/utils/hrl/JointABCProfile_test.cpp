#include <gtest/gtest.h>
#include "support/JointABCProfile.hpp"
#include "support/FixedBOracle.hpp"
#include <unsupported/Eigen/NonLinearOptimization>
#include <cmath>

namespace {
namespace m=second_stage_test::matched;
namespace p=m::joint_abc;
using Vector=Eigen::VectorXd;
using Matrix=Eigen::MatrixXd;
struct Sample
{
    m::unique_grid::Grid grid;
    std::vector<m::Atom> atoms{{{-.6,0,0},2,.42,-.3},{{.7,.2,0},1.5,.67,.4}};
    Vector y;
    Sample()
    {
        for (int x=-6;x<=6;++x) for (int z=-3;z<=3;++z) for (int y0=-4;y0<=4;++y0)
        {
            m::unique_grid::Voxel voxel; voxel.position={x*.4,y0*.4,z*.4}; voxel.index=grid.voxels.size();
            grid.voxels.push_back(voxel);
        }
        y.resize(static_cast<Eigen::Index>(grid.voxels.size()));
        for (Eigen::Index k=0;k<y.size();++k) y(k)=m::unique_grid::Direct(grid.voxels[static_cast<std::size_t>(k)].position,atoms,2.5);
    }
};
Vector Read(const boost::json::value & v)
{
    Vector out(static_cast<Eigen::Index>(v.as_array().size()));
    for (Eigen::Index k=0;k<out.size();++k) out(k)=boost::json::value_to<double>(v.at(static_cast<std::size_t>(k)));
    return out;
}
// Independent small direct joint solve, with numerical derivatives of A/C/log-B.
// This comparison uses an interior solution; boundary behavior is tested separately.
struct Direct
{
    const Sample & sample;
    int values() const {return static_cast<int>(sample.y.size());}
    int operator()(const Vector & x,Vector & residual) const
    {
        auto atoms=sample.atoms;
        for (std::size_t a=0;a<atoms.size();++a)
        {atoms[a].amplitude=x(static_cast<Eigen::Index>(3*a)); atoms[a].charge=x(static_cast<Eigen::Index>(3*a+1)); atoms[a].width=std::exp(x(static_cast<Eigen::Index>(3*a+2)));}
        residual.resize(sample.y.size());
        for (Eigen::Index k=0;k<residual.size();++k)
            residual(k)=m::unique_grid::Direct(sample.grid.voxels[static_cast<std::size_t>(k)].position,atoms,2.5)-sample.y(k);
        return 0;
    }
    int df(const Vector & x,Matrix & jacobian) const
    {
        jacobian.resize(sample.y.size(),x.size());
        for (Eigen::Index k=0;k<x.size();++k)
        {
            Vector plus=x,minus=x,rp,rm; plus(k)+=1e-5; minus(k)-=1e-5;
            (*this)(plus,rp); (*this)(minus,rm); jacobian.col(k)=(rp-rm)/2e-5;
        }
        return 0;
    }
};
}

TEST(JointABCProfileTest, LogWidthBasisIncludesCenterNearZeroAndFixedCutoff)
{
    for (double r:{0.,.5e-5,1e-5,1.1e-5,.4,2.5,2.50000001})
    {
        const double b=.53,h=1e-5;
        const auto value=m::EvaluateBasis(r*r,b,2.5),plus=m::EvaluateBasis(r*r,b*std::exp(h),2.5),minus=m::EvaluateBasis(r*r,b*std::exp(-h),2.5);
        EXPECT_NEAR(value.gaussian_log_width,(plus.gaussian-minus.gaussian)/(2*h),1e-8);
        EXPECT_NEAR(value.charge_log_width,(plus.charge-minus.charge)/(2*h),1e-8);
        if (r>2.5) EXPECT_EQ(value.gaussian+value.charge+value.gaussian_log_width+value.charge_log_width,0);
    }
}

TEST(JointABCProfileTest, FullJacobianAndEnvelopeGradientAtNonzeroResidual)
{
    Sample sample; const p::Domain domain(sample.grid,sample.atoms);
    for (Eigen::Index k=0;k<sample.y.size();++k) sample.y(k)+=.03*std::sin(static_cast<double>(k));
    const Vector eta=Eigen::Vector2d(.55,.51).array().log();
    const auto e=p::Evaluate(domain,sample.y,eta); ASSERT_TRUE(e.valid);
    const double scale=std::max(1.0,sample.y.norm()); const auto d=p::Differentiate(e,scale);
    ASSERT_TRUE(d.valid); EXPECT_GT((d.jacobian-d.projected).norm(),1e-3);
    for (Eigen::Index k=0;k<2;++k)
    {
        Vector plus=eta,minus=eta; plus(k)+=1e-5; minus(k)-=1e-5;
        const auto ep=p::Evaluate(domain,sample.y,plus,true),em=p::Evaluate(domain,sample.y,minus,true);
        ASSERT_TRUE(ep.valid && em.valid);
        const Vector finite=(ep.residual-em.residual)/(2e-5*scale);
        EXPECT_LT((finite-d.jacobian.col(k)).norm()/finite.norm(),1e-7);
        EXPECT_NEAR(e.gradient(k),(ep.residual.squaredNorm()-em.residual.squaredNorm())/(4e-5*scale*scale),1e-9);
    }
}

TEST(JointABCProfileTest, RecoveryAgreesWithDirectJointAndFixedB)
{
    Sample sample; const p::Domain domain(sample.grid,sample.atoms);
    const Vector widths=Eigen::Vector2d(.48,.59);
    const auto fit=p::Fit(domain,sample.y,widths);
    ASSERT_TRUE(fit.at("joint_qualified").as_bool())<<boost::json::serialize(fit);
    const auto beta=Read(fit.at("primary").at("beta")),b=Read(fit.at("primary").at("b"));
    EXPECT_LT((beta-Eigen::Vector4d(2,-.3,1.5,.4)).norm(),1e-10);
    EXPECT_LT((b-Eigen::Vector2d(.42,.67)).norm(),1e-10);
    EXPECT_EQ(fit.at("directional_evaluations"),12);
    const auto e=p::Evaluate(domain,sample.y,widths.array().log());
    const auto fixed=m::fixed_b::Fit(e.x,sample.y,m::joint_ac::SparseSpectrum(e.x,Vector::Ones(sample.y.size())));
    EXPECT_EQ((e.beta-Read(fixed.at("primary").at("beta"))).norm(),0);
    Direct direct{sample}; Eigen::LevenbergMarquardt<Direct> lm(direct);
    Vector joint(6); joint<<e.beta(0),e.beta(1),std::log(widths(0)),e.beta(2),e.beta(3),std::log(widths(1));
    lm.parameters.ftol=1e-14; lm.parameters.xtol=1e-12; lm.parameters.gtol=1e-12; lm.minimize(joint);
    EXPECT_NEAR(joint(0),beta(0),1e-10); EXPECT_NEAR(joint(3),beta(2),1e-10);
    EXPECT_NEAR(std::exp(joint(2)),b(0),1e-10); EXPECT_NEAR(std::exp(joint(5)),b(1),1e-10);
    EXPECT_GT(std::abs(b(0)-b(1)),.1); // Separate widths, not a shared B.
}

TEST(JointABCProfileTest, BoundaryAndChangingActiveFaceRebuildDerivative)
{
    Sample sample; sample.atoms.resize(1); sample.atoms[0].amplitude=.05; sample.atoms[0].charge=-1;
    for (Eigen::Index k=0;k<sample.y.size();++k) sample.y(k)=m::unique_grid::Direct(sample.grid.voxels[static_cast<std::size_t>(k)].position,sample.atoms,2.5);
    const p::Domain domain(sample.grid,sample.atoms); bool active=false,interior=false;
    for (double width:{.25,.42,.8})
    {
        const Vector eta=Vector::Constant(1,std::log(width)); const auto e=p::Evaluate(domain,sample.y,eta);
        ASSERT_TRUE(e.valid); active |= e.beta(0)==0; interior |= e.beta(0)>0;
        const double scale=std::max(1.0,sample.y.norm()); const auto d=p::Differentiate(e,scale); ASSERT_TRUE(d.valid);
        const auto plus=p::Evaluate(domain,sample.y,eta.array()+1e-5,true),minus=p::Evaluate(domain,sample.y,eta.array()-1e-5,true);
        ASSERT_EQ(plus.certificate.at("active_atoms"),minus.certificate.at("active_atoms"));
        EXPECT_LT(((plus.residual-minus.residual)/(2e-5*scale)-d.jacobian.col(0)).norm(),1e-8);
    }
    EXPECT_TRUE(active && interior);
}

TEST(JointABCProfileTest, RankFailureInvalidBAndUninformativeWidthCannotQualify)
{
    Sample sample; auto duplicate=sample.atoms; duplicate[1].position=duplicate[0].position;
    const p::Domain bad(sample.grid,duplicate),good(sample.grid,sample.atoms);
    const Vector eta=Vector::Constant(2,std::log(.5));
    const auto invalid=p::Evaluate(bad,sample.y,eta); EXPECT_FALSE(invalid.valid);
    EXPECT_FALSE(p::Differentiate(invalid,1).valid);
    EXPECT_FALSE(p::Fit(bad,sample.y,Vector::Constant(2,.5)).at("joint_qualified").as_bool());
    EXPECT_FALSE(p::Evaluate(good,sample.y,Vector::Constant(2,1000)).valid);
    EXPECT_FALSE(p::Evaluate(good,sample.y,Vector::Constant(2,-1000)).valid);
    const auto zero=p::Fit(good,Vector::Zero(sample.y.size()),Vector::Constant(2,.5));
    EXPECT_FALSE(zero.at("joint_qualified").as_bool());
    EXPECT_TRUE(zero.at("primary").at("kkt_passed").as_bool());
    EXPECT_EQ(zero.at("primary").at("rss"),0.0);
    EXPECT_EQ(zero.at("width_spectrum").at("rank"),0);
    // Geometry and observations are immutable across all evaluations.
    EXPECT_EQ(good.rows,sample.y.size()); EXPECT_EQ(good.atoms[0].size(),p::Domain(sample.grid,sample.atoms).atoms[0].size());
}
