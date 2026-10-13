#include <gtest/gtest.h>
#include "core/detail/joint_component/CompactSvd.hpp"
#include "core/detail/joint_component/Numerics.hpp"
#include <Eigen/SVD>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

namespace {
namespace n=rhbm_gem::core::joint_component;
using Matrix=Eigen::MatrixXd;
using Vector=Eigen::VectorXd;

struct RankFixture
{
    std::string name;
    Matrix design;
    std::vector<std::string> column_roles;
};

std::vector<RankFixture> ReferenceFixtures()
{
    Matrix ill(2,2); ill<<1,1,0,1e-7;
    Matrix duplicate(3,2); duplicate<<1,1,0,0,1,1;
    Matrix zero_column=Matrix::Zero(4,3); zero_column(0,0)=1; zero_column(1,1)=1;
    Matrix structural(3,3); structural<<1,1,0,0,0,1,0,0,0;
    Matrix numeric_deficient(3,2); numeric_deficient<<1,2,0,0,1,2;
    Matrix coupled(3,2); coupled<<1,1,0,1e-10,0,0;
    Matrix connected_full=Matrix::Zero(5,5);
    for(Eigen::Index k=0;k<5;++k)
    {
        connected_full(k,k)=1;
        if(k) connected_full(k,k-1)=.25;
    }
    Matrix connected_deficient(3,3); connected_deficient<<1,0,1,1,1,2,0,1,1;
    Matrix weak(3,2); weak<<1,1,0,1e-12,0,0;
    Matrix active_boundary(4,4); active_boundary<<0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,1;
    Matrix mixed(4,3); mixed<<1,0,1,0,1,0,0,0,1,0,0,0;
    return {
        {"full-rank-well-conditioned",Matrix::Identity(4,4),{}},
        {"full-rank-ill-conditioned",ill,{}},
        {"exact-duplicate-columns",duplicate,{}},
        {"zero-support-column",zero_column,{}},
        {"structural-rank-deficiency",structural,{}},
        {"structurally-full-numerically-deficient",numeric_deficient,{}},
        {"near-threshold",weak,{}},
        {"active-amplitude-boundary",active_boundary,{"active-A","free-C","free-A","free-C"}},
        {"mixed-target-halo",mixed,{"target","halo","target"}},
        {"strongly-coupled-A-C",coupled,{"A","C"}},
        {"connected-full-rank",connected_full,{}},
        {"connected-rank-deficient",connected_deficient,{}},
        {"weak-signal-cancellation-sensitive",weak,{"A","C"}},
    };
}

Matrix NormalizeColumns(Matrix x,Vector & norms)
{
    norms.resize(x.cols());
    for(Eigen::Index column=0;column<x.cols();++column)
    {
        norms(column)=x.col(column).norm();
        if(norms(column)>0) x.col(column)/=norms(column);
    }
    return x;
}

TEST(JointComponentBoundedRankTest, DenseSvdOracleCoversReferenceControls)
{
    for(const auto & fixture:ReferenceFixtures())
    {
        SCOPED_TRACE(fixture.name);
        const n::RankPolicy policy{fixture.design.rows(),fixture.design.cols(),0};
        Vector norms;
        const Matrix normalized=NormalizeColumns(fixture.design,norms);
        const auto reference=n::ComputeSpectrum(fixture.design,policy,fixture.design.cols(),true);
        ASSERT_TRUE(reference.available)<<reference.reason;
        Eigen::JacobiSVD<Matrix> svd(normalized,Eigen::ComputeThinU|Eigen::ComputeThinV);
        ASSERT_EQ(reference.singular_values.size(),svd.singularValues().size());
        EXPECT_LE((reference.singular_values-svd.singularValues()).lpNorm<Eigen::Infinity>(),
            2e-12*std::max(1.0,svd.singularValues()(0)));
        EXPECT_LE((reference.column_norms-norms).lpNorm<Eigen::Infinity>(),1e-14);
        const double expected_threshold=policy.Relative(fixture.design.cols())*svd.singularValues()(0);
        EXPECT_DOUBLE_EQ(reference.threshold,expected_threshold);
        EXPECT_EQ(reference.rank,(reference.singular_values.array()>reference.threshold).count());
        EXPECT_EQ(reference.rows,policy.rows);
        EXPECT_EQ(reference.columns,fixture.design.cols());
    }
}

TEST(JointComponentBoundedRankTest, FullDesignAndFreeFaceUseDifferentColumnDomains)
{
    const auto fixtures=ReferenceFixtures();
    const auto found=std::find_if(fixtures.begin(),fixtures.end(),[](const auto & f) {
        return f.name=="active-amplitude-boundary";
    });
    ASSERT_NE(found,fixtures.end());
    const auto full=n::ComputeSpectrum(found->design,{100,4,2},4,true);
    Matrix free_face=found->design(Eigen::all,Eigen::seq(1,3));
    const auto free=n::ComputeSpectrum(free_face,{100,4,2},3,true);
    ASSERT_TRUE(full.available && free.available);
    EXPECT_EQ(full.columns,4);
    EXPECT_EQ(free.columns,3);
    EXPECT_EQ(full.rank,3);
    EXPECT_EQ(free.rank,3);
}

TEST(JointComponentBoundedRankTest, StrictRankBoundaryUsesOriginalRowsAndStrictGreater)
{
    constexpr Eigen::Index original_rows=1000;
    const n::RankPolicy policy{original_rows,2,0};
    const double threshold=policy.Absolute(2,1.0);
    for(const auto & item:std::vector<std::pair<const char *,double>>{
        {"below",std::nextafter(threshold,0.0)},
        {"equal",threshold},
        {"above",std::nextafter(threshold,std::numeric_limits<double>::infinity())}})
    {
        SCOPED_TRACE(item.first);
        Matrix x=Matrix::Zero(2,2); x(0,0)=1; x(1,1)=item.second;
        const auto spectrum=n::ComputeSpectrum(x,policy,2,false);
        ASSERT_TRUE(spectrum.available)<<spectrum.reason;
        EXPECT_DOUBLE_EQ(spectrum.threshold,threshold);
        EXPECT_EQ(spectrum.rows,original_rows);
        EXPECT_EQ(spectrum.rank,item.second>threshold ? 2 : 1);
    }
}

} // namespace
