#include <gtest/gtest.h>
#include "core/sampling/MapInterpolation.hpp"
#include <cmath>
#include <functional>
#include <limits>
#include <memory>

namespace {
using namespace rhbm_gem;
using namespace rhbm_gem::core::detail;
using Position = std::array<double,3>;
MapObject Field(const std::array<int,3> & dims, const Position & spacing,
    const Position & origin, const std::function<double(Position)> & field)
{
    auto values = std::make_unique<double[]>(static_cast<std::size_t>(dims[0]*dims[1]*dims[2]));
    std::size_t index = 0;
    for (int z=0;z<dims[2];++z) for (int y=0;y<dims[1];++y) for (int x=0;x<dims[0];++x)
        values[index++] = field({origin[0]+x*spacing[0],origin[1]+y*spacing[1],origin[2]+z*spacing[2]});
    return MapObject{dims,spacing,origin,std::move(values)};
}
double Tricubic(const MapObject & map,const Position & p)
{
    return InterpolateTricubic(MakeTricubicStencil(map,p),
        [&](const auto & node){return map.GetMapValue(node[0],node[1],node[2]);});
}
constexpr std::array kernels{Tricubic,InterpolateTrilinear};
}

TEST(MapInterpolationTest, ConstantAndAffineFieldsWithAnisotropicSpacing)
{
    const auto constant = Field({8,8,8},{.3,.5,.7},{-1.,-2.,-3.},[](Position){return 7.25;});
    const auto affine = [](Position p){return 2.+3.*p[0]-2.*p[1]+.5*p[2];};
    const auto map = Field({8,8,8},{.3,.5,.7},{-1.,-2.,-3.},affine);
    for (auto kernel:kernels) for (Position p: {Position{-.4,-1.,-1.6},Position{-.31,-.72,-.93}})
    {
        EXPECT_NEAR(kernel(constant,p),7.25,1e-14);
        EXPECT_NEAR(kernel(map,p),affine(p),2e-14);
    }
}

TEST(MapInterpolationTest, NonlinearFieldHasIndependentMethodReferences)
{
    const auto square=[](Position p){return p[0]*p[0]+p[1]*p[1]+p[2]*p[2];};
    const auto map=Field({8,8,8},{.3,.5,.7},{-1.,-2.,-3.},square);
    const Position p{-.31,-.72,-.93};
    EXPECT_NEAR(Tricubic(map,p),square(p),2e-14);
    double linear=square(p);
    for (std::size_t k=0;k<3;++k)
    {
        const double g=(p[k]-map.GetOrigin()[k])/map.GetGridSpacing()[k],t=g-std::floor(g);
        linear+=map.GetGridSpacing()[k]*map.GetGridSpacing()[k]*t*(1.-t);
    }
    EXPECT_NEAR(InterpolateTrilinear(map,p),linear,2e-14);
}

TEST(MapInterpolationTest, AnalyticGaussianMatchesNodesAndTrilinearConvexWeights)
{
    const auto gaussian=[](Position p){return std::exp(-.5*(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]));};
    const auto map=Field({9,9,9},{.25,.25,.25},{-1.,-1.,-1.},gaussian);
    for (auto kernel:kernels) for (Position p:{Position{0.,0.,0.},Position{.25,-.5,.75}})
        EXPECT_DOUBLE_EQ(kernel(map,p),gaussian(p));
    const Position p{.1,.15,.2};
    double reference=0.;
    for (int z=0;z<2;++z) for (int y=0;y<2;++y) for (int x=0;x<2;++x)
        reference+=gaussian({x*.25,y*.25,z*.25})*(x?.4:.6)*(y?.6:.4)*(z?.8:.2);
    EXPECT_NEAR(InterpolateTrilinear(map,p),reference,1e-15);
    // Catmull-Rom is approximate on a Gaussian; preserve its numerical reference.
    EXPECT_NEAR(Tricubic(map,p),.96367925767246865,2e-15);
}

TEST(MapInterpolationTest, PreservesDistinctBoundarySemanticsAndSmallGrids)
{
    const auto map=Field({2,2,2},{1.,1.,1.},{0.,0.,0.},[](Position p){return 1.+p[0]+2.*p[1]+4.*p[2];});
    EXPECT_NEAR(Tricubic(map,{.3,.6,.8}),5.898,2e-14);
    EXPECT_NEAR(InterpolateTrilinear(map,{.3,.6,.8}),5.7,2e-14);
    EXPECT_DOUBLE_EQ(Tricubic(map,{-.5,0.,0.}),.9375);
    EXPECT_DOUBLE_EQ(InterpolateTrilinear(map,{-.5,0.,0.}),1.);
    EXPECT_DOUBLE_EQ(InterpolateTrilinear(map,{1.5,0.,0.}),2.);
    EXPECT_THROW(InterpolateTrilinear(map,{-.500001,0.,0.}),std::out_of_range);
    EXPECT_THROW(InterpolateTrilinear(map,{1.500001,0.,0.}),std::out_of_range);
    // MapObject's finite out-of-range sentinel is deliberately retained for cubic.
    EXPECT_NEAR(Tricubic(map,{1.6,0.,0.}),6.408,2e-14);
    const auto tiny=Field({1,1,1},{.2,.3,.4},{1.,2.,3.},[](Position){return 9.;});
    for (auto kernel:kernels)
        EXPECT_DOUBLE_EQ(kernel(tiny,{.9,2.15,3.2}),9.);
}

TEST(MapInterpolationTest, RejectsNonfinitePositionsBeforeIndexConversion)
{
    const auto map=Field({2,2,2},{1.,1.,1.},{0.,0.,0.},[](Position){return 1.;});
    for (auto kernel:kernels) for (std::size_t axis=0;axis<3;++axis)
        for (double bad:{std::numeric_limits<double>::quiet_NaN(),std::numeric_limits<double>::infinity(),-std::numeric_limits<double>::infinity()})
        {
            Position p{}; p[axis]=bad;
            EXPECT_THROW(kernel(map,p),std::out_of_range);
        }
}
