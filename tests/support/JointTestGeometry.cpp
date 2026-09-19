#include "support/JointTestGeometry.hpp"
#include "core/detail/joint_component/Numerics.hpp"
#include <rhbm_gem/data/object/MapObject.hpp>
#include <algorithm>
#include <cmath>
#include <numbers>
namespace second_stage_test::matched {
namespace {
std::array<double, 4> Coefficients(double t)
{
    const double t2{t*t}, t3{t2*t};
    return {-0.5*t3+t2-0.5*t, 1.5*t3-2.5*t2+1.0,
        -1.5*t3+2.0*t2+0.5*t, 0.5*t3-0.5*t2};
}
}
double SquareDistance(const Position & a, const Position & b)
{
    double out{}; for (std::size_t k=0;k<3;++k) out += (a[k]-b[k])*(a[k]-b[k]); return out;
}
Stencil MakeStencil(const rhbm_gem::MapObject & generation,
    const rhbm_gem::MapObject & sampling, const Position & p)
{
    if (generation.GetGridSize() != sampling.GetGridSize()) throw std::invalid_argument("Mismatched grid dimensions.");
    const auto index{sampling.GetIndexFromPosition(p)}, size{sampling.GetGridSize()};
    // Preserve the sampler's {-1,-1,-1} sentinel and subsequent clamping too:
    // float32 header coordinates can put an original boundary point just outside.
    const auto origin{sampling.GetOrigin()}, h{sampling.GetGridSpacing()};
    std::array<std::array<double,4>,3> c;
    Stencil out;
    for (std::size_t k=0;k<3;++k)
    {
        c[k] = Coefficients((p[k]-origin[k]-static_cast<double>(index[k])*h[k])/h[k]);
        out.boundary |= index[k]-1<0 || index[k]+2>=size[k];
    }
    std::size_t slot{};
    for (int x=0;x<4;++x) for (int y=0;y<4;++y) for (int z=0;z<4;++z)
    {
        const int ix{std::clamp(index[0]+x-1,0,size[0]-1)}, iy{std::clamp(index[1]+y-1,0,size[1]-1)},
            iz{std::clamp(index[2]+z-1,0,size[2]-1)};
        const auto global{static_cast<std::size_t>(ix)+static_cast<std::size_t>(size[0])*
            (static_cast<std::size_t>(iy)+static_cast<std::size_t>(size[1])*static_cast<std::size_t>(iz))};
        out.slots[slot++] = {global,generation.GetGridPosition(global),
            c[0][static_cast<std::size_t>(x)]*c[1][static_cast<std::size_t>(y)]*c[2][static_cast<std::size_t>(z)]};
    }
    return out;
}
Basis EvaluateBasis(double square,double width,double cutoff)
{
    const auto b=rhbm_gem::core::joint_component::EvaluateKernel(square,width,cutoff);
    return {b.gaussian,b.charge,b.gaussian_log_width,b.charge_log_width};
}
Design MakeDesign(const std::vector<Stencil> & stencils, const Position & center)
{
    Design out;
    for (const auto & stencil : stencils)
    {
        auto & row{out.emplace_back()}; row.reserve(64);
        for (const auto & slot : stencil.slots) row.push_back({SquareDistance(slot.position,center),slot.coefficient});
    }
    return out;
}
Eigen::MatrixXd EvaluateDesign(const Design & design, double width, double cutoff)
{
    Eigen::MatrixXd out{Eigen::MatrixXd::Zero(static_cast<Eigen::Index>(design.size()),4)};
    for (std::size_t i=0;i<design.size();++i) for (const auto & term : design[i])
    {
        const auto b{EvaluateBasis(term.square,width,cutoff)};
        out.row(static_cast<Eigen::Index>(i)) += term.coefficient*Eigen::RowVector4d{
            b.gaussian,b.charge,b.gaussian_log_width,b.charge_log_width};
    }
    return out;
}
double Predict(const Stencil & stencil, const std::vector<Atom> & atoms, double cutoff,
    bool quantized, double * absolute_sum, double * quantization_bound)
{
    double out{}, absolute{}, bound{};
    for (const auto & slot : stencil.slots)
    {
        double value{}, magnitudes{};
        for (const auto & atom : atoms)
        {
            const auto basis{EvaluateBasis(SquareDistance(slot.position,atom.position),atom.width,cutoff)};
            const double contribution{atom.amplitude*basis.gaussian+atom.charge*basis.charge};
            value += contribution; magnitudes += std::abs(contribution);
        }
        const double rounded{static_cast<double>(static_cast<float>(value))};
        out += slot.coefficient*(quantized ? rounded : value);
        absolute += std::abs(slot.coefficient)*magnitudes;
        bound += std::abs(slot.coefficient)*std::abs(rounded-value);
    }
    if (absolute_sum) *absolute_sum = absolute;
    if (quantization_bound) *quantization_bound = bound;
    return out;
}

namespace unique_grid {
// Independent scalar forward control; do not call the fitted kernel here.
double Direct(const Position & p,const std::vector<Atom> & atoms,double cutoff)
{
    double value{};
    for(const auto & atom:atoms)
    {
        const double square=SquareDistance(p,atom.position),b=atom.width;
        if(square>cutoff*cutoff) continue;
        const double gaussian=std::exp(-square/(2*b*b))/std::pow(2*std::numbers::pi*b*b,1.5);
        const double charge=square<1e-10 ? std::sqrt(2/std::numbers::pi)/b :
            std::sqrt(square)<=2.5 ? std::erf(std::sqrt(square)/(std::sqrt(2.)*b))/std::sqrt(square) : 0.;
        value+=atom.amplitude*gaussian+atom.charge*charge;
    }
    return value;
}
}
} // namespace second_stage_test::matched
