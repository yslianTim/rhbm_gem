#include "support/ObservationMatchedExperiment.hpp"
#include <rhbm_gem/data/object/MapObject.hpp>
#include <boost/math/tools/minima.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <stdexcept>

namespace second_stage_test::matched {
namespace {
constexpr double pi{std::numbers::pi};
constexpr double lower_width{0.1}, upper_width{2.0};
std::array<double, 4> Coefficients(double t)
{
    const double t2{t*t}, t3{t2*t};
    return {-0.5*t3+t2-0.5*t, 1.5*t3-2.5*t2+1.0,
        -1.5*t3+2.0*t2+0.5*t, 0.5*t3-0.5*t2};
}
boost::json::value Number(double x) { return std::isfinite(x) ? boost::json::value(x) : boost::json::value(nullptr); }
struct Endpoint
{
    double width{}, amplitude{}, charge{}, loss{std::numeric_limits<double>::infinity()};
    bool rank{}, amplitude_boundary{}, exhausted{};
};
struct Profile
{
    Endpoint best;
    std::vector<Endpoint> minima;
    boost::json::array scan;
    std::size_t evaluations{}, linear_solves{};
    bool exhausted{};
};
Endpoint Solve(const Design & design, const Eigen::VectorXd & y, double fixed_charge,
    bool fit_charge, double cutoff, double log_width, bool svd)
{
    Endpoint e; e.width = std::exp(log_width); e.charge = fixed_charge;
    const auto basis{EvaluateDesign(design,e.width,cutoff)};
    const Eigen::MatrixXd x{basis.leftCols(fit_charge ? 2 : 1)};
    Eigen::VectorXd rhs{y};
    if (!fit_charge) rhs -= fixed_charge*basis.col(1);
    const Eigen::JacobiSVD<Eigen::MatrixXd> decomposition(x,Eigen::ComputeThinU|Eigen::ComputeThinV);
    e.rank = decomposition.rank() == x.cols();
    if (!e.rank) return e;
    Eigen::VectorXd beta;
    if (svd) beta = decomposition.solve(rhs);
    else beta = x.colPivHouseholderQr().solve(rhs);
    e.amplitude = beta(0);
    if (fit_charge) e.charge = beta(1);
    if (e.amplitude <= 0.0)
    {
        e.amplitude_boundary = true; e.amplitude = 0.0;
        // Exact active-set solution for the single nonnegative coefficient.
        if (fit_charge) e.charge = basis.col(1).dot(y)/basis.col(1).squaredNorm();
    }
    e.loss = (y-e.amplitude*basis.col(0)-e.charge*basis.col(1)).squaredNorm()/static_cast<double>(y.size());
    return e;
}
Profile Search(const Design & design, const Eigen::VectorXd & y, double charge,
    bool fit_charge, double cutoff, int points, bool svd, int budget)
{
    Profile out;
    std::vector<Endpoint> grid;
    std::vector<double> coordinate;
    const double lo{std::log(lower_width)}, hi{std::log(upper_width)};
    auto evaluate = [&](double u) {
        ++out.evaluations; ++out.linear_solves;
        return Solve(design,y,charge,fit_charge,cutoff,u,svd);
    };
    for (int k=0;k<points;++k)
    {
        coordinate.push_back(lo+(hi-lo)*static_cast<double>(k)/static_cast<double>(points-1));
        grid.push_back(evaluate(coordinate.back()));
        out.scan.emplace_back(boost::json::array{coordinate.back(),Number(grid.back().loss)});
    }
    out.minima.push_back(grid.front()); out.minima.push_back(grid.back());
    for (int k=1;k<points-1;++k)
    {
        const auto i{static_cast<std::size_t>(k)};
        if (!std::isfinite(grid[i].loss) || grid[i].loss > grid[i-1].loss || grid[i].loss > grid[i+1].loss) continue;
        if (budget <= 0) { out.exhausted = true; out.minima.push_back(grid[i]); continue; }
        std::uintmax_t iterations{static_cast<std::uintmax_t>(budget)};
        const auto minimum{boost::math::tools::brent_find_minima(
            [&](double u) { const double loss{evaluate(u).loss};
                return std::isfinite(loss) ? loss : std::numeric_limits<double>::max(); },
            coordinate[i-1],coordinate[i+1],std::numeric_limits<double>::digits/2,iterations)};
        auto endpoint{evaluate(minimum.first)};
        endpoint.exhausted = iterations >= static_cast<std::uintmax_t>(budget);
        out.exhausted |= endpoint.exhausted; out.minima.push_back(endpoint);
    }
    for (const auto & candidate : out.minima) if (candidate.loss < out.best.loss) out.best = candidate;
    return out;
}
boost::json::object EndpointJSON(const Endpoint & e)
{
    return {{"A",Number(e.amplitude)},{"B",Number(e.width)},{"C",Number(e.charge)},
        {"loss",Number(e.loss)},{"rank",e.rank},{"amplitude_boundary",e.amplitude_boundary},
        {"budget_exhausted",e.exhausted}};
}
double Difference(const Endpoint & a, const Endpoint & b, bool fit_charge)
{
    double value{};
    for (const auto pair : {std::pair{a.amplitude,b.amplitude},std::pair{a.width,b.width},
        std::pair{fit_charge ? a.charge : 0.0,fit_charge ? b.charge : 0.0}})
        value = std::max(value,std::abs(pair.first-pair.second)/std::max({1.0,std::abs(pair.first),std::abs(pair.second)}));
    return value;
}
} // namespace

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
Basis EvaluateBasis(double square, double width, double cutoff)
{
    if (!(width>0.0) || !std::isfinite(width)) throw std::invalid_argument("Invalid matched width.");
    if (square > cutoff*cutoff) return {};
    const double r{std::sqrt(square)}, r2{r*r}, b2{width*width}, exponent{std::exp(-r2/(2.0*b2))};
    Basis out;
    out.gaussian = std::pow(2.0*pi*b2,-1.5)*exponent;
    out.gaussian_log_width = out.gaussian*(r2/b2-3.0);
    const double center{std::sqrt(2.0/pi)/width};
    if (r < 1e-5) {out.charge=center; out.charge_log_width=-center;}
    else if (r <= 2.5) {out.charge=std::erf(r/width/std::sqrt(2.0))/r; out.charge_log_width=-center*exponent;}
    return out;
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
boost::json::object Fit(const Design & design, const Eigen::VectorXd & y,
    double charge, bool fit_charge, double cutoff, int maximum_iterations)
{
    if (design.size()!=static_cast<std::size_t>(y.size()) || y.size()==0 || !y.allFinite())
        throw std::invalid_argument("Invalid matched fit data.");
    const auto primary{Search(design,y,charge,fit_charge,cutoff,129,false,maximum_iterations)};
    const auto reference{Search(design,y,charge,fit_charge,cutoff,257,true,maximum_iterations)};
    const auto & e{primary.best};
    auto out{EndpointJSON(e)};
    out["reference"] = EndpointJSON(reference.best);
    out["profile"] = primary.scan; out["reference_profile"] = reference.scan;
    out["evaluations"] = primary.evaluations; out["reference_evaluations"] = reference.evaluations;
    out["linear_solves"] = primary.linear_solves+reference.linear_solves;
    out["n"] = y.size();
    const double difference{Difference(e,reference.best,fit_charge)};
    out["reference_difference"] = Number(difference);
    bool ambiguous{};
    for (const auto * profile : {&primary,&reference})
    {
        boost::json::array minima;
        for (const auto & candidate : profile->minima)
        {
            minima.push_back(EndpointJSON(candidate));
            // Loss ties are measured relative to response energy, not parameter truth.
            if (std::isfinite(e.loss) && std::abs(candidate.loss-e.loss)<=64*std::numeric_limits<double>::epsilon()*
                std::max(1.0,y.squaredNorm()/static_cast<double>(y.size())) && Difference(candidate,e,fit_charge)>1e-6)
                ambiguous = true;
        }
        out[profile==&primary ? "minima" : "reference_minima"] = std::move(minima);
    }
    out["ambiguous"] = ambiguous;
    double stationarity{std::numeric_limits<double>::infinity()}, condition{stationarity}; bool rank{};
    if (e.rank && std::isfinite(e.loss))
    {
        const auto basis{EvaluateDesign(design,e.width,cutoff)};
        const Eigen::VectorXd residual{y-e.amplitude*basis.col(0)-e.charge*basis.col(1)};
        Eigen::MatrixXd jacobian(y.size(),fit_charge ? 3 : 2);
        jacobian.col(0) = e.amplitude*basis.col(0);
        jacobian.col(1) = e.amplitude*basis.col(2)+e.charge*basis.col(3);
        if (fit_charge) jacobian.col(2) = basis.col(1);
        stationarity=0.0;
        for (Eigen::Index k=0;k<jacobian.cols();++k)
        {
            const double norm{jacobian.col(k).norm()};
            if (norm>0.0) jacobian.col(k) /= norm;
            stationarity=std::max(stationarity,std::abs(jacobian.col(k).dot(residual))/std::max(1.0,y.norm()));
        }
        const Eigen::JacobiSVD<Eigen::MatrixXd> svd(jacobian);
        rank = svd.rank()==jacobian.cols();
        if (rank) condition=svd.singularValues()(0)/svd.singularValues()(jacobian.cols()-1);
    }
    const bool boundary{e.amplitude_boundary || e.width<=lower_width*(1+1e-6) || e.width>=upper_width*(1-1e-6)};
    const std::string reason{!rank ? "rank-deficient" : primary.exhausted || reference.exhausted ? "budget-exhausted" :
        boundary ? "parameter-boundary" : ambiguous ? "ambiguous-minima" : difference>1e-6 ? "reference-disagreement" :
        stationarity>1e-8 ? "nonstationary" : "qualified"};
    out["stationarity"] = Number(stationarity); out["condition"] = Number(condition);
    out["reason"] = reason; out["qualified"] = reason=="qualified";
    return out;
}
} // namespace second_stage_test::matched
