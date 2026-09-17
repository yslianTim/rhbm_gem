#include "support/AtomBlockGridComposite.hpp"
#include <algorithm>
#include <cmath>
#include <map>
#include <numeric>

namespace second_stage_test::matched::atom_block {
namespace {
namespace j=boost::json;
j::value Number(double v) {return std::isfinite(v) ? j::value(v) : j::value(nullptr);}
struct Group
{
    std::size_t count{}; double rss{},mass{},square{};
    void Add(double residual,double weight)
    {++count; rss+=residual*residual; mass+=weight; square+=weight*weight;}
    j::object Result(double total) const
    {
        return {{"rows",count},{"rmse",count ? Number(std::sqrt(rss/static_cast<double>(count))) : j::value(nullptr)},
            {"weight_share",total>0 ? Number(mass/total) : j::value(nullptr)},
            {"effective_n",square>0 ? Number(mass*mass/square) : j::value(nullptr)}};
    }
};
j::object Weights(std::vector<double> w)
{
    std::sort(w.begin(),w.end()); const double sum=std::accumulate(w.begin(),w.end(),0.0);
    const double square=std::inner_product(w.begin(),w.end(),w.begin(),0.0);
    const auto top=static_cast<std::size_t>(std::ceil(.01*static_cast<double>(w.size())));
    j::array quantiles;
    for (double p:{0.,.01,.1,.5,.9,.99,1.}) quantiles.emplace_back(j::array{p,w[static_cast<std::size_t>(p*static_cast<double>(w.size()-1))]});
    return {{"weight_quantiles",quantiles},{"effective_n",Number(sum*sum/square)},
        {"effective_fraction",Number(sum*sum/square/static_cast<double>(w.size()))},
        {"maximum_row_share",Number(w.back()/sum)},
        {"top_one_percent_share",Number(std::accumulate(w.end()-static_cast<std::ptrdiff_t>(top),w.end(),0.0)/sum)},
        {"underflow_count",static_cast<std::size_t>(std::count(w.begin(),w.end(),0.0))}};
}
}
j::object Diagnostics(const joint_ac::Evidence & e,const Eigen::VectorXd & residual,
    const Eigen::VectorXd & variances,const joint_ac::Blocks & blocks,const unique_grid::Grid & grid,
    const std::vector<Atom> & atoms)
{
    if (e.linear_weights.size()!=residual.size() || e.membership_weights.size()==0 || e.linear_weights.sum()<=0)
        return {{"available",false},{"reason",e.reason}};
    const double total=e.linear_weights.sum(); j::array diagnostics,coverage,radial;
    std::map<std::size_t,Group> cover; std::array<Group,5> distances{}; std::size_t membership{};
    for (std::size_t p=0;p<grid.voxels.size();++p)
        cover[grid.voxels[p].multiplicity].Add(residual(static_cast<Eigen::Index>(p)),e.linear_weights(static_cast<Eigen::Index>(p)));
    for (std::size_t i=0;i<blocks.size();++i)
    {
        const auto & b=blocks[i]; const auto bi=static_cast<Eigen::Index>(i); std::vector<double> w; w.reserve(b.rows.size());
        double mass{}; std::array<Group,5> local{};
        for (auto row:b.rows)
        {
            const double value=e.membership_weights(static_cast<Eigen::Index>(membership++)); w.push_back(value);
            const double weight=e.block_prefactors(bi)*value; mass+=weight;
            const double distance=std::sqrt(SquareDistance(grid.voxels[static_cast<std::size_t>(row)].position,atoms[b.owner].position));
            const auto bin=std::min(std::size_t{4},static_cast<std::size_t>(distance/.5));
            local[bin].Add(residual(row),weight); distances[bin].Add(residual(row),weight);
        }
        auto out=Weights(std::move(w)); out["owner"]=b.owner; out["alpha"]=b.alpha; out["rows"]=b.rows.size();
        out["variance"]=Number(variances(bi)); out["scale_equation"]=Number(e.scaled(e.scaled.size()-variances.size()+bi));
        out["denominator"]=Number(e.denominators(bi)); out["denominator_margin"]=Number(e.denominators(bi)/static_cast<double>(b.rows.size()));
        out["log_prefactor"]=Number(e.log_prefactors(bi)); out["scaled_prefactor"]=Number(e.block_prefactors(bi));
        out["linear_weight_share"]=Number(mass/total); j::array local_radial;
        for (std::size_t k=0;k<local.size();++k)
        {auto item=local[k].Result(mass); item["distance_bin"]=k; local_radial.push_back(std::move(item));}
        out["radial"]=std::move(local_radial); diagnostics.push_back(std::move(out));
    }
    for (const auto & [n,g]:cover) {auto item=g.Result(total); item["coverage"]=n; coverage.push_back(std::move(item));}
    for (std::size_t k=0;k<distances.size();++k)
    {auto item=distances[k].Result(total); item["distance_bin"]=k; radial.push_back(std::move(item));}
    return {{"available",true},{"blocks",diagnostics},{"coverage",coverage},{"membership_radial",radial},
        {"log_prefactor_offset",Number(e.log_prefactors.maxCoeff())},{"membership_count",membership}};
}
} // namespace second_stage_test::matched::atom_block
