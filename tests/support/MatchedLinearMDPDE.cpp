#include "support/MatchedJointAC.hpp"
#include "support/JointRuntimeJson.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <iostream>
#include <numbers>
#include <set>
#include <type_traits>
#include <Eigen/SparseQR>

namespace second_stage_test::matched::joint_ac {
template<class Matrix>
Evidence EvaluateImpl(const Matrix &, const Eigen::VectorXd &, const Eigen::VectorXd &,
    const Eigen::VectorXd &, const Blocks &, bool overlap=false);
double ObjectiveImpl(const Eigen::VectorXd &, const Eigen::VectorXd &, const Blocks &, bool);
Eigen::VectorXd BlockVariancesImpl(const Eigen::VectorXd &, const Blocks &, bool);
namespace {
namespace j = boost::json;
constexpr double eps{std::numeric_limits<double>::epsilon()};
j::value Number(double x) { return std::isfinite(x) ? j::value(x) : j::value(nullptr); }
j::array Vector(const Eigen::VectorXd & v)
{
    j::array out; for (double x : v) out.push_back(Number(x)); return out;
}
bool ValidBlocks(const Blocks & blocks, Eigen::Index rows, bool overlap=false)
{
    if (blocks.empty() || rows<=0) return false;
    std::vector<std::size_t> seen(static_cast<std::size_t>(rows)); std::set<std::size_t> owners;
    std::size_t block_id{};
    for (const auto & b:blocks)
    {
        ++block_id;
        if (!std::isfinite(b.alpha) || b.alpha<0 || b.rows.empty() || !owners.insert(b.owner).second) return false;
        for (auto p:b.rows)
        {
            if (p<0 || p>=rows) return false;
            auto & previous=seen[static_cast<std::size_t>(p)];
            if (previous && (!overlap || previous==block_id)) return false;
            previous=block_id;
        }
    }
    return std::all_of(seen.begin(),seen.end(),[](std::size_t v){return v!=0;});
}
double Difference(const Eigen::VectorXd & a, const Eigen::VectorXd & b)
{
    return ((a-b).array().abs()/(1.0+a.array().abs().max(b.array().abs()))).maxCoeff();
}
double ScaleDifference(const Eigen::VectorXd & a, const Eigen::VectorXd & b)
{
    return (a.array().log()-b.array().log()).abs().maxCoeff();
}
double BlockObjective(const Eigen::VectorXd & residual, double v, const Block & b)
{
    double sum{};
    if (b.alpha==0)
    {
        for (auto p:b.rows) sum+=residual(p)*residual(p);
        return .5*(std::log(2*std::numbers::pi*v)+sum/(static_cast<double>(b.rows.size())*v));
    }
    // Algebraically Q_DPD + 1/alpha, without subtracting two divergent constants.
    for (auto p:b.rows) sum+=std::expm1(-.5*b.alpha*residual(p)*residual(p)/v);
    const double logc{-.5*b.alpha*std::log(2*std::numbers::pi*v)};
    return std::exp(logc)*(std::expm1(-.5*std::log1p(b.alpha))-
        (1+b.alpha)/b.alpha*sum/static_cast<double>(b.rows.size()))-std::expm1(logc)/b.alpha;
}
struct Endpoint
{
    Eigen::VectorXd beta, variances;
    Evidence evidence;
    std::string stop;
    int solves{}, iterations{}, releases{}, backtracks{};
    j::array trace;
};
j::object Row(const Endpoint & e)
{
    return {{"beta",Vector(e.beta)},{"variances",Vector(e.variances)},
        {"objective",Number(e.evidence.objective)},{"stationarity",Number(e.evidence.stationarity)},
        {"denominators",Vector(e.evidence.denominators)}, {"stop",e.stop},
        {"failure_owner",e.evidence.failure_owner},{"scaled_equations",Vector(e.evidence.scaled)},
        {"iterations",e.iterations},{"linear_solves",e.solves},{"constraint_releases",e.releases},
        {"backtracks",e.backtracks}};
}
using Sparse = Eigen::SparseMatrix<double>;
bool Finite(const Eigen::MatrixXd & x) { return x.allFinite(); }
bool Finite(const Sparse & x)
{
    for (int k=0;k<x.outerSize();++k) for (Sparse::InnerIterator e(x,k);e;++e)
        if (!std::isfinite(e.value())) return false;
    return true;
}

Eigen::VectorXd ColumnNorms(const Eigen::MatrixXd & x) { return x.colwise().norm(); }
Eigen::VectorXd ColumnNorms(const Sparse & x)
{
    Eigen::VectorXd out=Eigen::VectorXd::Zero(x.cols());
    for (int k=0;k<x.outerSize();++k) for (Sparse::InnerIterator e(x,k);e;++e) out(k)+=e.value()*e.value();
    return out.cwiseSqrt();
}

Eigen::VectorXd RowSquares(const Eigen::MatrixXd & x) { return x.rowwise().squaredNorm(); }
Eigen::VectorXd RowSquares(const Sparse & x)
{
    Eigen::VectorXd out=Eigen::VectorXd::Zero(x.rows());
    for (int k=0;k<x.outerSize();++k) for (Sparse::InnerIterator e(x,k);e;++e) out(e.row())+=e.value()*e.value();
    return out;
}
double WeightedColumnSquare(const Eigen::MatrixXd & x,Eigen::Index k,const Eigen::VectorXd & w)
{ return (w.array()*x.col(k).array().square()).sum(); }
double WeightedColumnSquare(const Sparse & x,Eigen::Index k,const Eigen::VectorXd & w)
{
    double sum{}; for (Sparse::InnerIterator e(x,k);e;++e) sum+=w(e.row())*e.value()*e.value(); return sum;
}
boost::json::object DesignSpectrum(const Eigen::MatrixXd &,const Eigen::VectorXd &,bool);
boost::json::object DesignSpectrum(const Sparse &,const Eigen::VectorXd &,bool);
template<class Matrix>
Endpoint Iterate(const Matrix & x, const Eigen::VectorXd & y, Eigen::VectorXd beta,
    Eigen::VectorXd variances, const Blocks & blocks, int budget, double tolerance, const Eigen::SparseMatrix<double> * sparse_design, bool overlap)
{
    Endpoint out; out.beta=std::move(beta); out.variances=std::move(variances); out.stop="budget-exhausted";
    for (int iteration=0; iteration<=budget; ++iteration)
    {
        out.evidence=EvaluateImpl(x,y,out.beta,out.variances,blocks,overlap);
        if (x.cols()>2 && x.rows()>1000 && iteration%25==0)
            std::cout<<"joint composite iteration="<<iteration<<" tolerance="<<tolerance
                     <<" stationarity="<<out.evidence.stationarity<<std::endl;
        out.trace.push_back(Row(out));
        if (!out.evidence.valid) {out.stop=out.evidence.reason; break;}
        if (out.evidence.stationarity<=tolerance) {out.stop="stationary"; break;}
        if (iteration==budget) break;
        const auto next{WeightedSolve(x,y,out.evidence.linear_weights,false,false,sparse_design)};
        out.solves+=next.solves; out.releases+=next.releases;
        if (!next.valid) {out.stop=next.reason; break;}
        const Eigen::VectorXd residual{y-x*next.beta};
        Eigen::VectorXd proposed(out.variances.size());
        bool positive=true; std::size_t membership{};
        for (std::size_t i=0;i<blocks.size();++i)
        {
            double numerator{};
            for (auto p:blocks[i].rows)
            {
                const double w=overlap ? out.evidence.membership_weights(static_cast<Eigen::Index>(membership++)) : out.evidence.weights(p);
                numerator+=w*residual(p)*residual(p);
            }
            proposed(static_cast<Eigen::Index>(i))=numerator/out.evidence.denominators(static_cast<Eigen::Index>(i));
            if (!std::isfinite(proposed(static_cast<Eigen::Index>(i))) || proposed(static_cast<Eigen::Index>(i))<=0)
            {out.evidence.failure_owner=static_cast<int>(blocks[i].owner); positive=false; break;}
        }
        if (!positive) {out.stop="variance-boundary"; break;}
        bool accepted{};
        for (int backtrack=0; backtrack<=30; ++backtrack)
        {
            const double t{std::ldexp(1.0,-backtrack)};
            const Eigen::VectorXd trial{out.beta+t*(next.beta-out.beta)};
            const Eigen::VectorXd v{(out.variances.array().log()+t*(proposed.array().log()-out.variances.array().log())).exp()};
            const double objective{ObjectiveImpl(y-x*trial,v,blocks,overlap)};
            const double tau{32*eps*std::max(1.0,std::abs(out.evidence.objective))};
            if (std::isfinite(objective) && objective<=out.evidence.objective+tau)
            {
                if ((trial.array()==out.beta.array()).all() && (v.array()==out.variances.array()).all()) break;
                out.beta=trial; out.variances=v; out.backtracks+=backtrack; accepted=true; break;
            }
        }
        if (!accepted) {out.stop="stalled"; break;}
        ++out.iterations;
    }
    const int failure_owner{out.evidence.failure_owner};
    out.evidence=EvaluateImpl(x,y,out.beta,out.variances,blocks,overlap);
    if (failure_owner>=0) out.evidence.failure_owner=failure_owner;
    return out;
}
j::object Spectrum(const Eigen::MatrixXd & x, bool blocked_svd)
{
    Eigen::VectorXd s;
    if (blocked_svd)
    {
        const Eigen::JacobiSVD<Eigen::MatrixXd,Eigen::HouseholderQRPreconditioner> svd(x);
        s=svd.singularValues();
    }
    else
    {
        const Eigen::JacobiSVD<Eigen::MatrixXd> svd(x);
        s=svd.singularValues();
    }
    const double threshold{eps*static_cast<double>(std::max(x.rows(),x.cols()))};
    const int rank{static_cast<int>((s.array()>threshold*s(0)).count())};
    return {{"rank",rank},{"minimum_singular",Number(s.tail(1)(0))},
        {"condition",Number(s(0)/s.tail(1)(0))}};
}
j::object DesignSpectrum(const Eigen::MatrixXd & x,const Eigen::VectorXd & weights,bool blocked)
{
    Eigen::MatrixXd z=weights.cwiseSqrt().asDiagonal()*x;
    for (Eigen::Index k=0;k<z.cols();++k) {const double n=z.col(k).norm(); if (n>0) z.col(k)/=n;}
    return Spectrum(z,blocked);
}
j::object DesignSpectrum(const Sparse & x,const Eigen::VectorXd & weights,bool)
{return runtime_json::DesignSpectrum(rhbm_gem::core::joint_component::DesignSpectrum(x,weights));}
j::object WeightSummary(const Eigen::VectorXd & w)
{
    return {{"minimum",Number(w.minCoeff())},{"maximum",Number(w.maxCoeff())},
        {"mean",Number(w.mean())},{"underflow_count",static_cast<int>((w.array()==0).count())},
        {"below_0_1",static_cast<int>((w.array()<.1).count())},
        {"effective_n",Number(w.sum()*w.sum()/w.squaredNorm())}};
}
} // namespace

Eigen::VectorXd BlockVariancesImpl(const Eigen::VectorXd & r, const Blocks & blocks, bool overlap)
{
    if (!ValidBlocks(blocks,r.size(),overlap)) throw std::invalid_argument("Invalid observation blocks.");
    Eigen::VectorXd v(static_cast<Eigen::Index>(blocks.size()));
    for (std::size_t i=0;i<blocks.size();++i)
    {
        double sum{}; for (auto p:blocks[i].rows) sum+=r(p)*r(p);
        v(static_cast<Eigen::Index>(i))=sum/static_cast<double>(blocks[i].rows.size());
    }
    return v;
}
double ObjectiveImpl(const Eigen::VectorXd & residual, const Eigen::VectorXd & variances, const Blocks & blocks, bool overlap)
{
    if (!ValidBlocks(blocks,residual.size(),overlap) || !residual.allFinite() ||
        variances.size()!=static_cast<Eigen::Index>(blocks.size()) || !variances.allFinite() ||
        (variances.array()<=0).any()) return std::numeric_limits<double>::quiet_NaN();
    double sum{};
    for (std::size_t i=0;i<blocks.size();++i) sum+=BlockObjective(residual,variances(static_cast<Eigen::Index>(i)),blocks[i]);
    return sum/static_cast<double>(blocks.size());
}
template<class Matrix>
Evidence EvaluateImpl(const Matrix & x, const Eigen::VectorXd & y,
    const Eigen::VectorXd & beta, const Eigen::VectorXd & variances, const Blocks & blocks, bool overlap)
{
    Evidence out; out.reason="invalid-input"; out.objective=out.stationarity=std::numeric_limits<double>::quiet_NaN();
    if (x.rows()<=x.cols() || x.cols()==0 || x.cols()%2 || x.cols()!=beta.size() || x.rows()!=y.size() ||
        !ValidBlocks(blocks,y.size(),overlap) || variances.size()!=static_cast<Eigen::Index>(blocks.size())) return out;
    if (!Finite(x) || !y.allFinite() || !beta.allFinite()) {out.reason="nonfinite"; return out;}
    for (Eigen::Index k=0;k<beta.size();k+=2) if (beta(k)<0) {out.reason="infeasible-amplitude"; return out;}
    const Eigen::VectorXd r{y-x*beta};
    const Eigen::VectorXd row_squares=RowSquares(x);
    // A single exactly fitted observation block is already a scale boundary.
    for (std::size_t i=0;i<blocks.size();++i)
    {
        double rss{}, ys{}, xs{};
        for (auto p:blocks[i].rows) {rss+=r(p)*r(p); ys+=y(p)*y(p); xs+=row_squares(p);}
        const double bound{64*eps*(std::sqrt(ys)+std::sqrt(xs)*beta.norm())};
        if (std::sqrt(rss)<=bound)
        {out.reason="exact-fit-boundary"; out.failure_owner=static_cast<int>(blocks[i].owner); return out;}
        const double v{variances(static_cast<Eigen::Index>(i))};
        if (!std::isfinite(v) || v<=0)
        {out.reason=std::isfinite(v) ? "variance-boundary" : "nonfinite"; out.failure_owner=static_cast<int>(blocks[i].owner); return out;}
    }
    out.objective=ObjectiveImpl(r,variances,blocks,overlap);
    if (overlap)
    {
        std::size_t count{}; for (const auto & b:blocks) count+=b.rows.size();
        out.membership_weights.resize(static_cast<Eigen::Index>(count));
        out.block_prefactors.resize(variances.size());
        out.prefactors=Eigen::VectorXd::Zero(y.size());
        out.linear_weights=Eigen::VectorXd::Zero(y.size());
    }
    else {out.weights.resize(y.size()); out.prefactors.resize(y.size());}
    out.denominators.resize(variances.size()); out.scaled=Eigen::VectorXd::Zero(beta.size()+variances.size());
    Eigen::VectorXd logq(variances.size());
    for (std::size_t i=0;i<blocks.size();++i)
    {
        const double a{blocks[i].alpha},v{variances(static_cast<Eigen::Index>(i))};
        logq(static_cast<Eigen::Index>(i))=std::log1p(a)-std::log(static_cast<double>(blocks.size()))-
            std::log(static_cast<double>(blocks[i].rows.size()))-std::log(v)-.5*a*std::log(2*std::numbers::pi*v);
    }
    // One common multiplier preserves every relative block weight and KKT equation.
    const double offset{logq.maxCoeff()}; double qv{}; std::size_t membership{};
    if (overlap) out.log_prefactors=logq;
    for (std::size_t i=0;i<blocks.size();++i)
    {
        const auto bi{static_cast<Eigen::Index>(i)}; const auto & block{blocks[i]};
        const double a{block.alpha},v{variances(bi)},q{std::exp(logq(bi)-offset)};
        double sumw{},equation{};
        if (overlap) out.block_prefactors(bi)=q;
        for (auto p:block.rows)
        {
            const double t{r(p)*r(p)/v};
            const double w{a==0 ? 1 : std::exp(-.5*a*t)};
            if (overlap)
            {
                out.membership_weights(static_cast<Eigen::Index>(membership++))=w;
                out.prefactors(p)+=q; out.linear_weights(p)+=q*w;
            }
            else {out.weights(p)=w; out.prefactors(p)=q;}
            sumw+=w; equation+=w*(t-1); qv+=q*v;
        }
        out.denominators(bi)=sumw-static_cast<double>(block.rows.size())*a*std::pow(1+a,-1.5);
        out.scaled(beta.size()+bi)=equation/static_cast<double>(block.rows.size())+a*std::pow(1+a,-1.5);
    }
    if (!overlap) out.linear_weights=out.prefactors.array()*out.weights.array();
    for (std::size_t i=0;i<blocks.size();++i)
        if (!std::isfinite(out.denominators(static_cast<Eigen::Index>(i))) || out.denominators(static_cast<Eigen::Index>(i))<=0)
        {out.reason="invalid-denominator"; out.failure_owner=static_cast<int>(blocks[i].owner); return out;}
    out.scaled.head(beta.size())=x.transpose()*(out.linear_weights.array()*r.array()).matrix();
    for (Eigen::Index k=0;k<beta.size();++k)
    {
        const double norm{std::sqrt(WeightedColumnSquare(x,k,out.prefactors))};
        if (norm==0 || qv<=0) {out.reason="rank-deficient"; return out;}
        out.scaled(k)/=std::sqrt(qv)*norm;
        if (k%2==0 && beta(k)==0) out.scaled(k)=std::max(0.0,out.scaled(k));
    }
    out.stationarity=out.scaled.lpNorm<Eigen::Infinity>();
    out.valid=out.scaled.allFinite() && std::isfinite(out.objective);
    out.reason=out.valid ? "valid" : "nonfinite"; return out;
}

template<class Matrix>
j::object FitImpl(const Matrix & x, const Eigen::VectorXd & y,
    const Eigen::VectorXd & initial, const Blocks & blocks, int budget, int reference_budget, bool blocked_svd,
    const Eigen::SparseMatrix<double> * sparse_design, bool overlap=false)
{
    const auto start{std::chrono::steady_clock::now()};
    double verification_seconds{};
    if (initial.size()!=x.cols() || !initial.allFinite() || !ValidBlocks(blocks,y.size(),overlap) || budget<0 || reference_budget<0)
        throw std::invalid_argument("Invalid fixed-width composite MDPDE initialization.");
    for (Eigen::Index k=0;k<initial.size();k+=2) if (initial(k)<0)
        throw std::invalid_argument("Negative initial amplitude.");
    const auto ls{WeightedSolve(x,y,Eigen::VectorXd::Ones(y.size()),false,false,sparse_design)};
    j::array metadata;
    for (const auto & b:blocks) metadata.push_back(j::object{{"owner",b.owner},{"alpha",b.alpha},
        {"rows",b.rows.size()},{"lambda",1.0/static_cast<double>(blocks.size())}});
    j::object result{{"schema_version",2},{"blocks",metadata},{"rows",y.size()},{"columns",x.cols()},
        {"qualified",false},{"reason",ls.reason},{"initial_beta",Vector(initial)},
        {"beta",Vector(initial)},{"uncertainty",j::array{}},{"branches",j::array{}},
        {"initial_linear_solves",ls.solves}};
    if (!ls.valid) return result;
    const Eigen::VectorXd scales{ColumnNorms(x)};
    if constexpr (std::is_same_v<Matrix,Sparse>)
    {
        const auto check_start=std::chrono::steady_clock::now();
        result["design_spectrum"]=DesignSpectrum(x,Eigen::VectorXd::Ones(x.rows()),true);
        verification_seconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-check_start).count();
    }
    else
    {
        const Eigen::MatrixXd normalized{x*scales.cwiseInverse().asDiagonal()};
        result["design_spectrum"]=Spectrum(normalized,blocked_svd);
    }
    j::array branches; std::vector<Endpoint> primary, references; std::vector<Eigen::VectorXd> checked;
    std::vector<bool> qualified;
    for (int seed=0;seed<2;++seed)
    {
        const Eigen::VectorXd b{seed==1 ? ls.beta : initial};
        const Eigen::VectorXd v{BlockVariancesImpl(y-x*b,blocks,overlap)};
        auto endpoint{Iterate(x,y,b,v,blocks,budget,1e-8,sparse_design,overlap)};
        auto reference{endpoint.stop=="stationary" ?
            Iterate(x,y,endpoint.beta,endpoint.variances,blocks,reference_budget,1e-10,sparse_design,overlap) : endpoint};
        if (endpoint.stop!="stationary")
        {
            reference.trace.clear(); reference.solves=reference.iterations=reference.releases=reference.backtracks=0;
        }
        j::object branch{{"seed",seed==1 ? "constrained-ls" : "checkpoint"},
            {"initial_beta",Vector(b)},{"initial_variances",Vector(v)},
            {"primary",Row(endpoint)},{"reference",Row(reference)},
            {"trace",std::move(endpoint.trace)},{"reference_trace",std::move(reference.trace)}};
        bool pass{endpoint.stop=="stationary" && reference.stop=="stationary"};
        Eigen::VectorXd verified{reference.beta};
        if (pass)
        {
            const auto check_start=std::chrono::steady_clock::now();
            const auto fresh=EvaluateImpl(x,y,reference.beta,reference.variances,blocks,overlap);
            const auto check{WeightedSolve(x,y,fresh.linear_weights,true,blocked_svd)};
            if (check.valid) verified=check.beta;
            const double delta{std::max(Difference(endpoint.beta,reference.beta),Difference(reference.beta,verified))};
            const double vdelta{ScaleDifference(endpoint.variances,reference.variances)};
            branch["svd_beta"]=Vector(verified); branch["svd_reason"]=check.reason;
            branch["svd_linear_solves"]=check.solves; branch["coefficient_difference"]=Number(delta);
            branch["log_variance_difference"]=Number(vdelta);
            branch["weighted_spectrum"]=DesignSpectrum(x,reference.evidence.linear_weights,blocked_svd);
            verification_seconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-check_start).count();
            pass &= check.valid && delta<=1e-6 && vdelta<=1e-6 && branch.at("weighted_spectrum").at("rank").as_int64()==x.cols();
        }
        if (reference.evidence.linear_weights.size()==y.size())
        {
            const auto & e{reference.evidence}; branch["weights"]=WeightSummary(overlap ? e.linear_weights : e.weights);
            j::array diagnostics; std::size_t membership{};
            for (std::size_t i=0;i<blocks.size();++i)
            {
                const auto & block{blocks[i]}; Eigen::VectorXd w(static_cast<Eigen::Index>(block.rows.size())); double mass{};
                for (std::size_t k=0;k<block.rows.size();++k)
                {
                    const double value=overlap ? e.membership_weights(static_cast<Eigen::Index>(membership++)) : e.weights(block.rows[k]);
                    w(static_cast<Eigen::Index>(k))=value;
                    mass+=overlap ? e.block_prefactors(static_cast<Eigen::Index>(i))*value : e.linear_weights(block.rows[k]);
                }
                auto d{WeightSummary(w)};
                if (overlap) {d["log_prefactor"]=Number(e.log_prefactors(static_cast<Eigen::Index>(i))); d["scaled_prefactor"]=Number(e.block_prefactors(static_cast<Eigen::Index>(i)));}
                d["owner"]=block.owner; d["alpha"]=block.alpha;
                d["variance"]=Number(reference.variances(static_cast<Eigen::Index>(i)));
                d["scale_equation"]=Number(e.scaled(x.cols()+static_cast<Eigen::Index>(i)));
                d["denominator"]=Number(e.denominators(static_cast<Eigen::Index>(i)));
                d["denominator_margin"]=Number(e.denominators(static_cast<Eigen::Index>(i))/static_cast<double>(block.rows.size()));
                d["linear_weight_share"]=Number(mass/e.linear_weights.sum()); diagnostics.push_back(std::move(d));
            }
            branch["block_diagnostics"]=std::move(diagnostics);
        }
        branch["qualified"]=pass;
        primary.push_back(std::move(endpoint)); references.push_back(std::move(reference));
        checked.push_back(std::move(verified)); qualified.push_back(pass); branches.push_back(std::move(branch));
    }
    int chosen{-1};
    for (int k=0;k<2;++k) if (qualified[static_cast<std::size_t>(k)] &&
        (chosen<0 || primary[static_cast<std::size_t>(k)].evidence.objective<primary[static_cast<std::size_t>(chosen)].evidence.objective)) chosen=k;
    if (chosen<0) chosen=0;
    const auto index{static_cast<std::size_t>(chosen)}; const auto & end{primary[index]};
    result["qualified"]=static_cast<bool>(qualified[index]); result["reason"]=qualified[index] ? "qualified" :
        end.stop=="stationary" ? "reference-unverified" : end.stop;
    result["selected_seed"]=chosen; result["beta"]=Vector(end.beta); result["variances"]=Vector(end.variances);
    result["failure_owner"]=end.evidence.failure_owner;
    result["objective"]=Number(end.evidence.objective); result["stationarity"]=Number(end.evidence.stationarity);
    result["uncertainty"]=Vector((end.beta-references[index].beta).cwiseAbs()+
        (references[index].beta-checked[index]).cwiseAbs());
    result["branch_sensitive"]=qualified[0] && qualified[1] && (Difference(primary[0].beta,primary[1].beta)>1e-6 ||
        ScaleDifference(primary[0].variances,primary[1].variances)>1e-6);
    int active{}; for (Eigen::Index k=0;k<end.beta.size();k+=2) active+=end.beta(k)==0;
    result["active_amplitudes"]=active; result["branches"]=std::move(branches);
    if constexpr (std::is_same_v<Matrix,Sparse>) result["independent_validation_seconds"]=verification_seconds;
    result["seconds"]=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    return result;
}
Eigen::VectorXd BlockVariances(const Eigen::VectorXd & r,const Blocks & b) {return BlockVariancesImpl(r,b,false);}
double Objective(const Eigen::VectorXd & r,const Eigen::VectorXd & v,const Blocks & b) {return ObjectiveImpl(r,v,b,false);}
LinearResult WeightedSolve(const Eigen::MatrixXd & x,const Eigen::VectorXd & y,const Eigen::VectorXd & w,
    bool svd,bool blocked,const Sparse * cache) {return rhbm_gem::core::joint_component::SolveLinear(x,y,w,svd,blocked,cache);}
LinearResult WeightedSolve(const Sparse & x,const Eigen::VectorXd & y,const Eigen::VectorXd & w,
    bool svd,bool blocked,const Sparse * cache,const joint_abc::LinearPolicy * policy,const std::vector<LinearBlock> * blocks)
{return rhbm_gem::core::joint_component::SolveLinear(x,y,w,svd,blocked,cache,policy,blocks);}
Evidence Evaluate(const Eigen::MatrixXd & x,const Eigen::VectorXd & y,const Eigen::VectorXd & b,
    const Eigen::VectorXd & v,const Blocks & blocks) {return EvaluateImpl(x,y,b,v,blocks);}
Evidence Evaluate(const Sparse & x,const Eigen::VectorXd & y,const Eigen::VectorXd & b,
    const Eigen::VectorXd & v,const Blocks & blocks) {return EvaluateImpl(x,y,b,v,blocks);}
j::object Fit(const Eigen::MatrixXd & x,const Eigen::VectorXd & y,const Eigen::VectorXd & b,
    const Blocks & blocks,int budget,int reference,bool blocked,const Sparse * cache)
{return FitImpl(x,y,b,blocks,budget,reference,blocked,cache);}
j::object Fit(const Sparse & x,const Eigen::VectorXd & y,const Eigen::VectorXd & b,
    const Blocks & blocks,int budget,int reference,bool blocked,const Sparse * cache)
{return FitImpl(x,y,b,blocks,budget,reference,blocked,cache);}
Evidence EvaluateComposite(const Sparse & x,const Eigen::VectorXd & y,const Eigen::VectorXd & b,
    const Eigen::VectorXd & v,const Blocks & blocks) {return EvaluateImpl(x,y,b,v,blocks,true);}
j::object FitComposite(const Sparse & x,const Eigen::VectorXd & y,const Eigen::VectorXd & b,
    const Blocks & blocks,int budget,int reference)
{return FitImpl(x,y,b,blocks,budget,reference,true,&x,true);}
j::object SparseSpectrum(const Sparse & x,const Eigen::VectorXd & w) {return DesignSpectrum(x,w,true);}
} // namespace second_stage_test::matched::joint_ac
