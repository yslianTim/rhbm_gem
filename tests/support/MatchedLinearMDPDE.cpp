#include "support/MatchedJointAC.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <iostream>
#include <numbers>
#include <set>
#include <Eigen/SparseQR>

namespace second_stage_test::matched::joint_ac {
namespace {
namespace j = boost::json;
constexpr double eps{std::numeric_limits<double>::epsilon()};
j::value Number(double x) { return std::isfinite(x) ? j::value(x) : j::value(nullptr); }
j::array Vector(const Eigen::VectorXd & v)
{
    j::array out; for (double x : v) out.push_back(Number(x)); return out;
}
bool ValidBlocks(const Blocks & blocks, Eigen::Index rows)
{
    if (blocks.empty() || rows<=0) return false;
    std::vector<bool> seen(static_cast<std::size_t>(rows)); std::set<std::size_t> owners;
    for (const auto & b:blocks)
    {
        if (!std::isfinite(b.alpha) || b.alpha<0 || b.rows.empty() || !owners.insert(b.owner).second) return false;
        for (auto p:b.rows)
        {
            if (p<0 || p>=rows || seen[static_cast<std::size_t>(p)]) return false;
            seen[static_cast<std::size_t>(p)]=true;
        }
    }
    return std::all_of(seen.begin(),seen.end(),[](bool v){return v;});
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
Endpoint Iterate(const Eigen::MatrixXd & x, const Eigen::VectorXd & y, Eigen::VectorXd beta,
    Eigen::VectorXd variances, const Blocks & blocks, int budget, double tolerance, const Eigen::SparseMatrix<double> * sparse_design)
{
    Endpoint out; out.beta=std::move(beta); out.variances=std::move(variances); out.stop="budget-exhausted";
    for (int iteration=0; iteration<=budget; ++iteration)
    {
        out.evidence=Evaluate(x,y,out.beta,out.variances,blocks);
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
        bool positive=true;
        for (std::size_t i=0;i<blocks.size();++i)
        {
            double numerator{};
            for (auto p:blocks[i].rows) numerator+=out.evidence.weights(p)*residual(p)*residual(p);
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
            const double objective{Objective(y-x*trial,v,blocks)};
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
    out.evidence=Evaluate(x,y,out.beta,out.variances,blocks);
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
j::object WeightSummary(const Eigen::VectorXd & w)
{
    return {{"minimum",Number(w.minCoeff())},{"maximum",Number(w.maxCoeff())},
        {"mean",Number(w.mean())},{"underflow_count",static_cast<int>((w.array()==0).count())},
        {"below_0_1",static_cast<int>((w.array()<.1).count())},
        {"effective_n",Number(w.sum()*w.sum()/w.squaredNorm())}};
}
std::pair<Eigen::SparseMatrix<double>,Eigen::VectorXd> ReduceSparseRows(
    const Eigen::SparseMatrix<double> & a,const Eigen::VectorXd & rhs)
{
    // Orthogonal elimination within contiguous row tiles. Only structurally
    // present columns need local QR; all original rows and weights participate.
    // The discarded Q^T y tail contributes a coefficient-independent constant,
    // while objective and stationarity are always evaluated on the full design.
    const Eigen::SparseMatrix<double,Eigen::RowMajor> rows(a);
    std::vector<Eigen::Triplet<double>> entries; std::vector<double> response;
    constexpr Eigen::Index tile_size{1024};
    for (Eigen::Index first=0;first<a.rows();first+=tile_size)
    {
        const Eigen::Index count{std::min(tile_size,a.rows()-first)};
        std::vector<int> local_column(static_cast<std::size_t>(a.cols()),-1),columns;
        for (Eigen::Index row=first;row<first+count;++row)
            for (Eigen::SparseMatrix<double,Eigen::RowMajor>::InnerIterator e(rows,row);e;++e)
                local_column[static_cast<std::size_t>(e.col())]=0;
        for (Eigen::Index k=0;k<a.cols();++k) if (local_column[static_cast<std::size_t>(k)]==0)
        {local_column[static_cast<std::size_t>(k)]=static_cast<int>(columns.size()); columns.push_back(static_cast<int>(k));}
        if (columns.empty()) continue;
        Eigen::MatrixXd local{Eigen::MatrixXd::Zero(count,static_cast<Eigen::Index>(columns.size()))};
        for (Eigen::Index row=first;row<first+count;++row)
            for (Eigen::SparseMatrix<double,Eigen::RowMajor>::InnerIterator e(rows,row);e;++e)
                local(row-first,local_column[static_cast<std::size_t>(e.col())])=e.value();
        const Eigen::HouseholderQR<Eigen::MatrixXd> qr(local);
        const Eigen::VectorXd transformed{qr.householderQ().adjoint()*rhs.segment(first,count)};
        const Eigen::Index retained{std::min(count,local.cols())},offset{static_cast<Eigen::Index>(response.size())};
        for (Eigen::Index row=0;row<retained;++row)
        {
            response.push_back(transformed(row));
            for (Eigen::Index col=row;col<local.cols();++col)
                if (qr.matrixQR()(row,col)!=0) entries.emplace_back(offset+row,columns[static_cast<std::size_t>(col)],qr.matrixQR()(row,col));
        }
    }
    Eigen::SparseMatrix<double> reduced(static_cast<Eigen::Index>(response.size()),a.cols());
    reduced.setFromTriplets(entries.begin(),entries.end());
    Eigen::VectorXd transformed(static_cast<Eigen::Index>(response.size()));
    for (std::size_t row=0;row<response.size();++row) transformed(static_cast<Eigen::Index>(row))=response[row];
    return {std::move(reduced),std::move(transformed)};
}
} // namespace

Eigen::VectorXd BlockVariances(const Eigen::VectorXd & r, const Blocks & blocks)
{
    if (!ValidBlocks(blocks,r.size())) throw std::invalid_argument("Invalid observation blocks.");
    Eigen::VectorXd v(static_cast<Eigen::Index>(blocks.size()));
    for (std::size_t i=0;i<blocks.size();++i)
    {
        double sum{}; for (auto p:blocks[i].rows) sum+=r(p)*r(p);
        v(static_cast<Eigen::Index>(i))=sum/static_cast<double>(blocks[i].rows.size());
    }
    return v;
}
double Objective(const Eigen::VectorXd & residual, const Eigen::VectorXd & variances, const Blocks & blocks)
{
    if (!ValidBlocks(blocks,residual.size()) || !residual.allFinite() ||
        variances.size()!=static_cast<Eigen::Index>(blocks.size()) || !variances.allFinite() ||
        (variances.array()<=0).any()) return std::numeric_limits<double>::quiet_NaN();
    double sum{};
    for (std::size_t i=0;i<blocks.size();++i) sum+=BlockObjective(residual,variances(static_cast<Eigen::Index>(i)),blocks[i]);
    return sum/static_cast<double>(blocks.size());
}
LinearResult WeightedSolve(const Eigen::MatrixXd & x, const Eigen::VectorXd & y,
    const Eigen::VectorXd & weights, bool use_svd, bool blocked_svd, const Eigen::SparseMatrix<double> * sparse_design)
{
    LinearResult out; out.reason="invalid-input"; out.beta=Eigen::VectorXd::Zero(x.cols());
    if (x.cols()==0 || x.cols()%2!=0 || x.rows()<=x.cols() || y.size()!=x.rows() || weights.size()!=y.size() ||
        !x.allFinite() || !y.allFinite() || !weights.allFinite() || (weights.array()<0).any()) return out;
    const Eigen::VectorXd scales{x.colwise().norm()};
    if ((scales.array()==0).any()) {out.reason="rank-deficient"; return out;}
    const bool sparse{sparse_design && !use_svd};
    Eigen::MatrixXd z;
    Eigen::SparseMatrix<double> sparse_z;
    if (sparse)
    {
        if (sparse_design->rows()!=x.rows() || sparse_design->cols()!=x.cols()) return out;
        sparse_z=*sparse_design;
        for (int k=0;k<sparse_z.outerSize();++k)
            for (Eigen::SparseMatrix<double>::InnerIterator entry(sparse_z,k);entry;++entry)
                entry.valueRef()*=std::sqrt(weights(entry.row()))/scales(k);
    }
    else z=weights.cwiseSqrt().asDiagonal()*x*scales.cwiseInverse().asDiagonal();
    const double znorm{sparse ? sparse_z.norm() : z.norm()};
    const Eigen::VectorXd rhs{weights.cwiseSqrt().array()*y.array()};
    const double rank_threshold{eps*static_cast<double>(std::max(x.rows(),x.cols()))};
    Eigen::VectorXd b{Eigen::VectorXd::Zero(x.cols())};
    std::vector<bool> free(static_cast<std::size_t>(x.cols()),true);
    // Begin with the unconstrained problem. Block infeasible A coefficients and
    // release blocked ones when the dual gradient requires it (mixed NNLS).
    for (Eigen::Index iteration=0; iteration<20*x.cols()*x.cols(); ++iteration)
    {
        std::vector<Eigen::Index> columns;
        for (Eigen::Index k=0;k<x.cols();++k) if (free[static_cast<std::size_t>(k)]) columns.push_back(k);
        Eigen::MatrixXd a;
        if (!sparse)
        {
            a.resize(x.rows(),static_cast<Eigen::Index>(columns.size()));
            for (std::size_t k=0;k<columns.size();++k) a.col(static_cast<Eigen::Index>(k))=z.col(columns[k]);
        }
        Eigen::VectorXd solution;
        if (sparse)
        {
            Eigen::SparseMatrix<double> selected(x.rows(),static_cast<Eigen::Index>(columns.size()));
            selected.reserve(sparse_z.nonZeros()); double maximum_norm{};
            for (std::size_t k=0;k<columns.size();++k)
            {
                selected.startVec(static_cast<Eigen::Index>(k)); double squared{};
                for (Eigen::SparseMatrix<double>::InnerIterator entry(sparse_z,columns[k]);entry;++entry)
                {
                    selected.insertBack(entry.row(),static_cast<Eigen::Index>(k))=entry.value();
                    squared+=entry.value()*entry.value();
                }
                maximum_norm=std::max(maximum_norm,std::sqrt(squared));
            }
            selected.finalize();
            if (maximum_norm==0) {out.reason="rank-deficient"; return out;}
            Eigen::SparseQR<Eigen::SparseMatrix<double>,Eigen::COLAMDOrdering<int>> qr;
            // SparseQR uses an absolute pivot threshold. The scale is the
            // largest norm of the weighted, column-normalized design.
            const auto reduced{ReduceSparseRows(selected,rhs)};
            qr.setPivotThreshold(rank_threshold*maximum_norm); qr.compute(reduced.first);
            if (qr.info()!=Eigen::Success) {out.reason="nonfinite"; return out;}
            out.rank=static_cast<int>(qr.rank()); solution=qr.solve(reduced.second);
        }
        else if (use_svd && blocked_svd)
        {
            // Independent orthogonal reduction of the full weighted problem.
            // SVD(R) preserves its singular values and solves against Q^T y,
            // without materializing the tall left singular-vector matrix.
            const Eigen::HouseholderQR<Eigen::MatrixXd> reduction(a);
            const Eigen::MatrixXd r{reduction.matrixQR().topRows(a.cols()).triangularView<Eigen::Upper>()};
            const Eigen::VectorXd transformed{reduction.householderQ().adjoint()*rhs};
            Eigen::JacobiSVD<Eigen::MatrixXd> svd(r,Eigen::ComputeFullU|Eigen::ComputeFullV);
            svd.setThreshold(rank_threshold); out.rank=static_cast<int>(svd.rank()); solution=svd.solve(transformed.head(a.cols()));
        }
        else if (use_svd)
        {
            Eigen::JacobiSVD<Eigen::MatrixXd> svd(a,Eigen::ComputeThinU|Eigen::ComputeThinV);
            svd.setThreshold(rank_threshold); out.rank=static_cast<int>(svd.rank()); solution=svd.solve(rhs);
        }
        else if (x.cols()>2)
        {
            // Reduce the tall design with blocked orthogonal transformations,
            // then pivot its small R. This solves the same LS problem without
            // normal equations, sparsity truncation, or forming Q explicitly.
            const Eigen::HouseholderQR<Eigen::MatrixXd> reduction(a);
            const Eigen::MatrixXd r{reduction.matrixQR().topRows(a.cols()).triangularView<Eigen::Upper>()};
            const Eigen::VectorXd transformed{reduction.householderQ().adjoint()*rhs};
            Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(r); qr.setThreshold(rank_threshold);
            out.rank=static_cast<int>(qr.rank()); solution=qr.solve(transformed.head(a.cols()));
        }
        else
        {
            Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(a); qr.setThreshold(rank_threshold);
            out.rank=static_cast<int>(qr.rank()); solution=qr.solve(rhs);
        }
        ++out.solves;
        if (out.rank!=static_cast<int>(columns.size())) {out.reason="rank-deficient"; return out;}
        Eigen::VectorXd candidate{Eigen::VectorXd::Zero(x.cols())};
        for (std::size_t k=0;k<columns.size();++k) candidate(columns[k])=solution(static_cast<Eigen::Index>(k));
        if (!candidate.allFinite()) {out.reason="nonfinite"; return out;}
        double step{1}; Eigen::Index blocking{-1};
        for (Eigen::Index k=0;k<x.cols();k+=2) if (candidate(k)<0 && free[static_cast<std::size_t>(k)])
        {
            const double t{b(k)/(b(k)-candidate(k))};
            if (blocking<0 || t<step) {step=t; blocking=k;}
        }
        if (blocking>=0)
        {
            b+=step*(candidate-b); b(blocking)=0; free[static_cast<std::size_t>(blocking)]=false; continue;
        }
        b=candidate;
        Eigen::VectorXd gradient;
        if (sparse) gradient=sparse_z.transpose()*(rhs-sparse_z*b);
        else gradient=z.transpose()*(rhs-z*b);
        const double tolerance{128*eps*std::max(1.0,znorm*(rhs.norm()+znorm*b.norm()))};
        Eigen::Index release{-1}; double largest{tolerance};
        for (Eigen::Index k=0;k<x.cols();k+=2) if (!free[static_cast<std::size_t>(k)] && gradient(k)>largest)
        {release=k; largest=gradient(k);}
        if (release>=0) {free[static_cast<std::size_t>(release)]=true; ++out.releases; continue;}
        out.beta=b.cwiseQuotient(scales); out.valid=true; out.reason="solved"; return out;
    }
    out.reason="active-set-budget-exhausted"; return out;
}

Evidence Evaluate(const Eigen::MatrixXd & x, const Eigen::VectorXd & y,
    const Eigen::VectorXd & beta, const Eigen::VectorXd & variances, const Blocks & blocks)
{
    Evidence out; out.reason="invalid-input"; out.objective=out.stationarity=std::numeric_limits<double>::quiet_NaN();
    if (x.rows()<=x.cols() || x.cols()==0 || x.cols()%2 || x.cols()!=beta.size() || x.rows()!=y.size() ||
        !ValidBlocks(blocks,y.size()) || variances.size()!=static_cast<Eigen::Index>(blocks.size())) return out;
    if (!x.allFinite() || !y.allFinite() || !beta.allFinite()) {out.reason="nonfinite"; return out;}
    for (Eigen::Index k=0;k<beta.size();k+=2) if (beta(k)<0) {out.reason="infeasible-amplitude"; return out;}
    const Eigen::VectorXd r{y-x*beta};
    // A single exactly fitted observation block is already a scale boundary.
    for (std::size_t i=0;i<blocks.size();++i)
    {
        double rss{}, ys{}, xs{};
        for (auto p:blocks[i].rows) {rss+=r(p)*r(p); ys+=y(p)*y(p); xs+=x.row(p).squaredNorm();}
        const double bound{64*eps*(std::sqrt(ys)+std::sqrt(xs)*beta.norm())};
        if (std::sqrt(rss)<=bound)
        {out.reason="exact-fit-boundary"; out.failure_owner=static_cast<int>(blocks[i].owner); return out;}
        const double v{variances(static_cast<Eigen::Index>(i))};
        if (!std::isfinite(v) || v<=0)
        {out.reason=std::isfinite(v) ? "variance-boundary" : "nonfinite"; out.failure_owner=static_cast<int>(blocks[i].owner); return out;}
    }
    out.objective=Objective(r,variances,blocks);
    out.weights.resize(y.size()); out.prefactors.resize(y.size());
    out.denominators.resize(variances.size()); out.scaled=Eigen::VectorXd::Zero(beta.size()+variances.size());
    Eigen::VectorXd logq(variances.size());
    for (std::size_t i=0;i<blocks.size();++i)
    {
        const double a{blocks[i].alpha},v{variances(static_cast<Eigen::Index>(i))};
        logq(static_cast<Eigen::Index>(i))=std::log1p(a)-std::log(static_cast<double>(blocks.size()))-
            std::log(static_cast<double>(blocks[i].rows.size()))-std::log(v)-.5*a*std::log(2*std::numbers::pi*v);
    }
    // One common multiplier preserves every relative block weight and KKT equation.
    const double offset{logq.maxCoeff()}; double qv{};
    for (std::size_t i=0;i<blocks.size();++i)
    {
        const auto bi{static_cast<Eigen::Index>(i)}; const auto & block{blocks[i]};
        const double a{block.alpha},v{variances(bi)},q{std::exp(logq(bi)-offset)};
        double sumw{},equation{};
        for (auto p:block.rows)
        {
            const double t{r(p)*r(p)/v};
            const double w{a==0 ? 1 : std::exp(-.5*a*t)};
            out.weights(p)=w; out.prefactors(p)=q; sumw+=w; equation+=w*(t-1); qv+=q*v;
        }
        out.denominators(bi)=sumw-static_cast<double>(block.rows.size())*a*std::pow(1+a,-1.5);
        out.scaled(beta.size()+bi)=equation/static_cast<double>(block.rows.size())+a*std::pow(1+a,-1.5);
    }
    out.linear_weights=out.prefactors.array()*out.weights.array();
    for (std::size_t i=0;i<blocks.size();++i)
        if (!std::isfinite(out.denominators(static_cast<Eigen::Index>(i))) || out.denominators(static_cast<Eigen::Index>(i))<=0)
        {out.reason="invalid-denominator"; out.failure_owner=static_cast<int>(blocks[i].owner); return out;}
    out.scaled.head(beta.size())=x.transpose()*(out.linear_weights.array()*r.array()).matrix();
    for (Eigen::Index k=0;k<beta.size();++k)
    {
        const double norm{std::sqrt((out.prefactors.array()*x.col(k).array().square()).sum())};
        if (norm==0 || qv<=0) {out.reason="rank-deficient"; return out;}
        out.scaled(k)/=std::sqrt(qv)*norm;
        if (k%2==0 && beta(k)==0) out.scaled(k)=std::max(0.0,out.scaled(k));
    }
    out.stationarity=out.scaled.lpNorm<Eigen::Infinity>();
    out.valid=out.scaled.allFinite() && std::isfinite(out.objective);
    out.reason=out.valid ? "valid" : "nonfinite"; return out;
}

j::object Fit(const Eigen::MatrixXd & x, const Eigen::VectorXd & y,
    const Eigen::VectorXd & initial, const Blocks & blocks, int budget, int reference_budget, bool blocked_svd,
    const Eigen::SparseMatrix<double> * sparse_design)
{
    const auto start{std::chrono::steady_clock::now()};
    if (initial.size()!=x.cols() || !initial.allFinite() || !ValidBlocks(blocks,y.size()) || budget<0 || reference_budget<0)
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
    const Eigen::VectorXd scales{x.colwise().norm()};
    const Eigen::MatrixXd normalized{x*scales.cwiseInverse().asDiagonal()};
    result["design_spectrum"]=Spectrum(normalized,blocked_svd);
    j::array branches; std::vector<Endpoint> primary, references; std::vector<Eigen::VectorXd> checked;
    std::vector<bool> qualified;
    for (int seed=0;seed<2;++seed)
    {
        const Eigen::VectorXd b{seed==1 ? ls.beta : initial};
        const Eigen::VectorXd v{BlockVariances(y-x*b,blocks)};
        auto endpoint{Iterate(x,y,b,v,blocks,budget,1e-8,sparse_design)};
        auto reference{endpoint.stop=="stationary" ?
            Iterate(x,y,endpoint.beta,endpoint.variances,blocks,reference_budget,1e-10,sparse_design) : endpoint};
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
            const auto check{WeightedSolve(x,y,reference.evidence.linear_weights,true,blocked_svd)};
            if (check.valid) verified=check.beta;
            const double delta{std::max(Difference(endpoint.beta,reference.beta),Difference(reference.beta,verified))};
            const double vdelta{ScaleDifference(endpoint.variances,reference.variances)};
            branch["svd_beta"]=Vector(verified); branch["svd_reason"]=check.reason;
            branch["svd_linear_solves"]=check.solves; branch["coefficient_difference"]=Number(delta);
            branch["log_variance_difference"]=Number(vdelta);
            Eigen::MatrixXd weighted{reference.evidence.linear_weights.cwiseSqrt().asDiagonal()*x};
            for (Eigen::Index k=0;k<weighted.cols();++k)
            {
                const double norm{weighted.col(k).norm()}; if (norm>0) weighted.col(k)/=norm;
            }
            branch["weighted_spectrum"]=Spectrum(weighted,blocked_svd);
            pass &= check.valid && delta<=1e-6 && vdelta<=1e-6 && branch.at("weighted_spectrum").at("rank").as_int64()==x.cols();
        }
        if (reference.evidence.weights.size()==y.size())
        {
            const auto & e{reference.evidence}; branch["weights"]=WeightSummary(e.weights);
            j::array diagnostics;
            for (std::size_t i=0;i<blocks.size();++i)
            {
                const auto & block{blocks[i]}; Eigen::VectorXd w(static_cast<Eigen::Index>(block.rows.size())); double mass{};
                for (std::size_t k=0;k<block.rows.size();++k) {w(static_cast<Eigen::Index>(k))=e.weights(block.rows[k]); mass+=e.linear_weights(block.rows[k]);}
                auto d{WeightSummary(w)}; d["owner"]=block.owner; d["alpha"]=block.alpha;
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
    result["seconds"]=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    return result;
}
} // namespace second_stage_test::matched::joint_ac
