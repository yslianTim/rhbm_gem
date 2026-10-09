#include "Numerics.hpp"
#include "SparseFactor.hpp"
#include "CompactSvd.hpp"
#include <Eigen/QR>
#include <Eigen/SVD>
#include <algorithm>
#include <cmath>
#include <type_traits>

namespace rhbm_gem::core::joint_component {
namespace {
constexpr double eps=std::numeric_limits<double>::epsilon();
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
    Eigen::VectorXd out(x.cols());
    for(Eigen::Index k=0;k<x.cols();++k) out(k)=x.col(k).norm();
    return out;
}
}
std::pair<Eigen::MatrixXd,Eigen::VectorXd> ReferenceQR(const Sparse & x,
    const Eigen::VectorXd & weights,const Eigen::VectorXd & scales,VectorRef y)
{
    return SparseReferenceQR(x,weights,scales,y);
}
namespace {
struct BlockFace
{
    Eigen::VectorXd solution;
    int rank{}, factorizations{};
    bool valid{true};
};
// Only the free-face factorization changes. The caller owns the single global
// active-set trajectory, release tolerance and coefficient vector.
BlockFace SolveBlocks(const Sparse & x,VectorRef y,const Eigen::VectorXd & weights,
    const Eigen::VectorXd & scales,const std::vector<Eigen::Index> & free,
    const std::vector<LinearBlock> & blocks,bool reference,double relative)
{
    struct Factor
    {
        Sparse reduced;
        Eigen::VectorXd rhs;
        std::vector<Eigen::Index> positions;
        std::unique_ptr<Eigen::JacobiSVD<Eigen::MatrixXd>> svd;
    };
    BlockFace out; out.solution=Eigen::VectorXd::Zero(static_cast<Eigen::Index>(free.size()));
    std::vector<Factor> factors; double maximum{};
    for(const auto & block:blocks)
    {
        Factor f; std::vector<Eigen::Index> columns,row_map(static_cast<std::size_t>(x.rows()),-1);
        for(std::size_t k=0;k<block.rows.size();++k) row_map[static_cast<std::size_t>(block.rows[k])]=static_cast<Eigen::Index>(k);
        for(std::size_t k=0;k<free.size();++k) if(std::find(block.columns.begin(),block.columns.end(),free[k])!=block.columns.end())
        {columns.push_back(free[k]); f.positions.push_back(static_cast<Eigen::Index>(k));}
        const auto n=static_cast<Eigen::Index>(block.rows.size()),p=static_cast<Eigen::Index>(columns.size());
        if(!p) continue;
        if(n<p) {out.valid=false; return out;}
        Sparse selected(n,p); std::vector<Eigen::Triplet<double>> entries;
        Eigen::VectorXd response(n),w(n),norms(p);
        for(Eigen::Index r=0;r<n;++r) {const auto parent=block.rows[static_cast<std::size_t>(r)]; response(r)=y(parent); w(r)=weights(parent);}
        for(Eigen::Index k=0;k<p;++k)
        {
            const auto col=columns[static_cast<std::size_t>(k)]; norms(k)=scales(col);
            for(Sparse::InnerIterator entry(x,col);entry;++entry)
            {
                const auto row=row_map[static_cast<std::size_t>(entry.row())];
                if(row<0) throw std::invalid_argument("A linear block omits a contributor row.");
                entries.emplace_back(row,k,entry.value());
            }
        }
        selected.setFromTriplets(entries.begin(),entries.end());
        if(reference)
        {
            auto reduced=ReferenceQR(selected,w,norms,response); f.rhs=std::move(reduced.second);
            f.svd=std::make_unique<Eigen::JacobiSVD<Eigen::MatrixXd>>(reduced.first,Eigen::ComputeThinU|Eigen::ComputeThinV);
            maximum=std::max(maximum,f.svd->singularValues()(0));
        }
        else
        {
            for(Eigen::Index k=0;k<p;++k)
            {
                double squared{};
                for(Sparse::InnerIterator entry(selected,k);entry;++entry)
                {entry.valueRef()*=std::sqrt(w(entry.row()))/norms(k); squared+=entry.value()*entry.value();}
                maximum=std::max(maximum,std::sqrt(squared));
            }
            const Eigen::VectorXd rhs=w.cwiseSqrt().array()*response.array();
            f.reduced=std::move(selected); f.rhs=rhs;
        }
        factors.push_back(std::move(f));
    }
    const double absolute=relative*maximum;
    for(auto & f:factors)
    {
        Eigen::VectorXd solution;
        if(reference)
        {
            const double local=f.svd->singularValues()(0);
            f.svd->setThreshold(local>0 ? absolute/local : relative);
            out.rank+=static_cast<int>(f.svd->rank()); solution=f.svd->solve(f.rhs);
        }
        else
        {
            LinearWorkspace workspace;
            try {
                const auto factor=workspace.Factor(f.reduced,f.positions,absolute);
                out.rank+=factor->Rank(); solution=factor->LeastSquares(f.rhs);
            } catch(const std::runtime_error &) {out.valid=false; return out;}
        }
        ++out.factorizations;
        for(std::size_t k=0;k<f.positions.size();++k) out.solution(f.positions[k])=solution(static_cast<Eigen::Index>(k));
    }
    return out;
}
template<class Matrix>
LinearResult WeightedSolveImpl(const Matrix & x, VectorRef y,
    const Eigen::VectorXd & weights, bool use_svd, bool blocked_svd, const Eigen::SparseMatrix<double> * sparse_design,
    const LinearPolicy * policy=nullptr,const std::vector<LinearBlock> * blocks=nullptr,LinearWorkspace * workspace=nullptr)
{
    LinearWorkspace local_workspace;
    if(!workspace) workspace=&local_workspace;
    LinearResult out; out.reason="invalid-input"; out.beta=Eigen::VectorXd::Zero(x.cols());
    if (x.cols()==0 || x.cols()%2!=0 || x.rows()<=x.cols() || y.size()!=x.rows() || weights.size()!=y.size() ||
        !Finite(x) || !y.allFinite() || !weights.allFinite() || (weights.array()<0).any()) return out;
    const Eigen::VectorXd scales{ColumnNorms(x)};
    if ((scales.array()==0).any()) {out.reason="rank-deficient"; return out;}
    if constexpr (std::is_same_v<Matrix,Sparse>) sparse_design=&x;
    const bool sparse{sparse_design && !use_svd};
    Eigen::MatrixXd z;
    Eigen::SparseMatrix<double> sparse_z;
    if (sparse)
    {
        if (sparse_design->rows()!=x.rows() || sparse_design->cols()!=x.cols()) return out;
        sparse_z=*sparse_design;
        for (int k=0;k<sparse_z.outerSize();++k)
            for (Eigen::SparseMatrix<double>::InnerIterator entry(sparse_z,k);entry;++entry)
                entry.valueRef()=(entry.value()*std::sqrt(weights(entry.row())))/scales(k);
    }
    else if constexpr (!std::is_same_v<Matrix,Sparse>) z=weights.cwiseSqrt().asDiagonal()*x*scales.cwiseInverse().asDiagonal();
    Eigen::VectorXd rhs{weights.cwiseSqrt().array()*y.array()};
    const double original_rhs_norm=policy && policy->release_response_norm>=0 ? policy->release_response_norm : rhs.norm();
    if constexpr (std::is_same_v<Matrix,Sparse>) if (use_svd)
    {
        auto reduced=ReferenceQR(x,weights,scales,y); z=std::move(reduced.first); rhs=std::move(reduced.second);
    }
    const double znorm{sparse ? sparse_z.norm() : z.norm()};
    const double rank_threshold{policy ? policy->rank_relative : eps*static_cast<double>(std::max(x.rows(),x.cols()))};
    Eigen::VectorXd b{Eigen::VectorXd::Zero(x.cols())};
    std::vector<bool> free(static_cast<std::size_t>(x.cols()),true);
    // Begin with the unconstrained problem. Block infeasible A coefficients and
    // release blocked ones when the dual gradient requires it (mixed NNLS).
    for (Eigen::Index iteration=0; iteration<(policy ? policy->active_set_iteration_factor : 20)*x.cols()*x.cols(); ++iteration)
    {
        std::vector<Eigen::Index> columns;
        for (Eigen::Index k=0;k<x.cols();++k) if (free[static_cast<std::size_t>(k)]) columns.push_back(k);
        Eigen::MatrixXd a;
        if (!sparse)
        {
            a.resize(z.rows(),static_cast<Eigen::Index>(columns.size()));
            for (std::size_t k=0;k<columns.size();++k) a.col(static_cast<Eigen::Index>(k))=z.col(columns[k]);
        }
        Eigen::VectorXd solution;
        if (blocks)
        {
            if constexpr (std::is_same_v<Matrix,Sparse>)
            {
                const auto face=SolveBlocks(x,y,weights,scales,columns,*blocks,use_svd,rank_threshold);
                out.rank=face.rank; out.block_factorizations+=face.factorizations; solution=face.solution;
                if(!face.valid) {out.reason="rank-deficient-block"; return out;}
            }
            else throw std::invalid_argument("Partitioned linear solves require a sparse design.");
        }
        else if (sparse)
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
            try {
                out.factor=workspace->Factor(selected,columns,rank_threshold*maximum_norm);
                out.rank=out.factor->Rank();
                if(out.rank!=static_cast<int>(columns.size())) {out.reason="rank-deficient"; return out;}
                solution=out.factor->LeastSquares(rhs);
            } catch(const std::runtime_error &) {out.reason="sparse-factorization-failed"; return out;}
        }
        else if (use_svd && blocked_svd)
        {
            // Independent orthogonal reduction of the full weighted problem.
            // SVD(R) preserves its singular values and solves against Q^T y,
            // without materializing the tall left singular-vector matrix.
            Eigen::MatrixXd r; Eigen::VectorXd transformed;
            {
                ++SparseWorkForTesting().reference_compacts;
                WorkTimer timer(SparseWorkForTesting().reference_compact_seconds);
                const Eigen::HouseholderQR<Eigen::MatrixXd> reduction(a);
                r=reduction.matrixQR().topRows(a.cols()).triangularView<Eigen::Upper>();
                transformed=(reduction.householderQ().adjoint()*rhs).eval().head(a.cols());
            }
            const auto svd=CompactSvd(r,rank_threshold,-1,&transformed);
            if(!svd.valid) {out.reason="nonfinite"; return out;}
            out.rank=static_cast<int>(svd.rank); solution=svd.solution;
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
        const double tolerance{(policy ? policy->release_factor : 128)*eps*std::max(1.0,znorm*(original_rhs_norm+znorm*b.norm()))};
        Eigen::Index release{-1}; double largest{tolerance};
        for (Eigen::Index k=0;k<x.cols();k+=2) if (!free[static_cast<std::size_t>(k)] && gradient(k)>largest)
        {release=k; largest=gradient(k);}
        if (release>=0) {free[static_cast<std::size_t>(release)]=true; ++out.releases; continue;}
        out.beta=b.cwiseQuotient(scales); out.valid=true; out.reason="solved"; return out;
    }
    out.reason="active-set-budget-exhausted"; return out;
}


} // namespace
LinearResult SolveLinear(const Sparse & x,VectorRef y,const Vector & w,bool svd,bool blocked,const Sparse * cache,const LinearPolicy * policy,const std::vector<LinearBlock> * blocks,LinearWorkspace * workspace)
{return WeightedSolveImpl(x,y,w,svd,blocked,cache,policy,blocks,workspace);}
LinearResult SolveLinear(const Matrix & x,VectorRef y,const Vector & w,bool svd,bool blocked,const Sparse * cache)
{return WeightedSolveImpl(x,y,w,svd,blocked,cache);}
}
