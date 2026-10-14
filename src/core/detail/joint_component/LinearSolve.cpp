#include "Numerics.hpp"
#include "ResourceWork.hpp"
#include "CompactSvd.hpp"
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <Eigen/SparseQR>
#include <Eigen/QR>
#include <Eigen/SVD>
#include <algorithm>
#include <cmath>
#include <cstdint>
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
    Eigen::VectorXd out=Eigen::VectorXd::Zero(x.cols());
    for(int k=0;k<x.outerSize();++k) for(Sparse::InnerIterator e(x,k);e;++e) out(e.col())+=e.value()*e.value();
    out=out.cwiseSqrt();
    return out;
}
using RefinementFloat=boost::multiprecision::cpp_dec_float_50;
RefinementFloat ExactDouble(double value)
{
    if(value==0.0) return 0;
    int exponent{};
    const double fraction=std::frexp(std::abs(value),&exponent);
    const auto mantissa=static_cast<std::uint64_t>(std::ldexp(fraction,std::numeric_limits<double>::digits));
    RefinementFloat exact=boost::multiprecision::ldexp(RefinementFloat(mantissa),
        exponent-std::numeric_limits<double>::digits);
    if(std::signbit(value)) exact=-exact;
    return exact;
}
bool HighPrecisionFaceSolve(const Sparse & x,VectorRef y,const std::vector<Eigen::Index> & columns,
    const Eigen::VectorXd & scales,Vector & beta)
{
    constexpr Eigen::Index maximum_columns=32;
    constexpr std::uint64_t maximum_products=4'000'000;
    const auto count=static_cast<Eigen::Index>(columns.size());
    if(count==0 || count>maximum_columns || x.rows()<=0 ||
        static_cast<std::uint64_t>(x.rows())>maximum_products/static_cast<std::uint64_t>(count)/
            static_cast<std::uint64_t>(count)) return false;
    std::vector<std::vector<RefinementFloat>> normal(static_cast<std::size_t>(count),
        std::vector<RefinementFloat>(static_cast<std::size_t>(count)+1));
    std::vector<RefinementFloat> values(static_cast<std::size_t>(count));
    for(Eigen::Index row=0;row<x.rows();++row)
    {
        for(Eigen::Index column=0;column<count;++column)
        {
            const auto source=columns[static_cast<std::size_t>(column)];
            values[static_cast<std::size_t>(column)]=
                ExactDouble(x.coeff(row,source))/ExactDouble(scales(source));
        }
        for(Eigen::Index left=0;left<count;++left)
        {
            auto & equation=normal[static_cast<std::size_t>(left)];
            const auto & lhs=values[static_cast<std::size_t>(left)];
            equation[static_cast<std::size_t>(count)]+=lhs*ExactDouble(y(row));
            for(Eigen::Index right=0;right<count;++right)
                equation[static_cast<std::size_t>(right)]+=lhs*values[static_cast<std::size_t>(right)];
        }
    }
    for(Eigen::Index pivot=0;pivot<count;++pivot)
    {
        Eigen::Index selected=pivot;
        for(Eigen::Index row=pivot+1;row<count;++row)
            if(boost::multiprecision::abs(normal[static_cast<std::size_t>(row)][static_cast<std::size_t>(pivot)])>
                boost::multiprecision::abs(normal[static_cast<std::size_t>(selected)][static_cast<std::size_t>(pivot)]))
                selected=row;
        if(normal[static_cast<std::size_t>(selected)][static_cast<std::size_t>(pivot)]==0) return false;
        std::swap(normal[static_cast<std::size_t>(pivot)],normal[static_cast<std::size_t>(selected)]);
        const auto diagonal=normal[static_cast<std::size_t>(pivot)][static_cast<std::size_t>(pivot)];
        for(Eigen::Index row=pivot+1;row<count;++row)
        {
            auto & equation=normal[static_cast<std::size_t>(row)];
            const auto factor=equation[static_cast<std::size_t>(pivot)]/diagonal;
            for(Eigen::Index column=pivot+1;column<=count;++column)
                equation[static_cast<std::size_t>(column)]-=factor*
                    normal[static_cast<std::size_t>(pivot)][static_cast<std::size_t>(column)];
        }
    }
    std::vector<RefinementFloat> solution(static_cast<std::size_t>(count));
    for(Eigen::Index row=count;row-->0;)
    {
        RefinementFloat value=normal[static_cast<std::size_t>(row)][static_cast<std::size_t>(count)];
        for(Eigen::Index column=row+1;column<count;++column)
            value-=normal[static_cast<std::size_t>(row)][static_cast<std::size_t>(column)]*
                solution[static_cast<std::size_t>(column)];
        solution[static_cast<std::size_t>(row)]=value/
            normal[static_cast<std::size_t>(row)][static_cast<std::size_t>(row)];
    }
    beta=Vector::Zero(x.cols());
    for(Eigen::Index column=0;column<count;++column)
    {
        const auto source=columns[static_cast<std::size_t>(column)];
        beta(source)=(solution[static_cast<std::size_t>(column)]/ExactDouble(scales(source))).convert_to<double>();
    }
    return beta.allFinite();
}
}
std::pair<Eigen::MatrixXd,Eigen::VectorXd> ReferenceQR(const Sparse & x,
    const Eigen::VectorXd & weights,const Eigen::VectorXd & scales,VectorRef y)
{
    ++NumericsWorkForTesting().reference;
    const Eigen::SparseMatrix<double,Eigen::RowMajor> rows(x);
    Eigen::MatrixXd r(0,x.cols()); Eigen::VectorXd target(0);
    constexpr Eigen::Index tile=8192;
    for(Eigen::Index first=0;first<x.rows();first+=tile)
    {
        const Eigen::Index n=std::min(tile,x.rows()-first), prior=r.rows();
        Eigen::MatrixXd a=Eigen::MatrixXd::Zero(prior+n,x.cols());
        RecordDenseShape("reference-tile",prior+n,x.cols());
        Eigen::VectorXd rhs=Eigen::VectorXd::Zero(prior+n);
        if(prior) {a.topRows(prior)=r; rhs.head(prior)=target;}
        for(Eigen::Index row=first;row<first+n;++row)
        {
            const double w=std::sqrt(weights(row)); rhs(prior+row-first)=w*y(row);
            for(Eigen::SparseMatrix<double,Eigen::RowMajor>::InnerIterator e(rows,row);e;++e)
                a(prior+row-first,e.col())=w*e.value()/scales(e.col());
        }
        const Eigen::HouseholderQR<Eigen::MatrixXd> qr(a);
        const Eigen::VectorXd transformed=qr.householderQ().adjoint()*rhs;
        const Eigen::Index keep=std::min(a.rows(),a.cols());
        r=qr.matrixQR().topRows(keep).triangularView<Eigen::Upper>(); target=transformed.head(keep);
    }
    return {std::move(r),std::move(target)};
}
namespace {
std::pair<Eigen::SparseMatrix<double>,Eigen::VectorXd> ReduceSparseRows(
    const Eigen::SparseMatrix<double> & a,const Eigen::VectorXd & rhs)
{
    // Orthogonal elimination within contiguous row tiles. Only structurally
    // present columns need local QR; all original rows and weights participate.
    const Eigen::SparseMatrix<double,Eigen::RowMajor> rows(a);
    std::vector<Eigen::Triplet<double>> entries; std::vector<double> response;
    constexpr Eigen::Index tile_size{1024};
    for(Eigen::Index first=0;first<a.rows();first+=tile_size)
    {
        const Eigen::Index count{std::min(tile_size,a.rows()-first)};
        std::vector<int> local_column(static_cast<std::size_t>(a.cols()),-1),columns;
        for(Eigen::Index row=first;row<first+count;++row)
            for(Eigen::SparseMatrix<double,Eigen::RowMajor>::InnerIterator e(rows,row);e;++e)
                local_column[static_cast<std::size_t>(e.col())]=0;
        for(Eigen::Index k=0;k<a.cols();++k) if(local_column[static_cast<std::size_t>(k)]==0)
        {local_column[static_cast<std::size_t>(k)]=static_cast<int>(columns.size()); columns.push_back(static_cast<int>(k));}
        if(columns.empty()) continue;
        Eigen::MatrixXd local{Eigen::MatrixXd::Zero(count,static_cast<Eigen::Index>(columns.size()))};
        RecordDenseShape("ac-row-reduction",local.rows(),local.cols());
        for(Eigen::Index row=first;row<first+count;++row)
            for(Eigen::SparseMatrix<double,Eigen::RowMajor>::InnerIterator e(rows,row);e;++e)
                local(row-first,local_column[static_cast<std::size_t>(e.col())])=e.value();
        const Eigen::HouseholderQR<Eigen::MatrixXd> qr(local);
        const Eigen::VectorXd transformed{qr.householderQ().adjoint()*rhs.segment(first,count)};
        const Eigen::Index retained{std::min(count,local.cols())},offset{static_cast<Eigen::Index>(response.size())};
        for(Eigen::Index row=0;row<retained;++row)
        {
            response.push_back(transformed(row));
            for(Eigen::Index col=row;col<local.cols();++col)
                if(qr.matrixQR()(row,col)!=0) entries.emplace_back(offset+row,columns[static_cast<std::size_t>(col)],qr.matrixQR()(row,col));
        }
    }
    Eigen::SparseMatrix<double> reduced(static_cast<Eigen::Index>(response.size()),a.cols());
    reduced.setFromTriplets(entries.begin(),entries.end());
    Eigen::VectorXd transformed(static_cast<Eigen::Index>(response.size()));
    for(std::size_t row=0;row<response.size();++row) transformed(static_cast<Eigen::Index>(row))=response[row];
    return {std::move(reduced),std::move(transformed)};
}
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
            auto reduced=ReduceSparseRows(selected,rhs); f.reduced=std::move(reduced.first); f.rhs=std::move(reduced.second);
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
            Eigen::SparseQR<Sparse,Eigen::COLAMDOrdering<int>> qr;
            qr.setPivotThreshold(absolute); qr.compute(f.reduced);
            if(qr.info()!=Eigen::Success) {out.valid=false; return out;}
            out.rank+=static_cast<int>(qr.rank()); solution=qr.solve(f.rhs);
        }
        ++out.factorizations;
        for(std::size_t k=0;k<f.positions.size();++k) out.solution(f.positions[k])=solution(static_cast<Eigen::Index>(k));
    }
    return out;
}
template<class Matrix>
LinearResult WeightedSolveImpl(const Matrix & x, VectorRef y,
    const Eigen::VectorXd & weights, bool use_svd, bool blocked_svd, const Eigen::SparseMatrix<double> * sparse_design,
    const LinearPolicy * policy=nullptr,const std::vector<LinearBlock> * blocks=nullptr,double observation_scale_override=0)
{
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
            Eigen::SparseQR<Eigen::SparseMatrix<double>,Eigen::COLAMDOrdering<int>> qr;
            // SparseQR uses an absolute pivot threshold. The scale is the
            // largest norm of the weighted, column-normalized design.
            const auto reduced{ReduceSparseRows(selected,rhs)};
            qr.setPivotThreshold(rank_threshold*maximum_norm); qr.compute(reduced.first);
            if(qr.info()!=Eigen::Success) {out.reason="nonfinite"; return out;}
            out.rank=static_cast<int>(qr.rank()); solution=qr.solve(reduced.second);
        }
        else if (use_svd && blocked_svd)
        {
            // Independent orthogonal reduction of the full weighted problem.
            // SVD(R) preserves its singular values and solves against Q^T y,
            // without materializing the tall left singular-vector matrix.
            Eigen::MatrixXd r; Eigen::VectorXd transformed;
            {
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
        out.beta=b.cwiseQuotient(scales);
        if constexpr (std::is_same_v<Matrix,Sparse>)
        {
            constexpr Eigen::Index maximum_dense_elements=1'000'000; // Cap the dense face matrix at 8 MB.
            const bool unit_weights=(weights.array()==1.0).all();
            const bool bounded=!columns.empty() && x.rows()<=maximum_dense_elements/
                static_cast<Eigen::Index>(columns.size());
            const double observation_scale=observation_scale_override>0 ? observation_scale_override :
                std::max(1.0,original_rhs_norm);
            const Certificate initial=!use_svd && unit_weights && bounded ?
                CertifyLinear(x,y,out.beta,observation_scale) : Certificate{};
            if(initial.available && !initial.kkt_passed)
            {
                Eigen::MatrixXd face(x.rows(),static_cast<Eigen::Index>(columns.size()));
                Eigen::VectorXd theta(static_cast<Eigen::Index>(columns.size()));
                for(std::size_t k=0;k<columns.size();++k)
                {
                    const auto column=columns[k];
                    face.col(static_cast<Eigen::Index>(k)).setZero();
                    for(Sparse::InnerIterator entry(x,column);entry;++entry)
                        face(entry.row(),static_cast<Eigen::Index>(k))=entry.value()/scales(column);
                    theta(static_cast<Eigen::Index>(k))=b(column);
                }
                RecordDenseShape("ac-bounded-kkt-refinement",face.rows(),face.cols());
                Eigen::ColPivHouseholderQR<Eigen::MatrixXd> refinement_qr(face);
                refinement_qr.setThreshold(rank_threshold);
                if(refinement_qr.rank()==static_cast<Eigen::Index>(columns.size()))
                {
                    bool refinement_accepted=false;
                    for(int refinement=0;refinement<8;++refinement)
                    {
                        Eigen::VectorXd residual=y;
                        for(Eigen::Index row=0;row<face.rows();++row)
                        {
                            double prediction{},compensation{};
                            const auto accumulate=[&](double term) {
                                const double next=prediction+term;
                                compensation+=std::abs(prediction)>=std::abs(term) ?
                                    (prediction-next)+term : (term-next)+prediction;
                                prediction=next;
                            };
                            for(Eigen::Index column=0;column<face.cols();++column)
                            {
                                const double left=face(row,column),right=theta(column),term=left*right;
                                const double error=std::fma(left,right,-term);
                                accumulate(term); accumulate(error);
                            }
                            residual(row)-=prediction+compensation;
                        }
                        const Eigen::VectorXd correction=refinement_qr.solve(residual);
                        ++out.solves;
                        if(!correction.allFinite()) break;
                        theta+=correction;
                        Eigen::VectorXd refined_beta=Eigen::VectorXd::Zero(x.cols());
                        for(std::size_t k=0;k<columns.size();++k)
                            refined_beta(columns[k])=theta(static_cast<Eigen::Index>(k))/scales(columns[k]);
                        const Certificate refined=CertifyLinear(x,y,refined_beta,observation_scale);
                        if(refined.available && refined.kkt_passed && refined.objective<=initial.objective)
                        {out.beta=std::move(refined_beta); refinement_accepted=true; break;}
                    }
                    if(!refinement_accepted)
                    {
                        Eigen::VectorXd high_precision_beta;
                        if(HighPrecisionFaceSolve(x,y,columns,scales,high_precision_beta))
                        {
                            ++out.solves;
                            const Certificate refined=CertifyLinear(x,y,high_precision_beta,observation_scale);
                            if(refined.available && refined.kkt_passed && refined.objective<=initial.objective &&
                                refined.active_atoms==initial.active_atoms)
                                out.beta=std::move(high_precision_beta);
                        }
                    }
                }
            }
        }
        out.valid=true; out.reason="solved"; return out;
    }
    out.reason="active-set-budget-exhausted"; return out;
}


} // namespace
LinearResult SolveLinear(const Sparse & x,VectorRef y,const Vector & w,bool svd,bool blocked,const Sparse * cache,
    const LinearPolicy * policy,const std::vector<LinearBlock> * blocks,double observation_scale)
{return WeightedSolveImpl(x,y,w,svd,blocked,cache,policy,blocks,observation_scale);}
LinearResult SolveLinear(const Matrix & x,VectorRef y,const Vector & w,bool svd,bool blocked,const Sparse * cache,
    const LinearPolicy * policy,double observation_scale)
{return WeightedSolveImpl(x,y,w,svd,blocked,cache,policy,nullptr,observation_scale);}
}
