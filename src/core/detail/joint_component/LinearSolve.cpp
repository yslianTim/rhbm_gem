#include "Numerics.hpp"
#include "ResourceWork.hpp"
#include <Eigen/SparseQR>
#include <Eigen/QR>
#include <algorithm>
#include <chrono>
#include <cmath>

namespace rhbm_gem::core::joint_component {
namespace {
constexpr double eps=std::numeric_limits<double>::epsilon();
bool Finite(const Sparse & x)
{
    for (int k=0;k<x.outerSize();++k) for (Sparse::InnerIterator e(x,k);e;++e)
        if (!std::isfinite(e.value())) return false;
    return true;
}
Eigen::VectorXd ColumnNorms(const Sparse & x)
{
    Eigen::VectorXd out=Eigen::VectorXd::Zero(x.cols());
    for(int k=0;k<x.outerSize();++k) for(Sparse::InnerIterator e(x,k);e;++e) out(e.col())+=e.value()*e.value();
    out=out.cwiseSqrt();
    return out;
}
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
    const std::vector<LinearBlock> & blocks,double relative)
{
    struct Factor
    {
        Sparse reduced;
        Eigen::VectorXd rhs;
        std::vector<Eigen::Index> positions;
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
        for(Eigen::Index k=0;k<p;++k)
        {
            double squared{};
            for(Sparse::InnerIterator entry(selected,k);entry;++entry)
            {entry.valueRef()*=std::sqrt(w(entry.row()))/norms(k); squared+=entry.value()*entry.value();}
            maximum=std::max(maximum,std::sqrt(squared));
        }
        const Eigen::VectorXd rhs=w.cwiseSqrt().array()*response.array();
        auto reduced=ReduceSparseRows(selected,rhs); f.reduced=std::move(reduced.first); f.rhs=std::move(reduced.second);
        factors.push_back(std::move(f));
    }
    const double absolute=relative*maximum;
    for(auto & f:factors)
    {
        Eigen::SparseQR<Sparse,Eigen::COLAMDOrdering<int>> qr;
        qr.setPivotThreshold(absolute); qr.compute(f.reduced);
        if(qr.info()!=Eigen::Success) {out.valid=false; return out;}
        out.rank+=static_cast<int>(qr.rank()); const Eigen::VectorXd solution=qr.solve(f.rhs);
        ++out.factorizations;
        for(std::size_t k=0;k<f.positions.size();++k) out.solution(f.positions[k])=solution(static_cast<Eigen::Index>(k));
    }
    return out;
}
LinearResult WeightedSolveImpl(const Sparse & x, VectorRef y,const Eigen::VectorXd & weights,
    const LinearPolicy * policy,const std::vector<LinearBlock> * blocks,double observation_scale_override)
{
    LinearResult out; out.reason="invalid-input"; out.beta=Eigen::VectorXd::Zero(x.cols());
    if (x.cols()==0 || x.cols()%2!=0 || x.rows()<=x.cols() || y.size()!=x.rows() || weights.size()!=y.size() ||
        !Finite(x) || !y.allFinite() || !weights.allFinite() || (weights.array()<0).any()) return out;
    const Eigen::VectorXd scales{ColumnNorms(x)};
    if ((scales.array()==0).any()) {out.reason="rank-deficient"; return out;}
    Sparse sparse_z=x;
    for (int k=0;k<sparse_z.outerSize();++k)
        for (Sparse::InnerIterator entry(sparse_z,k);entry;++entry)
            entry.valueRef()=(entry.value()*std::sqrt(weights(entry.row())))/scales(k);
    Eigen::VectorXd rhs{weights.cwiseSqrt().array()*y.array()};
    const double original_rhs_norm=policy && policy->release_response_norm>=0 ? policy->release_response_norm : rhs.norm();
    const double znorm{sparse_z.norm()};
    const double rank_threshold{policy ? policy->rank_relative : eps*static_cast<double>(std::max(x.rows(),x.cols()))};
    Eigen::VectorXd b{Eigen::VectorXd::Zero(x.cols())};
    std::vector<bool> free(static_cast<std::size_t>(x.cols()),true);
    // Begin with the unconstrained problem. Block infeasible A coefficients and
    // release blocked ones when the dual gradient requires it (mixed NNLS).
    for (Eigen::Index iteration=0; iteration<(policy ? policy->active_set_iteration_factor : 20)*x.cols()*x.cols(); ++iteration)
    {
        ++out.active_set_iterations;
        std::vector<Eigen::Index> columns;
        for (Eigen::Index k=0;k<x.cols();++k) if (free[static_cast<std::size_t>(k)]) columns.push_back(k);
        Eigen::VectorXd solution;
        if (blocks)
        {
            const auto face=SolveBlocks(x,y,weights,scales,columns,*blocks,rank_threshold);
            out.rank=face.rank; out.block_factorizations+=face.factorizations; solution=face.solution;
            if(!face.valid) {out.reason="rank-deficient-block"; return out;}
        }
        else
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
        const Eigen::VectorXd gradient=sparse_z.transpose()*(rhs-sparse_z*b);
        const double tolerance{(policy ? policy->release_factor : 128)*eps*std::max(1.0,znorm*(original_rhs_norm+znorm*b.norm()))};
        Eigen::Index release{-1}; double largest{tolerance};
        for (Eigen::Index k=0;k<x.cols();k+=2) if (!free[static_cast<std::size_t>(k)] && gradient(k)>largest)
        {release=k; largest=gradient(k);}
        if (release>=0) {free[static_cast<std::size_t>(release)]=true; ++out.releases; continue;}
        out.beta=b.cwiseQuotient(scales);
        {
            constexpr Eigen::Index maximum_dense_elements=1'000'000; // Cap the dense face matrix at 8 MB.
            const bool unit_weights=(weights.array()==1.0).all();
            const bool bounded=!columns.empty() && x.rows()<=maximum_dense_elements/
                static_cast<Eigen::Index>(columns.size());
            const double observation_scale=observation_scale_override>0 ? observation_scale_override :
                std::max(1.0,original_rhs_norm);
            const Certificate initial=unit_weights && bounded ?
                CertifyLinear(x,y,out.beta,observation_scale) : Certificate{};
            if(initial.available)
            {
                ++out.refinement_eligibility_count;
                out.initial_projected_kkt=initial.projected_kkt;
                out.final_projected_kkt=initial.projected_kkt;
            }
            if(initial.available && !initial.kkt_passed)
            {
                ++out.refinement_trigger_count;
                out.refinement_before_objective=initial.objective;
                out.refinement_after_objective=initial.objective;
                const auto refinement_started=std::chrono::steady_clock::now();
                bool refinement_succeeded{};
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
                ++out.refinement_factorizations;
                if(refinement_qr.rank()==static_cast<Eigen::Index>(columns.size()))
                {
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
                        ++out.solves; ++out.refinement_iterations;
                        if(!correction.allFinite()) break;
                        theta+=correction;
                        Eigen::VectorXd refined_beta=Eigen::VectorXd::Zero(x.cols());
                        for(std::size_t k=0;k<columns.size();++k)
                            refined_beta(columns[k])=theta(static_cast<Eigen::Index>(k))/scales(columns[k]);
                        const Certificate refined=CertifyLinear(x,y,refined_beta,observation_scale);
                        if(refined.available && refined.kkt_passed && refined.objective<=initial.objective)
                        {
                            out.beta=std::move(refined_beta);
                            out.final_projected_kkt=refined.projected_kkt;
                            out.refinement_after_objective=refined.objective;
                            ++out.refinement_success_count;
                            refinement_succeeded=true;
                            break;
                        }
                    }
                }
                if(!refinement_succeeded) ++out.refinement_failure_count;
                out.refinement_seconds=std::chrono::duration<double>(
                    std::chrono::steady_clock::now()-refinement_started).count();
            }
        }
        out.valid=true; out.reason="solved"; return out;
    }
    out.reason="active-set-budget-exhausted"; return out;
}


} // namespace
LinearResult SolveLinear(const Sparse & x,VectorRef y,const Vector & w,
    const LinearPolicy * policy,const std::vector<LinearBlock> * blocks,double observation_scale)
{return WeightedSolveImpl(x,y,w,policy,blocks,observation_scale);}
}
