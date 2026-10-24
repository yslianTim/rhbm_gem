#include "Numerics.hpp"
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
    // present columns need local QR; all original rows and RHS entries participate.
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
LinearResult UnweightedSolveImpl(const Sparse & x,VectorRef y,
    const LinearPolicy * policy,double observation_scale_override)
{
    LinearResult out; out.reason="invalid-input"; out.beta=Eigen::VectorXd::Zero(x.cols());
    if (x.cols()==0 || x.cols()%2!=0 || x.rows()<=x.cols() || y.size()!=x.rows() ||
        !Finite(x) || !y.allFinite()) return out;
    const Eigen::VectorXd scales{ColumnNorms(x)};
    if ((scales.array()==0).any()) {out.reason="rank-deficient"; return out;}
    Sparse sparse_z=x;
    for (int k=0;k<sparse_z.outerSize();++k)
        for (Sparse::InnerIterator entry(sparse_z,k);entry;++entry)
            entry.valueRef()/=scales(k);
    Eigen::VectorXd rhs{y};
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
        // largest norm of the column-normalized design.
        const auto reduced{ReduceSparseRows(selected,rhs)};
        qr.setPivotThreshold(rank_threshold*maximum_norm); qr.compute(reduced.first);
        if(qr.info()!=Eigen::Success) {out.reason="nonfinite"; return out;}
        out.rank=static_cast<int>(qr.rank()); const Eigen::VectorXd solution=qr.solve(reduced.second);
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
            const bool bounded=!columns.empty() && x.rows()<=maximum_dense_elements/
                static_cast<Eigen::Index>(columns.size());
            const double observation_scale=observation_scale_override>0 ? observation_scale_override :
                std::max(1.0,original_rhs_norm);
            const Certificate initial=bounded ?
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
LinearResult SolveLinear(const Sparse & x,VectorRef y,const LinearPolicy * policy,double observation_scale)
{return UnweightedSolveImpl(x,y,policy,observation_scale);}
}
