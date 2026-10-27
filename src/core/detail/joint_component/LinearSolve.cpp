#include "Numerics.hpp"
#include "LinearSolve.hpp"
#include <Eigen/SparseQR>
#include <Eigen/QR>
#include <algorithm>
#include <chrono>
#include <cmath>

namespace rhbm_gem::core::joint_component {
namespace {
constexpr double eps=std::numeric_limits<double>::epsilon();
struct ScopeSeconds
{
    double & result;
    std::chrono::steady_clock::time_point started{std::chrono::steady_clock::now()};
    ~ScopeSeconds()
    {result+=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();}
};
bool Finite(const Sparse & x)
{
    for (int k=0;k<x.outerSize();++k) for (Sparse::InnerIterator e(x,k);e;++e)
        if (!std::isfinite(e.value())) return false;
    return true;
}
void ColumnNorms(const Sparse & x,Eigen::VectorXd & out)
{
    out.resize(x.cols()); out.setZero();
    for(int k=0;k<x.outerSize();++k) for(Sparse::InnerIterator e(x,k);e;++e) out(e.col())+=e.value()*e.value();
    out=out.cwiseSqrt();
}
}
namespace {
void ReduceSparseRows(const Eigen::SparseMatrix<double> & a,const Eigen::VectorXd & rhs,
    LinearSolveWorkspace & workspace,LinearResult & work)
{
    // Orthogonal elimination within contiguous row tiles. Only structurally
    // present columns need local QR; all original rows and RHS entries participate.
    const auto conversion_started=std::chrono::steady_clock::now();
    workspace.rows=a;
    work.row_major_conversion_seconds+=std::chrono::duration<double>(
        std::chrono::steady_clock::now()-conversion_started).count();
    workspace.entries.clear(); workspace.response.clear();
    constexpr Eigen::Index tile_size{1024};
    for(Eigen::Index first=0;first<a.rows();first+=tile_size)
    {
        const Eigen::Index count{std::min(tile_size,a.rows()-first)};
        const auto mapping_started=std::chrono::steady_clock::now();
        workspace.local_column.resize(static_cast<std::size_t>(a.cols()));
        std::fill(workspace.local_column.begin(),workspace.local_column.end(),-1);
        workspace.tile_columns.clear();
        for(Eigen::Index row=first;row<first+count;++row)
            for(LinearSolveWorkspace::RowMajorSparse::InnerIterator e(workspace.rows,row);e;++e)
                workspace.local_column[static_cast<std::size_t>(e.col())]=0;
        for(Eigen::Index k=0;k<a.cols();++k) if(workspace.local_column[static_cast<std::size_t>(k)]==0)
        {
            workspace.local_column[static_cast<std::size_t>(k)]=static_cast<int>(workspace.tile_columns.size());
            workspace.tile_columns.push_back(static_cast<int>(k));
        }
        ++work.tile_mapping_build_count;
        work.tile_mapping_build_seconds+=std::chrono::duration<double>(
            std::chrono::steady_clock::now()-mapping_started).count();
        if(workspace.tile_columns.empty()) continue;
        workspace.local.resize(count,static_cast<Eigen::Index>(workspace.tile_columns.size()));
        workspace.local.setZero();
        for(Eigen::Index row=first;row<first+count;++row)
            for(LinearSolveWorkspace::RowMajorSparse::InnerIterator e(workspace.rows,row);e;++e)
                workspace.local(row-first,workspace.local_column[static_cast<std::size_t>(e.col())])=e.value();
        const auto qr_started=std::chrono::steady_clock::now();
        const Eigen::HouseholderQR<Eigen::MatrixXd> qr(workspace.local);
        workspace.transformed.resize(count);
        workspace.transformed=qr.householderQ().adjoint()*rhs.segment(first,count);
        work.tiled_qr_seconds+=std::chrono::duration<double>(
            std::chrono::steady_clock::now()-qr_started).count();
        const auto assembly_started=std::chrono::steady_clock::now();
        const Eigen::Index retained{std::min(count,workspace.local.cols())},
            offset{static_cast<Eigen::Index>(workspace.response.size())};
        for(Eigen::Index row=0;row<retained;++row)
        {
            workspace.response.push_back(workspace.transformed(row));
            for(Eigen::Index col=row;col<workspace.local.cols();++col)
                if(qr.matrixQR()(row,col)!=0)
                    workspace.entries.emplace_back(offset+row,workspace.tile_columns[static_cast<std::size_t>(col)],
                        qr.matrixQR()(row,col));
        }
        work.reduced_sparse_assembly_seconds+=std::chrono::duration<double>(
            std::chrono::steady_clock::now()-assembly_started).count();
    }
    const auto assembly_started=std::chrono::steady_clock::now();
    workspace.reduced.resize(static_cast<Eigen::Index>(workspace.response.size()),a.cols());
    workspace.reduced.setFromTriplets(workspace.entries.begin(),workspace.entries.end());
    work.reduced_nonzeros+=static_cast<std::size_t>(workspace.reduced.nonZeros());
    workspace.reduced_rhs.resize(static_cast<Eigen::Index>(workspace.response.size()));
    for(std::size_t row=0;row<workspace.response.size();++row)
        workspace.reduced_rhs(static_cast<Eigen::Index>(row))=workspace.response[row];
    work.reduced_sparse_assembly_seconds+=std::chrono::duration<double>(
        std::chrono::steady_clock::now()-assembly_started).count();
}
LinearResult UnweightedSolveImpl(const Sparse & x,VectorRef y,
    const LinearPolicy * policy,double observation_scale_override,LinearSolveWorkspace * workspace_ptr)
{
    LinearResult out; out.reason="invalid-input"; out.beta=Eigen::VectorXd::Zero(x.cols());
    ScopeSeconds total_timer{out.solver_total_seconds};
    if (x.cols()==0 || x.cols()%2!=0 || x.rows()<=x.cols() || y.size()!=x.rows() ||
        !Finite(x) || !y.allFinite()) return out;
    LinearSolveWorkspace local_workspace;
    auto & workspace=workspace_ptr ? *workspace_ptr : local_workspace;
    ColumnNorms(x,workspace.scales);
    const auto & scales=workspace.scales;
    if ((scales.array()==0).any()) {out.reason="rank-deficient"; return out;}
    auto & sparse_z=workspace.normalized; sparse_z=x;
    for (int k=0;k<sparse_z.outerSize();++k)
        for (Sparse::InnerIterator entry(sparse_z,k);entry;++entry)
            entry.valueRef()/=scales(k);
    auto & rhs=workspace.rhs; rhs=y;
    const double original_rhs_norm=policy && policy->release_response_norm>=0 ? policy->release_response_norm : rhs.norm();
    const double znorm{sparse_z.norm()};
    const double rank_threshold{policy ? policy->rank_relative : eps*static_cast<double>(std::max(x.rows(),x.cols()))};
    auto & b=workspace.b; b.resize(x.cols()); b.setZero();
    auto & free=workspace.free; free.resize(static_cast<std::size_t>(x.cols()));
    std::fill(free.begin(),free.end(),true);
    auto & previous_columns=workspace.previous_columns;
    previous_columns.clear(); previous_columns.reserve(static_cast<std::size_t>(x.cols()));
    auto & columns=workspace.columns; columns.clear(); columns.reserve(static_cast<std::size_t>(x.cols()));
    auto & candidate=workspace.candidate; candidate.resize(x.cols());
    auto & gradient=workspace.gradient; gradient.resize(x.cols());
    auto & solution=workspace.solution; solution.resize(x.cols());
    auto & qr=workspace.qr;
    bool has_previous_face{};
    // Begin with the unconstrained problem. Block infeasible A coefficients and
    // release blocked ones when the dual gradient requires it (mixed NNLS).
    for (Eigen::Index iteration=0; iteration<(policy ? policy->active_set_iteration_factor : 20)*x.cols()*x.cols(); ++iteration)
    {
        ++out.active_set_iterations;
        columns.clear();
        for (Eigen::Index k=0;k<x.cols();++k) if (free[static_cast<std::size_t>(k)]) columns.push_back(k);
        if(has_previous_face)
        {
            if(columns==previous_columns) ++out.free_face_unchanged;
            else ++out.free_face_changes;
        }
        previous_columns=columns; has_previous_face=true;
        const auto structure_started=std::chrono::steady_clock::now();
        auto & selected=workspace.selected;
        selected.resize(x.rows(),static_cast<Eigen::Index>(columns.size()));
        selected.setZero();
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
        ++out.sparse_structure_build_count;
        out.selected_nonzeros+=static_cast<std::size_t>(selected.nonZeros());
        out.sparse_structure_build_seconds+=std::chrono::duration<double>(
            std::chrono::steady_clock::now()-structure_started).count();
        if (maximum_norm==0) {out.reason="rank-deficient"; return out;}
        // SparseQR uses an absolute pivot threshold. The scale is the
        // largest norm of the column-normalized design.
        ReduceSparseRows(selected,rhs,workspace,out);
        const auto sparse_qr_started=std::chrono::steady_clock::now();
        qr.setPivotThreshold(rank_threshold*maximum_norm); qr.compute(workspace.reduced);
        if(qr.info()!=Eigen::Success)
        {
            out.sparse_qr_seconds+=std::chrono::duration<double>(
                std::chrono::steady_clock::now()-sparse_qr_started).count();
            out.reason="nonfinite"; return out;
        }
        out.rank=static_cast<int>(qr.rank());
        solution.resize(static_cast<Eigen::Index>(columns.size()));
        solution=qr.solve(workspace.reduced_rhs);
        ++out.solves;
        out.sparse_qr_seconds+=std::chrono::duration<double>(
            std::chrono::steady_clock::now()-sparse_qr_started).count();
        if (out.rank!=static_cast<int>(columns.size())) {out.reason="rank-deficient"; return out;}
        candidate.setZero();
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
        gradient=sparse_z.transpose()*(rhs-sparse_z*b);
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
LinearResult SolveLinear(const Sparse & x,VectorRef y,const LinearPolicy * policy,double observation_scale,
    LinearSolveWorkspace * workspace)
{return UnweightedSolveImpl(x,y,policy,observation_scale,workspace);}
}
