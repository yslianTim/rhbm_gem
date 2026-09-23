#include "SparseFactor.hpp"
#include <chrono>
#include <stdexcept>
#ifdef RHBM_GEM_JOINT_SPQR
#include <SuiteSparseQR.hpp>
#endif

namespace rhbm_gem::core::joint_component {
SparseWork & SparseWorkForTesting() {static thread_local SparseWork work; return work;}
#ifdef RHBM_GEM_JOINT_SPQR
namespace {
using Clock=std::chrono::steady_clock;
double Seconds(Clock::time_point t) {return std::chrono::duration<double>(Clock::now()-t).count();}
using LongSparse=Eigen::SparseMatrix<double,Eigen::ColMajor,int64_t>;
cholmod_sparse View(LongSparse & a)
{
    cholmod_sparse v{}; v.nrow=static_cast<std::size_t>(a.rows()); v.ncol=static_cast<std::size_t>(a.cols());
    v.nzmax=static_cast<std::size_t>(a.nonZeros()); v.p=a.outerIndexPtr(); v.i=a.innerIndexPtr(); v.x=a.valuePtr();
    v.itype=CHOLMOD_LONG; v.xtype=CHOLMOD_REAL; v.dtype=CHOLMOD_DOUBLE; v.sorted=1; v.packed=1; return v;
}
cholmod_dense View(const Matrix & a)
{
    cholmod_dense v{}; v.nrow=static_cast<std::size_t>(a.rows()); v.ncol=static_cast<std::size_t>(a.cols());
    v.nzmax=static_cast<std::size_t>(a.size()); v.d=v.nrow; v.x=const_cast<double *>(a.data());
    v.xtype=CHOLMOD_REAL; v.dtype=CHOLMOD_DOUBLE; return v;
}
struct Dense
{
    cholmod_dense * p{}; cholmod_common * cc;
    ~Dense() {if(p) cholmod_l_free_dense(&p,cc);}
    Matrix Copy() const
    {
        if(!p) throw std::runtime_error("SPQR dense operation failed");
        return Eigen::Map<const Matrix,0,Eigen::OuterStride<>>(static_cast<const double *>(p->x),
            static_cast<Eigen::Index>(p->nrow),static_cast<Eigen::Index>(p->ncol),Eigen::OuterStride<>(static_cast<Eigen::Index>(p->d)));
    }
};
bool Pattern(const Sparse & a,const Sparse & b)
{
    return a.rows()==b.rows() && a.cols()==b.cols() && a.nonZeros()==b.nonZeros() &&
        std::equal(a.outerIndexPtr(),a.outerIndexPtr()+a.cols()+1,b.outerIndexPtr()) &&
        std::equal(a.innerIndexPtr(),a.innerIndexPtr()+a.nonZeros(),b.innerIndexPtr());
}
}
struct SparseFactorState
{
    cholmod_common cc{};
    SuiteSparseQR_factorization<double> * qr{};
    Sparse design;
    std::vector<Eigen::Index> columns;
    std::size_t generation{};
    const void * domain{}; const void * observations{};
    LinearPolicy policy{};
    double tolerance{};
    bool bound{},policy_present{};
    SparseFactorState() {cholmod_l_start(&cc); cc.SPQR_nthreads=1;}
    ~SparseFactorState() {if(qr) SuiteSparseQR_free<double>(&qr,&cc); cholmod_l_finish(&cc);}
    void Clear() {++generation; if(qr) SuiteSparseQR_free<double>(&qr,&cc);}
};
bool SparseBackendEnabled() {return true;}
LinearWorkspace::LinearWorkspace():state_(std::make_shared<SparseFactorState>()) {}
void LinearWorkspace::Bind(const void * domain,const void * observations,const LinearPolicy * policy)
{
    const auto p=policy ? *policy : LinearPolicy{};
    auto & s=*state_;
    if(s.bound && (s.domain!=domain || s.observations!=observations || s.policy_present!=(policy!=nullptr) ||
        s.policy.rank_relative!=p.rank_relative || s.policy.release_factor!=p.release_factor ||
        s.policy.release_response_norm!=p.release_response_norm || s.policy.active_set_iteration_factor!=p.active_set_iteration_factor)) s.Clear();
    s.domain=domain; s.observations=observations; s.policy=p; s.policy_present=policy!=nullptr; s.bound=true;
}
std::shared_ptr<FreeDesignFactor> LinearWorkspace::Factor(const Sparse & a,const std::vector<Eigen::Index> & columns,double tolerance)
{
    auto & s=*state_; ++s.generation;
    const auto conversion=Clock::now();
    LongSparse storage=a; storage.makeCompressed(); auto view=View(storage);
    SparseWorkForTesting().matrix_preparation_seconds+=Seconds(conversion);
    if(!s.qr || s.tolerance!=tolerance || s.columns!=columns || !Pattern(s.design,a))
    {
        s.Clear(); const auto start=Clock::now();
        s.qr=SuiteSparseQR_symbolic<double>(SPQR_ORDERING_COLAMD,true,&view,&s.cc);
        ++SparseWorkForTesting().symbolic; SparseWorkForTesting().symbolic_seconds+=Seconds(start);
        if(!s.qr) throw std::runtime_error("SPQR symbolic analysis failed");
    }
    else ++SparseWorkForTesting().symbolic_reuses;
    s.design=a; s.columns=columns; s.tolerance=tolerance;
    const auto start=Clock::now(); const auto ok=SuiteSparseQR_numeric<double>(tolerance,&view,s.qr,&s.cc);
    ++SparseWorkForTesting().numeric; SparseWorkForTesting().numeric_seconds+=Seconds(start);
    if(!ok) {s.Clear(); throw std::runtime_error("SPQR numeric factorization failed");}
    SparseWorkForTesting().factor_nonzeros=std::max(SparseWorkForTesting().factor_nonzeros,static_cast<std::size_t>(std::max<int64_t>(0,s.cc.SPQR_istat[0])));
    return std::shared_ptr<FreeDesignFactor>(new FreeDesignFactor(state_,s.generation));
}
FreeDesignFactor::FreeDesignFactor(std::shared_ptr<SparseFactorState> state,std::size_t generation):state_(std::move(state)),generation_(generation) {}
void FreeDesignFactor::Check() const
{if(!state_ || state_->generation!=generation_ || !state_->qr) throw std::logic_error("Expired free-design factor");}
bool FreeDesignFactor::Matches(const Sparse & a,const std::vector<Eigen::Index> & columns) const
{
    return state_ && state_->generation==generation_ && state_->qr && state_->columns==columns && Pattern(state_->design,a) &&
        std::equal(a.valuePtr(),a.valuePtr()+a.nonZeros(),state_->design.valuePtr());
}
int FreeDesignFactor::Rank() const {Check(); return static_cast<int>(state_->qr->rank);}
Matrix FreeDesignFactor::LeastSquares(const Matrix & rhs) const
{
    Check(); auto b=View(rhs); auto & s=*state_;
    Dense transformed{SuiteSparseQR_qmult<double>(SPQR_QTX,s.qr,&b,&s.cc),&s.cc};
    if(!transformed.p) throw std::runtime_error("SPQR Q transpose failed");
    Dense answer{SuiteSparseQR_solve<double>(SPQR_RETX_EQUALS_B,s.qr,transformed.p,&s.cc),&s.cc}; return answer.Copy();
}
Matrix FreeDesignFactor::NormalSolve(const Matrix & rhs) const
{
    Check(); auto b=View(rhs); auto & s=*state_;
    Dense adjoint{SuiteSparseQR_solve<double>(SPQR_RTX_EQUALS_ETB,s.qr,&b,&s.cc),&s.cc};
    if(!adjoint.p) throw std::runtime_error("SPQR adjoint triangular solve failed");
    Dense answer{SuiteSparseQR_solve<double>(SPQR_RETX_EQUALS_B,s.qr,adjoint.p,&s.cc),&s.cc}; return answer.Copy();
}
Matrix FreeDesignFactor::Compact() const
{
    Check(); auto & s=*state_; const auto p=s.design.cols(); Matrix compact(p,p);
    // Q^T X retains original column order; its singular values equal those of R.
    // Never materialize the observation-sized Q or all RHS columns together.
    for(Eigen::Index first=0;first<p;first+=16)
    {
        const auto count=std::min<Eigen::Index>(16,p-first); Matrix rhs(s.design.middleCols(first,count)); auto b=View(rhs);
        Dense transformed{SuiteSparseQR_qmult<double>(SPQR_QTX,s.qr,&b,&s.cc),&s.cc};
        compact.middleCols(first,count)=transformed.Copy().topRows(p);
    }
    return compact;
}
std::pair<Matrix,Vector> SparseReferenceQR(const Sparse & x,const Vector & weights,const Vector & scales,VectorRef y)
{
    ++SparseWorkForTesting().reference; WorkTimer timer(SparseWorkForTesting().reference_seconds);
    SparseFactorState state;
    LongSparse storage=x;
    for(Eigen::Index k=0;k<storage.outerSize();++k) for(LongSparse::InnerIterator e(storage,k);e;++e)
        e.valueRef()=(e.value()*std::sqrt(weights(e.row())))/scales(k);
    auto a=View(storage); const Matrix rhs=weights.cwiseSqrt().array()*y.array(); auto b=View(rhs);
    cholmod_sparse * r{}; int64_t * permutation{}; cholmod_dense * c{};
    const auto rank=SuiteSparseQR<double>(SPQR_ORDERING_COLAMD,SPQR_NO_TOL,static_cast<int64_t>(x.cols()),&a,&b,&c,&r,&permutation,&state.cc);
    Dense response{c,&state.cc}; Matrix result=Matrix::Zero(x.cols(),x.cols());
    if(rank>=0 && r)
    {
        const auto * outer=static_cast<const int64_t *>(r->p),* inner=static_cast<const int64_t *>(r->i);
        const auto * values=static_cast<const double *>(r->x);
        for(Eigen::Index k=0;k<x.cols();++k) for(auto i=outer[k];i<outer[k+1];++i)
            if(inner[i]<x.cols()) result(inner[i],permutation ? permutation[k] : k)=values[i];
    }
    if(r) cholmod_l_free_sparse(&r,&state.cc);
    if(permutation) cholmod_l_free(static_cast<std::size_t>(x.cols()),sizeof(int64_t),permutation,&state.cc);
    if(rank<0 || !c) throw std::runtime_error("Independent SPQR reduction failed");
    return {std::move(result),response.Copy().col(0)};
}
#else
struct SparseFactorState {};
bool SparseBackendEnabled() {return false;}
LinearWorkspace::LinearWorkspace()=default;
void LinearWorkspace::Bind(const void *,const void *,const LinearPolicy *) {}
std::shared_ptr<FreeDesignFactor> LinearWorkspace::Factor(const Sparse &,const std::vector<Eigen::Index> &,double)
{throw std::logic_error("SPQR backend is disabled");}
FreeDesignFactor::FreeDesignFactor(std::shared_ptr<SparseFactorState> s,std::size_t g):state_(std::move(s)),generation_(g) {}
void FreeDesignFactor::Check() const {throw std::logic_error("SPQR backend is disabled");}
bool FreeDesignFactor::Matches(const Sparse &,const std::vector<Eigen::Index> &) const {return false;}
int FreeDesignFactor::Rank() const {Check(); return 0;}
Matrix FreeDesignFactor::Compact() const {Check(); return {};}
Matrix FreeDesignFactor::LeastSquares(const Matrix &) const {Check(); return {};}
Matrix FreeDesignFactor::NormalSolve(const Matrix &) const {Check(); return {};}
std::pair<Matrix,Vector> SparseReferenceQR(const Sparse &,const Vector &,const Vector &,VectorRef)
{throw std::logic_error("SPQR backend is disabled");}
#endif
}
