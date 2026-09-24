#include "SparseFactor.hpp"
#include <chrono>
#include <stdexcept>
#include <Eigen/SparseQR>
#ifdef RHBM_GEM_JOINT_SPQR
#include <SuiteSparseQR.hpp>
#endif

namespace rhbm_gem::core::joint_component {
SparseWork & SparseWorkForTesting() {static thread_local SparseWork work; return work;}
ResourceWork & ResourceWorkForTesting() {static thread_local ResourceWork work; return work;}
void RecordDenseShape(const char * role,Eigen::Index rows,Eigen::Index columns)
{
    auto & w=ResourceWorkForTesting(); if(!w.enabled) return;
    const auto bytes=static_cast<std::size_t>(rows)*static_cast<std::size_t>(columns)*sizeof(double);
    for(auto & record:w.dense_shapes) if(record.phase==w.phase && record.role==role)
    {
        ++record.observations;
        if(bytes>record.maximum_bytes) {record.rows=rows; record.columns=columns; record.maximum_bytes=bytes;}
        return;
    }
    w.dense_shapes.push_back({w.phase,role,rows,columns,1,bytes});
}
ResourcePhase::ResourcePhase(const char * name):enabled_(ResourceWorkForTesting().enabled),name_(name)
{
    if(enabled_) {auto & w=ResourceWorkForTesting(); previous_=w.phase; w.phase=name; started_=std::chrono::steady_clock::now();}
}
ResourcePhase::~ResourcePhase()
{
    if(!enabled_) return;
    auto & w=ResourceWorkForTesting(); w.phase=previous_;
    const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-started_).count();
    for(auto & record:w.phases) if(record.phase==name_) {++record.calls; record.inclusive_seconds+=seconds; return;}
    w.phases.push_back({name_,1,seconds});
}
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
Matrix FreeDesignFactor::PseudoInverseTranspose(const Matrix & rhs) const
{
    Check(); auto & s=*state_;
    if(rhs.rows()!=s.design.cols()) throw std::invalid_argument("Invalid adjoint RHS size");
    auto b=View(rhs);
    Dense triangular{SuiteSparseQR_solve<double>(SPQR_RTX_EQUALS_ETB,s.qr,&b,&s.cc),&s.cc};
    Matrix padded=Matrix::Zero(s.design.rows(),rhs.cols());
    padded.topRows(s.design.cols())=triangular.Copy().topRows(s.design.cols());
    auto v=View(padded);
    Dense answer{SuiteSparseQR_qmult<double>(SPQR_QX,s.qr,&v,&s.cc),&s.cc};
    return answer.Copy();
}
Matrix FreeDesignFactor::ProjectComplement(const Matrix & rhs) const
{
    Check(); auto & s=*state_;
    if(rhs.rows()!=s.design.rows()) throw std::invalid_argument("Invalid projection RHS size");
    auto b=View(rhs);
    Dense transformed{SuiteSparseQR_qmult<double>(SPQR_QTX,s.qr,&b,&s.cc),&s.cc};
    Matrix tail=transformed.Copy(); tail.topRows(s.design.cols()).setZero(); auto v=View(tail);
    Dense answer{SuiteSparseQR_qmult<double>(SPQR_QX,s.qr,&v,&s.cc),&s.cc};
    return answer.Copy();
}
Matrix FreeDesignFactor::Compact() const
{
    Check(); auto & s=*state_; const auto p=s.design.cols(); Matrix compact(p,p);
    RecordDenseShape("free-design-compact",p,p);
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
    RecordDenseShape("reference-compact",x.cols(),x.cols());
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
// This factor is used only by the fixed-state operator. The production Eigen
// active-set solver retains its existing row reduction and rank policy.
struct SparseFactorState
{
    Eigen::SparseQR<Sparse,Eigen::COLAMDOrdering<int>> qr;
    Sparse design;
    std::vector<Eigen::Index> columns;
    std::size_t generation{};
    bool valid{};
};
bool SparseBackendEnabled() {return false;}
LinearWorkspace::LinearWorkspace()=default;
void LinearWorkspace::Bind(const void *,const void *,const LinearPolicy *) {}
std::shared_ptr<FreeDesignFactor> LinearWorkspace::Factor(const Sparse & a,const std::vector<Eigen::Index> & columns,double tolerance)
{
    if(!state_) state_=std::make_shared<SparseFactorState>();
    auto & s=*state_; ++s.generation; s.valid=false; s.design=a; s.columns=columns;
    s.qr.setPivotThreshold(tolerance);
    auto & work=SparseWorkForTesting();
    {WorkTimer timer(work.symbolic_seconds); s.qr.analyzePattern(a); ++work.symbolic;}
    {WorkTimer timer(work.numeric_seconds); s.qr.factorize(a); ++work.numeric;}
    if(s.qr.info()!=Eigen::Success) throw std::runtime_error("Eigen operator factorization failed");
    s.valid=true;
    work.factor_nonzeros=std::max(work.factor_nonzeros,static_cast<std::size_t>(s.qr.matrixR().nonZeros()));
    return std::shared_ptr<FreeDesignFactor>(new FreeDesignFactor(state_,s.generation));
}
FreeDesignFactor::FreeDesignFactor(std::shared_ptr<SparseFactorState> s,std::size_t g):state_(std::move(s)),generation_(g) {}
void FreeDesignFactor::Check() const
{if(!state_ || !state_->valid || generation_!=state_->generation) throw std::logic_error("Expired free-design factor");}
bool FreeDesignFactor::Matches(const Sparse & a,const std::vector<Eigen::Index> & columns) const
{
    return state_ && state_->valid && generation_==state_->generation && columns==state_->columns &&
        a.rows()==state_->design.rows() && a.cols()==state_->design.cols() && (a-state_->design).norm()==0;
}
int FreeDesignFactor::Rank() const {Check(); return static_cast<int>(state_->qr.rank());}
Matrix FreeDesignFactor::Compact() const
{
    Check(); const auto & s=*state_; const auto p=s.design.cols();
    RecordDenseShape("free-design-compact",p,p);
    // Z*P=Q*R. Restore the original column order without repeating N-by-p
    // orthogonal actions; this is still the same transient compact rank check.
    return Matrix(s.qr.matrixR().topRows(p))*s.qr.colsPermutation().transpose();
}
Matrix FreeDesignFactor::LeastSquares(const Matrix & rhs) const
{Check(); if(rhs.rows()!=state_->design.rows()) throw std::invalid_argument("Invalid least-squares RHS size"); return state_->qr.solve(rhs);}
Matrix FreeDesignFactor::PseudoInverseTranspose(const Matrix & rhs) const
{
    Check(); const auto & s=*state_; const auto p=s.design.cols();
    if(rhs.rows()!=p) throw std::invalid_argument("Invalid adjoint RHS size");
    Matrix padded=Matrix::Zero(s.design.rows(),rhs.cols());
    const Matrix permuted=s.qr.colsPermutation().transpose()*rhs;
    padded.topRows(p)=s.qr.matrixR().topLeftCorner(p,p).transpose().triangularView<Eigen::Lower>().solve(permuted);
    return s.qr.matrixQ()*padded;
}
Matrix FreeDesignFactor::ProjectComplement(const Matrix & rhs) const
{
    Check(); const auto & s=*state_;
    if(rhs.rows()!=s.design.rows()) throw std::invalid_argument("Invalid projection RHS size");
    Matrix tail=s.qr.matrixQ().adjoint()*rhs; tail.topRows(s.design.cols()).setZero();
    return s.qr.matrixQ()*tail;
}
Matrix FreeDesignFactor::NormalSolve(const Matrix & rhs) const
{return LeastSquares(PseudoInverseTranspose(rhs));}
std::pair<Matrix,Vector> SparseReferenceQR(const Sparse &,const Vector &,const Vector &,VectorRef)
{throw std::logic_error("SPQR backend is disabled");}
#endif
}
