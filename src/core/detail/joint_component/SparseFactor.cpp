#include "SparseFactor.hpp"
#include <algorithm>
#include <chrono>
#include <stdexcept>
#include <cstring>
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
#include <sys/resource.h>
#endif
#include <SuiteSparseQR.hpp>

namespace rhbm_gem::core::joint_component {
SparseWork & SparseWorkForTesting() {static thread_local SparseWork work; return work;}
ResourceWork & ResourceWorkForTesting() {static thread_local ResourceWork work; return work;}
ProfileEvaluationRole & ProfileEvaluationRoleForTesting()
{static thread_local auto role=ProfileEvaluationRole::Unspecified; return role;}
ProfileEvaluationRoleScopeForTesting::ProfileEvaluationRoleScopeForTesting(ProfileEvaluationRole role)
    :previous_(ProfileEvaluationRoleForTesting())
{ProfileEvaluationRoleForTesting()=role;}
ProfileEvaluationRoleScopeForTesting::~ProfileEvaluationRoleScopeForTesting()
{ProfileEvaluationRoleForTesting()=previous_;}
bool & LinearTelemetryEnabledForTesting()
{static thread_local bool enabled{}; return enabled;}
LinearTelemetryScopeForTesting::LinearTelemetryScopeForTesting(bool enabled)
    :previous_(LinearTelemetryEnabledForTesting())
{LinearTelemetryEnabledForTesting()=enabled;}
LinearTelemetryScopeForTesting::~LinearTelemetryScopeForTesting()
{LinearTelemetryEnabledForTesting()=previous_;}
namespace {
bool SameSparsePattern(const Sparse & a,const Sparse & b)
{
    return a.rows()==b.rows() && a.cols()==b.cols() && a.nonZeros()==b.nonZeros() &&
        std::equal(a.outerIndexPtr(),a.outerIndexPtr()+a.cols()+1,b.outerIndexPtr()) &&
        std::equal(a.innerIndexPtr(),a.innerIndexPtr()+a.nonZeros(),b.innerIndexPtr());
}
bool SameSparseValues(const Sparse & a,const Sparse & b)
{
    return SameSparsePattern(a,b) && std::equal(a.valuePtr(),a.valuePtr()+a.nonZeros(),b.valuePtr());
}
void RecordExactReuseOpportunity(SparseWork & work)
{
    ++work.numeric_factor_exact_reuse_opportunities;
    switch(ProfileEvaluationRoleForTesting())
    {
    case ProfileEvaluationRole::InitialProfile: ++work.initial_profile_exact_reuse_opportunities; break;
    case ProfileEvaluationRole::TrialProfile: ++work.trial_profile_exact_reuse_opportunities; break;
    case ProfileEvaluationRole::Reference: ++work.reference_exact_reuse_opportunities; break;
    case ProfileEvaluationRole::AcceptedEndpoint: ++work.accepted_endpoint_exact_reuse_opportunities; break;
    case ProfileEvaluationRole::Unspecified: break;
    }
}
}
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
FactorResidencyWork & FactorResidencyWorkForTesting() {static thread_local FactorResidencyWork work; return work;}
void ResetFactorResidencyWorkForTesting() {FactorResidencyWorkForTesting()=FactorResidencyWork{};}
std::string & FactorCreationRoleForTesting() {static thread_local std::string role; return role;}
FactorCreationRoleScopeForTesting::FactorCreationRoleScopeForTesting(std::string role)
    :previous_(FactorCreationRoleForTesting())
{FactorCreationRoleForTesting()=std::move(role);}
FactorCreationRoleScopeForTesting::~FactorCreationRoleScopeForTesting()
{FactorCreationRoleForTesting()=std::move(previous_);}
namespace {
std::string FactorStageForTesting()
{
    return ResourceWorkForTesting().phase;
}
std::string FactorSearchStageForTesting()
{
    return ResourceWorkForTesting().active_search_stage;
}
double FactorElapsedForTesting()
{
    return std::chrono::duration<double>(std::chrono::steady_clock::now()-FactorResidencyWorkForTesting().started).count();
}
std::size_t CurrentFactorBytesForTesting(const FactorResidencyWork & work)
{
    std::size_t bytes{};
    for(const auto & factor:work.factors) if(factor.alive) bytes+=factor.owned_factor_bytes;
    for(const auto & factor:work.constructing) bytes+=factor.owned_factor_bytes_estimate;
    return bytes;
}
std::size_t CurrentFactorCountForTesting(const FactorResidencyWork & work)
{
    return static_cast<std::size_t>(std::count_if(work.factors.begin(),work.factors.end(),
        [](const auto & factor){return factor.alive;}))+work.constructing.size();
}
void RecordFactorRssPeakForTesting(FactorResidencyWork & work)
{
    for(auto & factor:work.constructing) factor.owned_factor_bytes_estimate=factor.estimate_owned_bytes();
    const auto count=CurrentFactorCountForTesting(work), owned_bytes=CurrentFactorBytesForTesting(work);
    const std::string stage=FactorStageForTesting();
    if(count>work.maximum_concurrent_factor_count)
    {work.maximum_concurrent_factor_count=count; work.maximum_concurrent_factor_stage=stage;}
    if(owned_bytes>work.maximum_concurrent_owned_bytes)
    {work.maximum_concurrent_owned_bytes=owned_bytes; work.maximum_concurrent_owned_bytes_stage=stage;}
    rusage usage{}; if(getrusage(RUSAGE_SELF,&usage)!=0) return;
#ifdef __APPLE__
    const auto bytes=static_cast<std::size_t>(usage.ru_maxrss);
#else
    const auto bytes=static_cast<std::size_t>(usage.ru_maxrss)*1024;
#endif
    if(bytes<=work.peak_rss_bytes_observed) return;
    work.peak_rss_bytes_observed=bytes;
    work.peak_rss_factor_count=count;
    work.peak_rss_owned_bytes=owned_bytes;
    work.peak_rss_stage=stage;
    work.peak_rss_factor_ids.clear();
    work.peak_rss_factor_generations.clear();
    for(const auto & factor:work.factors) if(factor.alive)
    {
        work.peak_rss_factor_ids.push_back(factor.factor_id);
        work.peak_rss_factor_generations.push_back(std::to_string(factor.factor_id)+":"+
            std::to_string(factor.generation));
    }
    for(const auto & factor:work.constructing)
    {
        work.peak_rss_factor_ids.push_back(factor.factor_id);
        work.peak_rss_factor_generations.push_back(std::to_string(factor.factor_id)+":"+
            std::to_string(factor.generation));
    }
}
void StartFactorResidencyForTesting(std::size_t id,std::size_t generation,const char * kind,const std::string & role,
    Eigen::Index rows,Eigen::Index columns,std::size_t nonzeros,std::size_t r_nonzeros,std::size_t h_nonzeros,
    std::size_t owned_bytes)
{
    auto & work=FactorResidencyWorkForTesting();
    const std::string stage=FactorStageForTesting();
    work.factors.push_back({id,generation,kind,role.empty() ? kind : role,stage,{},FactorSearchStageForTesting(),{},rows,columns,nonzeros,
        r_nonzeros,h_nonzeros,owned_bytes,FactorElapsedForTesting(),0,true});
    const auto count=CurrentFactorCountForTesting(work), bytes=CurrentFactorBytesForTesting(work);
    if(count>work.maximum_concurrent_factor_count)
    {
        work.maximum_concurrent_factor_count=count;
        work.maximum_concurrent_factor_stage=stage;
    }
    if(bytes>work.maximum_concurrent_owned_bytes)
    {
        work.maximum_concurrent_owned_bytes=bytes;
        work.maximum_concurrent_owned_bytes_stage=stage;
    }
    RecordFactorRssPeakForTesting(work);
}
void EndFactorResidencyForTesting(std::size_t id,std::size_t generation)
{
    auto & work=FactorResidencyWorkForTesting();
    for(auto it=work.factors.rbegin();it!=work.factors.rend();++it)
    {
        if(it->factor_id!=id || it->generation!=generation || !it->alive) continue;
        RecordFactorRssPeakForTesting(work);
        it->alive=false; it->destroyed_at_stage=FactorStageForTesting();
        it->destroyed_search_stage=FactorSearchStageForTesting(); it->destroyed_seconds=FactorElapsedForTesting();
        return;
    }
}
std::size_t NewFactorIdForTesting()
{return FactorResidencyWorkForTesting().next_factor_id++;}
void BeginFactorConstructionForTesting(std::size_t id,std::size_t generation,const char * kind,const std::string & role,
    Eigen::Index rows,Eigen::Index columns,std::size_t nonzeros,std::function<std::size_t()> estimate)
{
    auto & work=FactorResidencyWorkForTesting();
    work.constructing.push_back({id,generation,kind,role.empty() ? kind : role,FactorStageForTesting(),rows,columns,
        nonzeros,estimate(),std::move(estimate)});
    RecordFactorRssPeakForTesting(work);
}
void EndFactorConstructionForTesting(std::size_t id,std::size_t generation)
{
    auto & constructing=FactorResidencyWorkForTesting().constructing;
    constructing.erase(std::remove_if(constructing.begin(),constructing.end(),[&](const auto & factor)
        {return factor.factor_id==id && factor.generation==generation;}),constructing.end());
}
std::size_t OwnedFactorBytesForTesting(const Sparse & design,
    const Eigen::SparseMatrix<double,Eigen::ColMajor,int64_t> * r,const cholmod_common & cc)
{
    const auto design_bytes=static_cast<std::size_t>(design.nonZeros())*(sizeof(double)+sizeof(int))+
        static_cast<std::size_t>(design.cols()+1)*sizeof(int);
    const auto r_bytes=r ? static_cast<std::size_t>(r->nonZeros())*(sizeof(double)+sizeof(int64_t))+
        static_cast<std::size_t>(r->cols()+1)*sizeof(int64_t) : 0;
    return design_bytes+r_bytes+cc.memory_inuse;
}
}
#endif
namespace {
ResourceStageObserverForTesting & ResourceStageObserver()
{static thread_local ResourceStageObserverForTesting observer{}; return observer;}
void * & ResourceStageObserverContext()
{static thread_local void * context{}; return context;}
SearchStageWork & SearchStage(ResourceWork & work,const char * name)
{
    const auto found=std::find_if(work.search_stages.begin(),work.search_stages.end(),
        [&](const auto & stage){return stage.stage==name;});
    if(found!=work.search_stages.end()) return *found;
    work.search_stages.push_back({name}); return work.search_stages.back();
}
void NotifyResourceStageObserver()
{if(const auto observer=ResourceStageObserver()) observer(ResourceWorkForTesting(),ResourceStageObserverContext());}
}
void SetResourceStageObserverForTesting(ResourceStageObserverForTesting observer,void * context)
{ResourceStageObserver()=observer; ResourceStageObserverContext()=context;}
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
void RecordSparseShape(const char * role,Eigen::Index rows,Eigen::Index columns,std::size_t nonzeros)
{
    auto & w=ResourceWorkForTesting(); if(!w.enabled) return;
    for(auto & record:w.search_stages) if(record.stage==w.active_search_stage)
    {
        record.rows=std::max(record.rows,rows); record.columns=std::max(record.columns,columns);
        record.nonzeros=std::max(record.nonzeros,nonzeros); break;
    }
    // Preserve the phase/role record separately from stage aggregation.
    const std::string phase=w.phase, name=role;
    auto & shapes=w.sparse_shapes;
    for(auto & record:shapes) if(record.phase==phase && record.role==name)
    {
        ++record.observations;
        if(nonzeros>record.maximum_nonzeros) {record.rows=rows; record.columns=columns; record.maximum_nonzeros=nonzeros;}
        return;
    }
    shapes.push_back({phase,name,rows,columns,1,nonzeros});
}
ResourcePhase::ResourcePhase(const char * name,bool search_stage,Eigen::Index rows,Eigen::Index columns,std::size_t nonzeros)
    :enabled_(ResourceWorkForTesting().enabled),name_(name)
{
    if(!enabled_) return;
    auto & w=ResourceWorkForTesting(); previous_=w.phase; w.phase=name;
    search_scope_=std::strcmp(name,"search")==0;
    if(search_scope_) ++w.search_depth;
    search_stage_=search_scope_ || (search_stage && w.search_depth>0);
    if(search_stage_)
    {
        previous_search_stage_=w.active_search_stage;
        w.active_search_stage=name;
        auto & stage=SearchStage(w,name); ++stage.calls;
        if(rows>=0) stage.rows=std::max(stage.rows,rows);
        if(columns>=0) stage.columns=std::max(stage.columns,columns);
        stage.nonzeros=std::max(stage.nonzeros,nonzeros);
        NotifyResourceStageObserver();
    }
    started_=std::chrono::steady_clock::now();
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    if(w.enabled) RecordFactorRssPeakForTesting(FactorResidencyWorkForTesting());
#endif
}
ResourcePhase::~ResourcePhase()
{
    if(!enabled_) return;
    auto & w=ResourceWorkForTesting();
    const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-started_).count();
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    RecordFactorRssPeakForTesting(FactorResidencyWorkForTesting());
#endif
    bool found=false;
    for(auto & record:w.phases) if(record.phase==name_) {++record.calls; record.inclusive_seconds+=seconds; found=true; break;}
    if(!found) w.phases.push_back({name_,1,seconds});
    if(search_stage_)
    {
        auto & stage=SearchStage(w,name_); ++stage.completed_calls; stage.seconds+=seconds;
        w.last_completed_search_stage=name_;
        if(std::find(w.completed_search_stages.begin(),w.completed_search_stages.end(),name_)==w.completed_search_stages.end())
            w.completed_search_stages.push_back(name_);
        w.active_search_stage=previous_search_stage_;
    }
    if(search_scope_ && w.search_depth>0) --w.search_depth;
    w.phase=previous_;
    if(search_stage_) NotifyResourceStageObserver();
}
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
cholmod_dense * MultiplyQ(int method,SuiteSparseQR_factorization<double> * qr,cholmod_dense * rhs,cholmod_common * cc)
{
    auto & work=SparseWorkForTesting(); ++work.q_actions; WorkTimer timer(work.q_seconds);
    return SuiteSparseQR_qmult<double>(method,qr,rhs,cc);
}
cholmod_dense * SolveR(int system,SuiteSparseQR_factorization<double> * qr,cholmod_dense * rhs,cholmod_common * cc)
{
    auto & work=SparseWorkForTesting(); ++work.triangular_solves; WorkTimer timer(work.triangular_seconds);
    return SuiteSparseQR_solve<double>(system,qr,rhs,cc);
}
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
    LongSparse r;
    cholmod_sparse * h{};
    cholmod_dense * tau{};
    int64_t * permutation{},* hpinv{};
    int fixed_rank{-1};
    Matrix Orthogonal(const Matrix & rhs,bool transpose)
    {
        auto b=View(rhs);
        if(fixed_rank<0) {Dense answer{MultiplyQ(transpose ? SPQR_QTX : SPQR_QX,qr,&b,&cc),&cc}; return answer.Copy();}
        auto & work=SparseWorkForTesting(); ++work.q_actions; WorkTimer timer(work.q_seconds);
        Dense answer{SuiteSparseQR_qmult<double>(transpose ? SPQR_QTX : SPQR_QX,h,tau,hpinv,&b,&cc),&cc};
        return answer.Copy();
    }
    Matrix Triangular(const Matrix & rhs,bool transpose)
    {
        const auto p=design.cols();
        auto & work=SparseWorkForTesting(); ++work.triangular_solves; WorkTimer timer(work.triangular_seconds);
        Matrix answer(p,rhs.cols());
        if(transpose)
        {
            for(Eigen::Index k=0;k<p;++k) answer.row(k)=rhs.row(permutation ? permutation[k] : k);
            return r.transpose().triangularView<Eigen::Lower>().solve(answer);
        }
        const Matrix solved=r.triangularView<Eigen::Upper>().solve(rhs.topRows(p));
        for(Eigen::Index k=0;k<p;++k) answer.row(permutation ? permutation[k] : k)=solved.row(k);
        return answer;
    }
    Sparse design;
    std::vector<Eigen::Index> columns;
    std::size_t generation{};
    const void * domain{}; const void * observations{};
    LinearPolicy policy{};
    double tolerance{};
    bool bound{},policy_present{},policy_mismatch_pending{};
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    std::size_t residency_id{},residency_generation{};
    bool residency_active{};
#endif
    SparseFactorState()
    {
        cholmod_l_start(&cc); cc.SPQR_nthreads=1;
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        residency_id=NewFactorIdForTesting();
#endif
    }
    ~SparseFactorState()
    {
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        if(residency_active) EndFactorResidencyForTesting(residency_id,residency_generation);
#endif
        if(qr) SuiteSparseQR_free<double>(&qr,&cc);
        if(h) cholmod_l_free_sparse(&h,&cc);
        if(tau) cholmod_l_free_dense(&tau,&cc);
        if(permutation) cholmod_l_free(static_cast<std::size_t>(design.cols()),sizeof(int64_t),permutation,&cc);
        if(hpinv) cholmod_l_free(static_cast<std::size_t>(design.rows()),sizeof(int64_t),hpinv,&cc);
        cholmod_l_finish(&cc);
    }
    void Clear()
    {
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        if(residency_active)
        {EndFactorResidencyForTesting(residency_id,residency_generation); residency_active=false;}
#endif
        ++generation; if(qr) SuiteSparseQR_free<double>(&qr,&cc);
    }
};
LinearWorkspace::LinearWorkspace():state_(std::make_shared<SparseFactorState>()) {}
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
void LinearWorkspace::HandoffForTesting() {state_=std::make_shared<SparseFactorState>();}
#endif
void LinearWorkspace::Bind(const void * domain,const void * observations,const LinearPolicy * policy)
{
    const auto p=policy ? *policy : LinearPolicy{};
    auto & s=*state_;
    s.policy_mismatch_pending=false;
    const bool policy_changed=s.bound && (s.policy_present!=(policy!=nullptr) ||
        s.policy.rank_relative!=p.rank_relative || s.policy.release_factor!=p.release_factor ||
        s.policy.release_response_norm!=p.release_response_norm || s.policy.active_set_iteration_factor!=p.active_set_iteration_factor);
    if(s.bound && (s.domain!=domain || s.observations!=observations || policy_changed))
    {
        s.policy_mismatch_pending=policy_changed && s.qr!=nullptr;
        s.Clear();
    }
    s.domain=domain; s.observations=observations; s.policy=p; s.policy_present=policy!=nullptr; s.bound=true;
}
std::shared_ptr<FreeDesignFactor> LinearWorkspace::Factor(const Sparse & a,const std::vector<Eigen::Index> & columns,double tolerance)
{
    if(copy_on_write_ && state_.use_count()>1)
        state_=std::make_shared<SparseFactorState>();
    auto & s=*state_;
    auto & work=SparseWorkForTesting();
    if(LinearTelemetryEnabledForTesting())
    {
        ++work.numeric_factor_requests;
        if(s.policy_mismatch_pending) ++work.numeric_factor_policy_mismatches;
        else if(s.qr)
        {
            if(s.columns!=columns) ++work.numeric_factor_column_mismatches;
            else if(s.tolerance!=tolerance) ++work.numeric_factor_policy_mismatches;
            else if(SameSparsePattern(s.design,a))
            {
                if(SameSparseValues(s.design,a)) RecordExactReuseOpportunity(work);
                else {++work.numeric_factor_pattern_only_matches; ++work.numeric_factor_value_mismatches;}
            }
            else ++work.numeric_factor_pattern_mismatches;
        }
    }
    s.policy_mismatch_pending=false;
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    if(s.residency_active)
    {EndFactorResidencyForTesting(s.residency_id,s.residency_generation); s.residency_active=false;}
#endif
    ++s.generation;
    LongSparse storage; cholmod_sparse view{};
    {
        ResourcePhase preparation("linear-factor-preparation",true,a.rows(),a.cols(),static_cast<std::size_t>(a.nonZeros()));
        const auto conversion=Clock::now(); storage=a; storage.makeCompressed(); view=View(storage);
        RecordSparseShape("free-design",a.rows(),a.cols(),static_cast<std::size_t>(a.nonZeros()));
        SparseWorkForTesting().matrix_preparation_seconds+=Seconds(conversion);
    }
    if(!s.qr || s.tolerance!=tolerance || s.columns!=columns || !Pattern(s.design,a))
    {
        s.Clear(); const auto start=Clock::now();
        ++work.symbolic; work.symbolic_rows=std::max(work.symbolic_rows,a.rows());
        work.symbolic_columns=std::max(work.symbolic_columns,a.cols());
        work.symbolic_input_nonzeros=std::max(work.symbolic_input_nonzeros,static_cast<std::size_t>(a.nonZeros()));
        {
            ResourcePhase linear_stage("linear-symbolic",true,a.rows(),a.cols(),static_cast<std::size_t>(a.nonZeros()));
            ResourcePhase stage("spqr-symbolic",true,a.rows(),a.cols(),static_cast<std::size_t>(a.nonZeros()));
            RecordSparseShape("free-design",a.rows(),a.cols(),static_cast<std::size_t>(a.nonZeros()));
            s.qr=SuiteSparseQR_symbolic<double>(SPQR_ORDERING_COLAMD,true,&view,&s.cc);
            work.symbolic_seconds+=Seconds(start);
        }
        if(!s.qr) throw std::runtime_error("SPQR symbolic analysis failed");
    }
    else ++SparseWorkForTesting().symbolic_reuses;
    s.design=a; s.columns=columns; s.tolerance=tolerance;
    const auto start=Clock::now(); int ok{};
    ++work.numeric; work.numeric_rows=std::max(work.numeric_rows,a.rows());
    work.numeric_columns=std::max(work.numeric_columns,a.cols());
    work.numeric_input_nonzeros=std::max(work.numeric_input_nonzeros,static_cast<std::size_t>(a.nonZeros()));
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    const auto role=FactorCreationRoleForTesting();
    BeginFactorConstructionForTesting(s.residency_id,s.generation,"workspace",role,a.rows(),a.cols(),
        static_cast<std::size_t>(a.nonZeros()),[state=&s]{return OwnedFactorBytesForTesting(state->design,nullptr,state->cc);});
#endif
    {
        ResourcePhase linear_stage("linear-numeric",true,a.rows(),a.cols(),static_cast<std::size_t>(a.nonZeros()));
        ResourcePhase stage("spqr-numeric",true,a.rows(),a.cols(),static_cast<std::size_t>(a.nonZeros()));
        RecordSparseShape("free-design",a.rows(),a.cols(),static_cast<std::size_t>(a.nonZeros()));
        ok=SuiteSparseQR_numeric<double>(tolerance,&view,s.qr,&s.cc);
        work.numeric_seconds+=Seconds(start);
    }
    if(!ok)
    {
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        EndFactorConstructionForTesting(s.residency_id,s.generation);
#endif
        s.Clear(); throw std::runtime_error("SPQR numeric factorization failed");
    }
    work.factor_nonzeros=std::max(work.factor_nonzeros,static_cast<std::size_t>(std::max<int64_t>(0,s.cc.SPQR_istat[0])));
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    EndFactorConstructionForTesting(s.residency_id,s.generation);
    const auto r_nonzeros=static_cast<std::size_t>(std::max<int64_t>(0,s.cc.SPQR_istat[0]));
    const auto h_nonzeros=s.qr->QRnum ? static_cast<std::size_t>(std::max<int64_t>(0,s.qr->QRnum->hisize)) : 0;
    StartFactorResidencyForTesting(s.residency_id,s.generation,"workspace",role,a.rows(),a.cols(),
        static_cast<std::size_t>(a.nonZeros()),r_nonzeros,h_nonzeros,OwnedFactorBytesForTesting(s.design,nullptr,s.cc));
    s.residency_generation=s.generation; s.residency_active=true;
#endif
    return std::shared_ptr<FreeDesignFactor>(new FreeDesignFactor(state_,s.generation));
}
std::shared_ptr<FreeDesignFactor> FreeDesignFactor::Fixed(const Sparse & a,const std::vector<Eigen::Index> & columns)
{
    auto state=std::make_shared<SparseFactorState>(); auto & s=*state;
    s.design=a; s.columns=columns;
    LongSparse storage=a; storage.makeCompressed(); auto view=View(storage);
    cholmod_sparse * r{};
    auto & work=SparseWorkForTesting();
    ++work.fixed_factorizations;
    work.fixed_factor_rows=std::max(work.fixed_factor_rows,a.rows());
    work.fixed_factor_columns=std::max(work.fixed_factor_columns,a.cols());
    work.fixed_factor_input_nonzeros=std::max(work.fixed_factor_input_nonzeros,static_cast<std::size_t>(a.nonZeros()));
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    BeginFactorConstructionForTesting(s.residency_id,0,"operator-fixed","operator-fixed",a.rows(),a.cols(),
        static_cast<std::size_t>(a.nonZeros()),[state=&s]{return OwnedFactorBytesForTesting(state->design,&state->r,state->cc);});
#endif
    {
        ResourcePhase stage("spqr-fixed-factor",true,a.rows(),a.cols(),static_cast<std::size_t>(a.nonZeros()));
        RecordSparseShape("free-design",a.rows(),a.cols(),static_cast<std::size_t>(a.nonZeros()));
        WorkTimer timer(work.fixed_factor_seconds);
        s.fixed_rank=static_cast<int>(SuiteSparseQR<double>(SPQR_ORDERING_COLAMD,0,static_cast<int64_t>(a.cols()),
            &view,&r,&s.permutation,&s.h,&s.hpinv,&s.tau,&s.cc));
    }
    if(r)
    {
        s.r=Eigen::Map<const LongSparse>(static_cast<Eigen::Index>(r->nrow),static_cast<Eigen::Index>(r->ncol),
            static_cast<Eigen::Index>(static_cast<int64_t *>(r->p)[r->ncol]),static_cast<int64_t *>(r->p),
            static_cast<int64_t *>(r->i),static_cast<double *>(r->x));
        cholmod_l_free_sparse(&r,&s.cc);
    }
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    EndFactorConstructionForTesting(s.residency_id,0);
#endif
    if(s.fixed_rank<0 || !s.h || !s.tau || s.r.rows()!=a.cols()) throw std::runtime_error("SPQR fixed factorization failed");
    work.factor_nonzeros=std::max(work.factor_nonzeros,static_cast<std::size_t>(s.r.nonZeros()));
    work.fixed_factor_nonzeros=std::max(work.fixed_factor_nonzeros,static_cast<std::size_t>(s.r.nonZeros()));
    // Owned exported arrays only, not an estimate of peak factorization scratch.
    const auto bytes=static_cast<std::size_t>(s.r.nonZeros())*(sizeof(double)+sizeof(int64_t))+
        static_cast<std::size_t>(s.r.cols()+1)*sizeof(int64_t)+s.h->nzmax*(sizeof(double)+sizeof(int64_t))+
        (s.h->ncol+1)*sizeof(int64_t)+s.tau->nzmax*sizeof(double)+
        (s.permutation ? static_cast<std::size_t>(a.cols())*sizeof(int64_t) : 0)+
        (s.hpinv ? static_cast<std::size_t>(a.rows())*sizeof(int64_t) : 0);
    work.factor_storage_bytes=std::max(work.factor_storage_bytes,bytes);
    work.fixed_factor_storage_bytes=std::max(work.fixed_factor_storage_bytes,bytes);
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    const auto h_nonzeros=static_cast<std::size_t>(static_cast<const int64_t *>(s.h->p)[s.h->ncol]);
    StartFactorResidencyForTesting(s.residency_id,0,"operator-fixed","operator-fixed",a.rows(),a.cols(),
        static_cast<std::size_t>(a.nonZeros()),static_cast<std::size_t>(s.r.nonZeros()),h_nonzeros,
        OwnedFactorBytesForTesting(s.design,&s.r,s.cc));
    s.residency_generation=0; s.residency_active=true;
#endif
    return std::shared_ptr<FreeDesignFactor>(new FreeDesignFactor(std::move(state),0));
}
FreeDesignFactor::FreeDesignFactor(std::shared_ptr<SparseFactorState> state,std::size_t generation):state_(std::move(state)),generation_(generation) {}
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
const Sparse & FreeDesignFactor::DesignForTesting() const {Check(); return state_->design;}
const std::vector<Eigen::Index> & FreeDesignFactor::ColumnsForTesting() const {Check(); return state_->columns;}
double FreeDesignFactor::ToleranceForTesting() const {Check(); return state_->tolerance;}
ProjectedTailTransformForTesting FreeDesignFactor::OrthogonalTransposeTailSparseForTesting(
    const Sparse & rhs,bool materialize_tail) const
{
    Check(); auto & s=*state_;
    if(rhs.rows()!=s.design.rows()) throw std::invalid_argument("Invalid sparse Q transpose RHS size");
    ProjectedTailTransformForTesting result;
    result.rows=rhs.rows(); result.columns=rhs.cols(); result.input_nonzeros=static_cast<std::size_t>(rhs.nonZeros());
    const auto p=s.design.cols();
    LongSparse storage=rhs; storage.makeCompressed(); auto view=View(storage);
    struct SparseOwner
    {
        cholmod_sparse * value{}; cholmod_common * common{};
        ~SparseOwner() {if(value) cholmod_l_free_sparse(&value,common);}
    } transformed{nullptr,&s.cc};
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    AssessmentStageTimerForTesting transform_stage("projected-tail-q-transform",rhs.rows(),rhs.cols());
#endif
    {
        ResourcePhase phase("spqr-q-transform-sparse",true,rhs.rows(),rhs.cols(),result.input_nonzeros);
        auto & work=SparseWorkForTesting(); ++work.q_actions; WorkTimer timer(work.q_seconds);
        const auto started=Clock::now();
        transformed.value=s.fixed_rank>=0 ? SuiteSparseQR_qmult<double>(SPQR_QTX,s.h,s.tau,s.hpinv,&view,&s.cc) :
            SuiteSparseQR_qmult<double>(SPQR_QTX,s.qr,&view,&s.cc);
        result.q_transform_seconds=Seconds(started);
    }
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    transform_stage.Finish();
#endif
    if(!transformed.value) throw std::runtime_error("SPQR sparse Q transpose failed");
    if(transformed.value->itype!=CHOLMOD_LONG) throw std::runtime_error("Unexpected SPQR sparse index type");
    result.tail_rows=static_cast<Eigen::Index>(transformed.value->nrow)-p;
    if(result.tail_rows<0 || transformed.value->nrow!=static_cast<std::size_t>(rhs.rows()) ||
        transformed.value->ncol!=static_cast<std::size_t>(rhs.cols()))
        throw std::runtime_error("Unexpected SPQR sparse Q transpose shape");
    const auto * outer=static_cast<const int64_t *>(transformed.value->p);
    const auto * inner=static_cast<const int64_t *>(transformed.value->i);
    result.transformed_nonzeros=static_cast<std::size_t>(outer[transformed.value->ncol]);
    result.transformed_storage_bytes=transformed.value->nzmax*(sizeof(double)+sizeof(int64_t))+
        (transformed.value->ncol+1)*sizeof(int64_t);
    if(!transformed.value->sorted) throw std::runtime_error("Unsorted SPQR sparse Q transpose output");
    const auto extraction_started=Clock::now();
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    AssessmentStageTimerForTesting extraction_stage("projected-tail-extract",rhs.rows(),rhs.cols());
#endif
    const auto * values=static_cast<const double *>(transformed.value->x);
    for(Eigen::Index column=0;column<result.columns;++column)
        for(auto item=outer[column];item<outer[column+1];++item)
            if(inner[item]>=p) ++result.tail_nonzeros;
    result.tail_storage_bytes=result.tail_nonzeros*(sizeof(double)+sizeof(int))+
        static_cast<std::size_t>(result.columns+1)*sizeof(int);
    if(materialize_tail)
    {
        result.tail.resize(result.tail_rows,result.columns);
        result.tail.reserve(static_cast<Eigen::Index>(result.tail_nonzeros));
        for(Eigen::Index column=0;column<result.columns;++column)
        {
            result.tail.startVec(static_cast<int>(column));
            for(auto item=outer[column];item<outer[column+1];++item)
                if(inner[item]>=p)
                    result.tail.insertBack(static_cast<int>(inner[item]-p),static_cast<int>(column))=values[item];
        }
        result.tail.finalize(); result.tail.makeCompressed();
    }
    result.tail_extract_seconds=Seconds(extraction_started);
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    extraction_stage.Finish();
#endif
    return result;
}
Matrix FreeDesignFactor::OrthogonalTransposeForTesting(const Matrix & rhs) const
{Check(); return state_->Orthogonal(rhs,true);}
ProjectedTailQrResultForTesting FreeDesignFactor::ProjectedTailQrForTesting(
    const Sparse & raw,VectorRef residual,double scale) const
{
    Check();
    ProjectedTailQrResultForTesting result;
    result.rows=raw.rows(); result.free_design_columns=state_->design.cols(); result.width_columns=raw.cols();
    const auto started=Clock::now();
    if(raw.rows()!=state_->design.rows() || residual.size()!=raw.rows() || raw.cols()<=0 ||
        !(scale>0) || !std::isfinite(scale))
    {result.reason="invalid-projected-tail-qr-input"; return result;}
    try
    {
        auto transformed=OrthogonalTransposeTailSparseForTesting(raw,true);
        result.q_transformed_nonzeros=transformed.transformed_nonzeros;
        result.q_transformed_storage_bytes=transformed.transformed_storage_bytes;
        result.tail_nonzeros=transformed.tail_nonzeros;
        result.tail_storage_bytes=transformed.tail_storage_bytes;
        result.tail_rows=transformed.tail_rows;
        result.tail_extract_seconds=transformed.tail_extract_seconds;
        result.q_transform_seconds=transformed.q_transform_seconds;
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        AssessmentStageTimerForTesting residual_transform_stage("projected-tail-residual-q-transform",raw.rows(),1);
#endif
        const auto residual_started=Clock::now();
        const Matrix qz_residual=OrthogonalTransposeForTesting(Matrix(residual/scale));
        result.q_transform_seconds+=Seconds(residual_started);
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        residual_transform_stage.Finish();
#endif
        const Matrix tail_response=qz_residual.bottomRows(result.tail_rows);

        auto tail_state=std::make_shared<SparseFactorState>();
        tail_state->design=std::move(transformed.tail);
        tail_state->columns.resize(static_cast<std::size_t>(raw.cols()));
        for(Eigen::Index column=0;column<raw.cols();++column)
            tail_state->columns[static_cast<std::size_t>(column)]=column;
        LongSparse tail_storage=tail_state->design; tail_storage.makeCompressed(); auto tail_view=View(tail_storage);
        struct ConstructionScope
        {
            std::size_t id{},generation{}; bool active{};
            ~ConstructionScope() {if(active) EndFactorConstructionForTesting(id,generation);}
            void Finish() {if(active) {EndFactorConstructionForTesting(id,generation); active=false;}}
        } construction;
        auto & work=SparseWorkForTesting();
        BeginFactorConstructionForTesting(tail_state->residency_id,tail_state->generation,
            "projected-tail","projected-tail",tail_state->design.rows(),tail_state->design.cols(),
            static_cast<std::size_t>(tail_state->design.nonZeros()),[state=tail_state.get()]{
                return OwnedFactorBytesForTesting(state->design,nullptr,state->cc);});
        construction={tail_state->residency_id,tail_state->generation,true};
        ++work.symbolic;
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        AssessmentStageTimerForTesting symbolic_stage("projected-tail-symbolic",tail_state->design.rows(),
            tail_state->design.cols());
#endif
        {
            ResourcePhase phase("spqr-tail-symbolic",true,tail_state->design.rows(),tail_state->design.cols(),
                static_cast<std::size_t>(tail_state->design.nonZeros()));
            const auto symbolic_started=Clock::now();
            RecordSparseShape("projected-tail",tail_state->design.rows(),tail_state->design.cols(),
                static_cast<std::size_t>(tail_state->design.nonZeros()));
            tail_state->qr=SuiteSparseQR_symbolic<double>(SPQR_ORDERING_COLAMD,false,&tail_view,&tail_state->cc);
            result.tail_symbolic_seconds=Seconds(symbolic_started);
            work.symbolic_seconds+=result.tail_symbolic_seconds;
        }
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        symbolic_stage.Finish();
#endif
        if(!tail_state->qr) {result.reason="projected-tail-symbolic-failed"; return result;}
        ++work.numeric;
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        AssessmentStageTimerForTesting numeric_stage("projected-tail-numeric",tail_state->design.rows(),
            tail_state->design.cols());
#endif
        {
            ResourcePhase phase("spqr-tail-numeric",true,tail_state->design.rows(),tail_state->design.cols(),
                static_cast<std::size_t>(tail_state->design.nonZeros()));
            const auto numeric_started=Clock::now();
            const auto success=SuiteSparseQR_numeric<double>(SPQR_NO_TOL,&tail_view,tail_state->qr,&tail_state->cc);
            result.tail_numeric_seconds=Seconds(numeric_started);
            work.numeric_seconds+=result.tail_numeric_seconds;
            if(!success) {result.reason="projected-tail-numeric-failed"; return result;}
        }
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        numeric_stage.Finish();
#endif
        construction.Finish();
        result.tail_factor_nonzeros=static_cast<std::size_t>(std::max<int64_t>(0,tail_state->cc.SPQR_istat[0]));
        result.tail_factor_storage_bytes=OwnedFactorBytesForTesting(tail_state->design,nullptr,tail_state->cc);
        const auto h_nonzeros=tail_state->qr->QRnum ?
            static_cast<std::size_t>(std::max<int64_t>(0,tail_state->qr->QRnum->hisize)) : 0;
        StartFactorResidencyForTesting(tail_state->residency_id,tail_state->generation,"projected-tail",
            "projected-tail",tail_state->design.rows(),tail_state->design.cols(),
            static_cast<std::size_t>(tail_state->design.nonZeros()),result.tail_factor_nonzeros,h_nonzeros,
            result.tail_factor_storage_bytes);
        tail_state->residency_generation=tail_state->generation; tail_state->residency_active=true;
        const auto tail_factor=std::shared_ptr<FreeDesignFactor>(
            new FreeDesignFactor(tail_state,tail_state->generation));
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        AssessmentStageTimerForTesting compact_stage("projected-tail-compact",raw.cols(),raw.cols());
#endif
        const auto compact_started=Clock::now();
        result.factor=tail_factor->Compact()/scale;
        result.tail_compact_seconds=Seconds(compact_started);
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        compact_stage.Finish();
        AssessmentStageTimerForTesting qmult_stage("projected-tail-qmult",raw.cols(),1);
#endif
        const auto qmult_started=Clock::now();
        result.response=tail_factor->OrthogonalTransposeForTesting(tail_response).topRows(raw.cols()).col(0);
        result.tail_qmult_seconds=Seconds(qmult_started);
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        qmult_stage.Finish();
#endif
        result.columns_restored=true;
        result.valid=result.factor.allFinite() && result.response.allFinite();
        result.reason=result.valid ? "projected-tail-qr" : "projected-tail-nonfinite";
    }
    catch(const std::exception & error)
    {result.reason=error.what();}
    result.seconds=Seconds(started);
    return result;
}
#endif
void FreeDesignFactor::Check() const
{if(!state_ || state_->generation!=generation_ || (!state_->qr && state_->fixed_rank<0)) throw std::logic_error("Expired free-design factor");}
bool FreeDesignFactor::Matches(const Sparse & a,const std::vector<Eigen::Index> & columns) const
{
    return state_ && state_->generation==generation_ && (state_->qr || state_->fixed_rank>=0) && state_->columns==columns && Pattern(state_->design,a) &&
        std::equal(a.valuePtr(),a.valuePtr()+a.nonZeros(),state_->design.valuePtr());
}
int FreeDesignFactor::Rank() const {Check(); return state_->fixed_rank>=0 ? state_->fixed_rank : static_cast<int>(state_->qr->rank);}
std::optional<RankFactorView> FreeDesignFactor::RankView() const
{
    Check(); const auto & s=*state_;
    if(s.fixed_rank<0) return std::nullopt;
    const auto n=s.design.rows(),p=s.design.cols();
    const auto * hp=static_cast<const int64_t *>(s.h->p);
    return RankFactorView{&s.design,n,p,static_cast<Eigen::Index>(s.h->ncol),
        {s.r.outerIndexPtr(),static_cast<std::size_t>(p+1)},
        {s.r.innerIndexPtr(),static_cast<std::size_t>(s.r.nonZeros())},
        {hp,s.h->ncol+1},{static_cast<const int64_t *>(s.h->i),static_cast<std::size_t>(hp[s.h->ncol])},
        {s.permutation,s.permutation ? static_cast<std::size_t>(p) : 0},
        {s.hpinv,s.hpinv ? static_cast<std::size_t>(n) : 0},
        {s.r.valuePtr(),static_cast<std::size_t>(s.r.nonZeros())},
        {static_cast<const double *>(s.h->x),static_cast<std::size_t>(hp[s.h->ncol])},
        {static_cast<const double *>(s.tau->x),s.h->ncol}};
}
Matrix FreeDesignFactor::LeastSquares(const Matrix & rhs) const
{
    Check(); auto b=View(rhs); auto & s=*state_;
    if(s.fixed_rank>=0) return s.Triangular(s.Orthogonal(rhs,true),false);
    Dense transformed{MultiplyQ(SPQR_QTX,s.qr,&b,&s.cc),&s.cc};
    if(!transformed.p) throw std::runtime_error("SPQR Q transpose failed");
    Dense answer{SolveR(SPQR_RETX_EQUALS_B,s.qr,transformed.p,&s.cc),&s.cc}; return answer.Copy();
}
Matrix FreeDesignFactor::NormalSolve(const Matrix & rhs) const
{
    Check(); auto b=View(rhs); auto & s=*state_;
    if(s.fixed_rank>=0) return s.Triangular(s.Triangular(rhs,true),false);
    Dense adjoint{SolveR(SPQR_RTX_EQUALS_ETB,s.qr,&b,&s.cc),&s.cc};
    if(!adjoint.p) throw std::runtime_error("SPQR adjoint triangular solve failed");
    Dense answer{SolveR(SPQR_RETX_EQUALS_B,s.qr,adjoint.p,&s.cc),&s.cc}; return answer.Copy();
}
Matrix FreeDesignFactor::PseudoInverseTranspose(const Matrix & rhs) const
{
    Check(); auto & s=*state_;
    if(rhs.rows()!=s.design.cols()) throw std::invalid_argument("Invalid adjoint RHS size");
    if(s.fixed_rank>=0)
    {
        Matrix padded=Matrix::Zero(s.design.rows(),rhs.cols());
        padded.topRows(s.design.cols())=s.Triangular(rhs,true);
        return s.Orthogonal(padded,false);
    }
    auto b=View(rhs);
    Dense triangular{SolveR(SPQR_RTX_EQUALS_ETB,s.qr,&b,&s.cc),&s.cc};
    Matrix padded=Matrix::Zero(s.design.rows(),rhs.cols());
    padded.topRows(s.design.cols())=triangular.Copy().topRows(s.design.cols());
    auto v=View(padded);
    Dense answer{MultiplyQ(SPQR_QX,s.qr,&v,&s.cc),&s.cc};
    return answer.Copy();
}
Matrix FreeDesignFactor::ProjectComplement(const Matrix & rhs) const
{
    Check(); auto & s=*state_;
    if(rhs.rows()!=s.design.rows()) throw std::invalid_argument("Invalid projection RHS size");
    if(s.fixed_rank>=0)
    {
        Matrix tail=s.Orthogonal(rhs,true); tail.topRows(s.design.cols()).setZero();
        return s.Orthogonal(tail,false);
    }
    auto b=View(rhs);
    Dense transformed{MultiplyQ(SPQR_QTX,s.qr,&b,&s.cc),&s.cc};
    Matrix tail=transformed.Copy(); tail.topRows(s.design.cols()).setZero(); auto v=View(tail);
    Dense answer{MultiplyQ(SPQR_QX,s.qr,&v,&s.cc),&s.cc};
    return answer.Copy();
}
Matrix FreeDesignFactor::Compact() const
{
    Check(); auto & work=SparseWorkForTesting(); ++work.compact_extractions; WorkTimer timer(work.compact_seconds);
    auto & s=*state_; const auto p=s.design.cols(); Matrix compact(p,p);
    RecordDenseShape("free-design-compact",p,p);
    if(s.fixed_rank>=0)
    {
        compact.setZero();
        for(Eigen::Index k=0;k<p;++k) for(LongSparse::InnerIterator e(s.r,k);e;++e)
            compact(e.row(),s.permutation ? s.permutation[k] : k)=e.value();
        return compact;
    }
    // Q^T X retains original column order; its singular values equal those of R.
    // Never materialize the observation-sized Q or all RHS columns together.
    for(Eigen::Index first=0;first<p;first+=16)
    {
        const auto count=std::min<Eigen::Index>(16,p-first); Matrix rhs(s.design.middleCols(first,count)); auto b=View(rhs);
        Dense transformed{MultiplyQ(SPQR_QTX,s.qr,&b,&s.cc),&s.cc};
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
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
StructuredProjectedQrResultForTesting StructuredProjectedQrForTesting(
    const Sparse & design,const Sparse & raw,VectorRef residual,double scale)
{
    StructuredProjectedQrResultForTesting result;
    if(design.rows()!=raw.rows() || residual.size()!=design.rows() || design.cols()<=0 || raw.cols()<=0 ||
        !(scale>0) || !std::isfinite(scale))
    {result.reason="invalid-structured-projected-qr-input"; return result;}
    using Triplet=Eigen::Triplet<double,int64_t>;
    std::vector<Triplet> entries;
    entries.reserve(static_cast<std::size_t>(design.nonZeros()+raw.nonZeros()));
    for(Eigen::Index column=0;column<design.outerSize();++column)
        for(Sparse::InnerIterator value(design,column);value;++value)
            entries.emplace_back(value.row(),column,value.value());
    for(Eigen::Index column=0;column<raw.outerSize();++column)
        for(Sparse::InnerIterator value(raw,column);value;++value)
            entries.emplace_back(value.row(),design.cols()+column,value.value());
    LongSparse combined(design.rows(),design.cols()+raw.cols());
    combined.setFromTriplets(entries.begin(),entries.end()); combined.makeCompressed();
    result.sparse_rows=combined.rows(); result.sparse_columns=combined.cols();
    result.sparse_nonzeros=static_cast<std::size_t>(combined.nonZeros());
    auto & work=SparseWorkForTesting();
    work.symbolic_rows=std::max(work.symbolic_rows,combined.rows());
    work.symbolic_columns=std::max(work.symbolic_columns,combined.cols());
    work.symbolic_input_nonzeros=std::max(work.symbolic_input_nonzeros,result.sparse_nonzeros);
    work.numeric_rows=std::max(work.numeric_rows,combined.rows());
    work.numeric_columns=std::max(work.numeric_columns,combined.cols());
    work.numeric_input_nonzeros=std::max(work.numeric_input_nonzeros,result.sparse_nonzeros);
    struct Owner
    {
        cholmod_common common{};
        SuiteSparseQR_factorization<double> * factor{};
        cholmod_sparse * transformed{};
        cholmod_dense * response{};
        bool started{};
        Owner():started(cholmod_l_start(&common)!=0) {common.SPQR_nthreads=1;}
        ~Owner()
        {
            if(transformed) cholmod_l_free_sparse(&transformed,&common);
            if(response) cholmod_l_free_dense(&response,&common);
            if(factor) SuiteSparseQR_free<double>(&factor,&common);
            if(started) cholmod_l_finish(&common);
        }
    } owner;
    if(!owner.started) {result.reason="spqr-common-init-failed"; return result;}
    auto matrix=View(combined);
    ++work.symbolic;
    {
        ResourcePhase phase("spqr-projected-symbolic",true,combined.rows(),combined.cols(),result.sparse_nonzeros);
        const auto started=Clock::now();
        owner.factor=SuiteSparseQR_symbolic<double>(SPQR_ORDERING_FIXED,false,&matrix,&owner.common);
        result.symbolic_seconds=Seconds(started); work.symbolic_seconds+=result.symbolic_seconds;
    }
    if(!owner.factor) {result.reason="spqr-projected-symbolic-failed"; return result;}
    const auto is_identity=[](const int64_t * permutation,Eigen::Index size) {
        if(!permutation) return true;
        for(Eigen::Index k=0;k<size;++k) if(permutation[k]!=k) return false;
        return true;
    };
    result.input_order_preserved=is_identity(owner.factor->QRsym->Qfill,combined.cols()) &&
        is_identity(owner.factor->Q1fill,combined.cols());
    if(!result.input_order_preserved)
    {result.reason="spqr-fixed-order-permutation"; return result;}
    ++work.numeric;
    {
        ResourcePhase phase("spqr-projected-numeric",true,combined.rows(),combined.cols(),result.sparse_nonzeros);
        const auto started=Clock::now();
        const int success=SuiteSparseQR_numeric<double>(SPQR_NO_TOL,&matrix,owner.factor,&owner.common);
        result.numeric_seconds=Seconds(started); work.numeric_seconds+=result.numeric_seconds;
        if(!success) {result.reason="spqr-projected-numeric-failed"; return result;}
    }
    if(owner.factor->QRsym->do_rank_detection)
    {result.reason="spqr-projected-rank-detection-enabled"; return result;}
    const Eigen::Index p=design.cols(),m=raw.cols(),k=p+m;
    Matrix rhs=residual/scale; auto dense_rhs=View(rhs);
    {
        ResourcePhase phase("spqr-projected-qmult",true,combined.rows(),combined.cols(),result.sparse_nonzeros);
        const auto started=Clock::now();
        owner.response=SuiteSparseQR_qmult<double>(SPQR_QTX,owner.factor,&dense_rhs,&owner.common);
        owner.transformed=SuiteSparseQR_qmult<double>(SPQR_QTX,owner.factor,&matrix,&owner.common);
        result.orthogonal_seconds=Seconds(started);
    }
    if(!owner.response || !owner.transformed)
    {result.reason="spqr-projected-qmult-failed"; return result;}
    if(owner.response->nrow<static_cast<std::size_t>(k) || owner.transformed->nrow<static_cast<std::size_t>(k) ||
        owner.transformed->ncol!=static_cast<std::size_t>(k))
    {result.reason="spqr-projected-factor-shape"; return result;}
    result.factor=Matrix::Zero(m,m); result.response.resize(m);
    result.factor_rows=m; result.factor_columns=m;
    result.maximum_dense_bytes=std::max({static_cast<std::size_t>(rhs.size())*sizeof(double),
        static_cast<std::size_t>(m)*static_cast<std::size_t>(m)*sizeof(double),
        static_cast<std::size_t>(k)*sizeof(double)});
    const auto * outer=static_cast<const int64_t *>(owner.transformed->p);
    const auto * inner=static_cast<const int64_t *>(owner.transformed->i);
    const auto * values=static_cast<const double *>(owner.transformed->x);
    for(Eigen::Index column=p;column<k;++column)
        for(int64_t item=outer[column];item<outer[column+1];++item)
        {
            const auto row=static_cast<Eigen::Index>(inner[item]);
            if(row>=p && row<k) result.factor(row-p,column-p)=values[item]/scale;
        }
    const auto * transformed_response=static_cast<const double *>(owner.response->x);
    for(Eigen::Index column=0;column<m;++column)
        result.response(column)=transformed_response[p+column];
    if(!result.factor.allFinite() || !result.response.allFinite())
    {result.reason="spqr-projected-nonfinite"; return result;}
    result.valid=true; result.reason="structured-fixed-order-qr";
    return result;
}
#endif
}
