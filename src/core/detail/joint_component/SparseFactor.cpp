#include "SparseFactor.hpp"
#include <algorithm>
#include <chrono>
#include <stdexcept>
#include <cstring>
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
#include <sys/resource.h>
#endif
#include <Eigen/SparseQR>

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

#ifdef RHBM_GEM_TEST_INSTRUMENTATION
StructuredProjectedQrResultForTesting StructuredProjectedQrForTesting(
    const Sparse&,
    const Sparse&,
    VectorRef,
    const double)
{
    StructuredProjectedQrResultForTesting result;
    result.reason =
        "structured projected QR is not part of the EIGEN production path";
    return result;
}
#endif

// Qualified EIGEN sparse factorization used by the active-set linear solve.
struct SparseFactorState
{
    Eigen::SparseQR<Sparse,Eigen::COLAMDOrdering<int>> qr;
    Sparse design;
    std::vector<Eigen::Index> columns;
    double tolerance{};
    std::size_t generation{};
    bool valid{};
};

LinearWorkspace::LinearWorkspace()=default;
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
void LinearWorkspace::HandoffForTesting() {state_=std::make_shared<SparseFactorState>();}
#endif
void LinearWorkspace::Bind(const void *,const void *,const LinearPolicy *) {}
std::shared_ptr<FreeDesignFactor> LinearWorkspace::Factor(const Sparse & a,const std::vector<Eigen::Index> & columns,double tolerance)
{
    if(!state_) state_=std::make_shared<SparseFactorState>();
    if(copy_on_write_ && state_.use_count()>1) state_=std::make_shared<SparseFactorState>();
    auto & s=*state_; auto & work=SparseWorkForTesting();
    if(LinearTelemetryEnabledForTesting())
    {
        ++work.numeric_factor_requests;
        if(s.valid)
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
    ++s.generation; s.valid=false; s.design=a; s.columns=columns; s.tolerance=tolerance;
    s.qr.setPivotThreshold(tolerance);
    {
        ResourcePhase phase("linear-symbolic",true,a.rows(),a.cols(),static_cast<std::size_t>(a.nonZeros()));
        WorkTimer timer(work.symbolic_seconds); s.qr.analyzePattern(a); ++work.symbolic;
    }
    if(s.qr.info()!=Eigen::Success) throw std::runtime_error("Eigen sparse symbolic analysis failed");
    {
        ResourcePhase phase("linear-numeric",true,a.rows(),a.cols(),static_cast<std::size_t>(a.nonZeros()));
        WorkTimer timer(work.numeric_seconds); s.qr.factorize(a); ++work.numeric;
    }
    if(s.qr.info()!=Eigen::Success) throw std::runtime_error("Eigen sparse numeric factorization failed");
    s.valid=true;
    work.factor_nonzeros=std::max(work.factor_nonzeros,static_cast<std::size_t>(s.qr.matrixR().nonZeros()));
    return std::shared_ptr<FreeDesignFactor>(new FreeDesignFactor(state_,s.generation));
}
std::shared_ptr<FreeDesignFactor> FreeDesignFactor::Fixed(const Sparse & a,const std::vector<Eigen::Index> & columns)
{LinearWorkspace workspace; return workspace.Factor(a,columns,0);}
FreeDesignFactor::FreeDesignFactor(std::shared_ptr<SparseFactorState> state,std::size_t generation)
    :state_(std::move(state)),generation_(generation) {}
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
const Sparse & FreeDesignFactor::DesignForTesting() const {Check(); return state_->design;}
const std::vector<Eigen::Index> & FreeDesignFactor::ColumnsForTesting() const {Check(); return state_->columns;}
double FreeDesignFactor::ToleranceForTesting() const {Check(); return state_->tolerance;}
ProjectedTailTransformForTesting FreeDesignFactor::OrthogonalTransposeTailSparseForTesting(const Sparse &,bool) const
{throw std::runtime_error("projected-tail sparse transform is not part of the EIGEN production path");}
ProjectedTailQrResultForTesting FreeDesignFactor::ProjectedTailQrForTesting(const Sparse &,VectorRef,double) const
{
    ProjectedTailQrResultForTesting result;
    result.reason="projected-tail QR is not part of the EIGEN production path";
    return result;
}
Matrix FreeDesignFactor::OrthogonalTransposeForTesting(const Matrix & rhs) const
{Check(); return state_->qr.matrixQ().adjoint()*rhs;}
#endif
void FreeDesignFactor::Check() const
{if(!state_ || !state_->valid || generation_!=state_->generation) throw std::logic_error("Expired free-design factor");}
bool FreeDesignFactor::Matches(const Sparse & a,const std::vector<Eigen::Index> & columns) const
{
    return state_ && state_->valid && generation_==state_->generation && columns==state_->columns &&
        a.rows()==state_->design.rows() && a.cols()==state_->design.cols() && (a-state_->design).norm()==0;
}
int FreeDesignFactor::Rank() const {Check(); return static_cast<int>(state_->qr.rank());}
std::optional<RankFactorView> FreeDesignFactor::RankView() const {Check(); return std::nullopt;}
Matrix FreeDesignFactor::Compact() const
{
    Check(); auto & work=SparseWorkForTesting(); ++work.compact_extractions; WorkTimer timer(work.compact_seconds);
    const auto & s=*state_; const auto p=s.design.cols();
    RecordDenseShape("free-design-compact",p,p);
    return Matrix(s.qr.matrixR().topRows(p))*s.qr.colsPermutation().transpose();
}
Matrix FreeDesignFactor::LeastSquares(const Matrix & rhs) const
{
    Check(); if(rhs.rows()!=state_->design.rows()) throw std::invalid_argument("Invalid least-squares RHS size");
    auto & work=SparseWorkForTesting(); ++work.q_actions; ++work.triangular_solves;
    WorkTimer timer(work.least_squares_seconds); return state_->qr.solve(rhs);
}
Matrix FreeDesignFactor::PseudoInverseTranspose(const Matrix & rhs) const
{
    Check(); const auto & s=*state_; const auto p=s.design.cols();
    if(rhs.rows()!=p) throw std::invalid_argument("Invalid adjoint RHS size");
    Matrix padded=Matrix::Zero(s.design.rows(),rhs.cols());
    const Matrix permuted=s.qr.colsPermutation().transpose()*rhs;
    auto & work=SparseWorkForTesting(); ++work.triangular_solves;
    {
        WorkTimer timer(work.triangular_seconds);
        padded.topRows(p)=s.qr.matrixR().topLeftCorner(p,p).transpose().triangularView<Eigen::Lower>().solve(permuted);
    }
    ++work.q_actions; WorkTimer timer(work.q_seconds); return s.qr.matrixQ()*padded;
}
Matrix FreeDesignFactor::ProjectComplement(const Matrix & rhs) const
{
    Check(); const auto & s=*state_;
    if(rhs.rows()!=s.design.rows()) throw std::invalid_argument("Invalid projection RHS size");
    auto & work=SparseWorkForTesting(); work.q_actions+=2; WorkTimer timer(work.q_seconds);
    Matrix tail=s.qr.matrixQ().adjoint()*rhs; tail.topRows(s.design.cols()).setZero();
    return s.qr.matrixQ()*tail;
}
Matrix FreeDesignFactor::NormalSolve(const Matrix & rhs) const
{
    Check(); const auto & s=*state_; const auto p=s.design.cols();
    if(rhs.rows()!=p) throw std::invalid_argument("Invalid normal RHS size");
    auto & work=SparseWorkForTesting(); work.triangular_solves+=2; WorkTimer timer(work.triangular_seconds);
    const Matrix permuted=s.qr.colsPermutation().transpose()*rhs;
    const Matrix adjoint=s.qr.matrixR().topLeftCorner(p,p).transpose().triangularView<Eigen::Lower>().solve(permuted);
    const Matrix solved=s.qr.matrixR().topLeftCorner(p,p).triangularView<Eigen::Upper>().solve(adjoint);
    return s.qr.colsPermutation()*solved;
}
std::pair<Matrix,Vector> SparseReferenceQR(const Sparse & x,const Vector & weights,const Vector & scales,VectorRef y)
{return ReferenceQR(x,weights,scales,y);}
}
