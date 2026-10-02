#include "FreeDesignRank.hpp"
#include <algorithm>
#include <bit>
#include <cfenv>
#include <array>
#include <map>
#include <numeric>
#include <unordered_map>

namespace rhbm_gem::core::joint_component {
namespace {
// Each elementary operation is rounded outwards. Zero identities avoid widening
// exact structural zeros. No fast-math or changed floating-point environment is allowed.
struct Stop {const char * reason;};
constexpr double inf=std::numeric_limits<double>::infinity();
double Up(double x) {const auto y=std::nextafter(x,inf); if(!std::isfinite(y)) throw Stop{"rank-bound-overflow"}; return y;}
double Down(double x) {const auto y=std::nextafter(x,-inf); if(!std::isfinite(y)) throw Stop{"rank-bound-overflow"}; return y;}
template<class C> decltype(auto) At(C & values,Eigen::Index index) {return values[static_cast<std::size_t>(index)];}
struct Interval {double lo{},hi{}; Interval()=default; Interval(double x):lo(x),hi(x) {} Interval(double l,double h):lo(l),hi(h) {}};
Interval Add(Interval a,Interval b)
{
    if(a.lo==0 && a.hi==0) return b;
    if(b.lo==0 && b.hi==0) return a;
    return {Down(a.lo+b.lo),Up(a.hi+b.hi)};
}
Interval Neg(Interval a) {return {-a.hi,-a.lo};}
Interval Mul(Interval a,Interval b)
{
    if((a.lo==0 && a.hi==0) || (b.lo==0 && b.hi==0)) return {};
    const double values[]{a.lo*b.lo,a.lo*b.hi,a.hi*b.lo,a.hi*b.hi};
    return {Down(*std::min_element(values,values+4)),Up(*std::max_element(values,values+4))};
}
double AbsLower(Interval x) {return x.lo<=0 && x.hi>=0 ? 0 : std::min(std::abs(x.lo),std::abs(x.hi));}
double AbsUpper(Interval x) {return std::max(std::abs(x.lo),std::abs(x.hi));}
double SqrtUp(double x) {return x==0 ? 0 : Up(std::sqrt(x));}
double SqrtDown(double x) {return x<=0 ? 0 : std::max(0.,Down(std::sqrt(x)));}
double SumUp(double a,double b) {return b==0 ? a : a==0 ? b : Up(a+b);}
double ProductUp(double a,double b) {return a==0 || b==0 ? 0 : Up(a*b);}
struct Audit
{
    const RankBudget & budget; FreeDesignRankResult & result;
    std::chrono::steady_clock::time_point start{std::chrono::steady_clock::now()};
    void Stage(FreeDesignRankWorkStage stage) {result.work_stage=stage;}
    void Charge(std::size_t entries)
    {
        if(result.entries>budget.entries || entries>budget.entries-result.entries) throw Stop{"rank-work-budget"};
        result.entries+=entries;
        if(result.estimated_total_entries && result.entries<=*result.estimated_total_entries)
            result.estimated_remaining_entries=*result.estimated_total_entries-result.entries;
        if(std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()>=budget.seconds) throw Stop{"rank-time-budget"};
    }
    ~Audit() {result.seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();}
};
void Finite(double value) {if(!std::isfinite(value)) throw Stop{"rank-bound-overflow"};}
std::optional<std::size_t> CheckedAdd(std::size_t a,std::size_t b)
{
    if(b>std::numeric_limits<std::size_t>::max()-a) return std::nullopt;
    return a+b;
}
std::optional<std::size_t> CheckedMultiply(std::size_t a,std::size_t b)
{
    if(a && b>std::numeric_limits<std::size_t>::max()/a) return std::nullopt;
    return a*b;
}
std::optional<std::size_t> WorkspaceSize(std::size_t rows,std::size_t columns,std::size_t nonzeros,bool local_witness)
{
    const auto row_bytes=CheckedMultiply(rows,3*sizeof(Interval));
    const auto column_bytes=CheckedMultiply(columns,192);
    if(!row_bytes || !column_bytes) return std::nullopt;
    auto total=CheckedAdd(*row_bytes,*column_bytes);
    if(!total || !local_witness) return total;
    const auto census_rows=CheckedMultiply(rows,4*sizeof(std::size_t));
    const auto census_design=CheckedMultiply(nonzeros,sizeof(Eigen::Index));
    const auto census_columns=CheckedMultiply(columns,sizeof(std::size_t));
    if(!census_rows || !census_design || !census_columns) return std::nullopt;
    total=CheckedAdd(*total,*census_rows);
    if(!total) return std::nullopt;
    total=CheckedAdd(*total,*census_design);
    if(!total) return std::nullopt;
    return CheckedAdd(*total,*census_columns);
}
FreeDesignLocalWitness DiagnoseLocalSupport(const Sparse & design,double threshold_upper,Audit * audit)
{
    FreeDesignLocalWitness out;
    out.total_columns=static_cast<std::size_t>(design.cols());
    out.threshold_upper=threshold_upper;
    if(!std::isfinite(threshold_upper)) {out.reason="rank-threshold-unavailable"; return out;}
    const auto nonzeros=static_cast<std::size_t>(design.nonZeros());
    if(audit)
    {
        audit->Stage(FreeDesignRankWorkStage::LocalWitness);
        const auto scans=CheckedMultiply(nonzeros,2);
        if(!scans) throw Stop{"rank-work-budget"};
        audit->Charge(*scans);
    }
    std::map<std::vector<Eigen::Index>,std::vector<Eigen::Index>> support_groups;
    for(Eigen::Index col=0;col<design.cols();++col)
    {
        if(audit) audit->Charge(0);
        std::vector<Eigen::Index> support;
        for(Sparse::InnerIterator value(design,col);value;++value) if(value.value()!=0) support.push_back(value.row());
        support_groups[std::move(support)].push_back(col);
    }
    out.groups=support_groups.size();
    const auto missing_group=std::numeric_limits<std::size_t>::max();
    std::vector<std::size_t> group_for_column(static_cast<std::size_t>(design.cols()));
    std::size_t group_id{};
    for(const auto & [support,columns]:support_groups)
    {
        (void)support;
        for(const auto col:columns) group_for_column[static_cast<std::size_t>(col)]=group_id;
        ++group_id;
    }
    std::vector<std::size_t> row_group(static_cast<std::size_t>(design.rows()),missing_group);
    std::vector<bool> row_conflict(static_cast<std::size_t>(design.rows()));
    for(Eigen::Index col=0;col<design.cols();++col)
    {
        if(audit) audit->Charge(0);
        for(Sparse::InnerIterator value(design,col);value;++value) if(value.value()!=0)
        {
            const auto row=static_cast<std::size_t>(value.row());
            const auto current=group_for_column[static_cast<std::size_t>(col)];
            if(row_group[row]==missing_group) row_group[row]=current;
            else if(row_group[row]!=current) row_conflict[row]=true;
        }
    }

    group_id=0;
    bool unsupported_size=false,missing_rows=false,weak_lower=false,bound_unavailable=false;
    for(const auto & [support,columns]:support_groups)
    {
        if(audit) audit->Charge(0);
        out.max_group_size=std::max(out.max_group_size,columns.size());
        std::vector<Eigen::Index> exclusive;
        for(const auto row:support)
        {
            const auto index=static_cast<std::size_t>(row);
            if(!row_conflict[index] && row_group[index]==group_id) exclusive.push_back(row);
        }
        out.exclusive_rows+=exclusive.size();
        if(columns.size()>2) {unsupported_size=true; ++group_id; continue;}
        if(exclusive.size()<(columns.size()==1 ? 1u : 2u)) {missing_rows=true; ++group_id; continue;}
        std::optional<double> best;
        if(columns.size()==1)
        {
            if(audit) audit->Charge(exclusive.size());
            for(const auto row:exclusive)
            {
                const auto lower=CertifiedMagnitudeLowerBound(design.coeff(row,columns[0]));
                if(lower && (!best || *lower>*best)) best=lower;
            }
        }
        else
        {
            std::vector<std::array<double,2>> row_values;
            row_values.reserve(exclusive.size());
            if(audit) audit->Charge(2*exclusive.size());
            for(const auto row:exclusive)
                row_values.push_back({design.coeff(row,columns[0]),design.coeff(row,columns[1])});
            std::array<std::size_t,3> anchors{};
            for(std::size_t k=0;k<row_values.size();++k)
            {
                if(std::abs(row_values[k][0])>std::abs(row_values[anchors[0]][0])) anchors[0]=k;
                if(std::abs(row_values[k][1])>std::abs(row_values[anchors[1]][1])) anchors[1]=k;
                const auto magnitude=[](const auto & row) {return std::max(std::abs(row[0]),std::abs(row[1]));};
                if(magnitude(row_values[k])>magnitude(row_values[anchors[2]])) anchors[2]=k;
            }
            for(std::size_t a=0;a<anchors.size();++a)
            {
                if(std::find(anchors.begin(),anchors.begin()+static_cast<std::ptrdiff_t>(a),anchors[a])!=
                    anchors.begin()+static_cast<std::ptrdiff_t>(a)) continue;
                for(std::size_t k=0;k<row_values.size();++k) if(k!=anchors[a])
                {
                    if(audit) audit->Charge(4);
                    const auto lower=CertifiedSmallestSingularLowerBound2x2(
                        row_values[anchors[a]][0],row_values[anchors[a]][1],row_values[k][0],row_values[k][1]);
                    if(lower && (!best || *lower>*best)) best=lower;
                }
            }
        }
        if(!best) bound_unavailable=true;
        else
        {
            out.minimum_lower=out.minimum_lower ? std::min(*out.minimum_lower,*best) : *best;
            if(*best>threshold_upper) out.covered_columns+=columns.size();
            else weak_lower=true;
        }
        ++group_id;
    }
    out.coverage_fraction=out.total_columns ?
        static_cast<double>(out.covered_columns)/static_cast<double>(out.total_columns) : 0.;
    out.would_certify=out.total_columns>0 && out.covered_columns==out.total_columns && out.exclusive_rows_disjoint;
    if(out.would_certify) out.reason="local-support-full-rank-witness";
    else if(unsupported_size) out.reason="unsupported-for-local-witness";
    else if(missing_rows) out.reason="missing-exclusive-rows";
    else if(bound_unavailable) out.reason="local-bound-unavailable";
    else if(weak_lower) out.reason="local-lower-not-above-threshold";
    else out.reason="partial-column-coverage";
    return out;
}
}
std::string_view FreeDesignRankWorkStageName(FreeDesignRankWorkStage stage)
{
    switch(stage)
    {
    case FreeDesignRankWorkStage::None: return "none";
    case FreeDesignRankWorkStage::StructuralScan: return "structural-scan";
    case FreeDesignRankWorkStage::DuplicateCheck: return "duplicate-check";
    case FreeDesignRankWorkStage::LocalWitness: return "local-witness";
    case FreeDesignRankWorkStage::FactorInspection: return "factor-inspection";
    case FreeDesignRankWorkStage::WeakDirection: return "weak-direction";
    case FreeDesignRankWorkStage::InverseBound: return "inverse-bound";
    case FreeDesignRankWorkStage::OrthogonalBound: return "orthogonal-bound";
    case FreeDesignRankWorkStage::Reconstruction: return "reconstruction";
    }
    return "none";
}
std::string_view FreeDesignRankCertificateName(FreeDesignRankCertificate certificate)
{
    switch(certificate)
    {
    case FreeDesignRankCertificate::None: return "none";
    case FreeDesignRankCertificate::Structural: return "structural";
    case FreeDesignRankCertificate::WeakDirection: return "weak-direction";
    case FreeDesignRankCertificate::LocalSupport: return "local-support";
    case FreeDesignRankCertificate::SpqrReconstruction: return "spqr-reconstruction";
    case FreeDesignRankCertificate::DenseOracle: return "dense-oracle";
    }
    return "none";
}
std::optional<double> CertifiedMagnitudeLowerBound(double value)
{
    if(!std::isfinite(value)) return std::nullopt;
    try {return AbsLower(Interval(value));}
    catch(const Stop &) {return std::nullopt;}
}
std::optional<double> CertifiedSmallestSingularLowerBound2x2(double a,double b,double c,double d)
{
    if(!std::isfinite(a) || !std::isfinite(b) || !std::isfinite(c) || !std::isfinite(d)) return std::nullopt;
    try
    {
        const auto determinant=Add(Mul(Interval(a),Interval(d)),Neg(Mul(Interval(b),Interval(c))));
        const double determinant_lower=AbsLower(determinant);
        Interval frobenius_squared;
        for(const auto value:{a,b,c,d}) frobenius_squared=Add(frobenius_squared,Mul(Interval(value),Interval(value)));
        const double frobenius_upper=SqrtUp(frobenius_squared.hi);
        if(!std::isfinite(frobenius_upper)) return std::nullopt;
        if(frobenius_upper==0 || determinant_lower==0) return 0.;
        return std::max(0.,Down(determinant_lower/frobenius_upper));
    }
    catch(const Stop &) {return std::nullopt;}
}
FreeDesignLocalWitness DiagnoseFreeDesignLocalWitnesses(const Sparse & design,double threshold_upper)
{
    return DiagnoseLocalSupport(design,threshold_upper,nullptr);
}
FreeDesignRankResult EvaluateFreeDesignRank(const Sparse & z,const FreeDesignFactor * factor,const RankRequest & request,const RankBudget & budget)
{
    FreeDesignRankResult out; out.rank_upper=std::min(z.rows(),z.cols());
    // Keep the timer inside a scope so its destructor runs before the return copy.
    {
    Audit audit{budget,out};
    try {
#ifdef __FAST_MATH__
        throw Stop{"rank-unsupported-arithmetic"};
#endif
        if(!std::numeric_limits<double>::is_iec559 || std::fegetround()!=FE_TONEAREST) throw Stop{"rank-unsupported-arithmetic"};
        const auto n=z.rows(),p=z.cols();
        if(n<=0 || p<=0 || request.policy.rows<0 || request.columns<=0 || !std::isfinite(request.absolute) ||
            !std::isfinite(budget.seconds) || budget.seconds<0) throw Stop{"rank-invalid-input"};
        out.design_nonzeros=static_cast<std::size_t>(z.nonZeros());
        const auto base_workspace=WorkspaceSize(static_cast<std::size_t>(n),static_cast<std::size_t>(p),out.design_nonzeros,false);
        out.workspace_bytes=base_workspace.value_or(std::numeric_limits<std::size_t>::max());
        if(out.workspace_bytes>budget.workspace_bytes) throw Stop{"rank-memory-budget"};
        RecordDenseShape("rank-observation-vector",n,1); RecordDenseShape("rank-free-vector",p,1);
        audit.Stage(FreeDesignRankWorkStage::StructuralScan);
        audit.Charge(static_cast<std::size_t>(z.nonZeros()));
        std::vector<double> row_sums(static_cast<std::size_t>(n));
        double frobenius{},one{},maximum_lower{};
        bool structural=n<p; Eigen::Index structural_upper=std::min(n,p);
        std::unordered_map<std::uint64_t,std::vector<Eigen::Index>> columns;
        for(Eigen::Index col=0;col<p;++col)
        {
            audit.Stage(FreeDesignRankWorkStage::StructuralScan);
            audit.Charge(0);
            Interval square; double sum{}; std::uint64_t hash=1469598103934665603ULL;
            bool zero=true;
            for(Sparse::InnerIterator v(z,col);v;++v)
            {
                if(!std::isfinite(v.value())) throw Stop{"rank-nonfinite-design"};
                square=Add(square,Mul(Interval(v.value()),Interval(v.value())));
                sum=SumUp(sum,std::abs(v.value()));
                At(row_sums,v.row())=SumUp(At(row_sums,v.row()),std::abs(v.value()));
                if(v.value()!=0) {zero=false; hash=(hash^static_cast<std::uint64_t>(v.row()))*1099511628211ULL; hash=(hash^std::bit_cast<std::uint64_t>(v.value()))*1099511628211ULL;}
            }
            maximum_lower=std::max(maximum_lower,SqrtDown(square.lo));
            frobenius=SumUp(frobenius,std::max(0.,square.hi)); one=std::max(one,sum);
            if(zero) {structural=true; structural_upper=std::min(structural_upper,p-1);}
            for(const auto previous:columns[hash])
            {
                audit.Stage(FreeDesignRankWorkStage::DuplicateCheck);
                audit.Charge(static_cast<std::size_t>(z.col(col).nonZeros()+z.col(previous).nonZeros()));
                {
                    // Compare entries exactly: a norm of the difference may underflow.
                    Sparse::InnerIterator a(z,col),b(z,previous); bool same=true;
                    while(a || b)
                    {
                        while(a && a.value()==0) ++a; while(b && b.value()==0) ++b;
                        if(!a || !b) {same=!a && !b; break;}
                        if(a.row()!=b.row() || a.value()!=b.value()) {same=false; break;}
                        ++a; ++b;
                    }
                    if(same) {structural=true; structural_upper=std::min(structural_upper,p-1); break;}
                }
            }
            columns[hash].push_back(col);
        }
        out.maximum_lower=maximum_lower;
        out.maximum_upper=std::min(SqrtUp(frobenius),SqrtUp(ProductUp(one,*std::max_element(row_sums.begin(),row_sums.end()))));
        Finite(out.maximum_upper);
        const double relative=request.policy.Relative(request.columns);
        if(!std::isfinite(relative) || relative<0) throw Stop{"rank-invalid-policy"};
        if(request.absolute>=0 && out.maximum_upper>0)
        {
            // Dense SVD computes (absolute / sigma_max) * sigma_max. Enclose
            // its two normal rounded operations, rather than assuming exact cancellation.
            if(request.absolute>0 && (out.maximum_lower==0 || request.absolute/out.maximum_upper<std::numeric_limits<double>::min()))
                throw Stop{"rank-subnormal-cutoff-unresolved"};
            out.threshold_lower=out.threshold_upper=request.absolute;
            if(request.absolute>0) for(int k=0;k<4;++k)
            {out.threshold_lower=std::max(0.,Down(out.threshold_lower)); out.threshold_upper=Up(out.threshold_upper);}
        }
        else
        {
            out.threshold_lower=std::max(0.,Mul(Interval(relative),Interval(out.maximum_lower)).lo);
            out.threshold_upper=ProductUp(relative,out.maximum_upper);
        }
        const double cutoff=request.boundary==RankBoundary::SvdNative ?
            std::max(out.threshold_upper,std::numeric_limits<double>::min()) : out.threshold_upper;
        if(structural)
        {
            if(n>=p && out.maximum_upper>0 && out.threshold_lower==0) throw Stop{"rank-boundary-unresolved"};
            out.status=FreeDesignRankStatus::Deficient; out.certificate=FreeDesignRankCertificate::Structural; out.reason="rank-structural-deficiency";
            out.rank_upper=out.maximum_upper==0 ? 0 : structural_upper;
            out.witness_upper=0;
        }
        else
        {
            audit.Stage(FreeDesignRankWorkStage::LocalWitness);
            const auto witness_workspace=WorkspaceSize(static_cast<std::size_t>(n),static_cast<std::size_t>(p),out.design_nonzeros,true);
            out.workspace_bytes=witness_workspace.value_or(std::numeric_limits<std::size_t>::max());
            if(out.workspace_bytes>budget.workspace_bytes) throw Stop{"rank-memory-budget"};
            const auto local_started=std::chrono::steady_clock::now();
            out.local_witness=DiagnoseLocalSupport(z,cutoff,&audit);
            out.local_witness_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-local_started).count();
            if(out.local_witness.would_certify)
            {
                out.status=FreeDesignRankStatus::FullRank; out.certificate=FreeDesignRankCertificate::LocalSupport;
                out.rank_lower=out.rank_upper=p; out.reason="rank-verified-full";
                out.minimum_lower=out.local_witness.minimum_lower.value_or(0.);
                out.estimated_total_entries=out.entries; out.estimated_remaining_entries=0;
            }
        }
        if(!structural && out.status==FreeDesignRankStatus::Unavailable)
        {
            audit.Stage(FreeDesignRankWorkStage::FactorInspection);
            if(!factor) throw Stop{"rank-factor-unavailable"};
            const auto view=factor->RankView(); if(!view) throw Stop{"rank-backend-unavailable"};
            const auto & f=*view;
            out.factor_r_nonzeros=f.r_values.size(); out.reflector_nonzeros=f.h_values.size();
            out.reflector_count=static_cast<std::size_t>(f.reflectors);
            if(f.rows!=n || f.columns!=p || !f.design || f.design->rows()!=n || f.design->cols()!=p ||
                f.design->nonZeros()!=z.nonZeros()) throw Stop{"rank-factor-mismatch"};
            // An exact entry comparison prevents an underflowed difference from matching.
            for(Eigen::Index col=0;col<p;++col)
            {
                Sparse::InnerIterator a(z,col),b(*f.design,col);
                for(;a && b;++a,++b) if(a.row()!=b.row() || a.value()!=b.value()) throw Stop{"rank-factor-mismatch"};
                if(a || b) throw Stop{"rank-factor-mismatch"};
            }
            const auto original=[&](Eigen::Index col) {return f.permutation.empty() ? col : static_cast<Eigen::Index>(At(f.permutation,col));};
            Vector diagonal=Vector::Zero(p);
            for(Eigen::Index col=0;col<p;++col) for(auto k=At(f.r_outer,col);k<At(f.r_outer,col+1);++k)
            {
                Finite(At(f.r_values,k));
                if(At(f.r_inner,k)>col && At(f.r_values,k)!=0) throw Stop{"rank-nontriangular-factor"};
                if(At(f.r_inner,k)==col) diagonal(col)=At(f.r_values,k);
            }
            audit.Charge(f.r_values.size());
            audit.Stage(FreeDesignRankWorkStage::WeakDirection);
            std::vector<Eigen::Index> pivots(static_cast<std::size_t>(p)); std::iota(pivots.begin(),pivots.end(),0);
            std::partial_sort(pivots.begin(),pivots.begin()+std::min<Eigen::Index>(3,p),pivots.end(),[&](auto a,auto b) {return std::abs(diagonal(a))<std::abs(diagonal(b)) || (std::abs(diagonal(a))==std::abs(diagonal(b)) && a<b);});
            for(Eigen::Index attempt=0;attempt<std::min<Eigen::Index>(3,p);++attempt)
            {
                const auto pivot=At(pivots,attempt); Vector v=Vector::Zero(p); v(pivot)=1;
                for(auto k=At(f.r_outer,pivot);k<At(f.r_outer,pivot+1);++k) if(At(f.r_inner,k)<pivot) v(At(f.r_inner,k))=-At(f.r_values,k);
                for(Eigen::Index col=pivot;col-->0;)
                {
                    if(diagonal(col)==0) {v(col)=unavailable; break;}
                    v(col)/=diagonal(col);
                    for(auto k=At(f.r_outer,col);k<At(f.r_outer,col+1);++k) if(At(f.r_inner,k)<col) v(At(f.r_inner,k))-=At(f.r_values,k)*v(col);
                }
                audit.Charge(f.r_values.size()); if(!v.allFinite()) continue;
                v/=v.cwiseAbs().maxCoeff();
                std::vector<Interval> action(static_cast<std::size_t>(n)); Interval norm;
                for(Eigen::Index col=0;col<p;++col)
                {
                    norm=Add(norm,Mul(Interval(v(col)),Interval(v(col))));
                    for(Sparse::InnerIterator zc(z,original(col));zc;++zc) At(action,zc.row())=Add(At(action,zc.row()),Mul(Interval(zc.value()),Interval(v(col))));
                }
                audit.Charge(static_cast<std::size_t>(z.nonZeros())); double sum{};
                for(const auto value:action) sum=SumUp(sum,ProductUp(AbsUpper(value),AbsUpper(value)));
                const double lower=SqrtDown(norm.lo);
                const double upper=lower>0 ? Up(SqrtUp(sum)/lower) : inf;
                out.witness_upper=std::isfinite(out.witness_upper) ? std::min(out.witness_upper,upper) : upper;
                if(upper<out.threshold_lower)
                {
                    out.status=FreeDesignRankStatus::Deficient; out.certificate=FreeDesignRankCertificate::WeakDirection;
                    out.rank_upper=p-1; out.reason="rank-verified-weak-direction"; break;
                }
            }
            if(out.status==FreeDesignRankStatus::Unavailable)
            {
                // |R^-1| <= M(R)^-1. Two nonnegative triangular solves bound
                // the infinity and one norms; no inverse or normal matrix exists.
                audit.Stage(FreeDesignRankWorkStage::InverseBound);
                std::vector<double> x(static_cast<std::size_t>(p),1),y(static_cast<std::size_t>(p),1);
                for(Eigen::Index col=p;col-->0;)
                {
                    if(diagonal(col)==0) throw Stop{"rank-boundary-unresolved"};
                    At(x,col)=Up(At(x,col)/std::abs(diagonal(col))); Finite(At(x,col));
                    for(auto k=At(f.r_outer,col);k<At(f.r_outer,col+1);++k) if(At(f.r_inner,k)<col)
                        At(x,At(f.r_inner,k))=SumUp(At(x,At(f.r_inner,k)),ProductUp(std::abs(At(f.r_values,k)),At(x,col)));
                }
                for(Eigen::Index col=0;col<p;++col)
                {
                    for(auto k=At(f.r_outer,col);k<At(f.r_outer,col+1);++k) if(At(f.r_inner,k)<col)
                        At(y,col)=SumUp(At(y,col),ProductUp(std::abs(At(f.r_values,k)),At(y,At(f.r_inner,k))));
                    At(y,col)=Up(At(y,col)/std::abs(diagonal(col))); Finite(At(y,col));
                }
                audit.Charge(2*f.r_values.size());
                const double inverse_upper=SqrtUp(ProductUp(*std::max_element(x.begin(),x.end()),*std::max_element(y.begin(),y.end())));
                Finite(inverse_upper);
                if(Down(1/inverse_upper)<=out.threshold_upper) throw Stop{"rank-bound-too-wide"};
                audit.Stage(FreeDesignRankWorkStage::OrthogonalBound);
                double qlower=1;
                for(Eigen::Index h=0;h<f.reflectors;++h)
                {
                    Interval norm;
                    for(auto k=At(f.h_outer,h);k<At(f.h_outer,h+1);++k) norm=Add(norm,Mul(Interval(At(f.h_values,k)),Interval(At(f.h_values,k))));
                    const auto eigenvalue=Add(Interval(1),Neg(Mul(Interval(At(f.tau,h)),norm)));
                    qlower=std::max(0.,Down(qlower*std::min(1.,AbsLower(eigenvalue))));
                }
                audit.Charge(f.h_values.size()); out.orthogonal_minimum=qlower;
                // Q = P_H' H_0 ... H_k. Reconstruct each R column with interval
                // Householder actions, retaining only one observation vector.
                audit.Stage(FreeDesignRankWorkStage::Reconstruction);
                const auto twice_reflector_nonzeros=CheckedMultiply(2,out.reflector_nonzeros);
                const auto per_column=twice_reflector_nonzeros ? CheckedAdd(*twice_reflector_nonzeros,static_cast<std::size_t>(n)) : std::nullopt;
                const auto reconstruction=per_column ? CheckedMultiply(static_cast<std::size_t>(p),*per_column) : std::nullopt;
                const auto total=reconstruction ? CheckedAdd(out.entries,*reconstruction) : std::nullopt;
                if(total)
                {
                    out.estimated_reconstruction_entries=*reconstruction;
                    out.estimated_total_entries=*total;
                    out.estimated_remaining_entries=*reconstruction;
                }
                double error_squared{};
                std::vector<Interval> column(static_cast<std::size_t>(n));
                for(Eigen::Index col=0;col<p;++col)
                {
                    std::fill(column.begin(),column.end(),Interval{});
                    for(auto k=At(f.r_outer,col);k<At(f.r_outer,col+1);++k) At(column,At(f.r_inner,k))=Interval(At(f.r_values,k));
                    for(Eigen::Index h=f.reflectors;h-->0;)
                    {
                        Interval dot;
                        for(auto k=At(f.h_outer,h);k<At(f.h_outer,h+1);++k) dot=Add(dot,Mul(Interval(At(f.h_values,k)),At(column,At(f.h_inner,k))));
                        const auto scaled=Mul(Interval(At(f.tau,h)),dot);
                        for(auto k=At(f.h_outer,h);k<At(f.h_outer,h+1);++k) At(column,At(f.h_inner,k))=Add(At(column,At(f.h_inner,k)),Neg(Mul(Interval(At(f.h_values,k)),scaled)));
                        audit.Charge(static_cast<std::size_t>(2*(At(f.h_outer,h+1)-At(f.h_outer,h))));
                    }
                    for(Sparse::InnerIterator zc(z,original(col));zc;++zc)
                    {
                        const auto row=f.row_permutation.empty() ? zc.row() : At(f.row_permutation,zc.row());
                        At(column,row)=Add(At(column,row),Interval(-zc.value()));
                    }
                    for(const auto value:column) error_squared=SumUp(error_squared,ProductUp(AbsUpper(value),AbsUpper(value)));
                    audit.Charge(static_cast<std::size_t>(n)); Finite(error_squared);
                }
                out.reconstruction_error=SqrtUp(error_squared);
                out.minimum_lower=std::max(0.,Down(Down(qlower/inverse_upper)-out.reconstruction_error));
                if(out.minimum_lower>cutoff)
                {
                    out.status=FreeDesignRankStatus::FullRank; out.certificate=FreeDesignRankCertificate::SpqrReconstruction;
                    out.rank_lower=out.rank_upper=p; out.reason="rank-verified-full";
                }
                else out.reason="rank-boundary-unresolved";
            }
        }
        audit.Charge(0); // Structural and weak-direction early decisions also obey the deadline.
    } catch(const Stop & stop) {out.status=FreeDesignRankStatus::Unavailable; out.reason=stop.reason;}
      catch(const std::exception &) {out.status=FreeDesignRankStatus::Unavailable; out.reason="rank-factor-verification-failed";}
    }
    return out;
}
}
