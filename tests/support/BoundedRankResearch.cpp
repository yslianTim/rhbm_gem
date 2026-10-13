#include "BoundedRankResearch.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <new>
#include <numeric>

namespace second_stage_test::bounded_rank {
namespace {
void AddCompensated(double value,double & sum,double & compensation)
{
    const double adjusted=value-compensation;
    const double next=sum+adjusted;
    compensation=(next-sum)-adjusted;
    sum=next;
}

std::size_t SparseBytes(const runtime::Sparse & x)
{
    std::size_t values{},outer{},total{};
    if(x.nonZeros()<0 || !runtime::CheckedMultiplyBytes(static_cast<std::size_t>(x.nonZeros()),
        sizeof(double)+sizeof(Eigen::Index),values) ||
        !runtime::CheckedAddBytes(static_cast<std::size_t>(x.outerSize()),1,outer) ||
        !runtime::CheckedMultiplyBytes(outer,sizeof(Eigen::Index),outer) ||
        !runtime::CheckedAddBytes(values,outer,total)) return std::numeric_limits<std::size_t>::max();
    return total;
}
}

StreamingRankEstimate ScalarGershgorin(const runtime::Sparse & x,const runtime::RankPolicy & policy,
    std::size_t resident_bytes,std::size_t budget)
{
    StreamingRankEstimate out; out.rows=x.rows(); out.columns=x.cols(); out.nonzero_entries=x.nonZeros();
    const auto started=std::chrono::steady_clock::now();
    out.current_rss_before_bytes=runtime::CurrentRssBytes();
    const auto finish=[&]() {
        out.seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
        out.current_rss_after_bytes=runtime::CurrentRssBytes();
        out.process_peak_rss_bytes=runtime::ProcessPeakRssBytes();
    };
    if(x.rows()<0 || x.cols()<=0 || x.nonZeros()<0 || policy.rows<0)
    {out.reason="invalid-design-shape"; finish(); return out;}
    if(!resident_bytes) resident_bytes=SparseBytes(x);
    if(resident_bytes==std::numeric_limits<std::size_t>::max())
    {out.reason="streaming-memory-overflow"; finish(); return out;}
    runtime::MemoryShape shape; shape.component_observations=static_cast<std::size_t>(x.rows());
    shape.component_design_columns=static_cast<std::size_t>(x.cols()); shape.resident_bytes=resident_bytes;
    const auto memory=runtime::PlanMemory(shape,runtime::MemoryStage::StreamingSpectralBounds,budget);
    out.estimated_peak_bytes=memory.estimated_peak_bytes;
    if(memory.status!=runtime::MemoryStatus::Pass)
    {out.reason=memory.status==runtime::MemoryStatus::Overflow ? "streaming-memory-overflow" : "streaming-memory-limit"; finish(); return out;}

    try
    {
        const auto rows=static_cast<std::size_t>(x.rows()),columns=static_cast<std::size_t>(x.cols());
        std::vector<double> norm_squared(columns),norm_compensation(columns),radius(columns),radius_compensation(columns),
            row_sum(rows),row_compensation(rows);
        std::vector<Eigen::Index> column_entries(columns),single_row(columns,-1);
        std::vector<std::uint8_t> single_unit(columns,1);
        for(Eigen::Index column=0;column<x.outerSize();++column)
            for(runtime::Sparse::InnerIterator entry(x,column);entry;++entry)
            {
                if(entry.value()==0) continue;
                const auto c=static_cast<std::size_t>(column),row=static_cast<std::size_t>(entry.row());
                if(column_entries[c]==0) {single_row[c]=entry.row(); single_unit[c]=std::abs(entry.value())==1.0;}
                else single_unit[c]=0;
                ++column_entries[c];
                const double square=entry.value()*entry.value();
                if(std::isfinite(square)) AddCompensated(square,norm_squared[c],norm_compensation[c]);
                else norm_squared[c]=std::numeric_limits<double>::infinity();
            }

        const double relative=policy.Relative(x.cols());
        out.relative_threshold=relative;
        bool singleton=true;
        for(std::size_t column=0;column<columns;++column)
            singleton=singleton && column_entries[column]==1 && single_unit[column]!=0;
        if(singleton)
        {
            std::sort(single_row.begin(),single_row.end());
            singleton=std::adjacent_find(single_row.begin(),single_row.end())==single_row.end();
        }
        if(singleton && relative<1.0)
        {
            out.available=true; out.lambda_min_lower=1.0; out.lambda_max_upper=1.0; out.bound_gap=0.0;
            out.threshold_upper=relative; out.gershgorin_suggests_full_rank=true;
            out.certified_full_rank=true; out.reason="exact-disjoint-unit-columns"; finish(); return out;
        }

        for(std::size_t column=0;column<columns;++column)
        {
            if(column_entries[column]==0) {out.reason="zero-column-needs-structural-evidence"; finish(); return out;}
            norm_squared[column]=std::sqrt(norm_squared[column]);
            if(!std::isfinite(norm_squared[column]) || norm_squared[column]<=0)
            {out.reason="column-norm-unavailable"; finish(); return out;}
        }
        for(Eigen::Index column=0;column<x.outerSize();++column)
        {
            const auto c=static_cast<std::size_t>(column);
            for(runtime::Sparse::InnerIterator entry(x,column);entry;++entry)
                if(entry.value()!=0)
                    AddCompensated(std::abs(entry.value())/norm_squared[c],row_sum[static_cast<std::size_t>(entry.row())],
                        row_compensation[static_cast<std::size_t>(entry.row())]);
        }
        for(Eigen::Index column=0;column<x.outerSize();++column)
        {
            const auto c=static_cast<std::size_t>(column);
            for(runtime::Sparse::InnerIterator entry(x,column);entry;++entry)
            {
                if(entry.value()==0) continue;
                const auto row=static_cast<std::size_t>(entry.row());
                const double value=std::abs(entry.value())/norm_squared[c];
                const double other=std::max(0.0,row_sum[row]-value);
                AddCompensated(value*other,radius[c],radius_compensation[c]);
            }
        }
        if(!std::all_of(row_sum.begin(),row_sum.end(),[](double v){return std::isfinite(v);}) ||
            !std::all_of(radius.begin(),radius.end(),[](double v){return std::isfinite(v);}))
        {out.reason="nonfinite-streaming-bound"; finish(); return out;}

        out.available=true;
        out.lambda_min_lower=1.0-radius.front();
        out.lambda_max_upper=1.0+radius.front();
        for(std::size_t column=1;column<columns;++column)
        {
            out.lambda_min_lower=std::min(out.lambda_min_lower,1.0-radius[column]);
            out.lambda_max_upper=std::max(out.lambda_max_upper,1.0+radius[column]);
        }
        out.bound_gap=out.lambda_max_upper-out.lambda_min_lower;
        out.threshold_upper=relative*std::sqrt(std::max(0.0,out.lambda_max_upper));
        out.gershgorin_suggests_full_rank=out.lambda_min_lower>relative*relative*out.lambda_max_upper;
        out.certified_full_rank=false;
        out.reason="experimental-floating-point-estimate-no-certified-enclosure";
        finish(); return out;
    }
    catch(const std::bad_alloc &)
    {out.reason="streaming-allocation-failed"; finish(); return out;}
}
}
