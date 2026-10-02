#pragma once
#include "core/detail/joint_component/FreeDesignRank.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <optional>
#include <string>
#include <vector>

namespace second_stage_test {
namespace n=rhbm_gem::core::joint_component;
struct LocalRankWitnessCensus
{
    std::size_t groups{},covered_columns{},total_columns{},exclusive_rows{},max_group_size{};
    double coverage_fraction{};
    std::optional<double> minimum_lower;
    double threshold_upper{};
    bool exclusive_rows_disjoint{true},would_certify{};
    std::string reason{"partial-column-coverage"};
};

// Groups are induced only by identical sparse column supports. Selected rows
// contain nonzeros from one group, so their union gives a block-separable row
// submatrix without relying on atom or topology identities. For selected rows
// S, Z'Z = S'S + R'R; hence sigma_min(Z) >= sigma_min(S).
inline LocalRankWitnessCensus DiagnoseLocalRankWitnesses(const n::Sparse & design,double threshold_upper)
{
    LocalRankWitnessCensus out; out.total_columns=static_cast<std::size_t>(design.cols());
    out.threshold_upper=threshold_upper;
    if(!std::isfinite(threshold_upper)) {out.reason="rank-threshold-unavailable"; return out;}
    std::vector<std::vector<Eigen::Index>> column_rows(static_cast<std::size_t>(design.cols()));
    std::vector<std::vector<Eigen::Index>> row_columns(static_cast<std::size_t>(design.rows()));
    for(Eigen::Index col=0;col<design.cols();++col)
        for(n::Sparse::InnerIterator value(design,col);value;++value) if(value.value()!=0)
        {
            column_rows[static_cast<std::size_t>(col)].push_back(value.row());
            row_columns[static_cast<std::size_t>(value.row())].push_back(col);
        }
    std::map<std::vector<Eigen::Index>,std::vector<Eigen::Index>> support_groups;
    for(Eigen::Index col=0;col<design.cols();++col)
        support_groups[column_rows[static_cast<std::size_t>(col)]].push_back(col);
    out.groups=support_groups.size();
    std::vector<std::size_t> group_for_column(static_cast<std::size_t>(design.cols()));
    std::size_t group_id{};
    for(const auto & [support,columns]:support_groups)
    {
        (void)support;
        for(const auto col:columns) group_for_column[static_cast<std::size_t>(col)]=group_id;
        ++group_id;
    }
    group_id=0;
    bool unsupported_size=false,missing_rows=false,weak_lower=false,bound_unavailable=false;
    for(const auto & [support,columns]:support_groups)
    {
        out.max_group_size=std::max(out.max_group_size,columns.size());
        std::vector<Eigen::Index> exclusive;
        for(const auto row:support)
        {
            const auto & members=row_columns[static_cast<std::size_t>(row)];
            if(!members.empty() && std::all_of(members.begin(),members.end(),[&](Eigen::Index col) {
                return group_for_column[static_cast<std::size_t>(col)]==group_id;
            })) exclusive.push_back(row);
        }
        out.exclusive_rows+=exclusive.size();
        if(columns.size()>2) {unsupported_size=true; ++group_id; continue;}
        if(exclusive.size()<(columns.size()==1 ? 1u : 2u)) {missing_rows=true; ++group_id; continue;}
        std::optional<double> best;
        if(columns.size()==1)
        {
            for(const auto row:exclusive)
            {
                const auto lower=n::CertifiedMagnitudeLowerBound(design.coeff(row,columns[0]));
                if(lower && (!best || *lower>*best)) best=lower;
            }
        }
        else
        {
            std::vector<std::array<double,2>> row_values; row_values.reserve(exclusive.size());
            for(const auto row:exclusive) row_values.push_back({design.coeff(row,columns[0]),design.coeff(row,columns[1])});
            for(std::size_t i=0;i<exclusive.size();++i) for(std::size_t j=i+1;j<exclusive.size();++j)
            {
                const auto lower=n::CertifiedSmallestSingularLowerBound2x2(
                    row_values[i][0],row_values[i][1],row_values[j][0],row_values[j][1]);
                if(lower && (!best || *lower>*best)) best=lower;
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
