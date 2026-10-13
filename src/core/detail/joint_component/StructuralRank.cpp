#include "Numerics.hpp"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numeric>
#include <new>
#include <utility>

namespace rhbm_gem::core::joint_component {
namespace {
struct Edge {Eigen::Index row{},column{};};

StructuralRankEvidence Unavailable(Eigen::Index rows,Eigen::Index columns,const char * reason)
{
    StructuralRankEvidence out; out.rows=rows; out.columns=columns; out.reason=reason; return out;
}

template<class Visit>
StructuralRankEvidence Analyze(Eigen::Index rows,Eigen::Index columns,std::size_t nonzeros,
    std::size_t resident_bytes,std::size_t budget,bool attempt_matching,Visit visit)
{
    auto unavailable_result=Unavailable(rows,columns,"structural-rank-unavailable");
    if(rows<0 || columns<=0) {unavailable_result.reason="invalid-design-shape"; return unavailable_result;}
    MemoryShape shape; shape.component_observations=static_cast<std::size_t>(rows);
    shape.component_design_columns=static_cast<std::size_t>(columns);
    shape.component_structural_nonzeros=nonzeros;
    shape.component_structural_supported_rows=std::min(shape.component_observations,nonzeros);
    shape.resident_bytes=resident_bytes;
    const auto base_memory=PlanMemory(shape,MemoryStage::StructuralRank,budget);
    unavailable_result.memory_estimate_bytes=base_memory.estimated_peak_bytes;
    if(base_memory.status!=MemoryStatus::Pass)
    {
        unavailable_result.reason=base_memory.status==MemoryStatus::Overflow ?
            "structural-rank-memory-overflow" : "structural-rank-memory-limit";
        return unavailable_result;
    }

    try
    {
        std::vector<Edge> edges; edges.reserve(nonzeros);
        std::vector<std::uint8_t> nonzero_columns(static_cast<std::size_t>(columns));
        visit([&](Eigen::Index row,Eigen::Index column) {
            if(row<0 || row>=rows || column<0 || column>=columns) return;
            edges.push_back({row,column}); nonzero_columns[static_cast<std::size_t>(column)]=1;
        });
        std::sort(edges.begin(),edges.end(),[](const Edge & a,const Edge & b) {
            return a.row<b.row || (a.row==b.row && a.column<b.column);
        });
        edges.erase(std::unique(edges.begin(),edges.end(),[](const Edge & a,const Edge & b) {
            return a.row==b.row && a.column==b.column;
        }),edges.end());

        StructuralRankEvidence out; out.available=true; out.rows=rows; out.columns=columns;
        out.nonzero_entries=static_cast<Eigen::Index>(edges.size());
        std::vector<Eigen::Index> supported_rows; supported_rows.reserve(std::min<std::size_t>(edges.size(),shape.component_observations));
        Eigen::Index last_row=-1;
        for(const auto & edge:edges) if(edge.row!=last_row)
        {supported_rows.push_back(edge.row); last_row=edge.row;}
        out.supported_rows=static_cast<Eigen::Index>(supported_rows.size());
        out.supported_columns=static_cast<Eigen::Index>(std::count(nonzero_columns.begin(),nonzero_columns.end(),std::uint8_t{1}));
        out.rank_upper_bound=std::min(out.supported_rows,out.supported_columns);
        out.memory_estimate_bytes=base_memory.estimated_peak_bytes;
        if(out.supported_rows<columns)
        {
            out.proven_deficient=true; out.reason="supported-row-count-below-design-columns"; return out;
        }
        if(out.supported_columns<columns)
        {
            out.proven_deficient=true; out.reason="zero-design-column"; return out;
        }
        if(!attempt_matching) {out.reason="matching-not-requested"; return out;}

        shape.component_structural_nonzeros=edges.size();
        shape.component_structural_supported_rows=supported_rows.size();
        shape.component_structural_matching=true;
        const auto matching_memory=PlanMemory(shape,MemoryStage::StructuralRank,budget);
        out.memory_estimate_bytes=matching_memory.estimated_peak_bytes;
        if(matching_memory.status!=MemoryStatus::Pass)
        {
            out.reason=matching_memory.status==MemoryStatus::Overflow ?
                "structural-matching-memory-overflow" : "structural-matching-memory-limit";
            return out;
        }

        for(auto & edge:edges)
            edge.row=static_cast<Eigen::Index>(std::lower_bound(supported_rows.begin(),supported_rows.end(),edge.row)-supported_rows.begin());
        std::sort(edges.begin(),edges.end(),[](const Edge & a,const Edge & b) {
            return a.column<b.column || (a.column==b.column && a.row<b.row);
        });
        std::vector<std::size_t> offsets(static_cast<std::size_t>(columns)+1);
        std::size_t edge_index{};
        for(Eigen::Index column=0;column<columns;++column)
        {
            offsets[static_cast<std::size_t>(column)]=edge_index;
            while(edge_index<edges.size() && edges[edge_index].column==column) ++edge_index;
        }
        offsets.back()=edge_index;

        std::vector<Eigen::Index> column_order(static_cast<std::size_t>(columns));
        std::iota(column_order.begin(),column_order.end(),Eigen::Index{0});
        const auto pattern_less=[&](Eigen::Index left,Eigen::Index right) {
            auto a=offsets[static_cast<std::size_t>(left)],b=offsets[static_cast<std::size_t>(right)];
            const auto a_end=offsets[static_cast<std::size_t>(left)+1],b_end=offsets[static_cast<std::size_t>(right)+1];
            while(a<a_end && b<b_end)
            {
                if(edges[a].row!=edges[b].row) return edges[a].row<edges[b].row;
                ++a; ++b;
            }
            return a==a_end && b!=b_end;
        };
        std::sort(column_order.begin(),column_order.end(),pattern_less);
        for(std::size_t first=0;first<column_order.size();)
        {
            std::size_t end=first+1;
            while(end<column_order.size() && !pattern_less(column_order[first],column_order[end]) &&
                !pattern_less(column_order[end],column_order[first])) ++end;
            if(end-first>1) ++out.identical_support_pattern_groups;
            first=end;
        }

        out.matching_available=true;
        const Eigen::Index infinity=std::numeric_limits<Eigen::Index>::max()/4;
        const auto row_count=static_cast<std::size_t>(out.supported_rows);
        std::vector<Eigen::Index> left_match(static_cast<std::size_t>(columns),-1),right_match(row_count,-1),distance(static_cast<std::size_t>(columns),-1);
        std::vector<std::size_t> cursor(static_cast<std::size_t>(columns));
        std::vector<Eigen::Index> queue; queue.reserve(static_cast<std::size_t>(columns));
        std::vector<Eigen::Index> stack; stack.reserve(static_cast<std::size_t>(columns));
        std::vector<Eigen::Index> via(static_cast<std::size_t>(columns),-1);
        Eigen::Index matching{};
        for(;;)
        {
            queue.clear(); std::fill(distance.begin(),distance.end(),-1);
            for(Eigen::Index column=0;column<columns;++column) if(left_match[static_cast<std::size_t>(column)]<0)
            {distance[static_cast<std::size_t>(column)]=0; queue.push_back(column);}
            Eigen::Index shortest=infinity;
            for(std::size_t head=0;head<queue.size();++head)
            {
                const auto column=queue[head]; const auto depth=distance[static_cast<std::size_t>(column)];
                if(depth>=shortest) continue;
                for(auto k=offsets[static_cast<std::size_t>(column)];k<offsets[static_cast<std::size_t>(column)+1];++k)
                {
                    const auto row=edges[k].row;
                    const auto mate=right_match[static_cast<std::size_t>(row)];
                    if(mate<0) shortest=depth+1;
                    else if(distance[static_cast<std::size_t>(mate)]<0)
                    {distance[static_cast<std::size_t>(mate)]=depth+1; queue.push_back(mate);}
                }
            }
            if(shortest==infinity) break;
            for(Eigen::Index column=0;column<columns;++column)
                cursor[static_cast<std::size_t>(column)]=offsets[static_cast<std::size_t>(column)];
            for(Eigen::Index root=0;root<columns;++root)
            {
                if(left_match[static_cast<std::size_t>(root)]>=0 || distance[static_cast<std::size_t>(root)]!=0) continue;
                stack.clear(); stack.push_back(root); bool augmented=false;
                while(!stack.empty() && !augmented)
                {
                    const auto column=stack.back(); const auto index=static_cast<std::size_t>(column);
                    bool descended=false;
                    while(cursor[index]<offsets[index+1])
                    {
                        const auto row=edges[cursor[index]++].row;
                        const auto mate=right_match[static_cast<std::size_t>(row)];
                        if(mate<0)
                        {
                            if(distance[index]+1!=shortest) continue;
                            auto assigned=row;
                            for(std::size_t depth=stack.size();depth>0;--depth)
                            {
                                const auto current=stack[depth-1];
                                left_match[static_cast<std::size_t>(current)]=assigned;
                                right_match[static_cast<std::size_t>(assigned)]=current;
                                if(depth>1) assigned=via[depth-1];
                            }
                            ++matching; augmented=true; break;
                        }
                        if(distance[static_cast<std::size_t>(mate)]==distance[index]+1 &&
                            distance[static_cast<std::size_t>(mate)]<shortest)
                        {
                            via[stack.size()]=row; stack.push_back(mate); descended=true; break;
                        }
                    }
                    if(augmented || descended) continue;
                    distance[index]=-1; stack.pop_back();
                }
            }
        }
        out.matching_size=matching;
        out.rank_upper_bound=std::min(out.rank_upper_bound,matching);
        out.proven_deficient=matching<columns;
        out.structurally_full_rank=matching==columns;
        out.reason=out.proven_deficient ? "maximum-matching-below-design-columns" : "structurally-full-matching";
        return out;
    }
    catch(const std::bad_alloc &)
    {
        unavailable_result.reason="structural-rank-allocation-failed"; return unavailable_result;
    }
}
}

StructuralRankEvidence AssessStructuralRank(const Sparse & x,std::size_t resident_bytes,std::size_t budget)
{
    if(x.rows()<0 || x.cols()<=0 || x.nonZeros()<0) return Unavailable(x.rows(),x.cols(),"invalid-design-shape");
    std::size_t nonzeros{};
    for(Eigen::Index column=0;column<x.outerSize();++column)
        for(Sparse::InnerIterator entry(x,column);entry;++entry) if(entry.value()!=0)
        {
            if(nonzeros==std::numeric_limits<std::size_t>::max())
                return Unavailable(x.rows(),x.cols(),"structural-nonzero-count-overflow");
            ++nonzeros;
        }
    return Analyze(x.rows(),x.cols(),nonzeros,resident_bytes,budget,true,[&](auto emit) {
        for(Eigen::Index column=0;column<x.outerSize();++column)
            for(Sparse::InnerIterator entry(x,column);entry;++entry)
                if(entry.value()!=0) emit(entry.row(),column);
    });
}

StructuralRankEvidence AssessDesignStructuralRank(const Domain & domain,VectorRef eta,std::size_t resident_bytes,
    std::size_t budget,bool attempt_matching)
{
    if(eta.size()>std::numeric_limits<Eigen::Index>::max()/2)
        return Unavailable(domain.rows,0,"invalid-design-shape");
    if(domain.rows<=0 || eta.size()<=0 || eta.size()!=static_cast<Eigen::Index>(domain.atoms.size()))
        return Unavailable(domain.rows,2*eta.size(),"invalid-design-shape");
    if(!eta.allFinite()) return Unavailable(domain.rows,2*eta.size(),"invalid-widths");
    std::size_t nonzeros{};
    for(Eigen::Index atom=0;atom<eta.size();++atom)
    {
        const double width=std::exp(eta(atom));
        if(!std::isfinite(width) || width<=0) return Unavailable(domain.rows,2*eta.size(),"invalid-widths");
        for(const auto & support:domain.atoms[static_cast<std::size_t>(atom)])
        {
            if(support.row<0 || support.row>=domain.rows || !std::isfinite(support.square) || support.square<0)
                return Unavailable(domain.rows,2*eta.size(),"invalid-support-row");
            const auto basis=EvaluateKernel(support.square,width,2.5);
            if(basis.gaussian!=0 && nonzeros==std::numeric_limits<std::size_t>::max())
                return Unavailable(domain.rows,2*eta.size(),"structural-nonzero-count-overflow");
            if(basis.gaussian!=0) ++nonzeros;
            if(basis.charge!=0 && nonzeros==std::numeric_limits<std::size_t>::max())
                return Unavailable(domain.rows,2*eta.size(),"structural-nonzero-count-overflow");
            if(basis.charge!=0) ++nonzeros;
        }
    }
    return Analyze(domain.rows,2*eta.size(),nonzeros,resident_bytes,budget,attempt_matching,[&](auto emit) {
        for(Eigen::Index atom=0;atom<eta.size();++atom)
        {
            const double width=std::exp(eta(atom));
            for(const auto & support:domain.atoms[static_cast<std::size_t>(atom)])
            {
                const auto basis=EvaluateKernel(support.square,width,2.5);
                if(basis.gaussian!=0) emit(support.row,2*atom);
                if(basis.charge!=0) emit(support.row,2*atom+1);
            }
        }
    });
}

} // namespace rhbm_gem::core::joint_component
