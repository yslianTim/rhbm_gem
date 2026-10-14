#include "StructuredProjectionResearch.hpp"
#include <Eigen/OrderingMethods>
#include <Eigen/SparseCholesky>
#include <algorithm>
#include <new>

namespace second_stage_test::projected_width {
namespace {
bool Add(std::size_t a,std::size_t b,std::size_t & out)
{return runtime::CheckedAddBytes(a,b,out);}
bool Product(std::size_t a,std::size_t b,std::size_t & out)
{return runtime::CheckedMultiplyBytes(a,b,out);}
bool AddProduct(std::size_t & total,std::size_t count,std::size_t width)
{
    std::size_t bytes{}; return Product(count,width,bytes) && Add(total,bytes,total);
}
void GraphMetrics(const runtime::Sparse & gram,StructuredProjectionResult & out)
{
    const auto columns=gram.cols();
    std::vector<std::vector<Eigen::Index>> adjacency(static_cast<std::size_t>(columns));
    for(Eigen::Index column=0;column<gram.outerSize();++column)
        for(runtime::Sparse::InnerIterator entry(gram,column);entry;++entry)
            if(entry.row()!=column && entry.value()!=0)
                adjacency[static_cast<std::size_t>(column)].push_back(entry.row());

    std::vector<std::uint8_t> visited(static_cast<std::size_t>(columns));
    for(Eigen::Index root=0;root<columns;++root)
    {
        const auto index=static_cast<std::size_t>(root);
        if(visited[index]) continue;
        ++out.exact_blocks; visited[index]=1;
        std::vector<Eigen::Index> pending{root};
        while(!pending.empty())
        {
            const auto current=pending.back(); pending.pop_back();
            for(const auto next:adjacency[static_cast<std::size_t>(current)])
            {
                auto & mark=visited[static_cast<std::size_t>(next)];
                if(!mark) {mark=1; pending.push_back(next);}
            }
        }
    }

    std::vector<Eigen::Index> discovery(static_cast<std::size_t>(columns),-1),
        low(static_cast<std::size_t>(columns),-1),parent(static_cast<std::size_t>(columns),-1);
    std::vector<std::uint8_t> articulation(static_cast<std::size_t>(columns));
    Eigen::Index clock{};
    const auto visit=[&](auto && self,Eigen::Index current)->void {
        const auto index=static_cast<std::size_t>(current);
        discovery[index]=low[index]=clock++;
        Eigen::Index children{};
        for(const auto next:adjacency[index])
        {
            const auto next_index=static_cast<std::size_t>(next);
            if(discovery[next_index]<0)
            {
                parent[next_index]=current; ++children; self(self,next);
                low[index]=std::min(low[index],low[next_index]);
                if(parent[index]<0 && children>1) articulation[index]=1;
                if(parent[index]>=0 && low[next_index]>=discovery[index]) articulation[index]=1;
            }
            else if(next!=parent[index]) low[index]=std::min(low[index],discovery[next_index]);
        }
    };
    for(Eigen::Index column=0;column<columns;++column)
        if(discovery[static_cast<std::size_t>(column)]<0) visit(visit,column);
    out.articulation_columns=std::count(articulation.begin(),articulation.end(),std::uint8_t{1});
    Eigen::Index best=columns;
    for(Eigen::Index removed=0;removed<columns;++removed)
    {
        if(!articulation[static_cast<std::size_t>(removed)]) continue;
        std::vector<std::uint8_t> seen(static_cast<std::size_t>(columns));
        Eigen::Index maximum{};
        for(Eigen::Index root=0;root<columns;++root)
        {
            if(root==removed || seen[static_cast<std::size_t>(root)]) continue;
            Eigen::Index size{}; seen[static_cast<std::size_t>(root)]=1;
            std::vector<Eigen::Index> pending{root};
            while(!pending.empty())
            {
                const auto current=pending.back(); pending.pop_back(); ++size;
                for(const auto next:adjacency[static_cast<std::size_t>(current)])
                {
                    if(next==removed) continue;
                    auto & mark=seen[static_cast<std::size_t>(next)];
                    if(!mark) {mark=1; pending.push_back(next);}
                }
            }
            maximum=std::max(maximum,size);
        }
        if(maximum<best) {best=maximum; out.separator_columns=1; out.interface_dimension=1;}
    }
    out.maximum_interior_columns=out.separator_columns ? best : columns;
}
}

StructuredProjectionPreflight PreflightStructuredProjection(std::size_t rows,std::size_t columns,
    std::size_t nonzeros,std::size_t resident_bytes,std::size_t budget)
{
    StructuredProjectionPreflight out; out.budget_bytes=budget;
    std::size_t input_bytes{},outer{},square_count{},factor_bound{},workspace{},vectors{},estimate{};
    if(!Add(columns,1,outer) ||
        !Product(nonzeros,sizeof(double)+sizeof(Eigen::Index),input_bytes) ||
        !AddProduct(input_bytes,outer,sizeof(Eigen::Index)) ||
        !Product(columns,columns,square_count) ||
        !Product(square_count,sizeof(double)+sizeof(Eigen::Index),factor_bound) ||
        !Product(factor_bound,4,workspace) ||
        !AddProduct(vectors,rows,5*sizeof(double)) ||
        !AddProduct(vectors,columns,16*sizeof(double)+8*sizeof(Eigen::Index)) ||
        !Add(resident_bytes ? resident_bytes : input_bytes,input_bytes,estimate) ||
        !Add(estimate,workspace,estimate) || !Add(estimate,vectors,estimate))
    {out.status="ResourceLimit"; out.reason="memory-estimate-overflow"; return out;}
    out.estimated_peak_bytes=estimate;
    if(estimate>budget)
    {out.status="ResourceLimit"; out.reason="estimated-factor-peak-exceeds-planning-budget"; return out;}
    out.status="Pass"; out.reason="worst-case-sparse-gram-and-fill-bound";
    return out;
}

StructuredProjectionResult ProjectWithSparseNormalEquations(const runtime::Sparse & design,
    const Eigen::Ref<const Vector> & observation_direction,std::size_t resident_bytes,std::size_t budget)
{
    StructuredProjectionResult out; out.rows=design.rows(); out.columns=design.cols();
    out.input_nonzeros=design.nonZeros(); out.current_rss_before_bytes=runtime::CurrentRssBytes();
    const auto started=std::chrono::steady_clock::now();
    const auto finish=[&] {
        out.seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
        out.current_rss_after_bytes=runtime::CurrentRssBytes();
        out.process_peak_rss_bytes=runtime::ProcessPeakRssBytes();
    };
    if(design.rows()<=0 || design.cols()<=0 || observation_direction.size()!=design.rows() ||
        !observation_direction.allFinite())
    {out.reason="invalid-shape-or-direction"; finish(); return out;}
    const auto preflight=PreflightStructuredProjection(static_cast<std::size_t>(design.rows()),
        static_cast<std::size_t>(design.cols()),static_cast<std::size_t>(design.nonZeros()),resident_bytes,budget);
    out.estimated_peak_bytes=preflight.estimated_peak_bytes;
    if(preflight.status!="Pass")
    {out.status=preflight.status; out.reason=preflight.reason; finish(); return out;}
    out.preflight_passed=true;

    try
    {
        runtime::Sparse normalized=design;
        Vector norms(design.cols());
        for(Eigen::Index column=0;column<design.cols();++column)
        {
            const double norm=design.col(column).norm();
            if(!(norm>0) || !std::isfinite(norm))
            {out.reason="zero-or-nonfinite-design-column"; finish(); return out;}
            norms(column)=norm;
            for(runtime::Sparse::InnerIterator entry(normalized,column);entry;++entry) entry.valueRef()/=norm;
        }
        const runtime::Sparse gram=normalized.transpose()*normalized;
        out.gram_nonzeros=gram.nonZeros();
        GraphMetrics(gram,out);

        Eigen::SimplicialLDLT<runtime::Sparse,Eigen::Lower,
            Eigen::AMDOrdering<runtime::Sparse::StorageIndex>> factor;
        factor.analyzePattern(gram);
        if(factor.info()!=Eigen::Success)
        {out.reason="symbolic-elimination-failed"; finish(); return out;}
        const auto & lower=factor.matrixL().nestedExpression();
        std::size_t symbolic_factor_nonzeros=static_cast<std::size_t>(lower.nonZeros())+
            static_cast<std::size_t>(design.cols());
        out.factor_nonzeros=static_cast<Eigen::Index>(symbolic_factor_nonzeros);
        for(Eigen::Index column=0;column<lower.outerSize();++column)
        {
            Eigen::Index count{};
            for(runtime::Sparse::InnerIterator entry(lower,column);entry;++entry) ++count;
            out.maximum_factor_size=std::max(out.maximum_factor_size,count+1);
        }
        factor.factorize(gram);
        if(factor.info()!=Eigen::Success)
        {out.reason="sparse-normal-equation-factorization-failed"; finish(); return out;}
        const auto diagonal=factor.vectorD();
        if(diagonal.size()==0 || !diagonal.allFinite() || diagonal.minCoeff()<=0)
        {out.reason="sparse-normal-equation-breakdown"; finish(); return out;}
        out.diagonal_condition_estimate=diagonal.maxCoeff()/diagonal.minCoeff();
        const Vector coefficients=factor.solve(normalized.transpose()*observation_direction);
        if(factor.info()!=Eigen::Success || !coefficients.allFinite())
        {out.reason="projection-solve-failed"; finish(); return out;}
        out.projection=normalized*coefficients;
        const Vector residual=observation_direction-out.projection;
        out.projection_residual=(normalized.transpose()*residual).norm();
        const Vector second_coefficients=factor.solve(normalized.transpose()*out.projection);
        if(factor.info()!=Eigen::Success || !second_coefficients.allFinite())
        {out.reason="projection-idempotence-solve-failed"; finish(); return out;}
        out.idempotence_error=(normalized*second_coefficients-out.projection).norm();
        out.valid=true; out.status="Completed"; out.reason="experimental-sparse-normal-equation-projection";
    }
    catch(const std::bad_alloc &)
    {out.status="ResourceLimit"; out.reason="structured-projection-allocation-failed";}
    finish(); return out;
}
}
