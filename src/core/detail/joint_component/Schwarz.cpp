#include "OperatorSearch.hpp"
#include <algorithm>
#include <deque>
#include <numeric>

namespace rhbm_gem::core::joint_component {
namespace {
constexpr double eps=std::numeric_limits<double>::epsilon();
void Budget(std::size_t bytes,std::size_t limit)
{if(bytes>limit) throw std::runtime_error("preconditioner-resource-limit");}
std::size_t TopologyBytes(const PreconditionerPartition & p)
{
    std::size_t bytes=sizeof(p)+8*(p.layout.full_atoms.size()+p.layout.informative_rows.size())+32*p.atom_blocks.size();
    for(const auto & block:p.blocks) bytes+=sizeof(block)+block.id.size()+8*(block.core_atoms.size()+block.overlap_atoms.size()+block.informative_rows.size());
    for(const auto & membership:p.atom_blocks) bytes+=8*membership.size();
    return bytes;
}
SolverBlockMapping CheckedWidthMapping(const PreconditionerPartition & p,const PreconditionerLimits & limits)
{
    std::size_t coordinates{};
    for(const auto & block:p.blocks)
    {
        const auto size=block.core_atoms.size()+block.overlap_atoms.size();
        if(size>limits.block_atoms) throw std::runtime_error("preconditioner-block-limit");
        coordinates+=size;
    }
    Budget(TopologyBytes(p)+64*coordinates,limits.storage_bytes);
    Budget(32*p.problem->atom_ids.size()+32*coordinates,limits.scratch_bytes);
    return WidthMapping(p);
}
std::vector<Indices> Coordinates(const PreconditionerPartition & p,const Indices & positions)
{
    std::vector<Indices> out;
    for(const auto & block:p.blocks)
    {
        Indices local;
        for(const auto * atoms:{&block.core_atoms,&block.overlap_atoms}) for(auto a:*atoms) local.push_back(positions.at(static_cast<std::size_t>(a)));
        out.push_back(std::move(local));
    }
    return out;
}
double Infinity(const Matrix & x) {return x.cwiseAbs().rowwise().sum().maxCoeff();}
void RecordRegularization(RegularizationRecord record,std::size_t held,const PreconditionerLimits & limits)
{
    if(!ResourceWorkForTesting().enabled) return;
    auto & records=SearchWorkForTesting().regularizations;
    Budget(held+2*(records.size()+1)*sizeof(RegularizationRecord),limits.storage_bytes);
    records.push_back(record);
}
}
std::shared_ptr<const PreconditionerPartition> BuildPreconditionerPartition(
    std::shared_ptr<const JointProblemInput> input,const JointParameterLayout & layout,const PreconditionerLimits & limits)
{
    ResourcePhase phase("preconditioner-partition"); auto & work=SearchWorkForTesting(); WorkTimer timer(work.partition_seconds);
    if(!input || limits.core_atoms==0 || limits.core_atoms>limits.block_atoms) throw std::invalid_argument("Invalid partition policy");
    const auto n=input->observations.size(),m=input->atom_ids.size();
    std::size_t memberships{};
    for(auto a:layout.full_atoms) memberships+=input->support.at(a).size();
    const auto scratch=24*(n+1)+24*(m+1)+16*memberships;
    Budget(scratch,limits.scratch_bytes); work.scratch_bytes=std::max(work.scratch_bytes,scratch);
    std::vector<bool> informative(n,false),assigned(m,false);
    for(auto r:layout.informative_rows) informative.at(r)=true;
    std::vector<std::size_t> offsets(n+1),order=layout.full_atoms;
    std::sort(order.begin(),order.end(),[&](auto a,auto b){return input->atom_ids.at(a)<input->atom_ids.at(b);});
    auto less=[&](Eigen::Index a,Eigen::Index b){return input->atom_ids.at(static_cast<std::size_t>(a))<input->atom_ids.at(static_cast<std::size_t>(b));};
    for(auto a:order) for(const auto & s:input->support.at(a)) if(informative.at(s.row)) ++offsets[s.row+1];
    std::partial_sum(offsets.begin(),offsets.end(),offsets.begin());
    Indices incidence(offsets.back()); auto cursor=offsets;
    for(auto a:order) for(const auto & s:input->support[a]) if(informative[s.row]) incidence[cursor[s.row]++]=static_cast<Eigen::Index>(a);
    std::vector<std::size_t> marks(m); std::size_t generation{};
    auto neighbors=[&](std::size_t a) {
        Indices out; ++generation;
        for(const auto & s:input->support[a]) if(informative[s.row])
            for(auto k=offsets[s.row];k<offsets[s.row+1];++k)
            {const auto b=static_cast<std::size_t>(incidence[k]); if(marks[b]!=generation) {marks[b]=generation; out.push_back(static_cast<Eigen::Index>(b));}}
        std::sort(out.begin(),out.end(),less); out.erase(std::unique(out.begin(),out.end()),out.end()); return out;
    };
    std::vector<PreconditionerBlock> blocks;
    std::size_t bytes=32*m+8*(layout.full_atoms.size()+layout.informative_rows.size());
    for(auto seed:order) if(!assigned[seed])
    {
        PreconditionerBlock block; block.id=input->atom_ids[seed];
        std::deque<Eigen::Index> queue{static_cast<Eigen::Index>(seed)}; std::vector<bool> queued(m,false); queued[seed]=true;
        while(!queue.empty() && block.core_atoms.size()<limits.core_atoms)
        {
            const auto a=queue.front(); queue.pop_front(); if(assigned[static_cast<std::size_t>(a)]) continue;
            assigned[static_cast<std::size_t>(a)]=true; block.core_atoms.push_back(a);
            for(auto b:neighbors(static_cast<std::size_t>(a))) if(!assigned[static_cast<std::size_t>(b)] && !queued[static_cast<std::size_t>(b)])
            {queue.push_back(b); queued[static_cast<std::size_t>(b)]=true;}
        }
        std::sort(block.core_atoms.begin(),block.core_atoms.end(),less);
        std::vector<bool> included(m,false); for(auto a:block.core_atoms) included[static_cast<std::size_t>(a)]=true;
        for(auto a:block.core_atoms) for(auto b:neighbors(static_cast<std::size_t>(a))) if(!included[static_cast<std::size_t>(b)])
        {
            included[static_cast<std::size_t>(b)]=true; block.overlap_atoms.push_back(b);
            if(block.core_atoms.size()+block.overlap_atoms.size()>limits.block_atoms) throw std::runtime_error("preconditioner-block-limit");
        }
        std::sort(block.overlap_atoms.begin(),block.overlap_atoms.end(),less);
        std::size_t row_bound{};
        for(const auto * atoms:{&block.core_atoms,&block.overlap_atoms}) for(auto a:*atoms) row_bound+=input->support[static_cast<std::size_t>(a)].size();
        Budget(bytes+8*row_bound+64*limits.block_atoms+sizeof(block)+block.id.size(),limits.storage_bytes);
        block.informative_rows.reserve(row_bound);
        for(const auto * atoms:{&block.core_atoms,&block.overlap_atoms}) for(auto a:*atoms)
            for(const auto & s:input->support[static_cast<std::size_t>(a)]) if(informative[s.row]) block.informative_rows.push_back(static_cast<Eigen::Index>(s.row));
        std::sort(block.informative_rows.begin(),block.informative_rows.end());
        block.informative_rows.erase(std::unique(block.informative_rows.begin(),block.informative_rows.end()),block.informative_rows.end());
        block.informative_rows.shrink_to_fit();
        bytes+=sizeof(block)+block.id.size()+32*(block.core_atoms.size()+block.overlap_atoms.size())+8*block.informative_rows.size();
        Budget(bytes,limits.storage_bytes);
        work.maximum_block_atoms=std::max(work.maximum_block_atoms,block.core_atoms.size()+block.overlap_atoms.size());
        blocks.push_back(std::move(block));
    }
    auto result=std::make_shared<const PreconditionerPartition>(std::move(input),layout,std::move(blocks));
    work.topology_bytes=TopologyBytes(*result); Budget(work.topology_bytes,limits.storage_bytes); return result;
}
std::shared_ptr<const PreconditionerPartition> SearchPartition(const Domain & domain,const EvaluationContext & context,const PreconditionerLimits & limits)
{
    auto snapshot=domain.atoms.Snapshot(); JointParameterLayout layout;
    // Hand-built numerical fixtures predate immutable problem identities.
    if(snapshot->atom_ids.size()!=domain.atoms.size() && snapshot->atom_ids.empty())
    {
        auto input=std::make_shared<JointProblemInput>(); input->atom_ids=static_cast<std::vector<std::string>>(context.atom_ids);
        input->observations.resize(static_cast<std::size_t>(domain.rows)); input->support.resize(domain.atoms.size());
        for(std::size_t a=0;a<domain.atoms.size();++a) for(const auto & s:domain.atoms[a]) input->support[a].push_back({static_cast<std::size_t>(s.row),s.square});
        snapshot=input; layout.full_atoms.resize(domain.atoms.size()); std::iota(layout.full_atoms.begin(),layout.full_atoms.end(),0);
        layout.informative_rows.resize(static_cast<std::size_t>(domain.rows)); std::iota(layout.informative_rows.begin(),layout.informative_rows.end(),0);
    }
    else
    {
        std::vector<bool> rows(snapshot->observations.size(),false);
        for(std::size_t a=0;a<domain.atoms.size();++a)
        {
            layout.full_atoms.push_back(domain.atoms.ParentAtom(a)); const auto support=domain.atoms[a];
            for(std::size_t k=0;k<support.size();++k) rows.at(static_cast<std::size_t>(support.ParentRow(k)))=true;
        }
        for(std::size_t r=0;r<rows.size();++r) if(rows[r]) layout.informative_rows.push_back(r);
    }
    return BuildPreconditionerPartition(std::move(snapshot),layout,limits);
}
SolverBlockMapping WidthMapping(const PreconditionerPartition & p)
{
    Indices positions(p.problem->atom_ids.size(),-1);
    for(std::size_t a=0;a<p.layout.full_atoms.size();++a) positions[p.layout.full_atoms[a]]=static_cast<Eigen::Index>(a);
    return {PreconditionerSpace::Width,static_cast<Eigen::Index>(p.layout.full_atoms.size()),Coordinates(p,positions)};
}
SolverBlockMapping FreeColumnMapping(const PreconditionerPartition & p,VectorRef beta)
{
    const auto widths=WidthMapping(p); if(beta.size()!=2*widths.dimension) throw std::invalid_argument("Invalid free face");
    Indices positions(static_cast<std::size_t>(beta.size()),-1); Eigen::Index count{};
    for(Eigen::Index k=0;k<beta.size();++k) if(k%2 || beta(k)>0) positions[static_cast<std::size_t>(k)]=count++;
    std::vector<Indices> coordinates;
    for(const auto & block:widths.blocks)
    {
        Indices columns; for(auto a:block.global) for(auto k:{2*a,2*a+1}) if(positions[static_cast<std::size_t>(k)]>=0) columns.push_back(positions[static_cast<std::size_t>(k)]);
        coordinates.push_back(std::move(columns));
    }
    return {PreconditionerSpace::FreeAC,count,coordinates};
}
Sparse RawWidthDerivative(const Evaluation & e)
{
    Sparse raw(e.derivative.rows(),e.eta.size()); raw.reserve(e.derivative.nonZeros());
    for(Eigen::Index a=0;a<e.eta.size();++a)
    {
        raw.startVec(a);
        // The two kernels have the same frozen sparse support, but explicit
        // sparse addition also handles test inputs with differing patterns.
        const Eigen::SparseVector<double> column=e.beta(2*a)*e.derivative.col(2*a)+e.beta(2*a+1)*e.derivative.col(2*a+1);
        for(Eigen::SparseVector<double>::InnerIterator it(column);it;++it) raw.insertBack(it.index(),a)=it.value();
    }
    raw.finalize(); return raw;
}
Vector WidthNorms(const Sparse & raw,double scale)
{
    if(!(scale>0) || !std::isfinite(scale)) throw std::invalid_argument("Invalid observation scale");
    Vector norms(raw.cols()); for(Eigen::Index a=0;a<raw.cols();++a) norms(a)=raw.col(a).norm()/scale;
    if(!norms.allFinite()) throw std::runtime_error("nonfinite-width-metric"); return norms;
}
Vector WidthMetric(VectorRef norms)
{
    if(norms.size()==0 || !norms.allFinite() || (norms.array()<0).any()) throw std::invalid_argument("Invalid width norms");
    const double maximum=norms.maxCoeff();
    return maximum==0 ? Vector::Ones(norms.size()).eval() : norms.array().max(std::sqrt(eps)*maximum).matrix().eval();
}
SchwarzModel::SchwarzModel(const PreconditionerPartition & partition,const Evaluation & e,double scale,
    const PreconditionerContext & context,const PreconditionerLimits & limits)
    :context_(context),mapping_(CheckedWidthMapping(partition,limits)),limits_(limits),bytes_(TopologyBytes(partition))
{
    ResourcePhase phase("preconditioner-local"); auto & work=SearchWorkForTesting(); WorkTimer timer(work.local_seconds); ++work.local_builds;
    build_=work.local_builds;
    if(!context.Valid() || context.space!=PreconditionerSpace::Width || context.metric.size()!=mapping_.dimension || e.eta.size()!=mapping_.dimension || !e.valid)
        throw std::invalid_argument("Invalid Schwarz linearization");
    Budget(bytes_,limits.storage_bytes);
    for(const auto & block:mapping_.blocks)
    {
        const auto b=block.global.size(); if(b>limits.block_atoms) throw std::runtime_error("preconditioner-block-limit");
        // Model plus a shifted LLT, mapping and metric. Sparse product temporaries
        // are bounded below before constructing any selected matrices.
        bytes_+=16*b*b+64*b; Budget(bytes_,limits.storage_bytes);
        std::size_t nnz{}; for(auto a:block.global) nnz+=static_cast<std::size_t>(e.x.col(2*a).nonZeros()+e.x.col(2*a+1).nonZeros()+e.derivative.col(2*a).nonZeros()+e.derivative.col(2*a+1).nonZeros());
        const std::size_t scratch=160*b*b+48*nnz+32*static_cast<std::size_t>(e.x.rows()+1);
        Budget(scratch,limits.scratch_bytes); work.scratch_bytes=std::max(work.scratch_bytes,scratch);
        Indices free; for(auto a:block.global) for(auto k:{2*a,2*a+1}) if(k%2 || e.beta(k)>0) free.push_back(k);
        Sparse z(e.x.rows(),static_cast<Eigen::Index>(free.size())),d(e.x.rows(),static_cast<Eigen::Index>(b));
        z.reserve(static_cast<Eigen::Index>(nnz)); d.reserve(static_cast<Eigen::Index>(nnz));
        for(std::size_t j=0;j<free.size();++j)
        {
            z.startVec(static_cast<Eigen::Index>(j)); const auto k=free[j]; const double norm=e.x.col(k).norm();
            if(!std::isfinite(norm)) throw std::runtime_error("preconditioner-nonfinite");
            if(norm>0) for(Sparse::InnerIterator it(e.x,k);it;++it) z.insertBack(it.row(),static_cast<Eigen::Index>(j))=it.value()/norm;
        }
        z.finalize();
        for(std::size_t j=0;j<b;++j)
        {
            d.startVec(static_cast<Eigen::Index>(j)); const auto a=block.global[j];
            const Eigen::SparseVector<double> column=(e.beta(2*a)*e.derivative.col(2*a)+e.beta(2*a+1)*e.derivative.col(2*a+1))/(scale*context.metric(a));
            for(Eigen::SparseVector<double>::InnerIterator it(column);it;++it) d.insertBack(it.index(),static_cast<Eigen::Index>(j))=it.value();
        }
        d.finalize();
        Matrix a=Matrix(z.transpose()*z),cross=Matrix(z.transpose()*d),c=Matrix(d.transpose()*d);
        RecordDenseShape("schwarz-ac-gram",a.rows(),a.cols()); RecordDenseShape("schwarz-cross",cross.rows(),cross.cols()); RecordDenseShape("schwarz-width-gram",c.rows(),c.cols());
        const double lambda=std::sqrt(eps)*std::max(1.,Infinity(a)); work.maximum_lambda=std::max(work.maximum_lambda,lambda);
        lambdas_.push_back(lambda);
        RecordRegularization({build_,0,schur_.size(),lambda,context.damping,0,1},bytes_,limits_);
        a.diagonal().array()+=lambda; Eigen::LLT<Matrix> factor(a);
        if(factor.info()!=Eigen::Success) throw std::runtime_error("preconditioner-ac-factorization");
        Matrix schur=c-cross.transpose()*factor.solve(cross); schur=(.5*(schur+schur.transpose())).eval();
        if(!schur.allFinite()) throw std::runtime_error("preconditioner-nonfinite");
        schur_.push_back(std::move(schur));
    }
    work.storage_bytes=std::max(work.storage_bytes,bytes_);
}
SchwarzPreconditioner::SchwarzPreconditioner(const SchwarzModel & model,const PreconditionerContext & context):model_(model),context_(context)
{
    ResourcePhase phase("preconditioner-factor"); auto & work=SearchWorkForTesting(); WorkTimer timer(work.factor_seconds); ++work.factor_builds;
    auto expected=context; expected.damping=model.Context().damping;
    if(!context.Valid() || !model.Context().Matches(expected)) throw std::logic_error("Stale preconditioner context");
    for(const auto & s:model.Matrices())
    {
        double tau=64*eps*std::max({1.,Infinity(s),context.damping}); Eigen::LLT<Matrix> factor;
        for(int retry=0;retry<6;++retry)
        {
            work.maximum_tau=std::max(work.maximum_tau,tau);
            RecordRegularization({model.Build(),work.factor_builds,factors_.size(),model.Lambdas()[factors_.size()],context.damping,tau,retry+1},model.Bytes(),model.Limits());
            Matrix shifted=s; shifted.diagonal().array()+=context.damping+tau;
            RecordDenseShape("schwarz-shifted-factor",s.rows(),s.cols()); factor.compute(shifted);
            if(factor.info()==Eigen::Success) break;
            if(retry==5) throw std::runtime_error("preconditioner-width-factorization");
            tau*=10;
        }
        work.maximum_tau=std::max(work.maximum_tau,tau); factors_.push_back(std::move(factor));
    }
}
Vector SchwarzPreconditioner::ApplyInverse(VectorRef rhs,const PreconditionerContext & context) const
{
    ResourcePhase phase("preconditioner-inverse"); auto & work=SearchWorkForTesting(); WorkTimer timer(work.inverse_seconds); ++work.inverse_actions;
    if(!context_.Matches(context)) throw std::logic_error("Stale preconditioner context");
    if(rhs.size()!=context_.metric.size() || !rhs.allFinite()) throw std::invalid_argument("Invalid preconditioner RHS");
    Vector out=Vector::Zero(rhs.size());
    for(std::size_t k=0;k<factors_.size();++k)
    {
        const auto & block=model_.Mapping().blocks[k]; Vector local=block.Restrict(rhs);
        for(std::size_t j=0;j<block.global.size();++j) local(static_cast<Eigen::Index>(j))/=context_.metric(block.global[j]);
        local=factors_[k].solve(local).eval();
        for(std::size_t j=0;j<block.global.size();++j) local(static_cast<Eigen::Index>(j))/=context_.metric(block.global[j]);
        block.Scatter(local,out);
    }
    if(!out.allFinite()) throw std::runtime_error("preconditioner-nonfinite"); return out;
}
}
