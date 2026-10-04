#include "OperatorSearch.hpp"
#include "StructuralPartition.hpp"
#include <algorithm>
#include <numeric>

namespace rhbm_gem::core::joint_component {
namespace {
constexpr double eps=std::numeric_limits<double>::epsilon();
void Budget(std::size_t bytes,std::size_t limit,const char * reason)
{if(bytes>limit) throw std::runtime_error(reason);}
std::size_t TopologyBytes(const PreconditionerPartition & p)
{
    std::size_t bytes=sizeof(p)+8*(p.layout.full_atoms.size()+p.layout.informative_rows.size())+32*p.atom_blocks.size();
    for(const auto & block:p.blocks) bytes+=sizeof(block)+block.id.size()+8*(block.core_atoms.size()+block.overlap_atoms.size()+block.informative_rows.size());
    for(const auto & membership:p.atom_blocks) bytes+=8*membership.size();
    return bytes;
}
SolverBlockMapping CheckedWidthMapping(const PreconditionerPartition & p,const SchwarzPolicy & policy)
{
    std::size_t coordinates{};
    for(const auto & block:p.blocks)
    {
        const auto size=block.core_atoms.size()+block.overlap_atoms.size();
        if(size>policy.max_block_atoms) throw std::runtime_error("preconditioner-block-limit");
        coordinates+=size;
    }
    Budget(TopologyBytes(p)+64*coordinates,policy.storage_bytes,"preconditioner-storage-limit");
    Budget(32*p.problem->atom_ids.size()+32*coordinates,policy.scratch_bytes,"preconditioner-scratch-limit");
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
void RecordRegularization(RegularizationRecord record,std::size_t held,const SchwarzPolicy & policy)
{
    if(!ResourceWorkForTesting().enabled) return;
    auto & records=SearchWorkForTesting().regularizations;
    Budget(held+2*(records.size()+1)*sizeof(RegularizationRecord),policy.storage_bytes,"preconditioner-storage-limit");
    records.push_back(record);
}
}
std::shared_ptr<const PreconditionerPartition> BuildPreconditionerPartition(
    std::shared_ptr<const JointProblemInput> input,const JointParameterLayout & layout,const SchwarzPolicy & policy)
{
    ResourcePhase phase("schwarz-partition",true,static_cast<Eigen::Index>(input ? input->observations.size() : 0),
        static_cast<Eigen::Index>(input ? input->atom_ids.size() : 0)); auto & work=SearchWorkForTesting(); WorkTimer timer(work.partition_seconds);
    if(!input || policy.core_atoms==0 || policy.core_atoms>policy.max_block_atoms) throw std::invalid_argument("Invalid partition policy");
    const auto n=input->observations.size(),m=input->atom_ids.size();
    std::size_t memberships{};
    for(auto a:layout.full_atoms) memberships+=input->support.at(a).size();
    const auto scratch=24*(n+1)+24*(m+1)+16*memberships+32*std::min(m,policy.max_block_atoms);
    Budget(scratch,policy.scratch_bytes,"preconditioner-scratch-limit"); work.scratch_bytes=std::max(work.scratch_bytes,scratch);
    StructuralTopology topology(*input,layout); auto structural=topology.BuildCores(policy.core_atoms,false);
    std::vector<bool> informative(n,false); for(auto row:layout.informative_rows) informative.at(row)=true;
    std::vector<PreconditionerBlock> blocks; blocks.reserve(structural.cores.size());
    std::size_t bytes=32*m+8*(layout.full_atoms.size()+layout.informative_rows.size());
    for(auto & core:structural.cores)
    {
        PreconditionerBlock block; block.id=std::move(core.id); block.core_atoms=std::move(core.atoms);
        try {block.overlap_atoms=topology.ContextAtoms(block.core_atoms,policy.overlap_hops,policy.max_block_atoms);}
        catch(const std::runtime_error & error)
        {
            if(std::string(error.what())=="structural-context-limit") throw std::runtime_error("preconditioner-block-limit");
            throw;
        }
        std::size_t row_bound{};
        for(const auto * atoms:{&block.core_atoms,&block.overlap_atoms}) for(auto a:*atoms) row_bound+=input->support[static_cast<std::size_t>(a)].size();
        Budget(bytes+8*row_bound+64*policy.max_block_atoms+sizeof(block)+block.id.size(),policy.storage_bytes,"preconditioner-storage-limit");
        block.informative_rows.reserve(row_bound);
        for(const auto * atoms:{&block.core_atoms,&block.overlap_atoms}) for(auto a:*atoms)
            for(const auto & s:input->support[static_cast<std::size_t>(a)])
                if(informative[s.row]) block.informative_rows.push_back(static_cast<Eigen::Index>(s.row));
        std::sort(block.informative_rows.begin(),block.informative_rows.end());
        block.informative_rows.erase(std::unique(block.informative_rows.begin(),block.informative_rows.end()),block.informative_rows.end());
        block.informative_rows.shrink_to_fit();
        bytes+=sizeof(block)+block.id.size()+32*(block.core_atoms.size()+block.overlap_atoms.size())+8*block.informative_rows.size();
        Budget(bytes,policy.storage_bytes,"preconditioner-storage-limit");
        work.maximum_block_atoms=std::max(work.maximum_block_atoms,block.core_atoms.size()+block.overlap_atoms.size());
        blocks.push_back(std::move(block));
    }
    auto result=std::make_shared<const PreconditionerPartition>(std::move(input),layout,std::move(blocks),policy);
    work.topology_bytes=TopologyBytes(*result); Budget(work.topology_bytes,policy.storage_bytes,"preconditioner-storage-limit"); return result;
}
std::shared_ptr<const PreconditionerPartition> SearchPartition(const Domain & domain,const EvaluationContext & context,const SchwarzPolicy & policy)
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
    return BuildPreconditionerPartition(std::move(snapshot),layout,policy);
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
    const PreconditionerContext & context)
    :context_(context),mapping_(CheckedWidthMapping(partition,partition.policy)),policy_(partition.policy),bytes_(TopologyBytes(partition))
{
    ResourcePhase phase("schwarz-local-build",true,e.x.rows(),e.beta.size()); auto & work=SearchWorkForTesting(); WorkTimer timer(work.local_seconds); ++work.local_builds;
    build_=work.local_builds;
    if(!context.Valid() || context.space!=PreconditionerSpace::Width || context.metric.size()!=mapping_.dimension || e.eta.size()!=mapping_.dimension || !e.valid)
        throw std::invalid_argument("Invalid Schwarz linearization");
    Budget(bytes_,policy_.storage_bytes,"preconditioner-storage-limit");
    for(const auto & block:mapping_.blocks)
    {
        const auto b=block.global.size(); if(b>policy_.max_block_atoms) throw std::runtime_error("preconditioner-block-limit");
        // Model plus a shifted LLT, mapping and metric. Sparse product temporaries
        // are bounded below before constructing any selected matrices.
        bytes_+=16*b*b+64*b; Budget(bytes_,policy_.storage_bytes,"preconditioner-storage-limit");
        std::size_t nnz{}; for(auto a:block.global) nnz+=static_cast<std::size_t>(e.x.col(2*a).nonZeros()+e.x.col(2*a+1).nonZeros()+e.derivative.col(2*a).nonZeros()+e.derivative.col(2*a+1).nonZeros());
        const std::size_t scratch=160*b*b+48*nnz+32*static_cast<std::size_t>(e.x.rows()+1);
        Budget(scratch,policy_.scratch_bytes,"preconditioner-scratch-limit"); work.scratch_bytes=std::max(work.scratch_bytes,scratch);
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
        RecordRegularization({build_,0,schur_.size(),lambda,context.damping,0,1},bytes_,policy_);
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
    ResourcePhase phase("schwarz-factor",true,static_cast<Eigen::Index>(model.Matrices().size()),
        static_cast<Eigen::Index>(model.Context().metric.size())); auto & work=SearchWorkForTesting(); WorkTimer timer(work.factor_seconds); ++work.factor_builds;
    auto expected=context; expected.damping=model.Context().damping;
    if(!context.Valid() || !model.Context().Matches(expected)) throw std::logic_error("Stale preconditioner context");
    for(const auto & s:model.Matrices())
    {
        double tau=64*eps*std::max({1.,Infinity(s),context.damping}); Eigen::LLT<Matrix> factor;
        for(int retry=0;retry<6;++retry)
        {
            work.maximum_tau=std::max(work.maximum_tau,tau);
            RecordRegularization({model.Build(),work.factor_builds,factors_.size(),model.Lambdas()[factors_.size()],context.damping,tau,retry+1},model.Bytes(),model.Policy());
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
