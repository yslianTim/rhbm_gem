#include "support/JointReplayExperiment.hpp"
#include <algorithm>
#include <bit>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>

namespace second_stage_test {
namespace {
using namespace joint;

Eigen::Index AffectedPosition(const Indices & affected_rows,Eigen::Index row)
{
    const auto found=std::lower_bound(affected_rows.begin(),affected_rows.end(),row);
    return found==affected_rows.end() || *found!=row ? -1 :
        static_cast<Eigen::Index>(found-affected_rows.begin());
}

Indices OrderedAtoms(const rhbm_gem::JointParameterLayout & layout,const Indices & first,const Indices * second=nullptr)
{
    std::vector<std::pair<std::size_t,Eigen::Index>> ordered;
    ordered.reserve(first.size()+(second ? second->size() : 0));
    const auto append=[&](Eigen::Index atom) {
        const auto found=std::find(layout.full_atoms.begin(),layout.full_atoms.end(),static_cast<std::size_t>(atom));
        if(found==layout.full_atoms.end()) throw std::invalid_argument("Replay atom is absent from the component layout.");
        ordered.emplace_back(static_cast<std::size_t>(found-layout.full_atoms.begin()),atom);
    };
    for(const auto atom:first) append(atom);
    if(second) for(const auto atom:*second) append(atom);
    std::sort(ordered.begin(),ordered.end());
    Indices result; result.reserve(ordered.size());
    for(const auto & [position,atom]:ordered)
    {
        (void)position;
        if(result.empty() || result.back()!=atom) result.push_back(atom);
    }
    return result;
}

std::size_t VectorBytes(const Vector & vector)
{return static_cast<std::size_t>(vector.size())*sizeof(double);}

void EvaluateObjective(const rhbm_gem::JointParameterLayout & layout,VectorRef observations,
    const BlockCoordinateState & accepted,const Indices & affected_rows,double scale,
    ReplayObjectiveStrategy strategy,ReplayCandidateResult & out)
{
    if(strategy==ReplayObjectiveStrategy::Incremental ||
        strategy==ReplayObjectiveStrategy::IncrementalScaled ||
        strategy==ReplayObjectiveStrategy::IncrementalKahan ||
        strategy==ReplayObjectiveStrategy::IncrementalNeumaier ||
        strategy==ReplayObjectiveStrategy::IncrementalKahanScaled ||
        strategy==ReplayObjectiveStrategy::IncrementalNeumaierScaled)
    {
        double sum{},correction{};
        const bool scaled=strategy==ReplayObjectiveStrategy::IncrementalScaled ||
            strategy==ReplayObjectiveStrategy::IncrementalKahanScaled ||
            strategy==ReplayObjectiveStrategy::IncrementalNeumaierScaled;
        const bool kahan=strategy==ReplayObjectiveStrategy::IncrementalKahan ||
            strategy==ReplayObjectiveStrategy::IncrementalKahanScaled;
        const bool neumaier=strategy==ReplayObjectiveStrategy::IncrementalNeumaier ||
            strategy==ReplayObjectiveStrategy::IncrementalNeumaierScaled;
        const double scale_squared=scale*scale;
        for(std::size_t k=0;k<affected_rows.size();++k)
        {
            const auto row=affected_rows[k];
            const double old_residual=accepted.residual(row),new_residual=out.residual(static_cast<Eigen::Index>(k));
            double term=new_residual*new_residual-old_residual*old_residual;
            if(scaled) term/=scale_squared;
            if(kahan)
            {
                const double adjusted=term-correction;
                const double next=sum+adjusted;
                correction=(next-sum)-adjusted;
                sum=next;
            }
            else if(neumaier)
            {
                const double next=sum+term;
                correction+=std::abs(sum)>=std::abs(term) ? (sum-next)+term : (term-next)+sum;
                sum=next;
            }
            else sum+=term;
        }
        sum+=neumaier ? correction : 0.0;
        out.objective=scaled ? accepted.objective+.5*sum :
            accepted.objective+.5*sum/scale_squared;
        return;
    }
    if(strategy==ReplayObjectiveStrategy::CachedExternally)
    {out.objective=accepted.objective; return;}
    double squared{};
    for(const auto row:layout.informative_rows)
    {
        const auto local=AffectedPosition(affected_rows,static_cast<Eigen::Index>(row));
        const double residual=local>=0 ? out.residual(local) : accepted.residual(static_cast<Eigen::Index>(row));
        squared+=residual*residual;
    }
    out.objective=.5*squared/(scale*scale);
    (void)observations;
}

double SparseValue(const std::vector<std::pair<std::size_t,double>> & values,
    std::size_t index,double fallback)
{
    const auto found=std::lower_bound(values.begin(),values.end(),index,
        [](const auto & value,std::size_t target){return value.first<target;});
    return found!=values.end() && found->first==index ? found->second : fallback;
}
}

std::size_t ReplayIncidence::StorageBytes() const
{
    std::size_t bytes=atoms.capacity()*sizeof(Eigen::Index)+rows.capacity()*sizeof(std::vector<ReplayContributor>);
    for(const auto & row:rows) bytes+=row.capacity()*sizeof(ReplayContributor);
    return bytes;
}

ReplayIncidence PrepareReplayIncidence(const joint::Domain & domain,const rhbm_gem::JointParameterLayout & layout,
    const joint::Indices & core_atoms,const joint::Indices & neighbor_atoms,const joint::Indices & affected_rows)
{
    ReplayIncidence out; out.atoms=OrderedAtoms(layout,core_atoms,&neighbor_atoms);
    out.rows.resize(affected_rows.size());
    for(std::size_t atom_position=0;atom_position<out.atoms.size();++atom_position)
    {
        const auto atom=static_cast<std::size_t>(out.atoms[atom_position]);
        for(const auto & support:domain.atoms[atom])
        {
            const auto row_position=AffectedPosition(affected_rows,support.row);
            if(row_position>=0)
                out.rows[static_cast<std::size_t>(row_position)].push_back({atom_position,support.square});
        }
    }
    return out;
}

bool VerifyReplayIncidenceCoverage(const joint::Domain & domain,const rhbm_gem::JointParameterLayout & layout,
    const joint::Indices & affected_rows,const ReplayIncidence & incidence)
{
    if(incidence.rows.size()!=affected_rows.size()) return false;
    for(std::size_t row_position=0;row_position<affected_rows.size();++row_position)
    {
        std::vector<ReplayContributor> expected;
        for(const auto atom_value:layout.full_atoms)
        {
            const auto atom=static_cast<std::size_t>(atom_value);
            for(const auto & support:domain.atoms[atom]) if(support.row==affected_rows[row_position])
                expected.push_back({static_cast<std::size_t>(std::find(incidence.atoms.begin(),incidence.atoms.end(),
                    static_cast<Eigen::Index>(atom))-incidence.atoms.begin()),support.square});
        }
        const auto & actual=incidence.rows[row_position];
        if(expected.size()!=actual.size()) return false;
        for(std::size_t k=0;k<expected.size();++k)
        {
            if(expected[k].atom_position>=incidence.atoms.size() || actual[k].atom_position>=incidence.atoms.size() ||
                incidence.atoms[expected[k].atom_position]!=incidence.atoms[actual[k].atom_position] ||
                std::bit_cast<std::uint64_t>(expected[k].square)!=std::bit_cast<std::uint64_t>(actual[k].square))
                return false;
        }
    }
    return true;
}

ReplayCandidateResult EvaluateReplayCandidate(const joint::Domain & domain,
    const rhbm_gem::JointParameterLayout & layout,joint::VectorRef observations,
    const joint::BlockCoordinateState & accepted,joint::VectorRef candidate_eta,
    joint::VectorRef candidate_beta,double scale,const joint::Indices & core_atoms,
    const joint::Indices & affected_rows,ReplayPredictionStrategy prediction_strategy,
    ReplayObjectiveStrategy objective_strategy,const ReplayIncidence * incidence)
{
    using namespace joint;
    ReplayCandidateResult out;
    out.prediction.resize(static_cast<Eigen::Index>(affected_rows.size()));
    out.residual.resize(static_cast<Eigen::Index>(affected_rows.size()));
    out.scratch_bytes=VectorBytes(out.prediction)+VectorBytes(out.residual);
    if(prediction_strategy==ReplayPredictionStrategy::AffectedRowExact)
    {
        if(!incidence || incidence->rows.size()!=affected_rows.size())
            throw std::invalid_argument("Affected-row replay needs complete local incidence.");
        Vector widths(static_cast<Eigen::Index>(incidence->atoms.size()));
        for(std::size_t k=0;k<incidence->atoms.size();++k)
            widths(static_cast<Eigen::Index>(k))=std::exp(candidate_eta(incidence->atoms[k]));
        out.scratch_bytes+=VectorBytes(widths)+incidence->StorageBytes();
        for(std::size_t row_position=0;row_position<affected_rows.size();++row_position)
        {
            double prediction{};
            for(const auto & contributor:incidence->rows[row_position])
            {
                const auto atom=incidence->atoms[contributor.atom_position];
                const auto basis=EvaluateKernel(contributor.square,widths(static_cast<Eigen::Index>(contributor.atom_position)),2.5);
                prediction+=candidate_beta(2*atom)*basis.gaussian+candidate_beta(2*atom+1)*basis.charge;
            }
            out.prediction(static_cast<Eigen::Index>(row_position))=prediction;
        }
    }
    else
    {
        const auto ordered_core=OrderedAtoms(layout,core_atoms);
        Vector old_widths(static_cast<Eigen::Index>(ordered_core.size()));
        Vector new_widths(static_cast<Eigen::Index>(ordered_core.size()));
        for(std::size_t k=0;k<ordered_core.size();++k)
        {
            old_widths(static_cast<Eigen::Index>(k))=std::exp(accepted.eta(ordered_core[k]));
            new_widths(static_cast<Eigen::Index>(k))=std::exp(candidate_eta(ordered_core[k]));
        }
        out.scratch_bytes+=VectorBytes(old_widths)+VectorBytes(new_widths);
        for(std::size_t k=0;k<affected_rows.size();++k)
            out.prediction(static_cast<Eigen::Index>(k))=accepted.prediction(affected_rows[k]);
        for(std::size_t atom_position=0;atom_position<ordered_core.size();++atom_position)
        {
            const auto atom=static_cast<std::size_t>(ordered_core[atom_position]);
            for(const auto & support:domain.atoms[atom])
            {
                const auto row_position=AffectedPosition(affected_rows,support.row);
                if(row_position<0) continue;
                const auto old_basis=EvaluateKernel(support.square,old_widths(static_cast<Eigen::Index>(atom_position)),2.5);
                const auto new_basis=EvaluateKernel(support.square,new_widths(static_cast<Eigen::Index>(atom_position)),2.5);
                const double old_contribution=accepted.beta(2*static_cast<Eigen::Index>(atom))*old_basis.gaussian+
                    accepted.beta(2*static_cast<Eigen::Index>(atom)+1)*old_basis.charge;
                const double new_contribution=candidate_beta(2*static_cast<Eigen::Index>(atom))*new_basis.gaussian+
                    candidate_beta(2*static_cast<Eigen::Index>(atom)+1)*new_basis.charge;
                auto & prediction=out.prediction(row_position);
                if(prediction_strategy==ReplayPredictionStrategy::DeltaAddDifference)
                    prediction+=new_contribution-old_contribution;
                else prediction=(prediction-old_contribution)+new_contribution;
            }
        }
    }
    for(std::size_t k=0;k<affected_rows.size();++k)
        out.residual(static_cast<Eigen::Index>(k))=out.prediction(static_cast<Eigen::Index>(k))-
            observations(affected_rows[k]);
    EvaluateObjective(layout,observations,accepted,affected_rows,scale,objective_strategy,out);
    return out;
}

TiledObjectiveCache::TiledObjectiveCache(const std::vector<std::size_t> & informative_rows,
    joint::VectorRef residuals,std::size_t tile_size,bool hierarchical)
    :informative_rows_(&informative_rows),tile_size_(tile_size),hierarchical_(hierarchical)
{
    if(tile_size_==0) throw std::invalid_argument("Objective tile size must be positive.");
    tile_count_=informative_rows.size()/tile_size_+
        (informative_rows.size()%tile_size_!=0 ? 1U : 0U);
    if(hierarchical_)
    {
        leaf_capacity_=1;
        while(leaf_capacity_<tile_count_)
        {
            if(leaf_capacity_>std::numeric_limits<std::size_t>::max()/2)
                throw std::overflow_error("Objective reduction tree size overflow.");
            leaf_capacity_*=2;
        }
        if(leaf_capacity_>std::numeric_limits<std::size_t>::max()/2)
            throw std::overflow_error("Objective reduction tree size overflow.");
        values_.assign(2*leaf_capacity_,0.0);
        for(std::size_t tile=0;tile<tile_count_;++tile)
        {
            const auto begin=tile*tile_size_;
            const auto end=begin+std::min(tile_size_,informative_rows.size()-begin);
            double squared{};
            for(std::size_t position=begin;position<end;++position)
            {
                const double residual=residuals(static_cast<Eigen::Index>(informative_rows[position]));
                squared+=residual*residual;
            }
            values_[leaf_capacity_+tile]=squared;
        }
        for(std::size_t node=leaf_capacity_;node-->1;)
            values_[node]=values_[2*node]+values_[2*node+1];
    }
    else
    {
        values_.assign(tile_count_,0.0);
        for(std::size_t tile=0;tile<tile_count_;++tile)
        {
            const auto begin=tile*tile_size_;
            const auto end=begin+std::min(tile_size_,informative_rows.size()-begin);
            for(std::size_t position=begin;position<end;++position)
            {
                const double residual=residuals(static_cast<Eigen::Index>(informative_rows[position]));
                values_[tile]+=residual*residual;
            }
        }
    }
}

ObjectiveCacheCandidate TiledObjectiveCache::Evaluate(joint::VectorRef accepted_residuals,
    const joint::Indices & affected_rows,joint::VectorRef candidate_residuals,double scale) const
{
    if(affected_rows.size()!=static_cast<std::size_t>(candidate_residuals.size()))
        throw std::invalid_argument("Objective candidate residual count mismatch.");
    const auto & informative_rows=*informative_rows_;
    std::vector<std::pair<std::size_t,double>> changed;
    changed.reserve(affected_rows.size());
    for(std::size_t k=0;k<affected_rows.size();++k)
    {
        const auto row=affected_rows[k];
        const auto found=std::lower_bound(informative_rows.begin(),informative_rows.end(),
            static_cast<std::size_t>(row));
        if(found!=informative_rows.end() && *found==static_cast<std::size_t>(row))
            changed.emplace_back(static_cast<std::size_t>(found-informative_rows.begin()),
                candidate_residuals(static_cast<Eigen::Index>(k)));
    }
    std::sort(changed.begin(),changed.end(),[](const auto & left,const auto & right) {
        return left.first<right.first;
    });
    std::vector<std::size_t> tiles;
    tiles.reserve(changed.size());
    for(const auto & [position,value]:changed)
    {
        (void)value;
        const auto tile=position/tile_size_;
        if(tiles.empty() || tiles.back()!=tile) tiles.push_back(tile);
    }
    std::vector<std::pair<std::size_t,double>> leaf_updates;
    leaf_updates.reserve(tiles.size());
    for(const auto tile:tiles)
    {
        const auto begin=tile*tile_size_;
        const auto end=begin+std::min(tile_size_,informative_rows.size()-begin);
        auto changed_it=std::lower_bound(changed.begin(),changed.end(),begin,
            [](const auto & value,std::size_t target){return value.first<target;});
        double squared{};
        for(std::size_t position=begin;position<end;++position)
        {
            double residual{};
            if(changed_it!=changed.end() && changed_it->first==position)
            {residual=changed_it->second; ++changed_it;}
            else residual=accepted_residuals(static_cast<Eigen::Index>(informative_rows[position]));
            squared+=residual*residual;
        }
        leaf_updates.emplace_back(tile,squared);
    }

    ObjectiveCacheCandidate out;
    out.tiles_recomputed=leaf_updates.size();
    if(!hierarchical_)
    {
        out.updates=leaf_updates;
        std::size_t update{};
        double total{};
        for(std::size_t tile=0;tile<tile_count_;++tile)
        {
            if(update<leaf_updates.size() && leaf_updates[update].first==tile)
                total+=leaf_updates[update++].second;
            else total+=values_[tile];
        }
        out.objective=.5*total/(scale*scale);
        return out;
    }

    std::vector<std::pair<std::size_t,double>> level;
    level.reserve(leaf_updates.size());
    for(const auto & [tile,value]:leaf_updates) level.emplace_back(leaf_capacity_+tile,value);
    out.updates=level;
    while(!level.empty() && level.front().first>1)
    {
        std::vector<std::pair<std::size_t,double>> parents;
        parents.reserve(level.size());
        for(const auto & [node,value]:level)
        {
            (void)value;
            const auto parent=node/2;
            if(parents.empty() || parents.back().first!=parent)
            {
                const auto left=2*parent,right=left+1;
                const double left_value=SparseValue(level,left,values_[left]);
                const double right_value=SparseValue(level,right,values_[right]);
                parents.emplace_back(parent,left_value+right_value);
            }
        }
        out.updates.insert(out.updates.end(),parents.begin(),parents.end());
        level=std::move(parents);
    }
    const double total=level.empty() ? values_[1] : level.front().second;
    out.objective=.5*total/(scale*scale);
    return out;
}

void TiledObjectiveCache::Commit(const ObjectiveCacheCandidate & candidate)
{
    for(const auto & [index,value]:candidate.updates) values_.at(index)=value;
}

std::size_t TiledObjectiveCache::StorageBytes() const
{return values_.capacity()*sizeof(double);}
}
