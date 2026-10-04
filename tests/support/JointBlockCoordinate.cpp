#include "support/JointBlockCoordinate.hpp"
#include <algorithm>
#include <cmath>

namespace second_stage_test::joint_block {
namespace {
n::Vector PredictRows(const rhbm_gem::core::JointProblemInput & input,const n::Indices & atoms,const n::Indices & rows,
    n::VectorRef eta,n::VectorRef beta)
{
    if(eta.size()!=static_cast<Eigen::Index>(input.atom_ids.size()) || beta.size()!=2*eta.size())
        throw std::invalid_argument("Invalid block state dimensions");
    std::vector<Eigen::Index> row_to_local(input.observations.size(),-1);
    for(std::size_t k=0;k<rows.size();++k) row_to_local.at(static_cast<std::size_t>(rows[k]))=static_cast<Eigen::Index>(k);
    n::Vector prediction=n::Vector::Zero(static_cast<Eigen::Index>(rows.size()));
    for(auto atom:atoms)
    {
        const auto a=static_cast<std::size_t>(atom); const double width=std::exp(eta(static_cast<Eigen::Index>(a)));
        for(const auto & support:input.support.at(a))
        {
            const auto row=row_to_local.at(support.row); if(row<0) continue;
            const auto basis=n::EvaluateKernel(support.squared_distance,width,2.5);
            prediction(row)+=beta(2*static_cast<Eigen::Index>(a))*basis.gaussian+
                beta(2*static_cast<Eigen::Index>(a)+1)*basis.charge;
        }
    }
    return prediction;
}
n::Vector PredictCoreRows(const rhbm_gem::core::JointProblemInput & input,const n::StructuralCore & core,
    const n::Indices & rows,n::VectorRef eta,n::VectorRef beta)
{
    if(eta.size()!=static_cast<Eigen::Index>(core.atoms.size()) || beta.size()!=2*eta.size())
        throw std::invalid_argument("Invalid local block dimensions");
    std::vector<Eigen::Index> row_to_local(input.observations.size(),-1);
    for(std::size_t k=0;k<rows.size();++k) row_to_local.at(static_cast<std::size_t>(rows[k]))=static_cast<Eigen::Index>(k);
    n::Vector prediction=n::Vector::Zero(static_cast<Eigen::Index>(rows.size()));
    for(std::size_t k=0;k<core.atoms.size();++k)
    {
        const auto atom=static_cast<std::size_t>(core.atoms[k]); const auto local=static_cast<Eigen::Index>(k);
        const double width=std::exp(eta(local));
        for(const auto & support:input.support.at(atom))
        {
            const auto row=row_to_local.at(support.row); if(row<0) continue;
            const auto basis=n::EvaluateKernel(support.squared_distance,width,2.5);
            prediction(row)+=beta(2*local)*basis.gaussian+beta(2*local+1)*basis.charge;
        }
    }
    return prediction;
}
double Objective(n::VectorRef residual,double scale)
{return .5*residual.squaredNorm()/(scale*scale);}
}
double ParentScale(n::VectorRef observations) {return std::max(1.0,observations.norm());}
BlockCoordinateState Replay(const rhbm_gem::core::JointProblemInput & input,const rhbm_gem::JointParameterLayout & layout,
    n::VectorRef observations,n::VectorRef eta,n::VectorRef beta,double parent_scale)
{
    if(observations.size()!=static_cast<Eigen::Index>(input.observations.size()) ||
        !(parent_scale>0) || !std::isfinite(parent_scale)) throw std::invalid_argument("Invalid parent block state");
    BlockCoordinateState state; state.observations=observations; state.eta=eta; state.beta=beta; state.scale=parent_scale;
    state.prediction=n::Vector::Zero(observations.size()); state.residual=n::Vector::Zero(observations.size());
    n::Indices atoms,rows; atoms.reserve(layout.full_atoms.size()); rows.reserve(layout.informative_rows.size());
    for(auto atom:layout.full_atoms) atoms.push_back(static_cast<Eigen::Index>(atom));
    for(auto row:layout.informative_rows) rows.push_back(static_cast<Eigen::Index>(row));
    const auto informative=PredictRows(input,atoms,rows,state.eta,state.beta);
    for(std::size_t k=0;k<layout.informative_rows.size();++k)
    {
        const auto row=static_cast<Eigen::Index>(layout.informative_rows[k]);
        state.prediction(row)=informative(static_cast<Eigen::Index>(k)); state.residual(row)=state.prediction(row)-state.observations(row);
    }
    state.objective=Objective(state.residual, state.scale);
    return state;
}
BlockCoordinateState Initialize(const rhbm_gem::core::JointProblemInput & input,const rhbm_gem::JointParameterLayout & layout,
    n::VectorRef observations,n::Vector eta,n::Vector beta,double parent_scale)
{return Replay(input,layout,observations,eta,beta,parent_scale);}
ConditionalBlockProblem BuildConditionalProblem(const rhbm_gem::core::JointProblemInput & input,
    const n::StructuralCore & core,const BlockCoordinateState & state)
{
    ConditionalBlockProblem problem; problem.rows=core.affected_rows; problem.scale=state.scale;
    n::Vector eta(static_cast<Eigen::Index>(core.atoms.size())),beta(2*static_cast<Eigen::Index>(core.atoms.size()));
    for(std::size_t k=0;k<core.atoms.size();++k)
    {
        const auto atom=static_cast<Eigen::Index>(core.atoms[k]); eta(static_cast<Eigen::Index>(k))=state.eta(atom);
        beta.segment<2>(2*static_cast<Eigen::Index>(k))=state.beta.segment<2>(2*atom);
    }
    problem.old_core_prediction=PredictCoreRows(input,core,problem.rows,eta,beta);
    problem.effective_response.resize(static_cast<Eigen::Index>(problem.rows.size()));
    for(std::size_t k=0;k<problem.rows.size();++k)
    {
        const auto row=problem.rows[k];
        problem.effective_response(static_cast<Eigen::Index>(k))=
            problem.old_core_prediction(static_cast<Eigen::Index>(k))-state.residual(row);
    }
    return problem;
}
double ConditionalObjective(const rhbm_gem::core::JointProblemInput & input,const n::StructuralCore & core,
    const ConditionalBlockProblem & problem,n::VectorRef eta,n::VectorRef beta)
{
    const auto prediction=PredictCoreRows(input,core,problem.rows,eta,beta);
    return Objective(prediction-problem.effective_response,problem.scale);
}
void ReplaceBlock(const rhbm_gem::core::JointProblemInput & input,const n::StructuralCore & core,
    n::VectorRef eta,n::VectorRef beta,BlockCoordinateState & state)
{
    if(eta.size()!=static_cast<Eigen::Index>(core.atoms.size()) || beta.size()!=2*eta.size() ||
        !eta.allFinite() || !beta.allFinite() || !((eta.array().exp()>0).all()) || !((eta.array().exp()).isFinite().all()))
        throw std::invalid_argument("Invalid conditional block state");
    for(Eigen::Index k=0;k<beta.size();k+=2) if(beta(k)<0) throw std::invalid_argument("Invalid conditional block state");
    n::Vector old_eta(eta.size()),old_beta(beta.size());
    for(std::size_t k=0;k<core.atoms.size();++k)
    {
        const auto atom=static_cast<Eigen::Index>(core.atoms[k]); old_eta(static_cast<Eigen::Index>(k))=state.eta(atom);
        old_beta.segment<2>(2*static_cast<Eigen::Index>(k))=state.beta.segment<2>(2*atom);
    }
    const auto old_core=PredictCoreRows(input,core,core.affected_rows,old_eta,old_beta);
    const auto next_core=PredictCoreRows(input,core,core.affected_rows,eta,beta);
    double old_squared{},new_squared{};
    for(std::size_t k=0;k<core.affected_rows.size();++k)
    {
        const auto row=core.affected_rows[k]; const double old=state.residual(row);
        const double change=next_core(static_cast<Eigen::Index>(k))-old_core(static_cast<Eigen::Index>(k));
        state.prediction(row)+=change; state.residual(row)+=change;
        old_squared+=old*old; new_squared+=std::pow(state.residual(row),2);
    }
    state.objective+=(new_squared-old_squared)/(2*state.scale*state.scale);
    for(std::size_t k=0;k<core.atoms.size();++k)
    {
        const auto atom=static_cast<std::size_t>(core.atoms[k]);
        state.eta(static_cast<Eigen::Index>(atom))=eta(static_cast<Eigen::Index>(k));
        state.beta.segment<2>(2*static_cast<Eigen::Index>(atom))=beta.segment<2>(2*static_cast<Eigen::Index>(k));
    }
}
}
