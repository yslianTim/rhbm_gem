#pragma once
#include "Preconditioner.hpp"
#include <functional>

namespace rhbm_gem::core::joint_component {
struct PreconditionerLimits
{
    std::size_t core_atoms{128},block_atoms{512};
    std::size_t storage_bytes{512ULL*1024*1024},scratch_bytes{256ULL*1024*1024};
};
struct RegularizationRecord
{
    std::size_t local_build{},factor_build{},block{};
    double lambda{},damping{},tau{};
    int attempt{};
};
struct SearchWork
{
    std::size_t linearizations{},pcg_solves{},pcg_iterations{},damping_trials{},local_builds{},factor_builds{},inverse_actions{};
    std::size_t topology_bytes{},storage_bytes{},scratch_bytes{},maximum_block_atoms{};
    double partition_seconds{},metric_seconds{},local_seconds{},factor_seconds{},inverse_seconds{},pcg_seconds{};
    double maximum_lambda{},maximum_tau{},last_relative_residual{};
    std::vector<RegularizationRecord> regularizations; // Resource-audit mode only.
};
SearchWork & SearchWorkForTesting();
std::shared_ptr<const PreconditionerPartition> BuildPreconditionerPartition(
    std::shared_ptr<const JointProblemInput>,const JointParameterLayout &,const PreconditionerLimits & = {});
std::shared_ptr<const PreconditionerPartition> SearchPartition(const Domain &,const EvaluationContext &,const PreconditionerLimits & = {});
SolverBlockMapping WidthMapping(const PreconditionerPartition &);
SolverBlockMapping FreeColumnMapping(const PreconditionerPartition &,VectorRef beta);
Sparse RawWidthDerivative(const Evaluation &);
Vector WidthNorms(const Sparse &,double);
Vector WidthMetric(VectorRef);
class SchwarzModel
{
    PreconditionerContext context_;
    SolverBlockMapping mapping_;
    std::vector<Matrix> schur_;
    std::vector<double> lambdas_;
    std::size_t build_{};
    PreconditionerLimits limits_;
    std::size_t bytes_{};
public:
    SchwarzModel(const PreconditionerPartition &,const Evaluation &,double,const PreconditionerContext &,const PreconditionerLimits & = {});
    const SolverBlockMapping & Mapping() const {return mapping_;}
    const std::vector<Matrix> & Matrices() const {return schur_;}
    const std::vector<double> & Lambdas() const {return lambdas_;}
    std::size_t Build() const {return build_;}
    const PreconditionerContext & Context() const {return context_;}
    const PreconditionerLimits & Limits() const {return limits_;}
    std::size_t Bytes() const {return bytes_;}
};
class SchwarzPreconditioner
{
    const SchwarzModel & model_;
    PreconditionerContext context_;
    std::vector<Eigen::LLT<Matrix>> factors_;
public:
    SchwarzPreconditioner(const SchwarzModel &,const PreconditionerContext &);
    Vector ApplyInverse(VectorRef,const PreconditionerContext &) const;
};
struct WidthStepResult
{
    Vector step;
    bool valid{};
    std::string reason;
    int iterations{};
    double relative_residual{unavailable},predicted{unavailable};
};
using VectorAction=std::function<Vector(VectorRef)>;
// Generic SPD kernel also permits an independent dense-system test oracle.
WidthStepResult SolvePcg(const VectorAction &,const VectorAction &,VectorRef,VectorRef,int=-1);
WidthStepResult WidthStepSolver(const ProfileJacobianOperator &,VectorRef,const PreconditionerContext &,const VectorAction &,int=-1);
SearchResult SearchOperatorProfile(const Domain &,VectorRef,const Vector &,const EvaluationContext &);
}
