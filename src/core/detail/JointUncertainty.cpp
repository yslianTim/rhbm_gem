#include "JointUncertainty.hpp"
#include "JointPostprocessing.hpp"
#include "joint_component/Problem.hpp"
#include "joint_component/Numerics.hpp"
#include "joint_component/TiledQR.hpp"
#include "joint_component/CompactSvd.hpp"
#include "joint_component/TargetEvidence.hpp"
#include <rhbm_gem/utils/hrl/LinearizationService.hpp>
#include <Eigen/SVD>
#include <cmath>
#include <limits>

namespace rhbm_gem::core::detail {
namespace {
StageUncertainty TargetUncertainty(const JointProblemInput & input,const JointAnalysisComponent & component,
    const joint_component::ComponentView & view,const std::vector<bool> & requested,
    std::map<std::size_t,Eigen::Matrix3d> & blocks,std::size_t original_rows,double observation_scale)
{
    StageUncertainty out; out.method="iid-ls-target-quotient"; out.status=EvidenceStatus::Unavailable;
    if(!component.state) {out.reason="missing-component-state"; return out;}
    const auto & state=*component.state;
    for(std::size_t a=0;a<state.b.size();++a) if(state.ac.at(2*a)==0)
    {out.reason="active-amplitude-boundary"; return out;}
    if(component.target_runtime_convergence!=JointCheckStatus::Passed)
    {out.reason="targets-not-converged"; return out;}
    joint_component::Vector y(component.rows.size());
    std::vector<std::vector<joint_component::Support>> support(component.atoms.size());
    for(std::size_t r=0;r<component.rows.size();++r) y(static_cast<Eigen::Index>(r))=input.observations.at(component.rows[r]);
    for(std::size_t a=0;a<component.atoms.size();++a) for(const auto & point:input.support.at(component.atoms[a]))
    {
        const auto row=view.LocalRow(static_cast<Eigen::Index>(point.row));
        if(row>=0) support[a].push_back({row,point.squared_distance});
    }
    const joint_component::Domain domain(y.size(),std::move(support));
    auto context=joint_component::CreateContext(y,static_cast<Eigen::Index>(component.atoms.size())); context.rank.rows=static_cast<Eigen::Index>(original_rows);
    context.scale=observation_scale;
    const auto geometry=joint_component::BuildTargetGeometry(domain,y,state,context);
    if(!geometry.valid || geometry.threshold_sensitive) {out.reason=geometry.reason; return out;}
    const auto & svd=geometry.svd; out.rank=static_cast<std::size_t>(svd.rank); out.rank_threshold=svd.threshold;
    if(component.rows.size()<=out.rank) {out.reason="insufficient-residual-degrees-of-freedom"; return out;}
    out.degrees_of_freedom=component.rows.size()-out.rank;
    const double variance=geometry.endpoint.residual.squaredNorm()/static_cast<double>(out.degrees_of_freedom);
    if(std::isfinite(variance)) out.residual_variance=variance;
    if(!(variance>0) || !std::isfinite(variance)) {out.reason="residual-variance-unavailable"; return out;}
    for(std::size_t a=0;a<component.atoms.size();++a) if(requested[component.atoms[a]])
    {
        const auto col=3*static_cast<Eigen::Index>(a);
        if(svd.right_vectors.block(col,svd.rank,3,svd.right_vectors.cols()-svd.rank).norm()>1e-10) continue;
        const Eigen::MatrixXd factor=geometry.scales.segment<3>(col).cwiseInverse().asDiagonal()*
            svd.right_vectors.block(col,0,3,svd.rank)*svd.singular_values.head(svd.rank).cwiseInverse().asDiagonal();
        const Eigen::Matrix3d covariance=variance*factor*factor.transpose();
        if(!covariance.allFinite()) {blocks.clear(); out.reason="nonfinite-covariance"; return out;}
        blocks.emplace(a,covariance);
    }
    out.status=EvidenceStatus::Available; return out;
}
StageUncertainty ComponentUncertainty(const JointProblemInput & input,
    const JointAnalysisComponent & component, const joint_component::ComponentView & view,
    const std::vector<bool> & requested, std::map<std::size_t, Eigen::Matrix3d> & blocks,std::size_t original_rows,double observation_scale)
{
    if(component.target_evidence) return TargetUncertainty(input,component,view,requested,blocks,original_rows,observation_scale);
    StageUncertainty out;
    out.method = "iid-ls-linearized";
    out.status = EvidenceStatus::Unavailable;
    if (!component.state) { out.reason = "missing-component-state"; return out; }
    if (component.runtime_convergence != JointCheckStatus::Passed)
    { out.reason = "component-not-converged"; return out; }
    const auto & state = *component.state;
    const auto atoms = component.atoms.size(), rows = component.rows.size();
    const auto columns = static_cast<Eigen::Index>(3 * atoms);
    for (std::size_t a = 0; a < atoms; ++a)
        if (!(state.ac.at(2 * a) > 0)) { out.reason = "active-amplitude-boundary"; return out; }
    if (rows <= 3 * atoms) { out.reason = "insufficient-residual-degrees-of-freedom"; return out; }
    out.degrees_of_freedom = rows - 3 * atoms;
    if (!std::equal(component.rows.begin(), component.rows.end(), view.rows.begin(), view.rows.end()) ||
        !std::equal(component.atoms.begin(), component.atoms.end(), view.atoms.begin(), view.atoms.end()))
        throw std::invalid_argument("Joint uncertainty component mapping mismatch.");
    using Membership = std::pair<std::size_t, const JointSupport *>;
    std::vector<std::vector<Membership>> tiles((rows + joint_component::derivative_tile_rows - 1) / joint_component::derivative_tile_rows);
    Eigen::VectorXd prediction = Eigen::VectorXd::Zero(static_cast<Eigen::Index>(rows));
    Eigen::VectorXd norms = Eigen::VectorXd::Zero(columns);
    const auto values = [&](std::size_t a, double squared) {
        const auto basis = joint_component::EvaluateKernel(squared, state.b.at(a), 2.5);
        return Eigen::Vector3d{basis.gaussian, basis.charge,
            state.ac[2 * a] * basis.gaussian_log_width + state.ac[2 * a + 1] * basis.charge_log_width};
    };
    for (std::size_t a = 0; a < atoms; ++a)
        for (const auto & support : input.support.at(component.atoms[a]))
        {
            const auto row = view.LocalRow(static_cast<Eigen::Index>(support.row));
            if (row < 0) continue; // Analytically profiled nuisance rows carry no remaining information.
            tiles[static_cast<std::size_t>(row / joint_component::derivative_tile_rows)].emplace_back(a, &support);
            const auto v = values(a, support.squared_distance);
            prediction(row) += state.ac[2 * a] * v(0) + state.ac[2 * a + 1] * v(1);
            for (Eigen::Index k = 0; k < 3; ++k)
            {
                const auto column = static_cast<Eigen::Index>(3 * a) + k;
                norms(column) = std::hypot(norms(column), v(k));
            }
        }
    if (!norms.allFinite() || (norms.array() <= 0).any())
    { out.reason = "rank-deficient-jacobian"; return out; }
    joint_component::TiledQR reduced(columns, 0);
    for (Eigen::Index first = 0; first < prediction.size(); first += joint_component::derivative_tile_rows)
    {
        const auto count = std::min(joint_component::derivative_tile_rows, prediction.size() - first);
        Eigen::MatrixXd tile = Eigen::MatrixXd::Zero(count, columns);
        for (const auto & [a, support] : tiles[static_cast<std::size_t>(first / joint_component::derivative_tile_rows)])
        {
            const auto row = view.LocalRow(static_cast<Eigen::Index>(support->row));
            const auto column = static_cast<Eigen::Index>(3 * a);
            tile.block<1, 3>(row - first, column) =
                (values(a, support->squared_distance).array() / norms.segment<3>(column).array()).matrix().transpose();
        }
        reduced.Append(tile, Eigen::MatrixXd(count, 0));
    }
    // Full-component right singular vectors retain the coupling to C and neighbors.
    const auto svd=joint_component::CompactSvd(reduced.r,
        std::numeric_limits<double>::epsilon()*static_cast<double>(std::max(original_rows,3*atoms)),
        -1,nullptr,joint_component::CompactSvdVectors::Right);
    if (!svd.valid)
    { out.reason = "jacobian-factorization-failed"; return out; }
    out.rank = static_cast<std::size_t>(svd.rank);
    out.rank_threshold = svd.threshold;
    if (svd.rank != columns) { out.reason = "rank-deficient-jacobian"; return out; }
    for (std::size_t row = 0; row < rows; ++row)
        prediction(static_cast<Eigen::Index>(row)) -= input.observations.at(component.rows[row]);
    const double variance = prediction.squaredNorm() / static_cast<double>(out.degrees_of_freedom);
    if (std::isfinite(variance)) out.residual_variance = variance;
    if (!(variance > 0) || !std::isfinite(variance))
    { out.reason = "residual-variance-unavailable"; return out; }
    for (std::size_t a = 0; a < atoms; ++a)
    {
        if (!requested[component.atoms[a]]) continue;
        const auto column = static_cast<Eigen::Index>(3 * a);
        const Eigen::MatrixXd factor = norms.segment<3>(column).cwiseInverse().asDiagonal() *
            svd.right_vectors.middleRows(column, 3) * svd.singular_values.cwiseInverse().asDiagonal();
        const Eigen::Matrix3d covariance = variance * factor * factor.transpose();
        if (!covariance.allFinite()) { out.reason = "nonfinite-covariance"; blocks.clear(); return out; }
        blocks.emplace(a, covariance);
    }
    out.status = EvidenceStatus::Available;
    return out;
}
}

std::map<int, StageUncertainty> ComputeJointUncertainty(const JointProblem & problem, const JointAnalysisResult & result,
    std::optional<std::span<const std::size_t>> outputs)
{
    joint_component::ResourcePhase phase("uncertainty");
    if (result.atom_ids != problem.Input().atom_ids || result.row_ids != problem.Input().row_ids)
        throw std::invalid_argument("Joint uncertainty snapshot identity mismatch.");
    const auto requested = JointOutputMask(problem.Input(), outputs);
    const auto & partition = JointProblemAccess::Get(problem).partition;
    std::map<int, StageUncertainty> output;
    for (const auto & component : result.components)
    {
        if (std::none_of(component.atoms.begin(), component.atoms.end(), [&](auto a) { return requested.at(a); })) continue;
        std::map<std::size_t, Eigen::Matrix3d> blocks;
        const auto index = partition.mappings->atom_component.at(component.atoms.at(0));
        auto profiled=component;
        auto view=partition.components.at(static_cast<std::size_t>(index));
        if(component.layout)
        {
            profiled.atoms=component.layout->full_atoms; profiled.rows=component.layout->informative_rows;
            auto mappings=std::make_shared<joint_component::PartitionMappings>();
            mappings->row_component.assign(problem.Input().observations.size(),-1);
            mappings->row_to_local.assign(problem.Input().observations.size(),-1);
            for(std::size_t r=0;r<profiled.rows.size();++r)
            {mappings->row_component[profiled.rows[r]]=0; mappings->row_to_local[profiled.rows[r]]=static_cast<Eigen::Index>(r);}
            view.rows.assign(profiled.rows.begin(),profiled.rows.end()); view.atoms.assign(profiled.atoms.begin(),profiled.atoms.end());
            view.mappings=std::move(mappings); view.component_index=0;
            for(const auto & group:component.layout->groups) for(auto atom:group.atoms) if(requested.at(atom))
            {StageUncertainty missing; missing.status=EvidenceStatus::Unavailable; missing.reason="observable-contribution-only"; output.emplace(std::stoi(result.atom_ids.at(atom)),missing);}
        }
        if(profiled.atoms.empty()) continue;
        const auto status = ComponentUncertainty(problem.Input(), profiled, view, requested, blocks,component.rows.size(),result.observation_scale);
        for (std::size_t a = 0; a < profiled.atoms.size(); ++a)
        {
            if (!requested[profiled.atoms[a]]) continue;
            auto uncertainty = status;
            if (status.status == EvidenceStatus::Available)
            {
                const auto found=blocks.find(a);
                if(found!=blocks.end()) uncertainty.covariance=found->second;
                else {uncertainty.status=EvidenceStatus::Unavailable; uncertainty.reason="nonunique-parameter-directions";}
            }
            output.emplace(std::stoi(result.atom_ids.at(profiled.atoms[a])), std::move(uncertainty));
        }
    }
    return output;
}

GroupParameterEvidence BuildJointParameterEvidence(const LocalStageEstimate & stage)
{
    GroupParameterEvidence evidence;
    evidence.atom_id = stage.source.atom_id;
    evidence.component_id = stage.source.component_id;
    evidence.source_id = stage.source.run_id;
    evidence.uncertainty_method = stage.uncertainty.method;
    evidence.status = EvidenceStatus::Ineligible;
    if (stage.source.role != FittingRole::Target) { evidence.reason = "not-a-target"; return evidence; }
    if (!stage.point) { evidence.reason = "missing-point-state"; return evidence; }
    const auto & point = *stage.point;
    if (!(point.GetAmplitude() > 0)) { evidence.reason = "nonpositive-amplitude"; return evidence; }
    evidence.estimate = linearization_service::EncodeGaussianToParameterVector(point);
    if (!evidence.estimate->allFinite()) { evidence.estimate.reset(); evidence.reason = "nonfinite-encoding"; return evidence; }
    if (stage.uncertainty.status != EvidenceStatus::Available || !stage.uncertainty.covariance)
    { evidence.status = EvidenceStatus::Unavailable; evidence.reason = stage.uncertainty.reason; return evidence; }
    Eigen::Matrix<double, 2, 3> derivative;
    derivative << 1 / point.GetAmplitude(), 0, -3,
        0, 0, -2 / (point.GetWidth() * point.GetWidth());
    const Eigen::Matrix2d covariance = derivative * *stage.uncertainty.covariance * derivative.transpose();
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> spectrum(covariance);
    if (!covariance.allFinite() || spectrum.info() != Eigen::Success ||
        spectrum.eigenvalues()(0) <= std::numeric_limits<double>::epsilon() * 2 * spectrum.eigenvalues()(1))
    { evidence.status = EvidenceStatus::Unavailable; evidence.reason = "singular-parameter-covariance"; return evidence; }
    evidence.covariance = covariance;
    evidence.status = EvidenceStatus::Available;
    return evidence;
}
}
