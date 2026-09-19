#pragma once
#include "support/JointABCProfile.hpp"

namespace second_stage_test::matched::joint_abc {
struct ComponentView
{
    std::string id;
    std::vector<Eigen::Index> atoms, rows, atom_to_local, row_to_local;
    Domain domain{0,{}};
};
struct ComponentPartition
{
    std::vector<ComponentView> components;
    std::vector<Eigen::Index> atom_component, row_component, constant_rows, unobserved_atoms;
};
ComponentPartition BuildPartition(const Domain &, const std::vector<std::string> & atom_ids);
Eigen::VectorXd Select(const Eigen::VectorXd &, const std::vector<Eigen::Index> &);
EvaluationContext ComponentContext(const EvaluationContext &, const ComponentView &, bool independent_search);
boost::json::object Census(const Domain &, const ComponentPartition &, const EvaluationContext &);
Evaluation EvaluatePartitioned(const Domain &, const Eigen::VectorXd &, const Eigen::VectorXd & eta,
    const ComponentPartition &, const EvaluationContext &, bool reference=false);
Evaluation ComponentEvaluation(const Evaluation &, const ComponentView &);
Differential DifferentiatePartitioned(const Evaluation &, const ComponentPartition &, const EvaluationContext &);
boost::json::object SameState(const Domain &, const Eigen::VectorXd &, const Eigen::VectorXd & eta,
    const Eigen::VectorXd & beta, const ComponentPartition &, const EvaluationContext &);
boost::json::object FitComponent(const ComponentView &, const Eigen::VectorXd & parent_y,
    const Eigen::VectorXd & parent_initial_b, const EvaluationContext &);
boost::json::object Assemble(const Domain &, const Eigen::VectorXd &, const ComponentPartition &,
    const EvaluationContext &, const boost::json::array & component_fits);
boost::json::object FitComponents(const Domain &, const Eigen::VectorXd &, const Eigen::VectorXd & initial_b,
    const ComponentPartition &, const EvaluationContext &);
} // namespace second_stage_test::matched::joint_abc
