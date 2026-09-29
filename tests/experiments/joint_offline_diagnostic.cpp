#include "support/JointPartialSelection.hpp"
#include "support/JointPrecisionAudit.hpp"
#include "support/JointComponentChecks.hpp"
#include "support/JointRuntimeJson.hpp"
#include "core/detail/joint_component/Problem.hpp"
#include "data/io/detail/JointResultJson.hpp"
#include <fstream>
#include <iostream>

namespace {
namespace j = boost::json;
namespace c = rhbm_gem::core;
namespace n = c::joint_component;
namespace p = second_stage_test::matched::joint_abc;
namespace r = second_stage_test::matched::runtime_json;
namespace audit = second_stage_test::matched::certification;
using V = Eigen::VectorXd;
using Path = std::filesystem::path;

void Write(const Path & path, const j::value & value)
{
    const auto temporary = Path(path.string() + ".tmp");
    std::ofstream out(temporary);
    out.exceptions(std::ios::failbit | std::ios::badbit);
    out << j::serialize(value) << '\n';
    out.close();
    std::filesystem::rename(temporary, path);
}

V Vector(const std::vector<double> & values)
{
    return Eigen::Map<const V>(values.data(), static_cast<Eigen::Index>(values.size()));
}

j::value Outcome(const c::JointFitResult & fit)
{
    j::parse_options options;
    options.numbers = j::number_precision::precise;
    auto result = j::parse(rhbm_gem::joint_result_io::Encode(c::CaptureJointAnalysisResult(fit)), {}, options);
    if (fit.assembled_state)
    {
        const auto & state = result.at("assembled_state");
        if (j::value_to<std::vector<double>>(state.at("ac")) != fit.assembled_state->ac ||
            j::value_to<std::vector<double>>(state.at("log_b")) != fit.assembled_state->log_b)
            throw std::runtime_error("Diagnostic JSON changed the runtime state.");
    }
    return result;
}

void Run(int start, const Path & output)
{
    auto fixture = joint_partial_test::Make("weak");
    const auto initialized = c::EstimateJointComponents(*fixture.map, *fixture.model);
    auto initial = initialized.initialization.b;
    if (start == 1) initial = {.5, .4};
    else if (start == 2) initial = {.6, .5};
    else if (start != 0) throw std::invalid_argument("Weak start must be 0, 1 or 2.");
    const auto fit = start == 0 ? initialized : c::FitJointComponents(*initialized.problem, initial);
    const auto & problem = *fit.problem;
    const auto & data = c::JointProblemAccess::Get(problem);
    const V y = data.y;
    const p::Domain domain(data.domain);
    j::array support;
    for (const auto & atom : problem.Input().support)
    {
        j::array entries;
        for (const auto & sample : atom) entries.push_back(j::array{sample.row, sample.squared_distance});
        support.push_back(std::move(entries));
    }
    j::object snapshot{{"observations", j::value_from(problem.Input().observations)},
        {"support", support}, {"atom_ids", j::value_from(problem.Input().atom_ids)},
        {"row_ids", j::value_from(problem.Input().row_ids)}, {"scale", problem.ObservationScale()},
        {"truth_a", j::value_from(fixture.a)}, {"truth_b", j::value_from(fixture.b)},
        {"truth_c", j::value_from(fixture.c)}};
    j::object out{{"start", start}, {"initial_b", j::value_from(initial)},
        {"production_initial_b", j::value_from(initialized.initialization.b)},
        {"snapshot", snapshot}, {"outcome", Outcome(fit)}};
    Write(output, out);
    if (!fit.assembled_state) return;

    const V eta = Vector(fit.assembled_state->log_b);
    const V beta = Vector(fit.assembled_state->ac);
    const auto evaluation = p::AtState(domain, y, eta, beta, data.context);
    const auto assessment = n::AssessProfile(domain, y, eta, data.context, &beta);
    out["actual_state_assessment"] = r::Assessment(assessment);
    out["same_state"] = p::SameState(domain, y, eta, beta, data.partition, data.context);
    const auto dense = p::DenseDifferentiate(evaluation, data.context.scale, &data.context);
    const auto tiled = p::MaterializeDerivative(evaluation, data.context.scale, &data.context);
    out["dense_tiled_jacobian_relative_difference"] =
        (dense.jacobian - tiled.jacobian).norm() / std::max(1e-300, dense.jacobian.norm());
    n::Matrix directions(eta.size(), 2);
    directions.col(0) = assessment.correction.normalized();
    directions.col(1) = assessment.weak_directions.col(0).normalized();
    out["directions"] = j::array{r::Values(directions.col(0)), r::Values(directions.col(1))};
    out["precision"] = audit::PrecisionAudit(domain, y, evaluation, directions, &data.context);
    Write(output, out);

    j::array scans;
    const auto base = n::EvaluateProfile(domain, y, eta, false, &data.context);
    for (int direction = 0; direction < 2; ++direction)
        for (int exponent = -6; exponent <= -1; ++exponent)
            for (double sign : {-1., 1.})
            {
                const double step = sign * std::pow(10., exponent);
                const V at = eta + step * directions.col(direction);
                const auto trial = n::EvaluateProfile(domain, y, at, false, &data.context);
                j::object row{{"direction", direction}, {"step", step}, {"eta", r::Values(at)},
                              {"endpoint", r::Endpoint(trial)}};
                if (trial.valid)
                {
                    row["trust"] = r::Trust(n::CheckTrust(domain, y, trial, data.context));
                    row["objective_change"] = (trial.residual.squaredNorm() - base.residual.squaredNorm()) /
                        (2 * data.context.scale * data.context.scale);
                    row["prediction_change_norm"] = (trial.residual - base.residual).norm();
                }
                if (exponent == -3 || exponent == -1)
                    row["precision_profile"] = audit::PrecisionProfileChange(domain, y, eta, at);
                scans.push_back(std::move(row));
                out["scans"] = scans;
                Write(output, out);
            }
    auto extended = data.context;
    extended.profile_budget = 1000;
    extended.update_budget = 500;
    out["diagnostic_restart"] = p::FitComponent(
        data.partition.components.front(), y, eta.array().exp(), extended);
    out["audit_complete"] = true;
    Write(output, out);
}
}

int main(int argc, char ** argv)
{
    try
    {
        Eigen::setNbThreads(1);
        if (argc != 4 || std::string(argv[1]) != "weak")
            throw std::invalid_argument("Usage: joint_offline_diagnostic weak START OUTPUT");
        Run(std::stoi(argv[2]), argv[3]);
        return 0;
    }
    catch (const std::exception & error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
