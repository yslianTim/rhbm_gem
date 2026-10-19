#include "core/detail/FirstStageInitialization.hpp"
#include "core/detail/joint_component/Numerics.hpp"
#include "core/detail/joint_component/ParameterStability.hpp"
#include "core/detail/joint_component/Problem.hpp"
#include "core/detail/joint_component/WidthCalibration.hpp"

#include <rhbm_gem/core/JointComponentEstimator.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <stdexcept>
#include <vector>

namespace {
namespace core = rhbm_gem::core;
namespace joint = rhbm_gem::core::joint_component;

constexpr std::size_t kAtomCount = 25;
constexpr double kTruthWidth = 5.0;
constexpr double kRegularWidth = 0.55;
constexpr double kAmplitude = 2.0;
constexpr double kCharge = 0.12;

struct SyntheticInput
{
    std::vector<double> truth_widths;
    rhbm_gem::ModelObject model;
    rhbm_gem::MapObject map;

    SyntheticInput() : model(MakeModel()), map(MakeMap(truth_widths)) {}

    static rhbm_gem::ModelObject MakeModel()
    {
        std::vector<std::unique_ptr<rhbm_gem::AtomObject>> atoms;
        atoms.reserve(kAtomCount);
        for (std::size_t index = 0; index < kAtomCount; ++index)
        {
            const int serial = static_cast<int>(index + 1);
            auto atom = std::make_unique<rhbm_gem::AtomObject>();
            atom->SetSerialID(serial);
            atom->SetSequenceID(serial);
            atom->SetComponentID("ALA");
            atom->SetAtomID("CA");
            atom->SetChainID("A");
            atom->SetIndicator(".");
            atom->SetComponentKey(1);
            atom->SetAtomKey(static_cast<AtomKey>(Spot::CA));
            atom->SetResidue(Residue::ALA);
            atom->SetElement(Element::CARBON);
            atom->SetSpot(Spot::CA);
            atom->SetStructure(Structure::FREE);
            atom->SetPosition(3.0 + 7.0 * static_cast<double>(index % 5),
                3.0 + 7.0 * static_cast<double>(index / 5), 3.0);
            atoms.push_back(std::move(atom));
        }
        rhbm_gem::ModelObject model(std::move(atoms));
        model.SelectAllAtoms();
        model.EditAnalysis().InitializeFromSelection();
        return model;
    }

    static rhbm_gem::MapObject MakeMap(std::vector<double> & truth_widths)
    {
        constexpr std::array<int, 3> grid{80, 80, 14};
        constexpr double spacing = 0.5;
        const std::array<double, 3> step{spacing, spacing, spacing};
        const std::array<double, 3> origin{0.0, 0.0, 0.0};
        const std::size_t voxel_count = static_cast<std::size_t>(grid[0]) * grid[1] * grid[2];
        auto values = std::make_unique<double[]>(voxel_count);
        std::fill_n(values.get(), voxel_count, 0.0);
        truth_widths.assign(kAtomCount, kRegularWidth);
        truth_widths.back() = kTruthWidth;
        for (std::size_t atom = 0; atom < kAtomCount; ++atom)
        {
            const double cx = 3.0 + 7.0 * static_cast<double>(atom % 5);
            const double cy = 3.0 + 7.0 * static_cast<double>(atom / 5);
            const double cz = 3.0;
            const int x0 = static_cast<int>(std::floor((cx - 2.5) / spacing));
            const int x1 = static_cast<int>(std::floor((cx + 2.5) / spacing));
            const int y0 = static_cast<int>(std::floor((cy - 2.5) / spacing));
            const int y1 = static_cast<int>(std::floor((cy + 2.5) / spacing));
            const int z0 = static_cast<int>(std::floor((cz - 2.5) / spacing));
            const int z1 = static_cast<int>(std::floor((cz + 2.5) / spacing));
            for (int z = z0; z <= z1; ++z)
                for (int y = y0; y <= y1; ++y)
                    for (int x = x0; x <= x1; ++x)
                    {
                        const double dx = spacing * x - cx;
                        const double dy = spacing * y - cy;
                        const double dz = spacing * z - cz;
                        const double square = dx * dx + dy * dy + dz * dz;
                        if (square > 6.25) continue;
                        const auto basis = joint::EvaluateKernel(square, truth_widths[atom], 2.5);
                        const std::size_t row = static_cast<std::size_t>(x) + grid[0] *
                            (static_cast<std::size_t>(y) + grid[1] * static_cast<std::size_t>(z));
                        const double noise = 1e-5 * std::sin(0.37 * x + 0.51 * y + 0.23 * z);
                        values[row] += kAmplitude * basis.gaussian + kCharge * basis.charge + noise;
                    }
        }
        return rhbm_gem::MapObject(grid, step, origin, std::move(values));
    }
};

const rhbm_gem::JointCheck * Check(const core::JointFitResult & fit, const char * name)
{
    const auto found = std::find_if(fit.evidence.begin(), fit.evidence.end(), [name](const auto & item) {
        return item.name == name;
    });
    return found == fit.evidence.end() ? nullptr : &*found;
}

double Rmse(const std::vector<double> & errors)
{
    long double sum = 0.0L;
    for (const double error : errors) sum += static_cast<long double>(error) * error;
    return std::sqrt(static_cast<double>(sum / static_cast<long double>(errors.size())));
}

void Run()
{
    SyntheticInput synthetic;
    const auto normalization = synthetic.map.MapValueArrayNormalization();
    if (!normalization) throw std::runtime_error("Synthetic map normalization failed.");
    auto problem = core::BuildPureHaloSupportExclusionProblem(synthetic.map, synthetic.model);
    const auto initialization = joint::InitializeContributors(synthetic.map, synthetic.model, problem);
    if (!initialization.valid)
        throw std::runtime_error("Synthetic FirstStage initialization failed: " + initialization.reason);
    const auto calibration = joint::CalibrateGlobalWidth(initialization, "production-bounded-synthetic");
    if (calibration.status != "Available" || !calibration.bmax4 || !calibration.p95)
        throw std::runtime_error("Synthetic qualified-donor calibration unavailable: " + calibration.unavailable_reason);
    if (calibration.qualified_donor_count < 20 || calibration.qualified_donor_coverage < 0.95)
        throw std::runtime_error("Synthetic FirstStage calibration failed its production preflight.");

    const double bound = *calibration.bmax4;
    auto fit = core::detail::FitJointComponentsInitializedBoundedWithProgress(problem, initialization, {});
    if (!fit.search_completed || fit.RuntimeConvergence() != rhbm_gem::JointCheckStatus::Passed || !fit.assembled_state)
        throw std::runtime_error("Synthetic bounded fit did not qualify.");
    const auto * replay = Check(fit, "independent-replay");
    const auto * width_kkt = Check(fit, "width-stationarity");
    const auto * bound_feasibility = Check(fit, "width-bound-feasibility");
    if (!replay || replay->status != rhbm_gem::JointCheckStatus::Passed ||
        !width_kkt || width_kkt->status != rhbm_gem::JointCheckStatus::Passed ||
        !bound_feasibility || bound_feasibility->status != rhbm_gem::JointCheckStatus::Passed)
        throw std::runtime_error("Synthetic independent replay, bound feasibility, or projected Width KKT failed.");

    std::size_t active{}, outside_truth{};
    double maximum_bound_violation{};
    double minimum_bound_slack=std::numeric_limits<double>::infinity();
    std::vector<double> amplitude_errors, width_errors, charge_errors,
        prediction_errors, truth_prediction_errors, truth_signal;
    const auto & problem_input = core::JointProblemAccess::Get(*fit.problem).input;
    std::vector<double> truth_prediction(problem_input->observations.size());
    for (std::size_t atom = 0; atom < fit.assembled_state->b.size(); ++atom)
    {
        const double estimate = fit.assembled_state->b[atom];
        minimum_bound_slack=std::min(minimum_bound_slack,bound-estimate);
        maximum_bound_violation = std::max(maximum_bound_violation, std::max(0.0, estimate - bound));
        if (estimate > bound + 1e-10 * std::max(1.0, bound))
            throw std::runtime_error("Synthetic estimate violated the global Width bound.");
        active += estimate >= bound - 1e-8 * std::max(1.0, bound);
        outside_truth += synthetic.truth_widths[atom] > bound;
        width_errors.push_back(estimate - synthetic.truth_widths[atom]);
        amplitude_errors.push_back(fit.assembled_state->ac[2 * atom] * *normalization - kAmplitude);
        charge_errors.push_back(fit.assembled_state->ac[2 * atom + 1] * *normalization - kCharge);
        for (const auto & support : problem_input->support[atom])
        {
            const auto basis = joint::EvaluateKernel(support.squared_distance,
                synthetic.truth_widths[atom], 2.5);
            truth_prediction[support.row] += kAmplitude * basis.gaussian + kCharge * basis.charge;
        }
    }
    if (!fit.prediction) throw std::runtime_error("Synthetic prediction unavailable.");
    const auto & observations = problem_input->observations;
    if (fit.prediction->size() != observations.size())
        throw std::runtime_error("Synthetic prediction row mapping mismatch.");
    for (std::size_t row = 0; row < observations.size(); ++row)
    {
        prediction_errors.push_back(fit.prediction->at(row) - observations[row]);
        truth_prediction_errors.push_back(truth_prediction[row] / *normalization - observations[row]);
        truth_signal.push_back(truth_prediction[row] / *normalization);
    }
    if (active == 0 || outside_truth == 0)
        throw std::runtime_error("Synthetic fixture did not activate the global bound against known truth.");

    const auto * ac_kkt = Check(fit, "conditional-ac-optimality");
    const rhbm_gem::JointCheck * active_set = nullptr;
    double active_set_kkt{};
    std::size_t active_set_solves{};
    for (const auto & component : fit.components)
        for (const auto & item : component.evidence)
            if (item.name == "active-set-lm-subproblems")
            {
                if (item.status != rhbm_gem::JointCheckStatus::Passed)
                    throw std::runtime_error("Synthetic active-set LM evidence failed.");
                active_set = &item;
                if (item.value) active_set_kkt = std::max(active_set_kkt, *item.value);
                const auto marker = item.reason.find("solves=");
                if (marker != std::string::npos)
                    active_set_solves += static_cast<std::size_t>(std::stoull(item.reason.substr(marker + 7)));
            }
    if (!active_set) throw std::runtime_error("Synthetic active-set LM evidence unavailable.");
    std::map<std::string,std::size_t> stability_status_counts;
    std::size_t stability_blocks{},stability_selected{},stability_selected_bound_active{},
        stability_selected_condition_bound_active{},stability_assessed{},
        stability_solves{},stability_evaluations{},conditioning_bound_active_blocks{};
    std::string selected_probe_ids,conditioning_bound_active_ids;
    double stability_cpu_seconds{},stability_wall_seconds{};
    std::size_t stability_memory_estimate{};
    for(const auto & component:fit.components)
    {
        conditioning_bound_active_blocks+=static_cast<std::size_t>(std::count_if(component.local_conditioning.begin(),
            component.local_conditioning.end(),[](const auto & block) {
                return block.active_width_bound_targets && *block.active_width_bound_targets>0;
            }));
        for(const auto & block:component.local_conditioning)
            if(block.active_width_bound_targets && *block.active_width_bound_targets>0)
            {
                if(!conditioning_bound_active_ids.empty()) conditioning_bound_active_ids+=',';
                conditioning_bound_active_ids+=block.block_id;
            }
        if(component.parameter_stability)
        {
            stability_wall_seconds+=component.parameter_stability->seconds;
            stability_cpu_seconds+=component.parameter_stability->cpu_seconds;
            if(component.parameter_stability->workspace_peak_memory_estimate_bytes)
                stability_memory_estimate=std::max(stability_memory_estimate,
                    *component.parameter_stability->workspace_peak_memory_estimate_bytes);
            stability_solves+=component.parameter_stability->perturbation_solves;
            stability_evaluations+=component.parameter_stability->perturbation_evaluations;
            for(std::size_t index=0;index<component.parameter_stability->blocks.size();++index)
            {
                const auto & block=component.parameter_stability->blocks[index];
                ++stability_blocks; ++stability_status_counts[block.status];
                stability_selected+=block.probe_selected ? 1u : 0u;
                if(block.probe_selected)
                {
                    if(!selected_probe_ids.empty()) selected_probe_ids+=',';
                    selected_probe_ids+=block.block_id;
                }
                stability_selected_bound_active+=(block.probe_selected && block.baseline_active_width_bound_targets &&
                    *block.baseline_active_width_bound_targets>0) ? 1u : 0u;
                const auto condition=std::find_if(component.local_conditioning.begin(),
                    component.local_conditioning.end(),[&](const auto & candidate) {
                        return candidate.block_id==block.block_id;
                    });
                stability_selected_condition_bound_active+=(block.probe_selected &&
                    condition!=component.local_conditioning.end() && condition->active_width_bound_targets &&
                    *condition->active_width_bound_targets>0) ? 1u : 0u;
                stability_assessed+=(block.status=="Stable" || block.status=="Sensitive" ||
                    block.status=="Unstable") ? 1u : 0u;
            }
        }
    }
    if(stability_selected>joint::MaximumPerturbedBlocksPerFit)
        throw std::runtime_error("Synthetic local stability exceeded its per-fit probe-block budget.");
    if(conditioning_bound_active_blocks>0 && stability_selected_condition_bound_active==0)
        throw std::runtime_error("Synthetic local stability did not select a Width-bound-active block.");
    std::cout << "targets=" << fit.assembled_state->b.size()
        << " donors=" << calibration.qualified_donor_count
        << " coverage=" << calibration.qualified_donor_coverage
        << " p95=" << *calibration.p95
        << " bmax=" << bound
        << " truth_outside_bound=" << outside_truth
        << " truth_coverage=" << 1.0 - static_cast<double>(outside_truth) /
            static_cast<double>(fit.assembled_state->b.size())
        << " active_estimates=" << active
        << " minimum_bound_slack=" << minimum_bound_slack
        << " maximum_bound_violation=" << maximum_bound_violation
        << " objective=" << fit.objective.value_or(std::numeric_limits<double>::quiet_NaN())
        << " ac_kkt=" << (ac_kkt && ac_kkt->value ? *ac_kkt->value : -1.0)
        << " projected_width_kkt=" << *width_kkt->value
        << " active_set_kkt=" << active_set_kkt
        << " active_set_solves=" << active_set_solves
        << " active_set_status=Passed"
        << " width_bound_feasibility=Passed"
        << " A_rmse=" << Rmse(amplitude_errors)
        << " B_rmse=" << Rmse(width_errors)
        << " C_rmse=" << Rmse(charge_errors)
        << " prediction_rmse=" << Rmse(prediction_errors)
        << " truth_prediction_rmse=" << Rmse(truth_prediction_errors)
        << " truth_signal_rmse=" << Rmse(truth_signal)
        << " prediction_error_pct_truth_signal=" << 100.0 * Rmse(prediction_errors) / Rmse(truth_signal)
        << " search_completed=" << fit.search_completed
        << " replay=Passed"
        << " truth_used_for_calibration=false"
        << " runtime_convergence=Passed"
        << " local_stability_blocks=" << stability_blocks
        << " local_stability_selected=" << stability_selected
        << " local_stability_selected_bound_active=" << stability_selected_bound_active
        << " local_stability_selected_condition_bound_active=" << stability_selected_condition_bound_active
        << " local_conditioning_bound_active_blocks=" << conditioning_bound_active_blocks
        << " local_conditioning_bound_active_ids=" << conditioning_bound_active_ids
        << " local_stability_selected_ids=" << selected_probe_ids
        << " local_stability_assessed=" << stability_assessed
        << " local_stability_stable=" << stability_status_counts["Stable"]
        << " local_stability_sensitive=" << stability_status_counts["Sensitive"]
        << " local_stability_unstable=" << stability_status_counts["Unstable"]
        << " local_stability_unavailable=" << stability_status_counts["Unavailable"]
        << " local_stability_not_assessed=" << stability_status_counts["NotAssessed"]
        << " local_stability_solves=" << stability_solves
        << " local_stability_evaluations=" << stability_evaluations
        << " local_stability_cpu_seconds=" << stability_cpu_seconds
        << " local_stability_wall_seconds=" << stability_wall_seconds
        << " local_stability_memory_estimate_bytes=" << stability_memory_estimate << '\n';
}
}

int main()
{
    try
    {
        Run();
        return 0;
    }
    catch (const std::exception & error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
