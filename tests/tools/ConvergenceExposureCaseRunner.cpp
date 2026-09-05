#include <CLI/CLI.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "core/detail/IterationProcess.hpp"
#include <rhbm_gem/core/GaussianEstimator.hpp>
#include <rhbm_gem/core/EstimatorTester.hpp>
#include <rhbm_gem/core/TestDataFactory.hpp>
#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>

namespace {

namespace rg = rhbm_gem;
namespace rt = rhbm_gem::core;
namespace rt_detail = rhbm_gem::core::detail;

struct Request
{
    std::string case_id;
    std::string family;
    std::string topology;
    int level{ 0 };
    int replica{ 0 };
    std::uint32_t seed{ 0 };
    int threads{ 1 };
    double noise_sigma{ 0.0 };
    double separation{ 1.0e-4 };
    double initial_perturbation{ 8.0e-4 };
    int active_target{ 16 };
    double largest_group_ratio{ 0.8 };
};

struct Scenario
{
    std::unique_ptr<rg::ModelObject> model;
    std::vector<rg::GaussianModel3D> truth;
};

double Distance(
    const std::array<double, 3> & lhs,
    const std::array<double, 3> & rhs)
{
    const auto dx{ lhs.at(0) - rhs.at(0) };
    const auto dy{ lhs.at(1) - rhs.at(1) };
    const auto dz{ lhs.at(2) - rhs.at(2) };
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

std::unique_ptr<rg::AtomObject> MakeAtom(
    int serial_id,
    Spot spot,
    Element element,
    ComponentKey component_key,
    const std::array<double, 3> & position)
{
    auto atom{ std::make_unique<rg::AtomObject>() };
    atom->SetSerialID(serial_id);
    atom->SetChainID("A");
    atom->SetSequenceID(serial_id);
    atom->SetComponentKey(component_key);
    atom->SetAtomKey(static_cast<AtomKey>(spot));
    atom->SetElement(element);
    atom->SetSpot(spot);
    atom->SetPosition(position);
    return atom;
}

rg::LocalGaussianResult MakeGaussianResult(const rg::GaussianModel3D & model)
{
    rg::LocalGaussianResult result;
    result.ols = rg::GaussianModel3DWithUncertainty{
        model, rg::GaussianModel3DUncertainty{} };
    result.mdpde = result.ols;
    return result;
}

LocalPotentialSampleList BuildSamples(
    const rg::AtomObject & target,
    const std::vector<rg::AtomObject *> & atom_list,
    const std::vector<rg::GaussianModel3D> & truth,
    double noise_sigma,
    std::mt19937 & generator)
{
    constexpr std::array<std::array<double, 3>, 6> directions{{
        { 1.0, 0.0, 0.0 }, { -1.0, 0.0, 0.0 },
        { 0.0, 1.0, 0.0 }, { 0.0, -1.0, 0.0 },
        { 0.0, 0.0, 1.0 }, { 0.0, 0.0, -1.0 }
    }};
    constexpr std::array<double, 4> radii{ 0.15, 0.35, 0.65, 0.95 };
    std::normal_distribution<double> noise{ 0.0, noise_sigma };
    LocalPotentialSampleList result;
    result.reserve(directions.size() * radii.size());
    for (const auto radius : radii)
    {
        for (const auto & direction : directions)
        {
            SamplingPoint point;
            point.distance = radius;
            for (std::size_t axis = 0; axis < point.position.size(); axis++)
            {
                point.position.at(axis) = target.GetPosition().at(axis) +
                    radius * direction.at(axis);
            }
            double response{ 0.0 };
            for (std::size_t index = 0; index < atom_list.size(); index++)
            {
                const auto distance{
                    Distance(point.position, atom_list.at(index)->GetPosition()) };
                if (distance <= 5.0)
                {
                    response += truth.at(index).ResponseAtDistance(distance);
                }
            }
            const auto scale{ std::max(1.0, std::abs(response)) };
            result.emplace_back(LocalPotentialSample{
                response + scale * noise(generator), point });
        }
    }
    return result;
}

rg::GaussianModel3D PerturbModel(
    const rg::GaussianModel3D & truth,
    double perturbation,
    std::size_t index)
{
    const auto sign{ index % 2 == 0 ? 1.0 : -1.0 };
    return rg::GaussianModel3D{
        truth.GetAmplitude() * std::exp(sign * perturbation),
        truth.GetWidth() * std::exp(-0.5 * sign * perturbation),
        truth.GetOffset() + sign * perturbation * truth.GetHeight()
    };
}

std::vector<std::array<double, 3>> BuildPositions(
    const Request & request,
    std::size_t atom_count)
{
    std::vector<std::array<double, 3>> positions;
    positions.reserve(atom_count);
    if (request.topology == "near-collinear-pair")
    {
        return {{ 0.0, 0.0, 0.0 },
            { request.separation, 0.0, 0.0 }};
    }
    if (request.topology == "dual-near-collinear")
    {
        return {{ 0.0, 0.0, 0.0 },
            { request.separation, 0.0, 0.0 },
            { 10.0, 0.0, 0.0 },
            { 10.0 + request.separation, 0.0, 0.0 }};
    }
    for (std::size_t index = 0; index < atom_count; index++)
    {
        const auto cluster{ index / 16 };
        const auto member{ index % 16 };
        const auto spacing{
            request.topology == "boundary-conflict" ? 0.45 : 0.40 };
        positions.push_back({
            static_cast<double>(cluster) * 10.0 +
                static_cast<double>(member) * spacing,
            static_cast<double>(index % 3) * 0.05,
            0.0 });
    }
    return positions;
}

std::size_t GetAtomCount(const Request & request)
{
    if (request.family == "natural")
    {
        if (request.topology == "unk-c") return 1;
        if (request.topology == "o-o") return 2;
        if (request.topology == "ca-c") return 5;
        return 4;
    }
    if (request.family == "stationarity")
    {
        if (request.topology == "near-collinear-pair") return 2;
        if (request.topology == "dual-near-collinear") return 4;
        if (request.topology == "eight-atom-chain") return 8;
        if (request.topology == "boundary-conflict") return 13;
        return 16;
    }
    if (request.topology == "unbalanced-shared-groups") return 104;
    if (request.topology == "combined-fixed-group-imbalance") return 128;
    return 100;
}

ComponentKey GetComponentKey(
    const Request & request,
    std::size_t atom_index,
    std::size_t atom_count)
{
    if (request.family != "population") return 1;
    if (request.topology == "unbalanced-shared-groups" ||
        request.topology == "combined-fixed-group-imbalance")
    {
        const auto largest_group_size{ static_cast<std::size_t>(std::clamp(
            request.largest_group_ratio, 0.0, 1.0) *
            static_cast<double>(atom_count)) };
        return atom_index < largest_group_size ? 1 :
            static_cast<ComponentKey>(atom_index - largest_group_size + 2);
    }
    return static_cast<ComponentKey>(atom_index / 8 + 1);
}

std::vector<Spot> BuildSpots(const Request & request, std::size_t atom_count)
{
    std::vector<Spot> spots(atom_count, Spot::C);
    if (request.family != "natural") return spots;
    if (request.topology == "unk-c") spots.front() = Spot::UNK;
    if (request.topology == "o-o") spots.front() = Spot::O;
    if (request.topology == "n-n") spots.front() = Spot::N;
    if (request.topology == "ca-c") spots.front() = Spot::CA;
    return spots;
}

Scenario BuildScenario(const Request & request)
{
    if (request.family == "natural")
    {
        const auto spot{
            request.topology == "o-o" ? Spot::O :
            request.topology == "n-n" ? Spot::N :
            request.topology == "c-c" ? Spot::C :
            request.topology == "ca-c" ? Spot::CA : Spot::UNK };
        const auto element{
            spot == Spot::O || request.topology == "unk-o" ? Element::OXYGEN :
            spot == Spot::N || request.topology == "unk-n" ? Element::NITROGEN :
            Element::CARBON };
        const auto charge{
            element == Element::OXYGEN ? -0.4 :
            element == Element::NITROGEN ? -0.1 : 0.3 };
        ElectricPotential potential_model;
        potential_model.SetModelChoice(0);
        potential_model.SetBlurringWidth(0.5);
        rt::PotentialModelScenario reference_scenario;
        reference_scenario.spot = Spot::UNK;
        reference_scenario.element = element;
        reference_scenario.charge = charge;
        reference_scenario.potential_model = potential_model;
        reference_scenario.data_error_sigma = 0.0;
        reference_scenario.replica_size = 1;
        reference_scenario.random_seed = 0;
        rt::FitOptions reference_options;
        reference_options.distance_min = 0.0;
        reference_options.distance_max = 1.0;
        reference_options.thread_size = request.threads;
        reference_options.quiet_mode = true;
        const auto reference_gaussian{
            rt::EstimateAtomicModelFullStageMean(
                rt::BuildPotentialModelTestData(reference_scenario),
                reference_options)
        };
        rt::PotentialModelScenario natural_scenario;
        natural_scenario.spot = spot;
        natural_scenario.element = element;
        natural_scenario.charge = charge;
        natural_scenario.gaus_true = reference_gaussian;
        natural_scenario.potential_model = potential_model;
        natural_scenario.data_error_sigma = request.noise_sigma;
        natural_scenario.replica_size = 1;
        natural_scenario.random_seed = request.seed;
        auto input{ rt::BuildPotentialModelTestData(natural_scenario) };
        std::vector<rg::GaussianModel3D> truth{
            natural_scenario.gaus_true };
        return Scenario{
            std::move(input.replica_model_objects.front()), std::move(truth) };
    }

    const auto atom_count{ GetAtomCount(request) };
    const auto positions{ BuildPositions(request, atom_count) };
    const auto spots{ BuildSpots(request, atom_count) };
    std::mt19937 generator{ request.seed };
    std::uniform_real_distribution<double> amplitude_scale{ 0.80, 1.20 };
    std::uniform_real_distribution<double> width{ 0.45, 0.65 };
    std::vector<rg::GaussianModel3D> truth;
    std::vector<std::unique_ptr<rg::AtomObject>> atoms;
    truth.reserve(atom_count);
    atoms.reserve(atom_count);
    for (std::size_t index = 0; index < atom_count; index++)
    {
        const auto component_key{ GetComponentKey(request, index, atom_count) };
        const auto amplitude{
            (6.0 + static_cast<double>(index % 3)) * amplitude_scale(generator) };
        const auto weak_peak{
            request.topology == "weak-peak-shared" && index + 1 == atom_count };
        truth.emplace_back(
            weak_peak ? 0.25 * amplitude : amplitude,
            width(generator),
            0.02 * static_cast<double>(static_cast<int>(component_key % 5) - 2));
        atoms.emplace_back(MakeAtom(
            static_cast<int>(index + 1),
            spots.at(index),
            spots.at(index) == Spot::O ? Element::OXYGEN :
                spots.at(index) == Spot::N ? Element::NITROGEN :
                Element::CARBON,
            component_key,
            positions.at(index)));
    }

    auto model{ std::make_unique<rg::ModelObject>(std::move(atoms)) };
    model->SelectAllAtoms();
    auto analysis{ model->EditAnalysis() };
    analysis.RebuildAtomGroupsFromSelection();
    analysis.InitializeGroupAlpha(0.0);
    std::vector<rg::AtomObject *> selected_atoms{ model->GetSelectedAtoms() };
    for (std::size_t index = 0; index < atom_count; index++)
    {
        auto * atom{ selected_atoms.at(index) };
        analysis.SetAtomLocalAlphaR(rg::FittingStage::Second, *atom, 0.0);
        auto perturbation{ request.initial_perturbation };
        if (request.family == "population")
        {
            const auto active_count{ std::min<std::size_t>(
                static_cast<std::size_t>(request.active_target), atom_count) };
            perturbation = index == 0 ? request.initial_perturbation : 0.0;
            if (index >= active_count) perturbation = 0.0;
        }
        analysis.SetAtomLocalGaussianResult(
            rg::FittingStage::Second,
            *atom,
            MakeGaussianResult(PerturbModel(
                truth.at(index), perturbation, index)));
        analysis.SetAtomLocalRawSamplingEntries(*atom, BuildSamples(
            *atom,
            selected_atoms,
            truth,
            request.noise_sigma,
            generator));
    }

    if (request.family == "population")
    {
        const auto active_count{ std::min<std::size_t>(
            static_cast<std::size_t>(request.active_target), atom_count) };
        for (std::size_t index = active_count; index < atom_count; index++)
        {
            auto entries{
                rg::AtomLocalPotentialView::For(*selected_atoms.at(index))
                    .GetRawSamplingEntries(false) };
            if (request.topology == "quarantine-dilution")
            {
                const auto target_position{ selected_atoms.at(index)->GetPosition() };
                entries.resize(256, entries.back());
                entries.front().response = 0.0;
                entries.front().point.distance = 0.0;
                entries.front().point.position = target_position;
                for (std::size_t sample_index = 1;
                    sample_index < entries.size(); sample_index++)
                {
                    auto & sample{ entries.at(sample_index) };
                    const auto response_scale{
                        0.5 + 0.5 * static_cast<double>(sample_index) /
                            static_cast<double>(entries.size()) };
                    sample.response = 1.0e4 * response_scale;
                    sample.point.position = target_position;
                    sample.point.position.at(0) += 100.0;
                    sample.point.distance = 100.0;
                }
            }
            else if (request.topology != "unbalanced-shared-groups")
            {
                entries.resize(request.topology == "offset-fixed-dilution" ? 2 : 1);
            }
            analysis.SetAtomLocalRawSamplingEntries(
                *selected_atoms.at(index), std::move(entries));
        }
    }
    return Scenario{ std::move(model), std::move(truth) };
}

void LogScenario(const Request & request, const Scenario & scenario)
{
    std::cout << "Convergence exposure scenario: schema=1"
        << ", case=" << request.case_id
        << ", family=" << request.family
        << ", topology=" << request.topology
        << ", level=" << request.level
        << ", replica=" << request.replica
        << ", seed=" << request.seed
        << ", atoms=" << scenario.model->GetSelectedAtoms().size()
        << ", truth-atoms=" << scenario.truth.size() << ".\n";
    for (std::size_t index = 0; index < scenario.truth.size(); index++)
    {
        const auto & model{ scenario.truth.at(index) };
        std::cout << "Convergence exposure truth: schema=1"
            << ", case=" << request.case_id
            << ", serial=" << index + 1
            << ", amplitude=" << model.GetAmplitude()
            << ", width=" << model.GetWidth()
            << ", offset=" << model.GetOffset() << ".\n";
    }
}

Request ParseRequest(int argc, char ** argv)
{
    Request request;
    CLI::App app{ "Run one build-gated convergence-exposure case" };
    app.add_option("--case-id", request.case_id)->required();
    app.add_option("--family", request.family)->required()
        ->check(CLI::IsMember({ "natural", "stationarity", "population" }));
    app.add_option("--topology", request.topology)->required();
    app.add_option("--level", request.level)->required()->check(CLI::Range(0, 4));
    app.add_option("--replica", request.replica)->required()->check(CLI::Range(0, 7));
    app.add_option("--seed", request.seed)->required();
    app.add_option("--threads", request.threads)->required()->check(CLI::PositiveNumber);
    app.add_option("--noise-sigma", request.noise_sigma)->check(CLI::NonNegativeNumber);
    app.add_option("--separation", request.separation)->check(CLI::PositiveNumber);
    app.add_option("--initial-perturbation", request.initial_perturbation)
        ->check(CLI::NonNegativeNumber);
    app.add_option("--active-target", request.active_target)->check(CLI::PositiveNumber);
    app.add_option("--largest-group-ratio", request.largest_group_ratio)
        ->check(CLI::Range(0.0, 1.0));
    app.parse(argc, argv);
    return request;
}

} // namespace

int main(int argc, char ** argv)
{
    try
    {
        const auto request{ ParseRequest(argc, argv) };
        auto scenario{ BuildScenario(request) };
        LogScenario(request, scenario);
        Logger::SetLogLevel(LogLevel::Debug);
        rt::FitOptions options;
        options.distance_min = 0.0;
        options.distance_max = 1.0;
        options.thread_size = request.threads;
        options.quiet_mode = false;
        if (request.family == "natural")
        {
            rt::RunPotentialFittingWorkflow(*scenario.model, options);
        }
        else
        {
            rt_detail::RunSecondStageIterations(*scenario.model, options);
        }
        return 0;
    }
    catch (const std::exception & error)
    {
        std::cerr << "Convergence exposure case failed: " << error.what() << "\n";
        return 1;
    }
}
