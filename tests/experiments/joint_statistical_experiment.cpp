#include "support/JointPartialSelection.hpp"
#include "support/JointComponentChecks.hpp"
#include "support/JointRuntimeJson.hpp"
#include "core/detail/joint_component/Problem.hpp"
#include "data/io/detail/JointResultJson.hpp"
#include <rhbm_gem/data/io/ModelMapFileIO.hpp>
#include <fstream>
#include <iostream>

namespace {
namespace j = boost::json;
namespace c = rhbm_gem::core;
using Path = std::filesystem::path;

j::value Read(const Path & path)
{
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Missing input: " + path.string());
    j::parse_options options;
    options.numbers = j::number_precision::precise;
    return j::parse(std::string(std::istreambuf_iterator<char>(in), {}), {}, options);
}

void Write(const Path & path, const j::value & value)
{
    const auto temporary = Path(path.string() + ".tmp");
    std::ofstream out(temporary);
    out.exceptions(std::ios::failbit | std::ios::badbit);
    out << j::serialize(value) << '\n';
    out.close();
    std::filesystem::rename(temporary, path);
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
            throw std::runtime_error("Experiment JSON changed the runtime state.");
    }
    return result;
}

j::object Census(const c::JointProblem & problem)
{
    const auto & data = c::JointProblemAccess::Get(problem);
    std::size_t memberships = 0;
    for (const auto & support : problem.Input().support) memberships += support.size();
    j::array sizes;
    for (const auto & component : data.partition.components) sizes.push_back(component.atoms.size());
    return {{"targets", problem.Input().selection_domain->target_indices.size()},
        {"contributors", problem.Input().atom_ids.size()}, {"rows", problem.Input().row_ids.size()},
        {"memberships", memberships}, {"component_sizes", sizes}};
}

std::unique_ptr<rhbm_gem::ModelObject> TwoAtoms(double distance)
{
    std::vector<std::unique_ptr<rhbm_gem::AtomObject>> atoms;
    for (int i = 0; i < 2; ++i)
    {
        auto atom = std::make_unique<rhbm_gem::AtomObject>();
        atom->SetSerialID(i + 1);
        atom->SetElement(Element::CARBON);
        atom->SetPosition(i * distance, 0, 0);
        atom->SetChainID("A");
        atom->SetComponentID("ALA");
        atom->SetAtomID(i ? "CB" : "CA");
        atoms.push_back(std::move(atom));
    }
    auto model = std::make_unique<rhbm_gem::ModelObject>(std::move(atoms));
    model->SelectAtoms([](const auto & atom) { return atom.GetSerialID() == 1; });
    return model;
}

void Run(const Path & input, const Path & output)
{
    const auto spec = Read(input);
    auto model = TwoAtoms(j::value_to<double>(spec.at("distance")));
    const auto values = j::value_to<std::vector<double>>(spec.at("values"));
    if (values.size() != 25 * 25 * 25) throw std::invalid_argument("Wrong statistical map size.");
    auto buffer = std::make_unique<double[]>(values.size());
    std::copy(values.begin(), values.end(), buffer.get());
    rhbm_gem::MapObject map({25, 25, 25}, {.5, .5, .5}, {-6, -6, -6}, std::move(buffer));
    const auto problem = c::BuildJointProblem(map, *model);
    const double distance = j::value_to<double>(spec.at("distance")) + j::value_to<double>(spec.at("shift"));
    const std::vector<second_stage_test::matched::Atom> truth{
        {{0, 0, 0}, 2, .5, .2}, {{distance, 0, 0}, 2.3, .55, .15}};
    const auto clean = j::value_to<std::vector<double>>(spec.at("clean"));
    double error = 0;
    for (std::size_t k = 0; k < values.size(); ++k)
    {
        const std::array<double, 3> position{.5 * static_cast<double>(k % 25) - 6,
            .5 * static_cast<double>((k / 25) % 25) - 6, .5 * static_cast<double>(k / 625) - 6};
        error = std::max(error, std::abs(clean.at(k) -
            second_stage_test::matched::unique_grid::Direct(position, truth, 2.5)));
    }
    if (error > 1e-12) throw std::runtime_error("Independent scalar generation disagreement.");
    j::object out{{"census", Census(problem)}, {"generator_max_error", error},
                  {"row_ids", j::value_from(problem.Input().row_ids)}};
    for (const auto & mode : spec.at("modes").as_array())
    {
        const auto fit = mode == "fixed" ? c::FitJointComponents(problem, {.6, .6}) :
                                            c::EstimateJointComponents(map, *model);
        j::object record{{"outcome", Outcome(fit)}, {"prediction", nullptr}};
        if (fit.prediction) record["prediction"] = j::value_from(*fit.prediction);
        out[std::string(mode.as_string())] = std::move(record);
    }
    Write(output, out);
}
}

int main(int argc, char ** argv)
{
    try
    {
        Eigen::setNbThreads(1);
        if (argc != 4 || std::string(argv[1]) != "stat")
            throw std::invalid_argument("Usage: joint_statistical_experiment stat INPUT OUTPUT");
        Run(argv[2], argv[3]);
        return 0;
    }
    catch (const std::exception & error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
