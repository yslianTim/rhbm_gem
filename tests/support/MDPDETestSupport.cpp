#include "support/MDPDETestSupport.hpp"
#include <cmath>
#include <fstream>
#include <limits>
#include <locale>
#include <stdexcept>

namespace second_stage_test {
namespace {
namespace json = boost::json;
json::value Number(double value) { return std::isfinite(value) ? json::value(value) : json::value(nullptr); }
json::array Vector(const Eigen::VectorXd & values)
{
    json::array result;
    for (double value : values) result.push_back(Number(value));
    return result;
}
Eigen::MatrixXd ReadMatrix(std::istream & input)
{
    int rows{}, columns{};
    if (!(input >> rows >> columns) || rows < 0 || rows > 1000000 || columns < 0 || columns > 1000)
        throw std::runtime_error("Invalid shape matrix dimensions.");
    Eigen::MatrixXd matrix(rows, columns);
    for (int row = 0; row < rows; ++row) for (int column = 0; column < columns; ++column)
    {
        std::string token;
        if (!(input >> token)) throw std::runtime_error("Truncated shape fixture.");
        matrix(row, column) = std::stod(token);
    }
    return matrix;
}
} // namespace

ShapeFixture ReadShapeFixture(const std::string & path)
{
    std::ifstream input(path); input.imbue(std::locale::classic());
    std::string kind; int version{}, status{};
    ShapeFixture fixture;
    if (!(input >> kind >> version) || kind != "shape" || version != 1)
        throw std::runtime_error("Expected shape version 1 fixture.");
    input >> fixture.alpha >> fixture.options.thread_size >> fixture.options.max_iterations >>
        fixture.options.tolerance >> fixture.options.data_weight_min;
    fixture.dataset.X = ReadMatrix(input); fixture.dataset.y = ReadMatrix(input);
    std::string variance;
    if (!(input >> status >> variance)) throw std::runtime_error("Missing expected shape result.");
    fixture.expected.status = static_cast<rhbm_gem::RHBMEstimationStatus>(status);
    fixture.expected.sigma_square = std::stod(variance);
    fixture.expected.beta_ols = ReadMatrix(input); fixture.expected.beta_mdpde = ReadMatrix(input);
    std::string extra;
    if (input >> extra) throw std::runtime_error("Unexpected shape fixture content.");
    return fixture;
}

boost::json::object EquationJSON(const MDPDEEquationEvidence & evidence)
{
    const double norm{ evidence.scaled.size() && evidence.scaled.allFinite()
        ? evidence.scaled.lpNorm<Eigen::Infinity>() : std::numeric_limits<double>::infinity() };
    return {{"valid", evidence.valid}, {"reason", evidence.reason}, {"raw_residual", Vector(evidence.raw)},
        {"scaled_residual", Vector(evidence.scaled)}, {"residual_inf", Number(norm)},
        {"denominator", Number(evidence.denominator)}, {"rank", evidence.rank},
        {"singular_values", Vector(evidence.singular_values)}, {"condition", Number(evidence.condition)},
        {"floor_count", evidence.floor_count}, {"equation_pass", evidence.valid && norm <= 1.0e-8},
        {"reference_pass", evidence.valid && norm <= 1.0e-10}};
}
} // namespace second_stage_test
