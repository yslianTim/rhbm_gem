#include "support/SolverFailureCapture.hpp"
#include "core/detail/second_stage/JointFitting.hpp"
#include <rhbm_gem/utils/hrl/RHBMHelper.hpp>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <mutex>
#include <set>
#include <stdexcept>

namespace second_stage_test {
namespace {
std::mutex capture_mutex;
std::set<std::string> captured;
template<class Writer> void Capture(const std::string & name, Writer writer) noexcept
{
    try
    {
        const auto * directory{ std::getenv("RHBM_TEST_SOLVER_CAPTURE_DIR") };
        if (!directory || !*directory) return;
        std::lock_guard lock(capture_mutex);
        const auto path{ std::filesystem::path(directory) / (name + ".txt") };
        if (captured.contains(path.string())) return;
        std::filesystem::create_directories(path.parent_path());
        std::ofstream out(path);
        out.exceptions(std::ios::failbit | std::ios::badbit);
        out.imbue(std::locale::classic()); out << std::setprecision(17);
        writer(out);
        out.close();
        captured.insert(path.string());
    }
    catch (...) { /* Capture failure cannot alter numerical execution. The harness checks artifacts. */ }
}
template<class Matrix> void WriteMatrix(std::ostream & out, const Matrix & matrix)
{
    out << matrix.rows() << ' ' << matrix.cols() << '\n';
    for (Eigen::Index row = 0; row < matrix.rows(); ++row)
    { for (Eigen::Index col = 0; col < matrix.cols(); ++col) out << matrix(row, col) << ' '; out << '\n'; }
}
double ReadNumber(std::istream & in)
{
    std::string text;
    if (!(in >> text)) throw std::runtime_error("Incomplete solver fixture.");
    return std::stod(text);
}
Eigen::MatrixXd ReadMatrix(std::istream & in)
{
    Eigen::Index rows, cols;
    if (!(in >> rows >> cols) || rows < 0 || cols < 0 || rows > 1000000 || cols > 1000)
        throw std::runtime_error("Invalid solver fixture dimensions.");
    Eigen::MatrixXd result(rows, cols);
    for (Eigen::Index row = 0; row < rows; ++row)
        for (Eigen::Index col = 0; col < cols; ++col) result(row, col) = ReadNumber(in);
    return result;
}
bool Equal(const Eigen::MatrixXd & left, const Eigen::MatrixXd & right)
{
    return left.rows() == right.rows() && left.cols() == right.cols() && (left.array() == right.array()).all();
}
}

void CaptureShapeFailure(const rhbm_gem::RHBMMemberDataset & dataset, double alpha,
    const rhbm_gem::RHBMExecutionOptions & options, const rhbm_gem::RHBMBetaEstimateResult & result) noexcept
{
    if (result.status == rhbm_gem::RHBMEstimationStatus::SUCCESS) return;
    Capture(std::string("shape-") + rhbm_gem::core::detail::LocalRefitStatusText(result.status), [&](auto & out) {
        out << "shape 1\n" << alpha << ' ' << options.thread_size << ' ' << options.max_iterations
            << ' ' << options.tolerance << ' ' << options.data_weight_min << '\n';
        WriteMatrix(out, dataset.X); WriteMatrix(out, dataset.y);
        out << static_cast<int>(result.status) << ' ' << result.sigma_square << '\n';
        WriteMatrix(out, result.beta_ols); WriteMatrix(out, result.beta_mdpde);
    });
}
void CaptureJointFailure(const rhbm_gem::algorithm::WeightedRidgeSystem & system,
    const rhbm_gem::core::detail::JointOffsetSolveResult & result) noexcept
{
    using namespace rhbm_gem::core::detail;
    if (result.status == JointOffsetSolveStatus::Converged) return;
    Capture(std::string("offset-") + JointOffsetSolveStatusText(result.status), [&](auto & out) {
        out << "offset 1\n";
        out << system.design_matrix.rows() << ' ' << system.design_matrix.cols() << ' ' << system.design_matrix.nonZeros() << '\n';
        for (Eigen::Index col = 0; col < system.design_matrix.outerSize(); ++col)
            for (Eigen::SparseMatrix<double>::InnerIterator entry(system.design_matrix, col); entry; ++entry)
                out << entry.row() << ' ' << entry.col() << ' ' << entry.value() << '\n';
        WriteMatrix(out, system.response);
        WriteMatrix(out, system.previous_parameter); WriteMatrix(out, system.ridge_diagonal);
        out << static_cast<int>(result.status) << '\n'; WriteMatrix(out, result.offset);
    });
}
bool ReplaySolverFailure(const std::string & path)
{
    using namespace rhbm_gem;
    std::ifstream in(path); in.imbue(std::locale::classic());
    std::string kind; int version, expected_status;
    if (!(in >> kind >> version) || version != 1) throw std::runtime_error("Unsupported solver fixture.");
    if (kind == "shape")
    {
        double alpha; RHBMExecutionOptions options;
        in >> alpha >> options.thread_size >> options.max_iterations >> options.tolerance >> options.data_weight_min;
        RHBMMemberDataset dataset; dataset.X = ReadMatrix(in); dataset.y = ReadMatrix(in);
        in >> expected_status; const auto variance{ ReadNumber(in) };
        const auto ols{ ReadMatrix(in) }, mdpde{ ReadMatrix(in) };
        const auto actual{ rhbm_helper::EstimateBetaMDPDE(alpha, dataset, options) };
        return static_cast<int>(actual.status) == expected_status && actual.sigma_square == variance &&
            Equal(actual.beta_ols, ols) && Equal(actual.beta_mdpde, mdpde);
    }
    if (kind == "offset")
    {
        algorithm::WeightedRidgeSystem system;
        Eigen::Index rows, cols, count;
        in >> rows >> cols >> count;
        if (!in || rows <= 0 || cols <= 0 || count < 0 || rows > 1000000 || cols > 1000 || count > rows * cols)
            throw std::runtime_error("Invalid sparse fixture.");
        std::vector<Eigen::Triplet<double>> entries;
        for (Eigen::Index i = 0; i < count; ++i)
        {
            Eigen::Index row, col; in >> row >> col;
            if (!in || row < 0 || row >= rows || col < 0 || col >= cols) throw std::runtime_error("Invalid sparse index.");
            entries.emplace_back(row, col, ReadNumber(in));
        }
        system.design_matrix.resize(rows, cols); system.design_matrix.setFromTriplets(entries.begin(), entries.end());
        system.response = ReadMatrix(in);
        system.previous_parameter = ReadMatrix(in); system.ridge_diagonal = ReadMatrix(in);
        in >> expected_status; const auto expected{ ReadMatrix(in) };
        algorithm::WeightedRidgeSolver solver;
        const auto actual{ core::detail::SolveJointOffsetSystem(system, solver) };
        return static_cast<int>(actual.status) == expected_status && Equal(actual.offset, expected);
    }
    throw std::runtime_error("Unknown solver fixture kind.");
}
}
