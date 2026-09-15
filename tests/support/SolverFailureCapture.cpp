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
#include <atomic>
#include <rhbm_gem/data/object/AtomObject.hpp>

namespace second_stage_test {
namespace {
std::mutex capture_mutex;
std::set<std::string> captured;
thread_local SolverCapturePhase phase{"proposal", 0};
thread_local boost::json::object member;
std::atomic<std::size_t> operator_count{};
boost::json::value Number(double x)
{ return std::isfinite(x) ? boost::json::value(x) : boost::json::value(nullptr); }
boost::json::array Model(const rhbm_gem::GaussianModel3D & m)
{ return {Number(m.GetAmplitude()), Number(m.GetWidth()), Number(m.GetOffset())}; }
void WriteMetadata(const std::filesystem::path & path)
{
    if (member.empty()) return;
    std::ofstream out(path.string() + ".json");
    out.exceptions(std::ios::failbit | std::ios::badbit);
    out << boost::json::serialize(member);
}
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
        WriteMetadata(path);
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

ScopedSolverCapturePhase::ScopedSolverCapturePhase(std::string name, std::size_t attempt)
    : previous{phase}
{ phase = {std::move(name), attempt ? attempt : previous.attempt}; }
ScopedSolverCapturePhase::~ScopedSolverCapturePhase() { phase = std::move(previous); }

SolverCaptureContext CaptureOperatorContext(const rhbm_gem::core::detail::SecondStageContext & context,
    const rhbm_gem::core::detail::FitState & state,
    const std::vector<rhbm_gem::core::detail::ClusterKey> & keys) noexcept
{
    try
    {
        const auto * directory{ std::getenv("RHBM_TEST_SOLVER_CAPTURE_DIR") };
        if (!directory || !*directory) return {};
        namespace j = boost::json;
        const auto id{ phase.name + "-" + std::to_string(phase.attempt) + "-" + std::to_string(++operator_count) };
        j::object data{{"schema_version", 1}, {"operator_id", id}, {"phase", phase.name},
            {"attempt", phase.attempt}, {"cluster_keys", j::value_from(keys)}};
        j::array models, atoms, background_models, background_response;
        for (std::size_t i = 0; i < state.size(); ++i)
        {
            models.push_back(Model(state[i].mdpde.GetModel()));
            const auto & a{ context.atom_list[i] }; const auto * atom{ a.atom };
            j::object entry{{"index", i}, {"alpha", a.alpha_r}};
            if (atom) entry["identity"] = j::object{{"serial_id", atom->GetSerialID()},
                {"chain_id", atom->GetChainID()}, {"sequence_id", atom->GetSequenceID()},
                {"component_id", atom->GetComponentID()}, {"atom_id", atom->GetAtomID()},
                {"alternate_indicator", atom->GetIndicator()}, {"position", j::value_from(atom->GetPosition())}};
            if (phase.name == "final" || phase.name == "recovery-current")
            {
                j::array samples;
                for (const auto & s : a.raw_sampling_entries) samples.push_back(j::object{
                    {"position", j::value_from(s.point.position)}, {"distance", s.point.distance},
                    {"selected", s.point.is_selected}, {"response", Number(s.response)}});
                entry["samples"] = std::move(samples);
            }
            atoms.push_back(std::move(entry));
        }
        if (context.frozen_background)
        {
            for (const auto & m : context.frozen_background->model_by_atom) background_models.push_back(Model(m));
            for (const auto & r : context.frozen_background->response_by_atom) background_response.push_back(j::value_from(r));
        }
        data["state"] = std::move(models); data["atoms"] = std::move(atoms);
        data["background_models"] = std::move(background_models);
        data["background_response"] = std::move(background_response);
        std::filesystem::create_directories(std::filesystem::path(directory) / "contexts");
        std::ofstream out(std::filesystem::path(directory) / "contexts" / (id + ".json"));
        out.exceptions(std::ios::failbit | std::ios::badbit); out << j::serialize(data); out.close();
        return std::make_shared<const j::object>(j::object{{"schema_version", 1},
            {"operator_id", id}, {"phase", phase.name}, {"attempt", phase.attempt},
            {"context_file", "contexts/" + id + ".json"}});
    }
    catch (...) { return {}; }
}

ScopedSolverCaptureMember::ScopedSolverCaptureMember(const SolverCaptureContext & context,
    const std::vector<std::size_t> & indices, const char * role) noexcept : previous{std::move(member)}
{
    try
    {
        member.clear();
        if (context) { member = *context; member["indices"] = boost::json::value_from(indices); member["role"] = role; }
    }
    catch (...) { member.clear(); }
}
ScopedSolverCaptureMember::~ScopedSolverCaptureMember() { member = std::move(previous); }

void CaptureShapeResponse(const std::vector<double> & response, const rhbm_gem::GaussianModel3D & model) noexcept
{
    try
    {
        if (!member.empty()) { member["adjusted_response"] = boost::json::value_from(response); member["offset_model"] = Model(model); }
    }
    catch (...) {}
}

void CaptureShapeFailure(const rhbm_gem::RHBMMemberDataset & dataset, double alpha,
    const rhbm_gem::RHBMExecutionOptions & options, const rhbm_gem::RHBMBetaEstimateResult & result) noexcept
{
    const bool contextual{ !member.empty() && (member.at("phase") == "final" || member.at("phase") == "recovery-current") };
    if (result.status == rhbm_gem::RHBMEstimationStatus::SUCCESS && !contextual) return;
    try
    {
        if (!member.empty()) member["diagnostics"] = boost::json::object{
            {"iterations", result.diagnostics.iterations},
            {"squared_beta_change", Number(result.diagnostics.squared_beta_change.value_or(0.0))},
            {"relative_variance_change", Number(result.diagnostics.relative_variance_change.value_or(0.0))}};
    }
    catch (...) {}
    std::string name{std::string("shape-") + rhbm_gem::core::detail::LocalRefitStatusText(result.status)};
    if (contextual) name += "-" + std::string(member.at("operator_id").as_string()) + "-" +
        std::to_string(boost::json::value_to<std::size_t>(member.at("indices").at(0)));
    Capture(name, [&](auto & out) {
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
