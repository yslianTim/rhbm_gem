#include <rhbm_gem/utils/algorithm/AndersonAcceleration.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace rhbm_gem::algorithm {
namespace {

bool IsFiniteVector(const Eigen::VectorXd & value)
{
    return value.allFinite();
}

void ValidateOptions(const AndersonAccelerationOptions & options)
{
    if (options.history_depth == 0)
    {
        throw std::invalid_argument("Anderson acceleration history depth must be positive.");
    }
    if (!std::isfinite(options.scale_floor) || options.scale_floor <= 0.0)
    {
        throw std::invalid_argument("Anderson acceleration scale floor must be positive and finite.");
    }
    if (!std::isfinite(options.coefficient_l1_limit) || options.coefficient_l1_limit < 1.0)
    {
        throw std::invalid_argument(
            "Anderson acceleration coefficient L1 limit must be finite and at least one.");
    }
    if (!std::isfinite(options.regularization) || options.regularization < 0.0)
    {
        throw std::invalid_argument("Anderson acceleration regularization must be finite and non-negative.");
    }
    if (!std::isfinite(options.coefficient_abs_limit) || options.coefficient_abs_limit < 1.0)
    {
        throw std::invalid_argument(
            "Anderson acceleration coefficient absolute limit must be finite and at least one.");
    }
}

bool HasFiniteConsistentShape(
    const std::vector<Eigen::VectorXd> & input_list,
    const std::vector<Eigen::VectorXd> & output_list)
{
    if (input_list.empty() || input_list.size() != output_list.size())
    {
        return false;
    }
    for (std::size_t i = 0; i < input_list.size(); i++)
    {
        if (input_list.at(i).size() == 0 ||
            input_list.at(i).size() != output_list.at(i).size() ||
            !IsFiniteVector(input_list.at(i)) ||
            !IsFiniteVector(output_list.at(i)))
        {
            return false;
        }
    }
    return true;
}

bool IsValidSample(
    const AndersonAccelerationHistory::Sample & sample,
    const std::vector<std::size_t> & active_index_list)
{
    if (!HasFiniteConsistentShape(sample.input_list, sample.output_list))
    {
        return false;
    }
    for (const auto active_index : active_index_list)
    {
        if (active_index >= sample.input_list.size())
        {
            return false;
        }
    }
    return true;
}

Eigen::MatrixXd BuildScaledResidualMatrix(
    const std::vector<AndersonAccelerationHistory::Sample> & sample_list,
    const std::vector<std::size_t> & active_index_list,
    double scale_floor)
{
    const auto active_size{ active_index_list.size() };
    const auto parameter_size{
        static_cast<std::size_t>(sample_list.front().input_list.at(active_index_list.front()).size())
    };
    Eigen::MatrixXd residual_matrix{
        Eigen::MatrixXd::Zero(
            static_cast<Eigen::Index>(active_size * parameter_size),
            static_cast<Eigen::Index>(sample_list.size()))
    };

    for (std::size_t sample_index = 0; sample_index < sample_list.size(); sample_index++)
    {
        const auto & sample{ sample_list.at(sample_index) };
        for (std::size_t active_position = 0; active_position < active_size; active_position++)
        {
            const auto state_index{ active_index_list.at(active_position) };
            const auto & input{ sample.input_list.at(state_index) };
            const auto & output{ sample.output_list.at(state_index) };
            if (static_cast<std::size_t>(input.size()) != parameter_size ||
                static_cast<std::size_t>(output.size()) != parameter_size)
            {
                return Eigen::MatrixXd{};
            }
            for (std::size_t parameter_index = 0; parameter_index < parameter_size; parameter_index++)
            {
                const auto scale{
                    std::max({
                        std::abs(input(static_cast<Eigen::Index>(parameter_index))),
                        std::abs(output(static_cast<Eigen::Index>(parameter_index))),
                        scale_floor
                    })
                };
                residual_matrix(
                    static_cast<Eigen::Index>(active_position * parameter_size + parameter_index),
                    static_cast<Eigen::Index>(sample_index)) =
                    (output(static_cast<Eigen::Index>(parameter_index)) -
                        input(static_cast<Eigen::Index>(parameter_index))) / scale;
            }
        }
    }
    return residual_matrix;
}

std::optional<Eigen::VectorXd> SolveConstrainedCoefficients(
    const Eigen::MatrixXd & residual_matrix,
    double regularization,
    double coefficient_l1_limit,
    double coefficient_abs_limit)
{
    if (residual_matrix.rows() == 0 || residual_matrix.cols() < 2 || !residual_matrix.allFinite())
    {
        return std::nullopt;
    }

    const auto sample_size{ residual_matrix.cols() };
    Eigen::MatrixXd system{
        Eigen::MatrixXd::Zero(sample_size + 1, sample_size + 1)
    };
    system.topLeftCorner(sample_size, sample_size) =
        residual_matrix.transpose() * residual_matrix;
    if (regularization > 0.0)
    {
        system.topLeftCorner(sample_size, sample_size).diagonal().array() += regularization;
    }
    system.block(0, sample_size, sample_size, 1).setOnes();
    system.block(sample_size, 0, 1, sample_size).setOnes();

    Eigen::VectorXd rhs{ Eigen::VectorXd::Zero(sample_size + 1) };
    rhs(sample_size) = 1.0;

    const Eigen::VectorXd solution{ system.fullPivLu().solve(rhs) };
    if (!solution.allFinite())
    {
        return std::nullopt;
    }

    const auto coefficients{ solution.head(sample_size).eval() };
    const auto residual{ (system * solution - rhs).norm() };
    if (!std::isfinite(residual) || residual > 1.0e-8 * std::max(rhs.norm(), 1.0))
    {
        return std::nullopt;
    }
    if (std::abs(coefficients.sum() - 1.0) > 1.0e-8)
    {
        return std::nullopt;
    }
    if (coefficients.array().abs().sum() > coefficient_l1_limit)
    {
        return std::nullopt;
    }
    if (coefficients.array().abs().maxCoeff() > coefficient_abs_limit)
    {
        return std::nullopt;
    }
    return coefficients;
}

std::optional<std::vector<Eigen::VectorXd>> CombineOutputs(
    const std::vector<AndersonAccelerationHistory::Sample> & sample_list,
    const Eigen::VectorXd & coefficients)
{
    const auto state_size{ sample_list.front().output_list.size() };
    std::vector<Eigen::VectorXd> candidate_list;
    candidate_list.reserve(state_size);
    for (std::size_t state_index = 0; state_index < state_size; state_index++)
    {
        Eigen::VectorXd candidate{
            Eigen::VectorXd::Zero(sample_list.front().output_list.at(state_index).size())
        };
        for (std::size_t sample_index = 0; sample_index < sample_list.size(); sample_index++)
        {
            const auto & output{ sample_list.at(sample_index).output_list.at(state_index) };
            if (output.size() != candidate.size())
            {
                return std::nullopt;
            }
            candidate += coefficients(static_cast<Eigen::Index>(sample_index)) * output;
        }
        if (!candidate.allFinite())
        {
            return std::nullopt;
        }
        candidate_list.emplace_back(std::move(candidate));
    }
    return candidate_list;
}

} // namespace

AndersonAccelerationHistory::AndersonAccelerationHistory(AndersonAccelerationOptions options)
    : m_options{ options }
{
    ValidateOptions(m_options);
}

void AndersonAccelerationHistory::Clear()
{
    m_active_index_list.clear();
    m_sample_list.clear();
}

bool AndersonAccelerationHistory::HasCompatibleActiveIndexList(
    const std::vector<std::size_t> & active_index_list) const
{
    return m_sample_list.empty() || m_active_index_list == active_index_list;
}

void AndersonAccelerationHistory::Commit(
    const std::vector<std::size_t> & active_index_list,
    const std::vector<Eigen::VectorXd> & input_list,
    const std::vector<Eigen::VectorXd> & output_list)
{
    Sample sample{ input_list, output_list };
    if (active_index_list.empty() || !IsValidSample(sample, active_index_list))
    {
        Clear();
        return;
    }
    if (!HasCompatibleActiveIndexList(active_index_list))
    {
        Clear();
    }
    m_active_index_list = active_index_list;
    m_sample_list.emplace_back(std::move(sample));
    if (m_sample_list.size() > m_options.history_depth)
    {
        m_sample_list.erase(m_sample_list.begin());
    }
}

std::optional<std::vector<Eigen::VectorXd>> AndersonAccelerationHistory::BuildCandidate(
    const std::vector<std::size_t> & active_index_list,
    const std::vector<Eigen::VectorXd> & input_list,
    const std::vector<Eigen::VectorXd> & output_list) const
{
    Sample current_sample{ input_list, output_list };
    if (active_index_list.empty() ||
        !HasCompatibleActiveIndexList(active_index_list) ||
        !IsValidSample(current_sample, active_index_list))
    {
        return std::nullopt;
    }

    std::vector<Sample> candidate_sample_list{ m_sample_list };
    candidate_sample_list.emplace_back(std::move(current_sample));
    if (candidate_sample_list.size() < 2)
    {
        return std::nullopt;
    }
    if (candidate_sample_list.size() > m_options.history_depth)
    {
        candidate_sample_list.erase(candidate_sample_list.begin());
    }
    for (const auto & sample : candidate_sample_list)
    {
        if (!IsValidSample(sample, active_index_list))
        {
            return std::nullopt;
        }
    }

    const auto residual_matrix{
        BuildScaledResidualMatrix(
            candidate_sample_list,
            active_index_list,
            m_options.scale_floor)
    };
    const auto coefficients{
        SolveConstrainedCoefficients(
            residual_matrix,
            m_options.regularization,
            m_options.coefficient_l1_limit,
            m_options.coefficient_abs_limit)
    };
    if (!coefficients.has_value())
    {
        return std::nullopt;
    }
    return CombineOutputs(candidate_sample_list, *coefficients);
}

} // namespace rhbm_gem::algorithm
