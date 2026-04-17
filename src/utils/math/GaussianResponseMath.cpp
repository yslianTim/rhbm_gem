#include <rhbm_gem/utils/math/GaussianResponseMath.hpp>

#include <rhbm_gem/utils/domain/Constants.hpp>

#include <cmath>
#include <stdexcept>

namespace rhbm_gem::GaussianResponseMath
{

double GetGaussianResponseAtDistance(double distance, double width, int dimension)
{
    if (width <= 0.0)
    {
        throw std::runtime_error("The gaus width should be positive value.");
    }

    const auto width_square{ width * width };
    const auto coefficient{
        1.0 / std::pow(Constants::two_pi * width_square, 0.5 * static_cast<double>(dimension))
    };
    return coefficient * std::exp(-0.5 * distance * distance / width_square);
}

double GetGaussianResponseAtPoint(
    const Eigen::VectorXd & point,
    const Eigen::VectorXd & center,
    double width)
{
    if (width <= 0.0)
    {
        throw std::runtime_error("The gaus width should be positive value.");
    }

    const auto width_square{ width * width };
    const auto coefficient{ 1.0 / std::pow(Constants::two_pi * width_square, 1.5) };
    return coefficient * std::exp(-0.5 * (point - center).squaredNorm() / width_square);
}

double GetGaussianResponseAtPointWithNeighborhood(
    const Eigen::VectorXd & point,
    const Eigen::VectorXd & center,
    const std::vector<Eigen::VectorXd> & neighbor_center_list,
    double width)
{
    if (width <= 0.0)
    {
        throw std::runtime_error("The gaus width should be positive value.");
    }

    const auto width_square{ width * width };
    const auto coefficient{ 1.0 / std::pow(Constants::two_pi * width_square, 1.5) };
    auto response{
        coefficient * std::exp(-0.5 * (point - center).squaredNorm() / width_square)
    };
    for (const auto & neighbor_center : neighbor_center_list)
    {
        response +=
            coefficient * std::exp(-0.5 * (point - neighbor_center).squaredNorm() / width_square);
    }
    return response;
}

} // namespace rhbm_gem::GaussianResponseMath
