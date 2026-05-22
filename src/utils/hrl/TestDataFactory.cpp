#include <rhbm_gem/utils/hrl/TestDataFactory.hpp>

#include <rhbm_gem/utils/domain/Constants.hpp>
#include <rhbm_gem/utils/domain/SampleFilter.hpp>
#include <rhbm_gem/utils/math/EigenHelper.hpp>
#include <rhbm_gem/utils/math/EigenValidation.hpp>
#include <rhbm_gem/utils/math/NumericValidation.hpp>
#include <rhbm_gem/utils/math/SphereSampler.hpp>

#include <array>
#include <cmath>
#include <cstddef>
#include <random>
#include <stdexcept>
#include <utility>
#include <vector>

namespace rhbm_gem::test_data_factory {

namespace {
constexpr double kDefaultSamplingDistanceMin{ 0.0 };
constexpr double kDefaultSamplingDistanceMax{ 1.0 };

struct AtomNeighborContribution
{
    Eigen::VectorXd center;
    double amplitude{ 0.0 };
};

std::mt19937 BuildReplicaGenerator(int replica_index, const std::optional<std::uint32_t> & random_seed)
{
    if (random_seed.has_value())
    {
        std::seed_seq seed_sequence{ *random_seed, static_cast<std::uint32_t>(replica_index) };
        return std::mt19937(seed_sequence);
    }

    std::random_device random_device;
    std::seed_seq seed_sequence{
        random_device(),
        random_device(),
        random_device(),
        static_cast<std::uint32_t>(replica_index)
    };
    return std::mt19937(seed_sequence);
}

double ComputeGaussianResponseAtPoint3D(
    const Eigen::VectorXd & point,
    const Eigen::VectorXd & center,
    double width)
{
    const auto width_square{ width * width };
    return 1.0 / std::pow(Constants::two_pi * width_square, 1.5) *
        std::exp(-0.5 * (point - center).squaredNorm() / width_square);
}

double ComputeGaussianResponseWithAtomNeighborhood3D(
    const Eigen::VectorXd & point,
    const Eigen::VectorXd & center,
    const std::vector<AtomNeighborContribution> & neighbor_list,
    double width)
{
    auto response{ ComputeGaussianResponseAtPoint3D(point, center, width) };
    for (const auto & neighbor : neighbor_list)
    {
        response += neighbor.amplitude *
            ComputeGaussianResponseAtPoint3D(point, neighbor.center, width);
    }
    return response;
}

AtomNeighborContribution MakeAtomNeighborContribution(
    const Eigen::Vector3d & direction,
    double distance,
    double amplitude)
{
    return AtomNeighborContribution{ distance * direction, amplitude };
}

std::vector<AtomNeighborContribution> BuildAtomNeighborList(const Spot & spot)
{
    std::vector<AtomNeighborContribution> neighbor_list;
    neighbor_list.reserve(4);

    if (spot == Spot::UNK)
    {
        return neighbor_list;
    }

    if (spot == Spot::O)
    {
        neighbor_list.emplace_back(MakeAtomNeighborContribution(
            Eigen::Vector3d{ 1.0, 0.0, 0.0 }, 1.23, 6.0 / 8.0));
        return neighbor_list;
    }

    if (spot == Spot::N)
    {
        neighbor_list.emplace_back(MakeAtomNeighborContribution(
            Eigen::Vector3d{ 1.0, 0.0, 0.0 }, 1.02, 1.0 / 7.0));
        neighbor_list.emplace_back(MakeAtomNeighborContribution(
            Eigen::Vector3d{ -0.5, std::sqrt(3) / 2.0, 0.0 }, 1.48, 6.0 / 7.0));
        neighbor_list.emplace_back(MakeAtomNeighborContribution(
            Eigen::Vector3d{ -0.5, -std::sqrt(3) / 2.0, 0.0 }, 1.48, 6.0 / 7.0));
        return neighbor_list;
    }

    if (spot == Spot::C)
    {
        neighbor_list.emplace_back(MakeAtomNeighborContribution(
            Eigen::Vector3d{ 1.0, 0.0, 0.0 }, 1.23, 8.0 / 6.0));
        neighbor_list.emplace_back(MakeAtomNeighborContribution(
            Eigen::Vector3d{ -0.5, std::sqrt(3) / 2.0, 0.0 }, 1.48, 7.0 / 6.0));
        neighbor_list.emplace_back(MakeAtomNeighborContribution(
            Eigen::Vector3d{ -0.5, -std::sqrt(3) / 2.0, 0.0 }, 1.54, 1.0));
        return neighbor_list;
    }

    if (spot == Spot::CA)
    {
        neighbor_list.emplace_back(MakeAtomNeighborContribution(
            Eigen::Vector3d{ 0.0, 0.0, 1.0 }, 1.06, 1.0 / 6.0));
        neighbor_list.emplace_back(MakeAtomNeighborContribution(
            Eigen::Vector3d{ 0.0, 2.0 * std::sqrt(2) / 3.0, -1.0 / 3.0 }, 1.48, 7.0 / 6.0));
        neighbor_list.emplace_back(MakeAtomNeighborContribution(
            Eigen::Vector3d{ -std::sqrt(6) / 3.0, -std::sqrt(2) / 3.0, -1.0 / 3.0 }, 1.54, 1.0));
        neighbor_list.emplace_back(MakeAtomNeighborContribution(
            Eigen::Vector3d{ std::sqrt(6) / 3.0, -std::sqrt(2) / 3.0, -1.0 / 3.0 }, 1.54, 1.0));
        return neighbor_list;
    }
    return neighbor_list;
}

LocalPotentialSampleList GenerateRadialSamples(
    size_t sample_count,
    const GaussianModel3D & model,
    std::mt19937 & generator)
{
    numeric_validation::RequirePositive(sample_count, "sample_count");
    GaussianModel3D::RequireFinitePositiveWidthModel(model);

    std::uniform_real_distribution<> dist_distance(kDefaultSamplingDistanceMin, kDefaultSamplingDistanceMax);
    LocalPotentialSampleList sample_list;
    sample_list.reserve(sample_count);
    for (size_t i = 0; i < sample_count; i++)
    {
        const auto distance{ dist_distance(generator) };
        const auto response{ model.ResponseAtDistance(distance) };
        sample_list.emplace_back(LocalPotentialSample{
            static_cast<float>(response),
            SamplingPoint{ static_cast<float>(distance) }
        });
    }

    return sample_list;
}

LocalPotentialSampleList GenerateAtomModelSamples(
    const GaussianModel3D & model,
    const Spot & spot,
    double reject_angle_deg)
{
    GaussianModel3D::RequireFinitePositiveWidthModel(model);

    const auto neighbor_list{ BuildAtomNeighborList(spot) };
    const Eigen::VectorXd atom_center{ Eigen::VectorXd::Zero(3) };

    const auto sample_point_list{
        sphere_sampler::GenerateSamplingPointList(
            { 0.0f, 0.0f, 0.0f },
            SphereSamplingMethod::FibonacciDeterministic)
    };
    std::vector<std::array<float, 3>> reject_position_list;
    reject_position_list.reserve(neighbor_list.size());
    for (const auto & neighbor : neighbor_list)
    {
        reject_position_list.emplace_back(std::array<float, 3>{
            static_cast<float>(neighbor.center(0)),
            static_cast<float>(neighbor.center(1)),
            static_cast<float>(neighbor.center(2))
        });
    }
    const auto filtered_sample_point_list{
        sample_filter::FilterSamplingPointList(
            sample_point_list,
            { 0.0f, 0.0f, 0.0f },
            reject_position_list,
            reject_angle_deg)
    };

    LocalPotentialSampleList sample_list;
    sample_list.reserve(filtered_sample_point_list.size());
    for (const auto & sampling_point : filtered_sample_point_list)
    {
        const auto response{
            model.GetAmplitude() * ComputeGaussianResponseWithAtomNeighborhood3D(
                eigen_helper::ToEigenVector(sampling_point.position),
                atom_center,
                neighbor_list,
                model.GetWidth()
            ) +
            model.GetIntercept()
        };
        sample_list.emplace_back(LocalPotentialSample{
            static_cast<float>(response),
            sampling_point
        });
    }

    if (reject_angle_deg == 0.0) return sample_list;

    return sample_filter::FilterLocalPotentialSampleList(std::move(sample_list));
}

LocalPotentialSampleList BuildGaussianSampling(
    size_t sampling_entry_size,
    const GaussianModel3D & model,
    double outlier_ratio,
    std::mt19937 & generator)
{
    auto sampling_entries{
        GenerateRadialSamples(sampling_entry_size, model, generator)
    };
    std::uniform_real_distribution<> dist_outlier(0.0, 1.0);
    const auto outlier_response{ 0.5 * model.Intensity() };
    for (auto & sampling_entry : sampling_entries)
    {
        if (dist_outlier(generator) < outlier_ratio)
        {
            sampling_entry.response = static_cast<float>(outlier_response);
        }
    }
    return sampling_entries;
}

LocalPotentialSampleList AddLogQuadraticNoise(
    LocalPotentialSampleList sampling_entries,
    const GaussianModel3D & model,
    double error_sigma,
    std::mt19937 & generator)
{
    std::normal_distribution<> dist_error(0.0, error_sigma * model.Intensity());
    for (auto & sampling_entry : sampling_entries)
    {
        sampling_entry.response =
            static_cast<float>(static_cast<double>(sampling_entry.response)
                * std::exp(dist_error(generator)));
    }
    return sampling_entries;
}

LocalPotentialSampleList AddNoise(
    LocalPotentialSampleList sampling_entries,
    const GaussianModel3D & model,
    double error_sigma,
    LocalGaussianFitModel fit_model,
    std::mt19937 & generator)
{
    if (fit_model != LocalGaussianFitModel::LogQuadratic)
    {
        throw std::invalid_argument("Unsupported local Gaussian fit model for beta test noise.");
    }
    return AddLogQuadraticNoise(std::move(sampling_entries), model, error_sigma, generator);
}

Eigen::MatrixXd BuildRandomGausParameters(
    int member_size,
    const Eigen::VectorXd & gaus_prior,
    const Eigen::VectorXd & gaus_sigma,
    const Eigen::VectorXd & outlier_prior,
    const Eigen::VectorXd & outlier_sigma,
    double outlier_ratio,
    std::mt19937 & generator)
{
    std::uniform_real_distribution<> dist_outlier(0.0, 1.0);
    std::vector<std::normal_distribution<>> dist_gaus_list;
    std::vector<std::normal_distribution<>> dist_outlier_list;
    for (int p = 0; p < GaussianModel3D::ParameterSize(); p++)
    {
        std::normal_distribution<> dist_gaus_par(gaus_prior(p), gaus_sigma(p));
        std::normal_distribution<> dist_outlier_par(outlier_prior(p), outlier_sigma(p));
        dist_gaus_list.emplace_back(dist_gaus_par);
        dist_outlier_list.emplace_back(dist_outlier_par);
    }

    Eigen::MatrixXd gaus_par_matrix{
        Eigen::MatrixXd::Zero(GaussianModel3D::ParameterSize(), member_size)
    };
    for (int i = 0; i < member_size; i++)
    {
        Eigen::VectorXd gaus_par{ Eigen::VectorXd::Zero(GaussianModel3D::ParameterSize()) };
        Eigen::VectorXd outlier_par{ Eigen::VectorXd::Zero(GaussianModel3D::ParameterSize()) };
        const bool outlier_flag{ dist_outlier(generator) < outlier_ratio };
        for (int p = 0; p < GaussianModel3D::ParameterSize(); p++)
        {
            gaus_par(p) = dist_gaus_list.at(static_cast<size_t>(p))(generator);
            outlier_par(p) = dist_outlier_list.at(static_cast<size_t>(p))(generator);
            if (outlier_flag)
            {
                gaus_par(p) = outlier_par(p);
            }
        }
        gaus_par_matrix.col(i) = gaus_par;
    }
    return gaus_par_matrix;
}

std::vector<LocalPotentialSampleList> BuildMuMemberSamplingEntries(
    const Eigen::MatrixXd & gaus_array,
    std::mt19937 & generator)
{
    static constexpr size_t kMuMemberSamplingEntrySize{ 10 };
    std::vector<LocalPotentialSampleList> member_sampling_entries;
    member_sampling_entries.reserve(static_cast<std::size_t>(gaus_array.cols()));
    for (Eigen::Index i = 0; i < gaus_array.cols(); i++)
    {
        member_sampling_entries.emplace_back(
            GenerateRadialSamples(
                kMuMemberSamplingEntrySize,
                GaussianModel3D::FromVector(gaus_array.col(i)),
                generator));
    }
    return member_sampling_entries;
}

test_data_factory::LocalTestInput MakeLocalTestInputShell(
    const GaussianModel3D & gaus_true,
    const std::vector<double> & requested_alpha_r_list,
    bool alpha_training,
    LocalGaussianFitModel local_fit_model,
    int replica_size)
{
    test_data_factory::LocalTestInput input;
    input.gaus_true = gaus_true;
    input.requested_alpha_r_list = requested_alpha_r_list;
    input.alpha_training = alpha_training;
    input.local_fit_model = local_fit_model;
    input.replica_sampling_entries.reserve(static_cast<size_t>(replica_size));
    return input;
}
} // namespace

LocalTestInput BuildLocalTestInput(const LocalScenario & scenario)
{
    numeric_validation::RequirePositive(scenario.sampling_entry_size, "sampling_entry_size");
    numeric_validation::RequirePositive(scenario.replica_size, "replica_size");
    GaussianModel3D::RequireFinitePositiveWidthModel(scenario.gaus_true, "scenario.gaus_true");

    auto input{
        MakeLocalTestInputShell(
            scenario.gaus_true,
            scenario.requested_alpha_r_list,
            scenario.alpha_training,
            scenario.local_fit_model,
            scenario.replica_size)
    };

    for (int i = 0; i < scenario.replica_size; i++)
    {
        auto generator{ BuildReplicaGenerator(i, scenario.random_seed) };
        auto sampling_entries{
            BuildGaussianSampling(
                static_cast<size_t>(scenario.sampling_entry_size),
                scenario.gaus_true,
                scenario.outlier_ratio,
                generator)
        };
        auto noisy_sampling_entries{
            AddNoise(
                std::move(sampling_entries),
                scenario.gaus_true,
                scenario.data_error_sigma,
                scenario.local_fit_model,
                generator)
        };
        input.replica_sampling_entries.emplace_back(std::move(noisy_sampling_entries));
    }

    return input;
}

GroupTestInput BuildGroupTestInput(const GroupScenario & scenario)
{
    numeric_validation::RequirePositive(scenario.member_size, "member_size");
    numeric_validation::RequirePositive(scenario.replica_size, "replica_size");
    eigen_validation::RequireVectorSize(
        scenario.gaus_prior,
        GaussianModel3D::ParameterSize(),
        "scenario.gaus_prior");
    eigen_validation::RequireVectorSize(
        scenario.gaus_sigma,
        GaussianModel3D::ParameterSize(),
        "scenario.gaus_sigma");
    eigen_validation::RequireVectorSize(
        scenario.outlier_prior,
        GaussianModel3D::ParameterSize(),
        "scenario.outlier_prior");
    eigen_validation::RequireVectorSize(
        scenario.outlier_sigma,
        GaussianModel3D::ParameterSize(),
        "scenario.outlier_sigma");

    GroupTestInput input;
    input.gaus_true = GaussianModel3D::FromVector(scenario.gaus_prior);
    input.requested_alpha_g_list = scenario.requested_alpha_g_list;
    input.alpha_training = scenario.alpha_training;
    input.replica_member_sampling_entries.reserve(static_cast<size_t>(scenario.replica_size));

    for (int i = 0; i < scenario.replica_size; i++)
    {
        auto generator{ BuildReplicaGenerator(i, scenario.random_seed) };
        const auto random_gaus_array{
            BuildRandomGausParameters(
                scenario.member_size,
                scenario.gaus_prior,
                scenario.gaus_sigma,
                scenario.outlier_prior,
                scenario.outlier_sigma,
                scenario.outlier_ratio,
                generator
            )
        };
        input.replica_member_sampling_entries.emplace_back(
            BuildMuMemberSamplingEntries(random_gaus_array, generator));
    }

    return input;
}

AtomModelTestInput BuildAtomModelTestInput(const AtomModelScenario & scenario)
{
    numeric_validation::RequirePositive(scenario.replica_size, "replica_size");
    GaussianModel3D::RequireFinitePositiveWidthModel(scenario.gaus_true, "scenario.gaus_true");

    AtomModelTestInput input;
    input.no_cut_input = MakeLocalTestInputShell(
        scenario.gaus_true,
        scenario.requested_alpha_r_list,
        scenario.alpha_training,
        scenario.local_fit_model,
        scenario.replica_size);
    input.cut_input = MakeLocalTestInputShell(
        scenario.gaus_true,
        scenario.requested_alpha_r_list,
        scenario.alpha_training,
        scenario.local_fit_model,
        scenario.replica_size);
    if (scenario.include_sampling_summary)
    {
        input.sampling_summaries.reserve(1);
        input.sampling_summaries.emplace_back(
            GenerateAtomModelSamples(
                scenario.gaus_true,
                scenario.spot,
                scenario.rejected_angle
            )
        );
    }

    for (int i = 0; i < scenario.replica_size; i++)
    {
        auto generator{ BuildReplicaGenerator(i, scenario.random_seed) };
        auto no_cut_sampling_entries{
            GenerateAtomModelSamples(scenario.gaus_true, scenario.spot, 0.0)
        };
        auto cut_sampling_entries{
            GenerateAtomModelSamples(scenario.gaus_true, scenario.spot, scenario.rejected_angle)
        };
        no_cut_sampling_entries = AddNoise(
            std::move(no_cut_sampling_entries),
            scenario.gaus_true,
            scenario.data_error_sigma,
            scenario.local_fit_model,
            generator);
        cut_sampling_entries = AddNoise(
            std::move(cut_sampling_entries),
            scenario.gaus_true,
            scenario.data_error_sigma,
            scenario.local_fit_model,
            generator);
        input.no_cut_input.replica_sampling_entries.emplace_back(std::move(no_cut_sampling_entries));
        input.cut_input.replica_sampling_entries.emplace_back(std::move(cut_sampling_entries));
    }

    return input;
}

} // namespace rhbm_gem::test_data_factory
