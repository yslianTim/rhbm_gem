#include <rhbm_gem/core/TestDataFactory.hpp>

#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <rhbm_gem/utils/domain/SampleFilter.hpp>
#include <rhbm_gem/utils/math/EigenHelper.hpp>
#include <rhbm_gem/utils/math/NumericValidation.hpp>
#include <rhbm_gem/utils/math/SphereSampler.hpp>

#include <array>
#include <cmath>
#include <cstddef>
#include <memory>
#include <random>
#include <utility>
#include <vector>

namespace rhbm_gem::core {

namespace {
constexpr double kDefaultSamplingDistanceMin{ 0.0 };
constexpr double kDefaultSamplingDistanceMax{ 1.0 };

struct AtomNeighborContribution
{
    Spot spot;
    Eigen::VectorXd center;
    double amplitude{ 0.0 };
};

struct AtomicModelContribution
{
    Spot spot;
    Eigen::VectorXd center;
    double amplitude_scale{ 1.0 };
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

double EvaluateZeroInterceptModelAtPoint3D(
    const GaussianModel3D & zero_intercept_model,
    const Eigen::VectorXd & point,
    const Eigen::VectorXd & center)
{
    const auto distance{ (point - center).norm() };
    return zero_intercept_model.ResponseAtDistance(distance);
}

double EvaluateModelOffsetAtDistance(const GaussianModel3D & model, double distance)
{
    return model.ResponseAtDistance(distance) -
        model.WithIntercept(0.0).ResponseAtDistance(distance);
}

double EvaluateAtomNeighborhoodModelResponse(
    const GaussianModel3D & model,
    const Eigen::VectorXd & point,
    const Eigen::VectorXd & center,
    const std::vector<AtomNeighborContribution> & neighbor_list,
    double sample_distance)
{
    const auto zero_intercept_model{ model.WithIntercept(0.0) };
    auto response{ EvaluateZeroInterceptModelAtPoint3D(zero_intercept_model, point, center) };
    for (const auto & neighbor : neighbor_list)
    {
        response += neighbor.amplitude *
            EvaluateZeroInterceptModelAtPoint3D(zero_intercept_model, point, neighbor.center);
    }
    return response + EvaluateModelOffsetAtDistance(model, sample_distance);
}

double EvaluateAtomicFieldModelResponse(
    const GaussianModel3D & model,
    const Eigen::VectorXd & point,
    const std::vector<AtomicModelContribution> & atom_field,
    double sample_distance)
{
    const auto zero_intercept_model{ model.WithIntercept(0.0) };
    double response{ 0.0 };
    for (const auto & atom : atom_field)
    {
        response += atom.amplitude_scale *
            EvaluateZeroInterceptModelAtPoint3D(zero_intercept_model, point, atom.center);
    }
    return response + EvaluateModelOffsetAtDistance(model, sample_distance);
}

AtomNeighborContribution MakeAtomNeighborContribution(
    Spot spot,
    const Eigen::Vector3d & direction,
    double distance,
    double amplitude)
{
    return AtomNeighborContribution{ spot, distance * direction, amplitude };
}

std::vector<AtomNeighborContribution> BuildAtomNeighborList(const Spot & spot)
{
    switch (spot)
    {
    case Spot::UNK:
        return {};
    case Spot::O:
        return {
            MakeAtomNeighborContribution(Spot::C,
                Eigen::Vector3d{ 1.0, 0.0, 0.0 }, 1.23, 6.0 / 8.0)
        };
    case Spot::N:
        return {
            MakeAtomNeighborContribution(Spot::H,
                Eigen::Vector3d{ 1.0, 0.0, 0.0 }, 1.02, 1.0 / 7.0),
            MakeAtomNeighborContribution(Spot::C,
                Eigen::Vector3d{ -0.5, std::sqrt(3) / 2.0, 0.0 }, 1.48, 6.0 / 7.0),
            MakeAtomNeighborContribution(Spot::CA,
                Eigen::Vector3d{ -0.5, -std::sqrt(3) / 2.0, 0.0 }, 1.48, 6.0 / 7.0)
        };
    case Spot::C:
        return {
            MakeAtomNeighborContribution(Spot::O,
                Eigen::Vector3d{ 1.0, 0.0, 0.0 }, 1.23, 8.0 / 6.0),
            MakeAtomNeighborContribution(Spot::N,
                Eigen::Vector3d{ -0.5, std::sqrt(3) / 2.0, 0.0 }, 1.48, 7.0 / 6.0),
            MakeAtomNeighborContribution(Spot::CA,
                Eigen::Vector3d{ -0.5, -std::sqrt(3) / 2.0, 0.0 }, 1.54, 1.0)
        };
    case Spot::CA:
        return {
            MakeAtomNeighborContribution(Spot::H,
                Eigen::Vector3d{ 0.0, 0.0, 1.0 }, 1.06, 1.0 / 6.0),
            MakeAtomNeighborContribution(Spot::N,
                Eigen::Vector3d{ 0.0, 2.0 * std::sqrt(2) / 3.0, -1.0 / 3.0 }, 1.48, 7.0 / 6.0),
            MakeAtomNeighborContribution(Spot::C,
                Eigen::Vector3d{ -std::sqrt(6) / 3.0, -std::sqrt(2) / 3.0, -1.0 / 3.0 }, 1.54, 1.0),
            MakeAtomNeighborContribution(Spot::CB,
                Eigen::Vector3d{ std::sqrt(6) / 3.0, -std::sqrt(2) / 3.0, -1.0 / 3.0 }, 1.54, 1.0)
        };
    default:
        return {};
    }
}

std::array<float, 3> ToArray3f(const Eigen::VectorXd & vector)
{
    return {
        static_cast<float>(vector(0)),
        static_cast<float>(vector(1)),
        static_cast<float>(vector(2))
    };
}

std::vector<AtomicModelContribution> BuildAtomicModelContributionList(const Spot & spot)
{
    std::vector<AtomicModelContribution> atom_field;
    const auto neighbor_list{ BuildAtomNeighborList(spot) };
    atom_field.reserve(neighbor_list.size() + 1);
    atom_field.emplace_back(AtomicModelContribution{ spot, Eigen::VectorXd::Zero(3), 1.0 });
    for (const auto & neighbor : neighbor_list)
    {
        atom_field.emplace_back(AtomicModelContribution{ neighbor.spot, neighbor.center, neighbor.amplitude });
    }
    return atom_field;
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

LocalPotentialSampleList GenerateAtomModelSamples(const GaussianModel3D & model, const Spot & spot)
{
    GaussianModel3D::RequireFinitePositiveWidthModel(model);
    const auto neighbor_list{ BuildAtomNeighborList(spot) };
    const Eigen::VectorXd atom_center{ Eigen::VectorXd::Zero(3) };

    auto sample_point_list{
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
    sample_filter::FilterSamplingPointList(
        sample_point_list,
        { 0.0f, 0.0f, 0.0f },
        reject_position_list);

    LocalPotentialSampleList sample_list;
    sample_list.reserve(sample_point_list.size());
    for (const auto & sampling_point : sample_point_list)
    {
        const auto point{ eigen_helper::ToEigenVector(sampling_point.position) };
        const auto response{
            EvaluateAtomNeighborhoodModelResponse(
                model,
                point,
                atom_center,
                neighbor_list,
                sampling_point.distance
            )
        };
        sample_list.emplace_back(LocalPotentialSample{
            static_cast<float>(response),
            sampling_point
        });
    }

    return sample_list;
}

LocalPotentialSampleList GenerateAtomicModelSamples(
    const GaussianModel3D & model,
    const AtomicModelContribution & local_atom,
    const std::vector<AtomicModelContribution> & atom_field)
{
    GaussianModel3D::RequireFinitePositiveWidthModel(model);

    const auto local_position{ ToArray3f(local_atom.center) };
    auto sample_point_list{
        sphere_sampler::GenerateSamplingPointList(
            local_position,
            SphereSamplingMethod::FibonacciDeterministic)
    };
    std::vector<std::array<float, 3>> reject_position_list;
    reject_position_list.reserve(atom_field.size());
    for (const auto & atom : atom_field)
    {
        reject_position_list.emplace_back(ToArray3f(atom.center));
    }
    sample_filter::FilterSamplingPointList(
        sample_point_list,
        local_position,
        reject_position_list);

    LocalPotentialSampleList sample_list;
    sample_list.reserve(sample_point_list.size());
    for (const auto & sampling_point : sample_point_list)
    {
        const auto point{ eigen_helper::ToEigenVector(sampling_point.position) };
        const auto response{
            EvaluateAtomicFieldModelResponse(
                model,
                point,
                atom_field,
                sampling_point.distance
            )
        };
        sample_list.emplace_back(LocalPotentialSample{
            static_cast<float>(response),
            sampling_point
        });
    }
    return sample_list;
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

LocalPotentialSampleList ApplyLogQuadraticNoise(
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

Eigen::MatrixXd BuildRandomGausParameters(
    int member_size,
    const GaussianParameterDistribution & inlier_distribution,
    const GaussianParameterDistribution & outlier_distribution,
    double outlier_ratio,
    std::mt19937 & generator)
{
    std::uniform_real_distribution<> dist_outlier(0.0, 1.0);
    Eigen::MatrixXd gaus_par_matrix{
        Eigen::MatrixXd::Zero(GaussianModel3D::ParameterSize(), member_size)
    };
    for (int i = 0; i < member_size; i++)
    {
        const auto & distribution{
            dist_outlier(generator) < outlier_ratio ? outlier_distribution : inlier_distribution
        };
        Eigen::VectorXd gaus_par{
            Eigen::VectorXd::Zero(GaussianModel3D::ParameterSize())
        };
        for (int p = 0; p < GaussianModel3D::ParameterSize(); p++)
        {
            std::normal_distribution<> dist_par(
                distribution.mean.GetModelParameter(p),
                distribution.sigma.GetModelParameter(p));
            gaus_par(p) = dist_par(generator);
        }
        gaus_par_matrix.col(i) = gaus_par;
    }
    return gaus_par_matrix;
}

std::vector<LocalPotentialSampleList> BuildMuMemberSamplingEntries(
    const Eigen::MatrixXd & gaus_array,
    size_t sampling_entry_size,
    std::mt19937 & generator)
{
    std::vector<LocalPotentialSampleList> member_sampling_entries;
    member_sampling_entries.reserve(static_cast<std::size_t>(gaus_array.cols()));
    for (Eigen::Index i = 0; i < gaus_array.cols(); i++)
    {
        member_sampling_entries.emplace_back(
            GenerateRadialSamples(
                sampling_entry_size,
                GaussianModel3D::FromVector(gaus_array.col(i)),
                generator));
    }
    return member_sampling_entries;
}

std::unique_ptr<AtomObject> MakeAtomicModelAtom(int serial_id, const AtomicModelContribution & atom)
{
    auto atom_object{ std::make_unique<AtomObject>() };
    atom_object->SetSerialID(serial_id);
    atom_object->SetComponentKey(1);
    atom_object->SetAtomKey(static_cast<AtomKey>(atom.spot));
    atom_object->SetSpot(atom.spot);
    atom_object->SetPosition(ToArray3f(atom.center));
    return atom_object;
}

std::unique_ptr<ModelObject> BuildAtomicModelObject(
    const GaussianModel3D & model,
    const std::vector<AtomicModelContribution> & atom_field,
    double error_sigma,
    std::mt19937 & generator)
{
    std::vector<std::unique_ptr<AtomObject>> atom_list;
    atom_list.reserve(atom_field.size());
    atom_list.emplace_back(MakeAtomicModelAtom(1, atom_field.front()));
    for (std::size_t i = 1; i < atom_field.size(); i++)
    {
        atom_list.emplace_back(MakeAtomicModelAtom(static_cast<int>(i + 1), atom_field.at(i)));
    }

    auto model_object{ std::make_unique<ModelObject>(std::move(atom_list)) };
    model_object->SelectAllAtoms();
    auto analysis{ model_object->EditAnalysis() };
    analysis.RebuildAtomGroupsFromSelection();
    for (std::size_t i = 0; i < atom_field.size(); i++)
    {
        auto sampling_entries{ GenerateAtomicModelSamples(model, atom_field.at(i), atom_field) };
        sampling_entries = ApplyLogQuadraticNoise(
            std::move(sampling_entries),
            model,
            error_sigma,
            generator);
        analysis.EnsureAtomLocalPotential(*model_object->GetSelectedAtoms().at(i))
            .SetSamplingEntries(std::move(sampling_entries));
    }
    return model_object;
}
} // namespace

LocalTestData BuildLocalTestData(const LocalScenario & scenario)
{
    numeric_validation::RequirePositive(scenario.sampling_entry_size, "sampling_entry_size");
    numeric_validation::RequirePositive(scenario.replica_size, "replica_size");
    GaussianModel3D::RequireFinitePositiveWidthModel(scenario.gaus_true, "scenario.gaus_true");

    LocalTestData input;
    input.gaus_true = scenario.gaus_true;
    input.replica_sampling_entries.reserve(static_cast<size_t>(scenario.replica_size));

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
            ApplyLogQuadraticNoise(
                std::move(sampling_entries),
                scenario.gaus_true,
                scenario.data_error_sigma,
                generator)
        };
        input.replica_sampling_entries.emplace_back(std::move(noisy_sampling_entries));
    }

    return input;
}

GroupTestData BuildGroupTestData(const GroupScenario & scenario)
{
    numeric_validation::RequirePositive(scenario.member_size, "member_size");
    numeric_validation::RequirePositive(scenario.sampling_entry_size, "sampling_entry_size");
    numeric_validation::RequirePositive(scenario.replica_size, "replica_size");
    GaussianModel3D::RequireFinitePositiveWidthModel(
        scenario.inlier_distribution.mean,
        "scenario.inlier_distribution.mean");
    GaussianModel3D::RequireFinitePositiveWidthModel(
        scenario.outlier_distribution.mean,
        "scenario.outlier_distribution.mean");
    GaussianModel3DUncertainty::RequireFiniteNonNegativeUncertainty(
        scenario.inlier_distribution.sigma,
        "scenario.inlier_distribution.sigma");
    GaussianModel3DUncertainty::RequireFiniteNonNegativeUncertainty(
        scenario.outlier_distribution.sigma,
        "scenario.outlier_distribution.sigma");

    GroupTestData input;
    input.gaus_true = scenario.inlier_distribution.mean;
    input.replica_member_sampling_entries.reserve(static_cast<size_t>(scenario.replica_size));

    for (int i = 0; i < scenario.replica_size; i++)
    {
        auto generator{ BuildReplicaGenerator(i, scenario.random_seed) };
        const auto random_gaus_array{
            BuildRandomGausParameters(
                scenario.member_size,
                scenario.inlier_distribution,
                scenario.outlier_distribution,
                scenario.outlier_ratio,
                generator
            )
        };
        input.replica_member_sampling_entries.emplace_back(
            BuildMuMemberSamplingEntries(
                random_gaus_array,
                static_cast<size_t>(scenario.sampling_entry_size),
                generator));
    }

    return input;
}

LocalTestData BuildLocalTestData(const AtomModelScenario & scenario)
{
    numeric_validation::RequirePositive(scenario.replica_size, "replica_size");
    GaussianModel3D::RequireFinitePositiveWidthModel(scenario.gaus_true, "scenario.gaus_true");

    LocalTestData input;
    input.gaus_true = scenario.gaus_true;
    input.replica_sampling_entries.reserve(static_cast<size_t>(scenario.replica_size));

    for (int i = 0; i < scenario.replica_size; i++)
    {
        auto generator{ BuildReplicaGenerator(i, scenario.random_seed) };
        auto sampling_entries{
            GenerateAtomModelSamples(scenario.gaus_true, scenario.spot)
        };
        sampling_entries = ApplyLogQuadraticNoise(
            std::move(sampling_entries),
            scenario.gaus_true,
            scenario.data_error_sigma,
            generator);
        input.replica_sampling_entries.emplace_back(std::move(sampling_entries));
    }

    return input;
}

AtomicModelTestData BuildAtomicModelTestData(const AtomModelScenario & scenario)
{
    numeric_validation::RequirePositive(scenario.replica_size, "replica_size");
    GaussianModel3D::RequireFinitePositiveWidthModel(scenario.gaus_true, "scenario.gaus_true");

    AtomicModelTestData input;
    input.gaus_true = scenario.gaus_true;
    input.replica_model_objects.reserve(static_cast<size_t>(scenario.replica_size));
    const auto atom_field{ BuildAtomicModelContributionList(scenario.spot) };

    for (int i = 0; i < scenario.replica_size; i++)
    {
        auto generator{ BuildReplicaGenerator(i, scenario.random_seed) };
        input.replica_model_objects.emplace_back(
            BuildAtomicModelObject(
                scenario.gaus_true,
                atom_field,
                scenario.data_error_sigma,
                generator));
    }

    return input;
}

} // namespace rhbm_gem::core
