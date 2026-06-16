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

struct AtomicModelContribution
{
    Spot spot;
    Element element;
    double charge{ 0.0 };
    Eigen::VectorXd center;
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

double EvaluatePotentialModelResponse(
    const ElectricPotential & model,
    const Eigen::VectorXd & sample_point,
    const std::vector<AtomicModelContribution> & atom_field)
{
    double response{ 0.0 };
    for (const auto & atom : atom_field)
    {
        auto distance{ (sample_point - atom.center).norm() };
        response += model.GetPotentialValue(atom.element, distance, atom.charge);
    }
    return response;
}

AtomicModelContribution MakeAtomNeighborContribution(
    Spot spot,
    Element element,
    double charge,
    const Eigen::Vector3d & direction,
    double distance)
{
    return AtomicModelContribution{ spot, element, charge, distance * direction };
}

std::vector<AtomicModelContribution> BuildAtomNeighborList(const Spot & spot)
{
    switch (spot)
    {
    case Spot::UNK:
        return {};
    case Spot::O:
        return {
            MakeAtomNeighborContribution(Spot::C, Element::CARBON, 0.0,
                Eigen::Vector3d{ 1.0, 0.0, 0.0 }, 1.23)
        };
    case Spot::N:
        return {
            MakeAtomNeighborContribution(Spot::H, Element::HYDROGEN, 0.0,
                Eigen::Vector3d{ 1.0, 0.0, 0.0 }, 1.02),
            MakeAtomNeighborContribution(Spot::C, Element::CARBON, 0.0,
                Eigen::Vector3d{ -0.5, std::sqrt(3) / 2.0, 0.0 }, 1.48),
            MakeAtomNeighborContribution(Spot::CA, Element::CARBON, 0.0,
                Eigen::Vector3d{ -0.5, -std::sqrt(3) / 2.0, 0.0 }, 1.48)
        };
    case Spot::C:
        return {
            MakeAtomNeighborContribution(Spot::O, Element::OXYGEN, 0.0,
                Eigen::Vector3d{ 1.0, 0.0, 0.0 }, 1.23),
            MakeAtomNeighborContribution(Spot::N, Element::NITROGEN, 0.0,
                Eigen::Vector3d{ -0.5, std::sqrt(3) / 2.0, 0.0 }, 1.48),
            MakeAtomNeighborContribution(Spot::CA, Element::CARBON, 0.0,
                Eigen::Vector3d{ -0.5, -std::sqrt(3) / 2.0, 0.0 }, 1.54)
        };
    case Spot::CA:
        return {
            MakeAtomNeighborContribution(Spot::H, Element::HYDROGEN, 0.0,
                Eigen::Vector3d{ 0.0, 0.0, 1.0 }, 1.06),
            MakeAtomNeighborContribution(Spot::N, Element::NITROGEN, 0.0,
                Eigen::Vector3d{ 0.0, 2.0 * std::sqrt(2) / 3.0, -1.0 / 3.0 }, 1.48),
            MakeAtomNeighborContribution(Spot::C, Element::CARBON, 0.0,
                Eigen::Vector3d{ -std::sqrt(6) / 3.0, -std::sqrt(2) / 3.0, -1.0 / 3.0 }, 1.54),
            MakeAtomNeighborContribution(Spot::CB, Element::CARBON, 0.0,
                Eigen::Vector3d{ std::sqrt(6) / 3.0, -std::sqrt(2) / 3.0, -1.0 / 3.0 }, 1.54)
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

LocalPotentialSampleList GenerateAtomicModelSampleList(
    const ElectricPotential & model,
    const AtomicModelContribution & local_atom,
    const std::vector<AtomicModelContribution> & atom_field)
{
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
    sample_filter::FilterSamplingPointList(sample_point_list, local_position, reject_position_list);

    LocalPotentialSampleList sample_list;
    sample_list.reserve(sample_point_list.size());
    for (const auto & sampling_point : sample_point_list)
    {
        const auto point{ eigen_helper::ToEigenVector(sampling_point.position) };
        const auto response{
            EvaluatePotentialModelResponse(model, point, atom_field)
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
    double model_response_max,
    double error_sigma,
    std::mt19937 & generator)
{
    std::normal_distribution<> dist_error(0.0, error_sigma * model_response_max);
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
    atom_object->SetElement(atom.element);
    atom_object->SetSpot(atom.spot);
    atom_object->SetPosition(ToArray3f(atom.center));
    return atom_object;
}

std::unique_ptr<ModelObject> BuildAtomicModelObject(
    const ElectricPotential & model,
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
        const auto & atom_object{ *model_object->GetSelectedAtoms().at(i) };
        auto sampling_entries{ GenerateAtomicModelSampleList(model, atom_field.at(i), atom_field) };
        auto model_response_max{ model.GetPotentialValue(atom_field.at(i).element, 0.0, atom_field.at(i).charge) };
        sampling_entries = ApplyLogQuadraticNoise(std::move(sampling_entries), model_response_max, error_sigma, generator);
        analysis.EnsureAtomLocalPotential(atom_object).SetSamplingEntries(std::move(sampling_entries));
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
                scenario.gaus_true.Intensity(),
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

AtomicModelTestData BuildPotentialModelTestData(const PotentialModelScenario & scenario)
{
    AtomicModelTestData input;
    input.gaus_true = scenario.gaus_true;
    input.replica_model_objects.reserve(static_cast<size_t>(scenario.replica_size));

    auto atom_field{ BuildAtomNeighborList(scenario.spot) };
    atom_field.insert(atom_field.begin(), AtomicModelContribution{
        scenario.spot,
        scenario.element,
        scenario.charge,
        Eigen::VectorXd::Zero(3) }
    );

    for (int i = 0; i < scenario.replica_size; i++)
    {
        auto generator{ BuildReplicaGenerator(i, scenario.random_seed) };
        input.replica_model_objects.emplace_back(
            BuildAtomicModelObject(
                scenario.potential_model,
                atom_field,
                scenario.data_error_sigma,
                generator));
    }

    return input;
}

} // namespace rhbm_gem::core
