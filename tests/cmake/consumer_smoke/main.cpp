#include <rhbm_gem/core/CommandSystem.hpp>
#include <rhbm_gem/core/JointComponentEstimator.hpp>
#include <rhbm_gem/core/GaussianEstimator.hpp>
#include <rhbm_gem/data/io/DataRepository.hpp>
#include <rhbm_gem/data/io/ModelMapFileIO.hpp>
#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/BondObject.hpp>
#include <rhbm_gem/data/object/ChemicalComponentEntry.hpp>
#include <rhbm_gem/utils/hrl/RHBMTypes.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <rhbm_gem/data/object/ModelAnalysisView.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/algorithm/KDTreeAlgorithm.hpp>
#include <rhbm_gem/utils/domain/FilePathHelper.hpp>

#include <array>
#include <string>

template <typename T>
constexpr bool kHeaderExposesCompleteType = sizeof(T) > 0;

int main()
{
    static_assert(kHeaderExposesCompleteType<rhbm_gem::DataRepository>);
    static_assert(kHeaderExposesCompleteType<rhbm_gem::AtomObject>);
    static_assert(kHeaderExposesCompleteType<rhbm_gem::BondObject>);
    static_assert(kHeaderExposesCompleteType<rhbm_gem::ChemicalComponentEntry>);
    static_assert(kHeaderExposesCompleteType<rhbm_gem::GaussianModel3D>);
    static_assert(kHeaderExposesCompleteType<rhbm_gem::GaussianModel3DWithUncertainty>);
    static_assert(kHeaderExposesCompleteType<rhbm_gem::MapObject>);
    static_assert(kHeaderExposesCompleteType<rhbm_gem::AtomLocalPotentialView>);
    static_assert(kHeaderExposesCompleteType<rhbm_gem::ModelAnalysisEditor>);
    static_assert(kHeaderExposesCompleteType<rhbm_gem::ModelObject>);
    static_assert(kHeaderExposesCompleteType<rhbm_gem::core::FitOptions>);
#ifdef RHBM_GEM_ENABLE_UMAP
    static_assert(kHeaderExposesCompleteType<rhbm_gem::core::UmapEmbeddingRequest>);
#endif

    rhbm_gem::core::JointProblemInput input;
    input.atom_ids={"unobserved"}; input.support.resize(1);
    const rhbm_gem::core::JointProblem problem(std::move(input));
    const auto joint=rhbm_gem::core::FitJointComponents(problem,{0.});
    if(joint.initialization.valid || joint.regular_certificate!=rhbm_gem::core::JointCheckStatus::NotRun) return 2;
    const auto saved=rhbm_gem::core::CaptureJointAnalysisResult(joint);
    rhbm_gem::ModelObject saved_model;
    saved_model.EditAnalysis().SetJointResult(saved);
    if(!saved_model.GetAnalysisView().GetJointResult()) return 5;
    const auto unobserved=rhbm_gem::core::FitJointComponents(problem,{.5});
    if(!unobserved.initialization.valid || unobserved.components.size()!=1 || unobserved.prediction) return 3;
    if(joint.RuntimeConvergence()!=rhbm_gem::core::JointCheckStatus::Unavailable ||
        unobserved.RuntimeConvergence()!=rhbm_gem::core::JointCheckStatus::Unavailable ||
        unobserved.components[0].RuntimeConvergence()!=rhbm_gem::core::JointCheckStatus::Unavailable) return 4;
    rhbm_gem::ModelObject model_object;
    (void)model_object;
    rhbm_gem::core::FitOptions fit_options;
    using PotentialFittingWorkflow = void (*)(
        rhbm_gem::ModelObject &,
        const rhbm_gem::core::FitOptions &);
    PotentialFittingWorkflow volatile workflow_entry{
        &rhbm_gem::core::RunPotentialFittingWorkflow
    };
    const rhbm_gem::GaussianModel3D estimate{ 4.0, 2.0 };
    const auto transformed_estimate{ estimate.ToTransformedCoordinates() };
    const auto transformed_round_trip{
        transformed_estimate.has_value() ?
            rhbm_gem::GaussianModel3D::FromTransformedCoordinates(
                *transformed_estimate) :
            std::nullopt
    };
    const std::array<int, 3> compile_only_sizes{
        static_cast<int>(sizeof(rhbm_gem::AtomObject)),
        static_cast<int>(sizeof(rhbm_gem::BondObject)),
        static_cast<int>(sizeof(rhbm_gem::ChemicalComponentEntry)) };

    const std::string extension{ rhbm_gem::path_helper::GetExtension("sample.cif") };
    const auto default_database_path{ rhbm_gem::core::GetDefaultDatabasePath() };
    return extension == ".cif"
            && compile_only_sizes.front() > 0
            && default_database_path.filename() == "database.sqlite"
            && fit_options.thread_size == 1
            && workflow_entry != nullptr
            && estimate.Intensity() > 0.0
            && transformed_round_trip.has_value() ? 0 : 1;
}
