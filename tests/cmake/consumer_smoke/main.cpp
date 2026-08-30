#include <rhbm_gem/core/CommandSystem.hpp>
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
            && estimate.Intensity() > 0.0 ? 0 : 1;
}
