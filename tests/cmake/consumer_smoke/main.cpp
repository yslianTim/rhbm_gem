#include <rhbm_gem/core/command/CommandApi.hpp>
#include <rhbm_gem/data/io/DataRepository.hpp>
#include <rhbm_gem/data/io/ModelMapFileIO.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/BondObject.hpp>
#include <rhbm_gem/data/object/ChemicalComponentEntry.hpp>
#include <rhbm_gem/data/object/GaussianStatistics.hpp>
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
    static_assert(kHeaderExposesCompleteType<rhbm_gem::GaussianEstimate>);
    static_assert(kHeaderExposesCompleteType<rhbm_gem::GaussianPosterior>);
    static_assert(kHeaderExposesCompleteType<rhbm_gem::MapObject>);
    static_assert(kHeaderExposesCompleteType<rhbm_gem::ModelAnalysisEditor>);
    static_assert(kHeaderExposesCompleteType<rhbm_gem::LocalPotentialView>);
    static_assert(kHeaderExposesCompleteType<rhbm_gem::ModelObject>);

    rhbm_gem::ModelObject model_object;
    (void)model_object;
    const rhbm_gem::GaussianEstimate estimate{ 4.0, 2.0 };
    const std::array<int, 3> compile_only_sizes{
        static_cast<int>(sizeof(rhbm_gem::AtomObject)),
        static_cast<int>(sizeof(rhbm_gem::BondObject)),
        static_cast<int>(sizeof(rhbm_gem::ChemicalComponentEntry)) };

    const std::string extension{ FilePathHelper::GetExtension("sample.cif") };
    const auto default_database_path{ rhbm_gem::GetDefaultDatabasePath() };
    return extension == ".cif"
            && compile_only_sizes.front() > 0
            && default_database_path.filename() == "database.sqlite"
            && estimate.Intensity() > 0.0 ? 0 : 1;
}
