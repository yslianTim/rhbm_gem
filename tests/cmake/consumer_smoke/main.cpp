#include <rhbm_gem/core/command/CommandApi.hpp>
#include <rhbm_gem/data/object/GaussianStatistics.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/domain/FilePathHelper.hpp>

#include <string>

int main()
{
    rhbm_gem::ModelObject model_object;
    (void)model_object;
    const rhbm_gem::GaussianEstimate estimate{ 4.0, 2.0 };

    const std::string extension{ FilePathHelper::GetExtension("sample.cif") };
    const auto default_database_path{ rhbm_gem::GetDefaultDatabasePath() };
    return extension == ".cif"
            && default_database_path.filename() == "database.sqlite"
            && estimate.Intensity() > 0.0 ? 0 : 1;
}
