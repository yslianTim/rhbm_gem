#include <rhbm_gem/core/command/CommandMetadata.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/domain/FilePathHelper.hpp>

#include <string>

int main()
{
    rhbm_gem::ModelObject model_object;
    (void)model_object;

    const std::string extension{ FilePathHelper::GetExtension("sample.cif") };
    const auto option{ rhbm_gem::CommonOption::OutputFolder };
    (void)option;
    return extension == ".cif" ? 0 : 1;
}
