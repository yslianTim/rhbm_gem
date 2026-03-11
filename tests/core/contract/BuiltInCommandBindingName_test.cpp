#include <gtest/gtest.h>

#include <algorithm>
#include <string_view>

#include "internal/BuiltInCommandCatalogInternal.hpp"
#include <rhbm_gem/core/command/HRLModelTestCommand.hpp>
#include <rhbm_gem/core/command/MapSimulationCommand.hpp>
#include <rhbm_gem/core/command/MapVisualizationCommand.hpp>
#include <rhbm_gem/core/command/PositionEstimationCommand.hpp>
#include <rhbm_gem/core/command/PotentialAnalysisCommand.hpp>
#include <rhbm_gem/core/command/PotentialDisplayCommand.hpp>
#include <rhbm_gem/core/command/ResultDumpCommand.hpp>

namespace rg = rhbm_gem;

namespace {

template <typename CommandType>
void ExpectBindingNameMatchesCatalog()
{
    const auto & catalog{ rg::BuiltInCommandCatalog() };
    const auto descriptor_iter{
        std::find_if(
            catalog.begin(),
            catalog.end(),
            [](const rg::CommandDescriptor & descriptor)
            {
                return descriptor.id == CommandType::kCommandId;
            })
    };
    ASSERT_NE(descriptor_iter, catalog.end());
    EXPECT_EQ(
        std::string_view{ descriptor_iter->python_binding_name },
        rg::BuiltInPythonBindingName(CommandType::kCommandId));
}

} // namespace

TEST(BuiltInCommandBindingNameTest, BindingNamesAreResolvedFromBuiltInCatalog)
{
    ExpectBindingNameMatchesCatalog<rg::PotentialAnalysisCommand>();
    ExpectBindingNameMatchesCatalog<rg::PotentialDisplayCommand>();
    ExpectBindingNameMatchesCatalog<rg::ResultDumpCommand>();
    ExpectBindingNameMatchesCatalog<rg::MapSimulationCommand>();
    ExpectBindingNameMatchesCatalog<rg::MapVisualizationCommand>();
    ExpectBindingNameMatchesCatalog<rg::PositionEstimationCommand>();
    ExpectBindingNameMatchesCatalog<rg::HRLModelTestCommand>();
}
