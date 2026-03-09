#include <gtest/gtest.h>

#include <algorithm>
#include <string_view>

#include "internal/BuiltInCommandBindingInternal.hpp"
#include "internal/BuiltInCommandCatalogInternal.hpp"
#include <rhbm_gem/core/HRLModelTestCommand.hpp>
#include <rhbm_gem/core/MapSimulationCommand.hpp>
#include <rhbm_gem/core/MapVisualizationCommand.hpp>
#include <rhbm_gem/core/PositionEstimationCommand.hpp>
#include <rhbm_gem/core/PotentialAnalysisCommand.hpp>
#include <rhbm_gem/core/PotentialDisplayCommand.hpp>
#include <rhbm_gem/core/ResultDumpCommand.hpp>

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
        rg::BuiltInCommandBindingName<CommandType>::value);
}

} // namespace

TEST(BuiltInCommandBindingNameTest, BindingNamesMatchBuiltInCatalogMetadata)
{
    ExpectBindingNameMatchesCatalog<rg::PotentialAnalysisCommand>();
    ExpectBindingNameMatchesCatalog<rg::PotentialDisplayCommand>();
    ExpectBindingNameMatchesCatalog<rg::ResultDumpCommand>();
    ExpectBindingNameMatchesCatalog<rg::MapSimulationCommand>();
    ExpectBindingNameMatchesCatalog<rg::MapVisualizationCommand>();
    ExpectBindingNameMatchesCatalog<rg::PositionEstimationCommand>();
    ExpectBindingNameMatchesCatalog<rg::HRLModelTestCommand>();
}
