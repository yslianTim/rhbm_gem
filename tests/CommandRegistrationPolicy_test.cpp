#include <gtest/gtest.h>

#include <CLI/CLI.hpp>

#include "BuiltInCommandCatalog.hpp"
#include "CommandBase.hpp"
#include "CommandRegistry.hpp"

namespace rg = rhbm_gem;

namespace {

class ExtensionCommand final : public rg::CommandBase
{
public:
    struct Options : public rg::CommandOptions
    {
    };

    void RegisterCLIOptionsExtend(CLI::App * /*command*/) override {}
    const rg::CommandOptions & GetOptions() const override { return m_options; }
    rg::CommandOptions & GetOptions() override { return m_options; }
    rg::CommandId GetCommandId() const override { return rg::CommandId::ModelTest; }

private:
    bool ExecuteImpl() override { return true; }

    Options m_options{};
};

} // namespace

TEST(CommandRegistrationPolicyTest, CommandDescriptorDirectRegistrationCapturesExplicitMetadata)
{
    const auto descriptor{
        rg::CommandDescriptor{
            rg::CommandId::ModelTest,
            "test_extension_command",
            "Test extension command",
            rg::MakeCommandSurface(
                rg::ToMask(rg::CommonOption::Verbose),
                rg::ToMask(rg::CommonOption::Database)),
            rg::DatabaseUsage::Optional,
            rg::BindingExposure::CliOnly,
            "",
            []() { return std::make_unique<ExtensionCommand>(); }
        }
    };

    EXPECT_EQ(descriptor.id, rg::CommandId::ModelTest);
    EXPECT_EQ(descriptor.name, "test_extension_command");
    EXPECT_EQ(descriptor.description, "Test extension command");
    EXPECT_EQ(descriptor.surface.common_options, rg::ToMask(rg::CommonOption::Verbose));
    EXPECT_EQ(
        descriptor.surface.deprecated_hidden_options,
        rg::ToMask(rg::CommonOption::Database));
    EXPECT_EQ(descriptor.database_usage, rg::DatabaseUsage::Optional);
    EXPECT_EQ(descriptor.binding_exposure, rg::BindingExposure::CliOnly);
    EXPECT_TRUE(descriptor.python_binding_name.empty());

    auto command{ descriptor.factory() };
    ASSERT_NE(command, nullptr);
    EXPECT_EQ(command->GetCommandId(), rg::CommandId::ModelTest);
}

TEST(CommandRegistrationPolicyTest, RegistryAcceptsDirectDescriptorRegistration)
{
    const auto descriptor{
        rg::CommandDescriptor{
            rg::CommandId::ModelTest,
            "test_extension_command_direct_registration",
            "Test direct registry registration",
            rg::MakeCommandSurface(
                rg::ToMask(rg::CommonOption::Verbose),
                rg::ToMask(rg::CommonOption::Database)),
            rg::DatabaseUsage::Optional,
            rg::BindingExposure::CliOnly,
            "",
            []() { return std::make_unique<ExtensionCommand>(); }
        }
    };

    auto & registry{ rg::CommandRegistry::Instance() };
    EXPECT_TRUE(registry.RegisterCommand(descriptor));
    EXPECT_FALSE(registry.RegisterCommand(descriptor));
}
