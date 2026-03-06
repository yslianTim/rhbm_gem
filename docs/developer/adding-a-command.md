# Adding a Built-in Command

This guide is the practical companion to [`architecture/command-architecture.md`](architecture/command-architecture.md). Use it when you need to add a new built-in command or refactor an existing command into the current command framework.

## 1. Before you add a new top-level command

Start by deciding whether you really need a new built-in command.

Add a new top-level command when:

- the workflow has a distinct entry point and output contract
- the option set is large enough that an extra enum mode would become hard to understand
- the command needs a different combination of shared options such as `--database` or `--folder`

Prefer extending an existing command when:

- the new behavior is just another mode of the same workflow
- the inputs and outputs are mostly the same
- the validation and data-loading path would be largely duplicated

If you are unsure, compare the shape of the change to the existing commands listed in [`architecture/command-architecture.md`](architecture/command-architecture.md).

## 2. Files you will usually touch

For a new built-in command, expect to update these areas in one change:

- `include/core/<YourCommand>.hpp`
- `src/core/<YourCommand>.cpp`
- `src/core/BuiltInCommandCatalog.cpp`
- `bindings/CoreBindings.cpp`
- `tests/<YourCommand>_test.cpp`
- `tests/CMakeLists.txt`
- relevant docs under `docs/developer/` or `docs/user/`

You should not add public registration headers under `include/core/`. Built-in catalog metadata and registration helpers are internal and stay under `src/core/`.

## 3. Minimal command shape

The normal pattern is:

1. define an options struct derived from `CommandOptions`
2. derive the command from `CommandWithProfileOptions<...>`
3. implement setters on top of `CommandBase` helper methods
4. register command-local CLI options in `RegisterCLIOptionsExtend(...)`
5. put cross-field validation in `ValidateOptions()`
6. clear transient runtime state in `ResetRuntimeState()`
7. keep `ExecuteImpl()` focused on orchestration

## 4. Header template

This is a schematic template. Replace placeholder names, option fields, and `CommandId` with real values.

```cpp
#pragma once

#include <filesystem>
#include <memory>
#include <string>

#include <CLI/CLI.hpp>

#include "CommandBase.hpp"

namespace rhbm_gem {

class ModelObject;
class MapObject;

struct ExampleCommandOptions : public CommandOptions
{
    std::filesystem::path input_file_path;
    std::string output_stem{"example"};
    double threshold{ 1.0 };
    bool enable_feature{ false };
};

class ExampleCommand
    : public CommandWithProfileOptions<
          ExampleCommandOptions,
          CommandId::Example,
          CommonOptionProfile::FileWorkflow>
{
public:
    using Options = ExampleCommandOptions;

    ExampleCommand();
    ~ExampleCommand() override;

    void SetInputFilePath(const std::filesystem::path & path);
    void SetOutputStem(const std::string & value);
    void SetThreshold(double value);
    void SetEnableFeature(bool value);

private:
    std::shared_ptr<ModelObject> m_model_object;
    std::shared_ptr<MapObject> m_map_object;

    void RegisterCLIOptionsExtend(CLI::App * cmd) override;
    void ValidateOptions() override;
    void ResetRuntimeState() override;
    bool ExecuteImpl() override;

    bool BuildDataObject();
    void RunWorkflow();
    std::filesystem::path BuildOutputFilePath() const;
};

} // namespace rhbm_gem
```

## 5. Source template

This is the smallest useful implementation shape. It shows where option registration, setter validation, and execution phases belong.

```cpp
#include "ExampleCommand.hpp"

#include "CommandDataLoaderInternal.hpp"
#include "CommandOptionBinding.hpp"
#include "DataObjectManager.hpp"
#include "Logger.hpp"
#include "ModelObject.hpp"
#include "ScopeTimer.hpp"

#include <stdexcept>

namespace {
constexpr std::string_view kInputFlags{ "-i,--input" };
constexpr std::string_view kInputOption{ "--input" };
constexpr std::string_view kOutputStemOption{ "--output-stem" };
constexpr std::string_view kThresholdOption{ "--threshold" };
constexpr std::string_view kEnableFeatureOption{ "--enable-feature" };
constexpr std::string_view kThresholdRangeIssue{ "--threshold-range" };
constexpr std::string_view kInputKey{ "input" };
} // namespace

namespace rhbm_gem {

ExampleCommand::ExampleCommand() = default;
ExampleCommand::~ExampleCommand() = default;

void ExampleCommand::RegisterCLIOptionsExtend(CLI::App * cmd)
{
    command_cli::AddPathOption(
        cmd,
        kInputFlags,
        [&](const std::filesystem::path & value) { SetInputFilePath(value); },
        "Input file path",
        std::nullopt,
        true);
    command_cli::AddStringOption(
        cmd,
        kOutputStemOption,
        [&](const std::string & value) { SetOutputStem(value); },
        "Output file stem",
        m_options.output_stem);
    command_cli::AddScalarOption<double>(
        cmd,
        kThresholdOption,
        [&](double value) { SetThreshold(value); },
        "Threshold value",
        m_options.threshold);
    command_cli::AddScalarOption<bool>(
        cmd,
        kEnableFeatureOption,
        [&](bool value) { SetEnableFeature(value); },
        "Enable optional workflow branch",
        m_options.enable_feature);
}

void ExampleCommand::SetInputFilePath(const std::filesystem::path & path)
{
    SetRequiredExistingPathOption(
        m_options.input_file_path,
        path,
        kInputOption,
        "Input file");
}

void ExampleCommand::SetOutputStem(const std::string & value)
{
    MutateOptions([&]() { m_options.output_stem = value; });
}

void ExampleCommand::SetThreshold(double value)
{
    SetFinitePositiveScalarOption(
        m_options.threshold,
        value,
        kThresholdOption,
        1.0,
        "Threshold must be a finite positive value.");
}

void ExampleCommand::SetEnableFeature(bool value)
{
    MutateOptions([&]() { m_options.enable_feature = value; });
}

void ExampleCommand::ValidateOptions()
{
    ResetPrepareIssues(kThresholdRangeIssue);
    if (m_options.output_stem.empty())
    {
        AddValidationError(
            kThresholdRangeIssue,
            "Output stem must not be empty.");
    }
}

void ExampleCommand::ResetRuntimeState()
{
    m_model_object.reset();
    m_map_object.reset();
}

bool ExampleCommand::ExecuteImpl()
{
    if (!BuildDataObject())
    {
        return false;
    }

    RunWorkflow();
    return true;
}

bool ExampleCommand::BuildDataObject()
{
    ScopeTimer timer("ExampleCommand::BuildDataObject");
    try
    {
        m_model_object = command_data_loader::ProcessModelFile(
            m_data_manager,
            m_options.input_file_path,
            kInputKey,
            "input file");
    }
    catch (const std::exception & ex)
    {
        Logger::Log(
            LogLevel::Error,
            "ExampleCommand::BuildDataObject : " + std::string(ex.what()));
        return false;
    }
    return true;
}

void ExampleCommand::RunWorkflow()
{
    ScopeTimer timer("ExampleCommand::RunWorkflow");

    const auto output_path{ BuildOutputFilePath() };
    (void)output_path;

    // Run command-specific orchestration here.
}

std::filesystem::path ExampleCommand::BuildOutputFilePath() const
{
    return BuildOutputPath(m_options.output_stem, ".txt");
}

} // namespace rhbm_gem
```

## 6. Setter rules

Use command setters as the first validation boundary.

Good patterns:

- `SetRequiredExistingPathOption(...)` for required file inputs
- `SetOptionalExistingPathOption(...)` for optional file inputs
- `SetFinitePositiveScalarOption(...)` for positive numeric values
- `SetFiniteNonNegativeScalarOption(...)` for zero-or-positive numeric values
- `SetValidatedEnumOption(...)` for enum-backed mode selection
- `MutateOptions(...)` for custom parsing or custom issue handling

Avoid:

- writing directly into `m_options` outside `MutateOptions(...)`
- deferring obvious single-field validation until `ExecuteImpl()`
- creating directories inside setters

## 7. Where validation belongs

Keep these responsibilities separate:

- setters: single-field validation, normalization, parse-phase issues
- `ValidateOptions()`: cross-field rules and mode-dependent rules
- `ExecuteImpl()`: runtime failures that can only be known after preparation

Typical `ValidateOptions()` checks:

- min/max ordering such as `--fit-min <= --fit-max`
- mode-dependent requirements such as "this option is required when mode X is selected"
- semantic conflicts between multiple options

## 8. Data loading patterns

Use the existing command-facing boundaries rather than building a parallel loading path.

For file-driven commands:

- prefer non-template loaders from `CommandDataLoaderInternal.hpp`, for example
  `command_data_loader::ProcessModelFile(...)` and
  `command_data_loader::ProcessMapFile(...)`
- use `m_data_manager` as the owner of processed objects

For database-backed commands:

1. call `RequireDatabaseManager()`
2. load or save objects through `m_data_manager`

Use `BuildOutputPath(...)` for command outputs when the command uses the
`CommonOptionProfile::FileWorkflow` or `CommonOptionProfile::DatabaseWorkflow` profile.

## 9. Register the command in the built-in catalog

Add the descriptor entry in [`src/core/BuiltInCommandCatalog.cpp`](../../src/core/BuiltInCommandCatalog.cpp).

Template:

```cpp
MakeBuiltInDescriptor<ExampleCommand>(
    "example_command",
    "Run example command",
    "ExampleCommand"),
```

Rules:

- the descriptor name is the CLI subcommand name
- `python_binding_name` is required for every built-in command
- the concrete command type remains the source of `CommandId`, `kCommonOptionProfile`,
  and `kCommonOptions`

Do not add a public `BuiltInCommandCatalog.hpp` header.

## 10. Add the Python binding

Add the binding in [`bindings/CoreBindings.cpp`](../../bindings/CoreBindings.cpp).

Typical pattern:

```cpp
auto example_command{
    BindBuiltInCommand<rhbm_gem::ExampleCommand>(m)
};
example_command
    .def(py::init<>())
    .def("Execute", &rhbm_gem::ExampleCommand::Execute)
    .def("SetInputFilePath", &rhbm_gem::ExampleCommand::SetInputFilePath)
    .def("SetOutputStem", &rhbm_gem::ExampleCommand::SetOutputStem)
    .def("SetThreshold", &rhbm_gem::ExampleCommand::SetThreshold)
    .def("SetEnableFeature", &rhbm_gem::ExampleCommand::SetEnableFeature);
BindCommandDiagnostics(example_command);
```

Expose only the supported command surface. Internal lifecycle hooks such as `RegisterCLIOptionsExtend(...)` and `ValidateOptions()` should not be bound.

## 11. Tests to add

For a new command, add at least:

- one happy-path test that executes the main workflow
- invalid-input tests for required paths and scalar / enum validation
- prepare-time validation tests for cross-field or mode-dependent rules
- output verification for generated files or persisted objects

Also update:

- `tests/CMakeLists.txt`
- any catalog or metadata tests affected by the new built-in

The most relevant existing test references are:

- `tests/CommandExecutionContract_test.cpp`
- `tests/PreparedStateInvalidation_test.cpp`
- `tests/CommandCommonOptions_test.cpp`
- `tests/CommandDescriptorShape_test.cpp`
- command-specific tests such as `tests/MapSimulationCommand_test.cpp`

## 12. Final checklist

Before you consider the command done, verify all of the following:

1. the command header is under `include/core/` and the implementation is under `src/core/`
2. the command derives from `CommandWithProfileOptions<...>`
3. setter mutations go through `MutateOptions(...)` or a helper built on top of it
4. shared option policy is declared on the concrete command type
5. the command is registered in `BuiltInCommandCatalog.cpp`
6. the Python binding is added in `bindings/CoreBindings.cpp`
7. tests are added and wired into `tests/CMakeLists.txt`
8. docs and examples are updated if the command is part of the supported workflow

If you are changing an existing command rather than adding a new one, use the same checklist but skip any item that is already structurally in place.
