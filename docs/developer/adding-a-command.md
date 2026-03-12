# Adding a Command

This guide describes how to add a command in this repository.

Related references:

- [`architecture/command-architecture.md`](architecture/command-architecture.md)
- [`development-guidelines.md`](development-guidelines.md)

## 1. Decide whether to add a top-level command

Add a new top-level command when:

- the workflow has a distinct entry contract
- option semantics are materially different from existing commands
- shared-option policy (`--database`, `--folder`) differs from existing command modes

Prefer extending an existing command when:

- behavior is a mode switch of the same workflow
- inputs/outputs are mostly identical
- validation and loading paths would otherwise be duplicated

## 2. Required change set

For a new command, the minimum manual change set is:

- `include/rhbm_gem/core/command/<YourCommand>.hpp`
- `src/core/command/<YourCommand>.cpp` (and optional `src/core/workflow/*Workflow.cpp` files)
- `include/rhbm_gem/core/command/CommandMetadata.hpp` (add new `CommandId`)
- `src/core/internal/CommandList.def` (register command metadata)
- `bindings/<YourFeature>Bindings.cpp` (add command-specific Python bindings)
- `tests/core/command/<YourCommand>_test.cpp`
- related contract tests under `tests/core/contract/`
- developer docs affected by command surface changes

Generated registration surfaces are derived from `CommandList.def`.
After changing commands, run:

```bash
python3 scripts/generate_command_artifacts.py
```

This refreshes:

- `src/core/internal/CommandCatalogIncludes.generated.inc`
- `src/core/internal/CommandCatalogEntries.generated.inc`
- `src/CommandSources.generated.cmake`
- `bindings/CommandBindingUnits.generated.cmake`
- `tests/cmake/CoreCommandTests.generated.cmake`
- generated sections in `docs/developer/architecture/command-architecture.md`

Scaffold helper:

```bash
python3 scripts/command_scaffold.py --name Example --profile FileWorkflow
```

The scaffold creates command/binding/test/doc skeleton files only.
Use `--wire` to apply registration/manifests automatically.

Note:

- scaffolded binding files do not add command-specific setters or `BindCommonCommandSetters(...)`; add them manually

Scaffold + auto-wire helper:

```bash
python3 scripts/command_scaffold.py --name Example --profile FileWorkflow --wire
```

`--wire` will also update:

- `src/core/internal/CommandList.def`
- generated command artifacts (catalog includes/entries, CMake source lists, docs blocks)

Strict wiring mode (recommended in CI/local guard runs):

```bash
python3 scripts/command_scaffold.py --name Example --profile FileWorkflow --wire --strict
```

`--wire --strict` will fail-fast when manifest wiring or generation fails.

## 3. Implement the command type

Standard command shape:

1. define `Options` derived from `CommandOptions`
2. derive from `CommandWithProfileOptions<OptionsT, CommandId::..., ...>`
3. implement setters using `CommandBase` helper APIs
4. register command-local CLI options in `RegisterCLIOptionsExtend(...)`
5. keep cross-field or mode-dependent checks in `ValidateOptions()`
6. clear transient state in `ResetRuntimeState()`
7. keep `ExecuteImpl()` orchestration-focused

Header skeleton:

```cpp
#pragma once

#include <filesystem>
#include <memory>
#include <string>

#include <rhbm_gem/core/command/CommandBase.hpp>

namespace CLI
{
    class App;
}

namespace rhbm_gem {

class ModelObject;

struct ExampleCommandOptions : public CommandOptions
{
    std::filesystem::path input_path;
    std::string output_stem{"example"};
    double threshold{1.0};
};

class ExampleCommand
    : public CommandWithProfileOptions<
          ExampleCommandOptions,
          CommandId::Example,
          CommonOptionProfile::FileWorkflow>
{
public:
    using Options = ExampleCommandOptions;

    explicit ExampleCommand();
    ~ExampleCommand() override = default;

    void SetInputPath(const std::filesystem::path & path);
    void SetOutputStem(const std::string & value);
    void SetThreshold(double value);

private:
    std::shared_ptr<ModelObject> m_model_object;

    void RegisterCLIOptionsExtend(CLI::App * cmd) override;
    void ValidateOptions() override;
    void ResetRuntimeState() override;
    bool ExecuteImpl() override;
    bool BuildDataObject();
};

} // namespace rhbm_gem
```

Source skeleton:

```cpp
#include <rhbm_gem/core/command/ExampleCommand.hpp>

#include "internal/CommandDataLoaderInternal.hpp"
#include <rhbm_gem/core/command/CommandOptionBinding.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/domain/ScopeTimer.hpp>

namespace {
constexpr std::string_view kInputFlags{ "-i,--input" };
constexpr std::string_view kInputOption{ "--input" };
constexpr std::string_view kOutputStemOption{ "--output-stem" };
constexpr std::string_view kThresholdOption{ "--threshold" };
constexpr std::string_view kOutputRuleIssue{ "--output-rule" };
constexpr std::string_view kInputKey{ "input" };
} // namespace

namespace rhbm_gem {

ExampleCommand::ExampleCommand() :
    CommandWithProfileOptions<
        ExampleCommandOptions,
        CommandId::Example,
        CommonOptionProfile::FileWorkflow>{},
    m_model_object{ nullptr }
{
}

void ExampleCommand::RegisterCLIOptionsExtend(CLI::App * cmd)
{
    command_cli::AddPathOption(
        cmd, kInputFlags,
        [&](const std::filesystem::path & value) { SetInputPath(value); },
        "Input file path",
        std::nullopt,
        true);
    command_cli::AddStringOption(
        cmd, kOutputStemOption,
        [&](const std::string & value) { SetOutputStem(value); },
        "Output file stem",
        m_options.output_stem);
    command_cli::AddScalarOption<double>(
        cmd, kThresholdOption,
        [&](double value) { SetThreshold(value); },
        "Threshold value",
        m_options.threshold);
}

void ExampleCommand::SetInputPath(const std::filesystem::path & path)
{
    SetRequiredExistingPathOption(m_options.input_path, path, kInputOption, "Input file");
}

void ExampleCommand::SetOutputStem(const std::string & value)
{
    MutateOptions([&]()
    {
        ResetParseIssues(kOutputStemOption);
        m_options.output_stem = value;
        if (m_options.output_stem.empty())
        {
            AddValidationError(
                kOutputStemOption,
                "Output stem must not be empty.",
                ValidationPhase::Parse);
        }
    });
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

void ExampleCommand::ValidateOptions()
{
    ResetPrepareIssues(kOutputRuleIssue);
    if (m_options.output_stem == "tmp")
    {
        AddValidationError(
            kOutputRuleIssue,
            "Output stem 'tmp' is reserved. Please choose another value.");
    }
}

void ExampleCommand::ResetRuntimeState()
{
    m_model_object.reset();
}

bool ExampleCommand::BuildDataObject()
{
    ScopeTimer timer("ExampleCommand::BuildDataObject");
    try
    {
        m_model_object = command_data_loader::ProcessModelFile(
            m_data_manager,
            m_options.input_path,
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

bool ExampleCommand::ExecuteImpl()
{
    if (!BuildDataObject())
    {
        return false;
    }

    const auto output_path{ BuildOutputPath(m_options.output_stem, ".txt") };
    Logger::Log(LogLevel::Info, "Output path: " + output_path.string());

    // command-specific workflow here
    return true;
}

} // namespace rhbm_gem
```

Implementation rules:

- use `command_cli::Add*Option(...)` helpers for CLI options
- use `MutateOptions(...)` or convenience setters from `CommandBase`
- do not write to `m_options` outside mutation helpers
- do not create directories in setters

Before pushing, run:

```bash
python3 scripts/check_command_sync.py
```

This is also enforced by the repository guard (`lint_repo`).

## 4. Register the command

The command manifest is `src/core/internal/CommandList.def`.
Add one entry:

```cpp
RHBM_GEM_COMMAND(
    ExampleCommand,
    "example_command",
    "Run example command")
```

Important:

- command order in `CommandList.def` affects CLI help order
- Python class name is derived from `COMMAND_TYPE` (for example `ExampleCommand`)

## 5. Add Python bindings

Python bindings are explicit/manual and split by command under `bindings/`.

Required actions:

1. create/update the command-specific binding source (for example `bindings/ExampleBindings.cpp`)
2. bind class with `BindCommandClass<rhbm_gem::ExampleCommand>(m)`
3. expose default constructor
4. expose supported setters and `Execute`
5. call `BindCommonCommandSetters(example_command)` so common options stay synchronized
6. call `BindCommandDiagnostics(example_command)`
7. implement command registration specialization `template <> void BindCommand<rhbm_gem::ExampleCommand>(py::module_ & module)`

Binding skeleton:

```cpp
template <>
void BindCommand<rhbm_gem::ExampleCommand>(py::module_ & module)
{
    auto example_command{ BindCommandClass<rhbm_gem::ExampleCommand>(module) };
    example_command
        .def(py::init<>())
        .def("Execute", &rhbm_gem::ExampleCommand::Execute)
        .def("SetInputPath", &rhbm_gem::ExampleCommand::SetInputPath)
        .def("SetOutputStem", &rhbm_gem::ExampleCommand::SetOutputStem)
        .def("SetThreshold", &rhbm_gem::ExampleCommand::SetThreshold);
    BindCommonCommandSetters(example_command);
    BindCommandDiagnostics(example_command);
}
```

Python callers construct commands directly:

```python
command = m.ExampleCommand()
```

`bindings/CoreBindings.cpp` registers all commands from `CommandList.def`; no manual
module entrypoint edits are required for new command bindings.

Do not expose internal hooks such as `RegisterCLIOptionsExtend(...)` or `ValidateOptions()`.

## 6. Build wiring

Add extra workflow files to `src/rhbm_gem_sources.cmake` when needed.
Command source/binding/test lists are generated from `CommandList.def`.

## 7. Testing updates

Add command tests and update build wiring:

- create `tests/core/command/<YourCommand>_test.cpp`
- run `python3 scripts/generate_command_artifacts.py`

Update contract tests that assert explicit command lists or per-command mappings:

- `tests/core/contract/CommandCatalog_test.cpp`
- `tests/core/contract/DocsSync_test.cpp`
- `tests/core/contract/PublicHeaderSurface_test.cpp`

When command membership or command metadata changes, update generated blocks in:

- `docs/developer/architecture/command-architecture.md`

Then run docs sync tests (`tests/core/contract/DocsSync_test.cpp`).

## 8. Documentation updates

When command surface changes, update:

- `docs/developer/architecture/command-architecture.md`
- user-facing docs/examples when applicable

## 9. Practical checklist

Before merging:

1. `CommandId` was added in `include/rhbm_gem/core/command/CommandMetadata.hpp`
2. command class compiles and any workflow source files are wired in `src/rhbm_gem_sources.cmake`
3. command macro entry exists in `src/core/internal/CommandList.def`
4. generated artifacts were refreshed (`python3 scripts/generate_command_artifacts.py`)
5. Python binding wiring is complete (`bindings/<YourFeature>Bindings.cpp`)
6. command tests and contract tests are updated
7. docs are updated to reflect final command surface

Recommended verification:

1. `ctest --output-on-failure -R "Command|DocsSync"`
