# Adding a Built-in Command

This guide describes the current workflow for adding a built-in command in this repository.
It is implementation-oriented and aligned with the codebase as it exists today.

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

For a new built-in command, update these areas in one change:

- `include/rhbm_gem/core/command/<YourCommand>.hpp`
- `src/core/command/<YourCommand>.cpp` (and optional `src/core/workflow/*Workflow.cpp` files)
- `src/core/CMakeLists.txt` (add source file)
- `include/rhbm_gem/core/command/CommandMetadata.hpp` (add new `CommandId`)
- `src/core/internal/BuiltInCommandList.def` (register built-in metadata)
- `src/core/command/BuiltInCommandCatalog.cpp` (add command header include)
- `bindings/*Bindings.cpp` (add Python bindings in a command-specific binding file)
- `tests/core/command/<YourCommand>_test.cpp`
- `tests/cmake/core_tests.cmake` (wire tests)
- related contract tests under `tests/core/contract/`
- developer docs affected by command surface changes

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

    ExampleCommand();
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

## 4. Register the built-in command

The built-in command manifest is `src/core/internal/BuiltInCommandList.def`.
Add one entry:

```cpp
RHBM_GEM_BUILTIN_COMMAND(
    ExampleCommand,
    "example_command",
    "Run example command",
    "ExampleCommand")
```

Also add required include(s):

- include `<rhbm_gem/core/command/ExampleCommand.hpp>` in `src/core/command/BuiltInCommandCatalog.cpp`
- include `<rhbm_gem/core/command/ExampleCommand.hpp>` in `bindings/CoreBindings.cpp`

Important:

- built-in order in `BuiltInCommandList.def` affects CLI help order
- `python_binding_name` in the macro must match Python class name exposed in bindings

## 5. Add Python bindings

Python bindings are explicit/manual and split by command under `bindings/`.

Required actions:

1. create/update the command-specific binding source (for example `bindings/ExampleCommandBindings.cpp`)
2. bind class with `BindBuiltInCommand<rhbm_gem::ExampleCommand>(m)`
3. expose constructor that requires `DataIoServices` context
4. expose supported setters and `Execute`
5. call `BindCommandDiagnostics(example_command)`

Binding skeleton:

```cpp
auto example_command{
    BindBuiltInCommand<rhbm_gem::ExampleCommand>(m)
};
example_command
    .def(py::init<const rhbm_gem::DataIoServices &>())
    .def("Execute", &rhbm_gem::ExampleCommand::Execute)
    .def("SetInputPath", &rhbm_gem::ExampleCommand::SetInputPath)
    .def("SetOutputStem", &rhbm_gem::ExampleCommand::SetOutputStem)
    .def("SetThreshold", &rhbm_gem::ExampleCommand::SetThreshold);
BindCommandDiagnostics(example_command);
```

Python callers must create one context and pass it into command constructors:

```python
services = m.DataIoServices()
command = m.ExampleCommand(services)
```

Also register the new binding function call in `bindings/CoreBindings.cpp` (module entrypoint).

Do not expose internal hooks such as `RegisterCLIOptionsExtend(...)` or `ValidateOptions()`.

## 6. Build wiring

Add command source files to `src/core/CMakeLists.txt` in `CORE_SOURCES`.
Use the `command/<YourCommand>.cpp` path entry.

If you add extra workflow files, include them in the same list.

## 7. Testing updates

Add command tests and update build wiring:

- create `tests/core/command/<YourCommand>_test.cpp`
- append it to `CORE_COMMAND_TEST_SOURCES` in `tests/cmake/core_tests.cmake`

Update contract tests that assert explicit built-in lists or per-command mappings:

- `tests/core/contract/BuiltInCommandCatalog_test.cpp`
- `tests/core/contract/BuiltInCommandBindingName_test.cpp`
- `tests/core/contract/CommandDescriptorShape_test.cpp`
- `tests/core/contract/CommandCommonOptions_test.cpp`

If command docs include generated built-in lists/matrices, update those sections and run docs sync tests.

## 8. Documentation updates

When command surface changes, update:

- `docs/developer/architecture/command-architecture.md`
- user-facing docs/examples when applicable

## 9. Practical checklist

Before merging:

1. `CommandId` was added in `include/rhbm_gem/core/command/CommandMetadata.hpp`
2. command class compiles and is added to `src/core/CMakeLists.txt`
3. built-in macro entry exists in `src/core/internal/BuiltInCommandList.def`
4. command headers are included where required (`BuiltInCommandCatalog.cpp`, `CoreBindings.cpp`)
5. Python class and binding methods are added in `bindings/CoreBindings.cpp`
6. command tests and contract tests are updated
7. docs are updated to reflect final command surface

Recommended verification:

1. `ctest --output-on-failure -R "Command|BuiltIn|DocsSync"`
