#pragma once

#include <filesystem>
#include <string>
#include <vector>
#include <Eigen/Dense>

#include "CommandBase.hpp"
#include <rhbm_gem/core/command/OptionEnumClass.hpp>

namespace CLI
{
    class App;
}

namespace rhbm_gem {

struct HRLModelTestRequest;

struct HRLModelTestCommandOptions : public CommandOptions
{
    TesterType tester_choice{ TesterType::BENCHMARK };
    double fit_range_min{ 0.0 };
    double fit_range_max{ 1.0 };
    double alpha_r{ 0.1 };
    double alpha_g{ 0.2 };
};

class HRLModelTestCommand
    : public CommandWithProfileOptions<
          HRLModelTestCommandOptions,
          CommandId::ModelTest>
{
public:
    using Options = HRLModelTestCommandOptions;

    explicit HRLModelTestCommand();
    ~HRLModelTestCommand() = default;
    void ApplyRequest(const HRLModelTestRequest & request);

private:
    void SetTesterChoice(TesterType value);
    void SetFitRangeMinimum(double value);
    void SetFitRangeMaximum(double value);
    void SetAlphaR(double value);
    void SetAlphaG(double value);
    void ValidateOptions() override;
    bool ExecuteImpl() override;
};

} // namespace rhbm_gem
