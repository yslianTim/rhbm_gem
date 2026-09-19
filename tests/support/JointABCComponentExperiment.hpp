#pragma once
#include <string>
namespace second_stage_test::matched::joint_abc {
void ComponentRegression(const std::string & dataset, const std::string & output);
void ComponentSameState(const std::string & dataset, const std::string & output);
void ComponentRun(const std::string & dataset, const std::string & output);
void ComponentAudit(const std::string & dataset, const std::string & run, const std::string & output,
    const std::string & case_name = "");
void ComponentRerun(const std::string & dataset, const std::string & case_name, const std::string & component,
    const std::string & context, const std::string & output);
} // namespace second_stage_test::matched::joint_abc
