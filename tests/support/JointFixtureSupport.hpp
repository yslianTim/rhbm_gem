#pragma once
#include "support/JointComponentChecks.hpp"
#include <filesystem>
namespace second_stage_test::matched::joint_abc {
struct FrozenFixture {
    Domain domain{0,{}};
    Eigen::VectorXd y64,y32;
    std::string hash,name;
    std::vector<std::string> ids;
};
FrozenFixture LoadFixture(const std::filesystem::path &);
void RunFrozenFixture(const std::filesystem::path &,const std::string &,const std::filesystem::path &);
}
