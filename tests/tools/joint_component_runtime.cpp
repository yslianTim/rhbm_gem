#include "support/JointComponentRuntime.hpp"
#include "support/JointFixtureSupport.hpp"
#include <iostream>
int main(int argc,char ** argv)
{
    namespace p=second_stage_test::matched::joint_abc;
    try {
        const std::string command=argc>1 ? argv[1] : "";
        if(command=="run" && argc==5) p::RunJointRuntime(argv[2],argv[3],argv[4]);
        else if(command=="physical" && argc==3) p::RunPhysicalJointRuntime(argv[2]);
        else if(command=="physical-inputs" && argc==3) p::RunPhysicalJointRuntime(argv[2],true);
        else if(command=="fixture" && argc==5) p::RunFrozenFixture(argv[2],argv[3],argv[4]);
        else throw std::invalid_argument("Usage: joint_component_runtime run MODEL MAP OUTPUT | physical|physical-inputs OUTPUT | fixture DATASET CASE OUTPUT_JSON");
        return 0;
    } catch(const std::exception & e) {std::cerr<<e.what()<<'\n'; return 1;}
}
