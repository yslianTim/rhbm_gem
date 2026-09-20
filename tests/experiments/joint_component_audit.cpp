#include "support/JointOfflineTools.hpp"
#include "support/JointFixtureSupport.hpp"
#include "support/JointDerivativeAudit.hpp"
#include <iostream>
#include <fstream>
int main(int argc,char ** argv)
{
    namespace p=second_stage_test::matched::joint_abc;
    namespace j=boost::json;
    try {
        Eigen::setNbThreads(1);
        if(argc==4 && std::string(argv[1])=="local-bundle") p::ComponentLocalBundleRerun(argv[2],argv[3]);
        else if(argc==5 && (std::string(argv[1])=="fixture" || std::string(argv[1])=="two-step-fixture")) {
            const std::filesystem::path path(argv[2]),output(argv[4]);
            if(std::filesystem::exists(output)) throw std::invalid_argument("Use a fresh audit output directory.");
            const auto in=p::LoadFixture(path); const std::string name(argv[3]);
            std::ifstream stream(path/"cases.json"); j::parse_options options; options.numbers=j::number_precision::precise;
            const auto cases=j::parse(std::string(std::istreambuf_iterator<char>(stream),{}),{},options);
            const auto & values=cases.at(name).at("initial_b").as_array(); Eigen::VectorXd initial(values.size());
            for(std::size_t k=0;k<values.size();++k) initial(static_cast<Eigen::Index>(k))=j::value_to<double>(values[k]);
            const auto & y=name.ends_with("double") ? in.y64 : in.y32;
            auto plan=p::RegisteredAudit(initial.size(),in.name,name); plan.cache_precision=true;
            auto context=p::MakeContext(y,initial.size(),in.hash,&plan); context.atom_ids=in.ids;
            const auto partition=p::BuildPartition(in.domain,in.ids);
            for(std::size_t k=0;k<partition.components.size();++k) {
                auto fit=p::FitComponent(partition.components[k],y,initial,context);
                fit["dataset"]=in.name; fit["case"]=name;
                if(std::string(argv[1])=="two-step-fixture")
                {
                    const auto & saved=fit.at("search_endpoint_eta").as_array(); Eigen::VectorXd eta(saved.size());
                    for(std::size_t a=0;a<saved.size();++a) eta(static_cast<Eigen::Index>(a))=j::value_to<double>(saved[a]);
                    auto assessed=second_stage_test::matched::certification::AssessWithDerivativeAudit(
                        partition.components[k].domain,p::Select(y,partition.components[k].rows),eta,
                        p::ComponentContext(context,partition.components[k],true));
                    j::object record;
                    for(const char * key:{"qualification_checks","qualification_failure"})
                        if(assessed.contains(key)) record[key]=assessed.at(key);
                    std::filesystem::create_directories(output);
                    std::ofstream record_stream(output/(std::to_string(k)+".json")); record_stream.exceptions(std::ios::failbit|std::ios::badbit);
                    record_stream<<j::serialize(record)<<'\n';
                }
                else p::WriteLocalAudit(partition.components[k],y,fit,context,output/std::to_string(k));
            }
        } else throw std::invalid_argument("Usage: joint_component_audit local-bundle INPUT OUTPUT | fixture DATASET CASE OUTPUT | two-step-fixture DATASET CASE OUTPUT");
        return 0;
    } catch(const std::exception & e) {std::cerr<<e.what()<<'\n'; return 1;}
}
