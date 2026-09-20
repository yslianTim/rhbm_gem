// Public-API-only benchmark, also buildable against the pre-change library.
#include <rhbm_gem/core/JointComponentEstimator.hpp>
#include <boost/json.hpp>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <sys/resource.h>
namespace {
namespace j=boost::json;
using Clock=std::chrono::steady_clock;
double Seconds(Clock::time_point t) {return std::chrono::duration<double>(Clock::now()-t).count();}
j::value Read(const std::filesystem::path & p)
{
    std::ifstream f(p); if(!f) throw std::runtime_error("Missing benchmark input.");
    j::parse_options options; options.numbers=j::number_precision::precise;
    return j::parse(std::string(std::istreambuf_iterator<char>(f),{}),{},options);
}
}
int main(int argc,char ** argv)
{
    try {
        if(argc!=3) throw std::invalid_argument("Usage: joint_component_benchmark DATASET CASE");
        const auto started=Clock::now(); const std::filesystem::path root(argv[1]); const std::string name(argv[2]);
        const auto dataset=Read(root/"dataset.json"),cases=Read(root/"cases.json");
        rhbm_gem::core::JointProblemInput input;
        for(const auto & a:dataset.at("atoms").as_array()) input.atom_ids.push_back(std::to_string(j::value_to<int>(a.at("serial_id"))));
        input.support.resize(input.atom_ids.size()); std::string line;
        std::ifstream voxels(root/"voxels.csv"); std::getline(voxels,line);
        while(std::getline(voxels,line))
        {
            std::istringstream row(line); std::string cell; std::vector<double> values;
            while(std::getline(row,cell,',')) values.push_back(std::stod(cell));
            input.row_ids.push_back(std::to_string(input.observations.size()));
            input.observations.push_back(values.at(name.ends_with("double") ? 7 : 8));
        }
        std::ifstream contributors(root/"contributors.csv"); std::getline(contributors,line);
        while(std::getline(contributors,line))
        {
            std::istringstream row(line); std::string r,a,s;
            std::getline(row,r,','); std::getline(row,a,','); std::getline(row,s,',');
            input.support.at(std::stoull(a)).push_back({std::stoull(r),std::stod(s)});
        }
        std::vector<double> initial; for(const auto & b:cases.at(name).at("initial_b").as_array()) initial.push_back(j::value_to<double>(b));
        const auto construction=Clock::now(); const rhbm_gem::core::JointProblem problem(std::move(input));
        const double construction_seconds=Seconds(construction);
        const auto fit=rhbm_gem::core::FitJointComponents(problem,initial);
        rusage usage{}; getrusage(RUSAGE_SELF,&usage);
#ifdef __APPLE__
        const auto peak=usage.ru_maxrss;
#else
        const auto peak=usage.ru_maxrss*1024;
#endif
        std::cout<<j::serialize(j::object{{"construction_seconds",construction_seconds},{"search_seconds",fit.costs.search_seconds},
            {"assessment_seconds",fit.costs.assessment_seconds},{"assembly_seconds",fit.costs.assembly_seconds},
            {"total_seconds",Seconds(started)},{"process_peak_rss_bytes",peak},
            {"available",fit.assembled_state.has_value()},{"objective",fit.objective ? j::value(*fit.objective) : j::value(nullptr)}})<<'\n';
        return fit.assembled_state ? 0 : 1;
    } catch(const std::exception & e) {std::cerr<<e.what()<<'\n'; return 1;}
}
