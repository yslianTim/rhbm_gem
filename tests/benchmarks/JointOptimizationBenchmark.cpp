#include "core/detail/joint_component/ActiveSetBoundedLM.hpp"
#include "core/detail/joint_component/Numerics.hpp"
#include "core/detail/joint_component/Problem.hpp"
#include "core/detail/joint_component/TiledDerivative.hpp"
#include "support/JointSyntheticWorkload.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
namespace joint=rhbm_gem::core::joint_component;
using Clock=std::chrono::steady_clock;
using Matrix=Eigen::MatrixXd;
using Vector=Eigen::VectorXd;

double Seconds(Clock::time_point started)
{return std::chrono::duration<double>(Clock::now()-started).count();}

void ReportSamples(const std::string & name,std::vector<double> samples)
{
    auto sorted=samples; std::sort(sorted.begin(),sorted.end());
    const double median=sorted[sorted.size()/2];
    std::cout<<name<<"_ms_per_call_median="<<median*1000.0
        <<" range_ms="<<sorted.front()*1000.0<<".."<<sorted.back()*1000.0<<" samples_ms=";
    for(std::size_t k=0;k<samples.size();++k)
        std::cout<<(k ? "," : "")<<samples[k]*1000.0;
    std::cout<<'\n';
}

int Run()
{
    constexpr int samples=5,repeats=5;
    Matrix j(4096,8);
    for(Eigen::Index row=0;row<j.rows();++row)
        for(Eigen::Index col=0;col<j.cols();++col)
        {
            const double phase=static_cast<double>(row+1)*static_cast<double>(col+1);
            j(row,col)=(std::sin(.00019*phase)+.4*std::cos(.00031*(phase+7.0)))/std::sqrt(4096.0);
        }
    Vector truth=Vector::Constant(j.cols(),.6),r=-(j*truth),d=Vector::Ones(j.cols());
    Vector upper=Vector::Constant(j.cols(),.15); upper(0)=.02;
    std::vector<double> bounded_times;
    joint::ActiveSetStep bounded;
    for(int sample=0;sample<samples+1;++sample)
    {
        const auto started=Clock::now();
        for(int repeat=0;repeat<repeats;++repeat)
        {
            bounded=joint::SolveBoundedTrustRegionStep(j,r,d,upper,.3,.5);
            if(!bounded.valid) throw std::runtime_error("bounded LM: "+bounded.status);
        }
        if(sample) bounded_times.push_back(Seconds(started)/repeats);
    }
    ReportSamples("bounded_lm",std::move(bounded_times));

    auto input=std::make_shared<rhbm_gem::core::JointProblemInput>(
        second_stage_test::SyntheticJointWorkload("chain",30));
    joint::Domain domain(input);
    const auto context=joint::CreateContext(input);
    Vector y=Eigen::Map<const Vector>(input->observations.data(),
        static_cast<Eigen::Index>(input->observations.size()));
    Vector eta=Vector::Constant(static_cast<Eigen::Index>(input->atom_ids.size()),std::log(.55));
    const auto evaluation=joint::EvaluateProfile(domain,y,eta,&context);
    if(!evaluation.valid) throw std::runtime_error("profile preparation: "+evaluation.reason);

    std::vector<double> prepare_times,reduce_times;
    joint::TiledDifferential prepared;
    joint::ReducedDifferential reduced;
    for(int sample=0;sample<samples+1;++sample)
    {
        const auto prepare_started=Clock::now();
        for(int repeat=0;repeat<repeats;++repeat)
        {
            prepared=joint::PrepareDerivative(evaluation,context.scale,&context);
            if(!prepared.valid) throw std::runtime_error("derivative preparation: "+prepared.reason);
        }
        if(sample) prepare_times.push_back(Seconds(prepare_started)/repeats);

        const auto reduce_started=Clock::now();
        for(int repeat=0;repeat<repeats;++repeat)
        {
            reduced=joint::ReduceDerivative(prepared,evaluation.residual);
            if(!reduced.valid) throw std::runtime_error("derivative reduction: "+reduced.reason);
        }
        if(sample) reduce_times.push_back(Seconds(reduce_started)/repeats);
    }
    ReportSamples("derivative_prepare",std::move(prepare_times));
    ReportSamples("derivative_reduce",std::move(reduce_times));
    std::cout<<"derivative_rows="<<evaluation.residual.size()
        <<" derivative_columns="<<evaluation.x.cols()<<'\n';
    return 0;
}
}

int main()
{
    try {return Run();}
    catch(const std::exception & error)
    {std::cerr<<"joint optimization benchmark failed: "<<error.what()<<'\n'; return 1;}
}
