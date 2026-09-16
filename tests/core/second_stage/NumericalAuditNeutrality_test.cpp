#include <gtest/gtest.h>
#include "support/SecondStageTestSupport.hpp"
#include "support/SecondStageNumericalProbe.hpp"
#include "support/EndpointRefinementExperiment.hpp"
#include "core/detail/second_stage/IterationProcess.hpp"
#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <filesystem>
#include <fstream>
#include <chrono>
#include <cstdlib>
#include <set>

TEST(SecondStageNumericalProbe, CapturesSmallProductionRuns)
{
    using namespace second_stage_test;
    namespace rg=rhbm_gem;
    const auto saved=Logger::GetLogLevel();
    for(int fixture=0;fixture<3;++fixture)
    for(int threads:{1,4})
    for(bool quiet:{false,true})
    {
        auto model=fixture==0 ? BuildJointPolishDefenseModel() :
            fixture==1 ? BuildIndependentOffsetDefenseModel() : BuildSeparatedRollbackDefenseModel();
        model->EditAnalysis().CopyLocalFittingStageResult(rg::FittingStage::Second,rg::FittingStage::First);
        auto options=MakeSecondStageOptions(); options.thread_size=threads; options.quiet_mode=quiet;
        Logger::SetLogLevel(LogLevel::Debug);
        testing::internal::CaptureStdout(); BeginNumericalCapture();
        rg::core::detail::RunSecondStageIterations(*model,options);
        const auto capture=EndNumericalCapture();
        const auto log=testing::internal::GetCapturedStdout();
        for(auto count:capture.work) EXPECT_GT(count,0U);
        std::ostringstream result; result << std::setprecision(17);
        result << "NUMERICAL " << fixture << ' ' << threads << ' ' << quiet << " WORK";
        for(auto count:capture.work) result << ' ' << count;
        result << " COMMITS"; for(const auto & commit:capture.commits) result << ' ' << std::quoted(commit);
        const auto summary=log.find("Second-Stage Local Fitting Summary");
        if(!quiet) { ASSERT_NE(summary,std::string::npos); result << " SUMMARY " << std::quoted(log.substr(summary,log.find("Second-Stage Local Fitting Performance",summary)-summary)); }
        result << " TERMINAL " << std::quoted(capture.terminal) << " BACKGROUNDS";
        for (const auto & background : capture.backgrounds)
            for (const auto & parameters : background) for (auto value : parameters) result << ' ' << value;
        result << " STATE";
        for(const auto * atom:model->GetSelectedAtoms())
        {
            const auto view=rg::AtomLocalPotentialView::For(*atom);
            const auto local=view.GetGaussianResult(rg::FittingStage::Second);
            for(const auto * fit:{&local.ols,&local.mdpde})
                for(int par=0;par<3;++par) result << ' ' << fit->GetModelParameter(par) << ' ' << fit->GetModelStandardDeviation(par);
            result << ' ' << view.GetAlphaR(rg::FittingStage::Second);
            for(const auto & sample:view.GetPeelingSamplingEntries(false))
            {
                result << ' ' << sample.response << ' ' << sample.point.distance << ' ' << sample.point.is_selected;
                for(auto value:sample.point.position) result << ' ' << value;
            }
        }
        std::istringstream audit_lines{log}; std::string audit_line;
        while (std::getline(audit_lines, audit_line))
            if (const auto marker=audit_line.find("Second-stage audit: schema="); marker!=std::string::npos)
                std::cout << audit_line.substr(marker) << '\n';
        // One physical output line for stable cross-build comparison.
        auto text=result.str(); for(auto & c:text) if(c=='\n'||c=='\r') c=' ';
        std::cout << text << '\n';
    }
    Logger::SetLogLevel(saved);
}

TEST(SecondStageNumericalProbe, EndpointComparisonIsReadOnlyAndReachesEveryWorker)
{
    using namespace second_stage_test;
    namespace rg=rhbm_gem;
    const auto root{std::filesystem::temp_directory_path()/
        ("rhbm-endpoint-"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()))};
    struct Environment
    {
        std::map<std::string,std::optional<std::string>> saved;
        void Set(const std::string & key,const std::string & value)
        {
            if (!saved.contains(key)) saved[key]=std::getenv(key.c_str()) ? std::optional<std::string>{std::getenv(key.c_str())} : std::nullopt;
            setenv(key.c_str(),value.c_str(),1);
        }
        ~Environment() { for (const auto & [key,value]:saved) { if(value) setenv(key.c_str(),value->c_str(),1); else unsetenv(key.c_str()); } }
    } env;
    const auto saved_level{Logger::GetLogLevel()}; Logger::SetLogLevel(LogLevel::Info);
    for (int threads : {1,4})
    {
        auto execute=[&](bool probe) {
            auto model=BuildSeparatedRollbackDefenseModel();
            model->EditAnalysis().CopyLocalFittingStageResult(rg::FittingStage::Second,rg::FittingStage::First);
            auto options=MakeSecondStageOptions(); options.thread_size=threads; options.quiet_mode=true;
            env.Set("RHBM_TEST_ENDPOINT_DIR",probe ? (root/std::to_string(threads)).string() : "");
            env.Set("RHBM_TEST_ENDPOINT_POLICY","legacy"); env.Set("RHBM_TEST_ENDPOINT_BUDGET","128");
            env.Set("RHBM_TEST_ENDPOINT_COMPARE",probe ? "1" : "0");
            BeginNumericalCapture(); rg::core::detail::RunSecondStageIterations(*model,options);
            auto numerical=EndNumericalCapture();
            std::vector<double> state;
            for (const auto * atom:model->GetSelectedAtoms())
            {
                const auto view=rg::AtomLocalPotentialView::For(*atom);
                const auto fit=view.GetGaussianResult(rg::FittingStage::Second);
                for(int k=0;k<3;++k) state.push_back(fit.mdpde.GetModelParameter(k));
                for(const auto & sample:view.GetPeelingSamplingEntries(false)) state.push_back(sample.response);
            }
            return std::pair{numerical,state};
        };
        const auto baseline{execute(false)}, comparison{execute(true)};
        EXPECT_EQ(baseline.first.work,comparison.first.work);
        EXPECT_EQ(baseline.first.commits,comparison.first.commits);
        EXPECT_EQ(baseline.first.terminal,comparison.first.terminal);
        EXPECT_EQ(baseline.first.backgrounds,comparison.first.backgrounds);
        EXPECT_EQ(baseline.second,comparison.second);
        const auto read=[&](const auto & path) { std::ifstream in(path); return boost::json::parse(std::string{std::istreambuf_iterator<char>(in),{}}); };
        const auto report{read(root/std::to_string(threads)/"final/comparison.json")};
        EXPECT_TRUE(report.at("input_unchanged").as_bool());
        for(const auto & variant:report.at("variants").as_array())
        {
            EXPECT_TRUE(variant.at("offset_solves_equal").as_bool());
            if(variant.at("policy")=="legacy") EXPECT_TRUE(variant.at("legacy_exact").as_bool());
            std::set<std::size_t> atoms;
            for(const auto & solve:variant.at("solves").as_array())
            {
                EXPECT_EQ(solve.at("policy"),variant.at("policy"));
                for(const auto & index:solve.at("indices").as_array()) atoms.insert(boost::json::value_to<std::size_t>(index));
            }
            EXPECT_EQ(atoms.size(),4u);
        }
        EXPECT_TRUE(read(root/std::to_string(threads)/"session.json").at("complete").as_bool());
    }
    Logger::SetLogLevel(saved_level);
    std::filesystem::remove_all(root);
}
