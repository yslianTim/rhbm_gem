#include <gtest/gtest.h>
#include "support/SecondStageTestSupport.hpp"
#include "support/SecondStageNumericalProbe.hpp"
#include "core/detail/second_stage/IterationProcess.hpp"
#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <iomanip>
#include <iostream>
#include <sstream>

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
