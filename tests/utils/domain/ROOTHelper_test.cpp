#include <gtest/gtest.h>

#include <stdexcept>

#include <rhbm_gem/utils/domain/ROOTHelper.hpp>
#include <rhbm_gem/utils/math/GaussianModel3D.hpp>

#ifdef HAVE_ROOT
#include <TF1.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TROOT.h>

TEST(ROOTHelperTest, CreateGaus3DFunctionIn1DIncludesIntercept)
{
    gROOT->SetBatch(kTRUE);
    const rhbm_gem::GaussianModel3D model{ 1.25, 0.55, -0.15 };
    auto function{
        rhbm_gem::root_helper::CreateGaus3DFunctionIn1D(
            "gaus", model.GetAmplitude(), model.GetWidth(), model.GetIntercept())
    };

    ASSERT_NE(function, nullptr);
    EXPECT_EQ(function->GetNpar(), 3);
    EXPECT_DOUBLE_EQ(function->GetParameter(2), model.GetIntercept());
    EXPECT_NEAR(function->Eval(0.75), model.ResponseAtDistance(0.75), 1e-12);
}

TEST(ROOTHelperTest, SetPadMarginInCanvasThrowsForNullPad)
{
    gROOT->SetBatch(kTRUE);
    EXPECT_THROW(
        rhbm_gem::root_helper::SetPadMarginInCanvas(nullptr, 0.1, 0.1, 0.1, 0.1),
        std::invalid_argument);
}

TEST(ROOTHelperTest, SetPaveTextMarginInCanvasThrowsForNullPave)
{
    gROOT->SetBatch(kTRUE);
    EXPECT_THROW(
        rhbm_gem::root_helper::SetPaveTextMarginInCanvas(nullptr, nullptr, 0.1, 0.1, 0.1, 0.1),
        std::invalid_argument);
}

TEST(ROOTHelperTest, SetPaveTextMarginInCanvasThrowsForNullPad)
{
    gROOT->SetBatch(kTRUE);
    auto pave{ rhbm_gem::root_helper::CreatePaveText(0.0, 0.0, 1.0, 1.0) };
    EXPECT_THROW(
        rhbm_gem::root_helper::SetPaveTextMarginInCanvas(nullptr, pave.get(), 0.1, 0.1, 0.1, 0.1),
        std::invalid_argument);
}

TEST(ROOTHelperTest, SetLegendMarginInCanvasThrowsForNullLegend)
{
    gROOT->SetBatch(kTRUE);
    EXPECT_THROW(
        rhbm_gem::root_helper::SetLegendMarginInCanvas(nullptr, nullptr, 0.1, 0.1, 0.1, 0.1),
        std::invalid_argument);
}

TEST(ROOTHelperTest, SetLegendMarginInCanvasThrowsForNullPad)
{
    gROOT->SetBatch(kTRUE);
    auto legend{ rhbm_gem::root_helper::CreateLegend(0.0, 0.0, 1.0, 1.0) };
    EXPECT_THROW(
        rhbm_gem::root_helper::SetLegendMarginInCanvas(nullptr, legend.get(), 0.1, 0.1, 0.1, 0.1),
        std::invalid_argument);
}

#endif
