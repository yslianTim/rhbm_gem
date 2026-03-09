#include <gtest/gtest.h>

#include <stdexcept>

#include <rhbm_gem/utils/domain/ROOTHelper.hpp>

#ifdef HAVE_ROOT
#include <TLegend.h>
#include <TPaveText.h>
#include <TROOT.h>

TEST(ROOTHelperTest, SetPadMarginInCanvasThrowsForNullPad)
{
    gROOT->SetBatch(kTRUE);
    EXPECT_THROW(
        ROOTHelper::SetPadMarginInCanvas(nullptr, 0.1, 0.1, 0.1, 0.1),
        std::invalid_argument);
}

TEST(ROOTHelperTest, SetPaveTextMarginInCanvasThrowsForNullPave)
{
    gROOT->SetBatch(kTRUE);
    EXPECT_THROW(
        ROOTHelper::SetPaveTextMarginInCanvas(nullptr, nullptr, 0.1, 0.1, 0.1, 0.1),
        std::invalid_argument);
}

TEST(ROOTHelperTest, SetPaveTextMarginInCanvasThrowsForNullPad)
{
    gROOT->SetBatch(kTRUE);
    auto pave{ ROOTHelper::CreatePaveText(0.0, 0.0, 1.0, 1.0) };
    EXPECT_THROW(
        ROOTHelper::SetPaveTextMarginInCanvas(nullptr, pave.get(), 0.1, 0.1, 0.1, 0.1),
        std::invalid_argument);
}

TEST(ROOTHelperTest, SetLegendMarginInCanvasThrowsForNullLegend)
{
    gROOT->SetBatch(kTRUE);
    EXPECT_THROW(
        ROOTHelper::SetLegendMarginInCanvas(nullptr, nullptr, 0.1, 0.1, 0.1, 0.1),
        std::invalid_argument);
}

TEST(ROOTHelperTest, SetLegendMarginInCanvasThrowsForNullPad)
{
    gROOT->SetBatch(kTRUE);
    auto legend{ ROOTHelper::CreateLegend(0.0, 0.0, 1.0, 1.0) };
    EXPECT_THROW(
        ROOTHelper::SetLegendMarginInCanvas(nullptr, legend.get(), 0.1, 0.1, 0.1, 0.1),
        std::invalid_argument);
}

#endif
