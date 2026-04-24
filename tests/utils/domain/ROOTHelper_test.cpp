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
