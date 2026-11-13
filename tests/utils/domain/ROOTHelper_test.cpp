#include <gtest/gtest.h>

#include <stdexcept>

#include <rhbm_gem/utils/domain/ROOTHelper.hpp>
#include <rhbm_gem/utils/math/GaussianModel3D.hpp>

#ifdef HAVE_ROOT
#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TDirectoryFile.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TROOT.h>

TEST(ROOTHelperTest, HistogramsAreOwnedByTheirCallers)
{
    TDirectoryFile directory("histogram_owner_test", "");
    TDirectory::TContext context(&directory);
    TH1D existing("repeated", "", 4, 0., 2.);
    const auto one{ rhbm_gem::root_helper::CreateHist1D("repeated", "", 4, 0., 2.) };
    ASSERT_EQ(one->GetDirectory(), nullptr);
    const auto two{ rhbm_gem::root_helper::CreateHist2D("repeated", "", 4, 0., 2., 4, -2., 2.) };
    ASSERT_EQ(two->GetDirectory(), nullptr);
    one->Fill(.25);
    two->Fill(.25, 1.25);
    EXPECT_DOUBLE_EQ(one->GetBinContent(1), 1.);
    EXPECT_DOUBLE_EQ(two->GetBinContent(1, 4), 1.);
    EXPECT_EQ(directory.GetList()->GetSize(), 1);
    EXPECT_EQ(directory.GetList()->FindObject("repeated"), &existing);
}

TEST(ROOTHelperTest, CreateGaus3DFunctionIn1DIncludesOffset)
{
    gROOT->SetBatch(kTRUE);
    const rhbm_gem::GaussianModel3D model{ 1.25, 0.55, -0.15 };
    auto function{
        rhbm_gem::root_helper::CreateGaus3DFunctionIn1D(
            "gaus", model.GetAmplitude(), model.GetWidth(), model.GetOffset())
    };

    ASSERT_NE(function, nullptr);
    EXPECT_EQ(function->GetNpar(), 3);
    EXPECT_DOUBLE_EQ(function->GetParameter(2), model.GetOffset());
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
