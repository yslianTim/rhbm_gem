#include <gtest/gtest.h>

#include <memory>
#include <type_traits>

#include <rhbm_gem/core/painter/AtomPainter.hpp>
#include <rhbm_gem/core/painter/ComparisonPainter.hpp>
#include <rhbm_gem/core/painter/DemoPainter.hpp>
#include <rhbm_gem/core/painter/GausPainter.hpp>
#include <rhbm_gem/core/painter/ModelPainter.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include "data/detail/LocalPotentialEntry.hpp"
#include "support/DataObjectTestSupport.hpp"

namespace rg = rhbm_gem;

TEST(DataObjectPainterIngestionTest, PainterTypedIngestionUsesConcreteModelApis)
{
    auto model{ data_test::MakeModelWithBond() };
    model->SetKeyTag("model");
    model->GetAtomList().front()->SetSelectedFlag(true);
    model->GetAtomList().front()->SetLocalPotentialEntry(
        std::make_unique<rg::LocalPotentialEntry>());

    rg::AtomPainter atom_painter;
    EXPECT_NO_THROW(atom_painter.AddModel(*model));

    rg::ModelPainter model_painter;
    EXPECT_NO_THROW(model_painter.AddModel(*model));

    rg::GausPainter gaus_painter;
    EXPECT_NO_THROW(gaus_painter.AddModel(*model));

    rg::ComparisonPainter comparison_painter;
    EXPECT_NO_THROW(comparison_painter.AddModel(*model));
    EXPECT_NO_THROW(comparison_painter.AddReferenceModel(*model, "ref"));

    rg::DemoPainter demo_painter;
    EXPECT_NO_THROW(demo_painter.AddModel(*model));
    EXPECT_NO_THROW(demo_painter.AddReferenceModel(*model, "ref"));
}

TEST(DataObjectPainterIngestionTest, PainterInterfacesAreConcreteAndTyped)
{
    static_assert(std::is_invocable_v<decltype(&rg::AtomPainter::AddModel), rg::AtomPainter &, rg::ModelObject &>);
    static_assert(std::is_invocable_v<decltype(&rg::ModelPainter::AddModel), rg::ModelPainter &, rg::ModelObject &>);
    static_assert(std::is_invocable_v<decltype(&rg::GausPainter::AddModel), rg::GausPainter &, rg::ModelObject &>);
    static_assert(std::is_invocable_v<decltype(&rg::ComparisonPainter::AddModel), rg::ComparisonPainter &, rg::ModelObject &>);
    static_assert(std::is_invocable_v<decltype(&rg::ComparisonPainter::AddReferenceModel), rg::ComparisonPainter &, rg::ModelObject &, std::string_view>);
    static_assert(std::is_invocable_v<decltype(&rg::DemoPainter::AddModel), rg::DemoPainter &, rg::ModelObject &>);
    static_assert(std::is_invocable_v<decltype(&rg::DemoPainter::AddReferenceModel), rg::DemoPainter &, rg::ModelObject &, std::string_view>);

    SUCCEED();
}
