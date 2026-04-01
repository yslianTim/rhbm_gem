#include <gtest/gtest.h>

#include <memory>
#include <stdexcept>

#include <rhbm_gem/core/painter/AtomPainter.hpp>
#include <rhbm_gem/core/painter/ComparisonPainter.hpp>
#include <rhbm_gem/core/painter/DemoPainter.hpp>
#include <rhbm_gem/core/painter/GausPainter.hpp>
#include <rhbm_gem/core/painter/ModelPainter.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/DataObjectDispatch.hpp>
#include <rhbm_gem/data/object/LocalPotentialEntry.hpp>
#include "support/DataObjectTestSupport.hpp"

namespace rg = rhbm_gem;

TEST(DataObjectDispatchTest, AsHelpersResolveModelAndMap)
{
    auto model{ data_test::MakeModelWithBond() };
    auto map{ data_test::MakeMapObject() };
    rg::AtomObject atom;

    EXPECT_EQ(rg::AsModelObject(*model), model.get());
    EXPECT_EQ(rg::AsMapObject(*model), nullptr);
    EXPECT_EQ(rg::AsMapObject(map), &map);
    EXPECT_EQ(rg::AsModelObject(map), nullptr);
    EXPECT_EQ(rg::AsModelObject(atom), nullptr);
    EXPECT_EQ(rg::AsMapObject(atom), nullptr);

    const rg::DataObjectBase & const_model_ref{ *model };
    const rg::DataObjectBase & const_map_ref{ map };
    EXPECT_EQ(rg::AsModelObject(const_model_ref), model.get());
    EXPECT_EQ(rg::AsMapObject(const_map_ref), &map);
}

TEST(DataObjectDispatchTest, ExpectHelpersResolveModelAndMap)
{
    auto model{ data_test::MakeModelWithBond() };
    auto map{ data_test::MakeMapObject() };

    const auto & model_ref{ rg::ExpectModelObject(*model, "dispatch-test-model") };
    const auto & map_ref{ rg::ExpectMapObject(map, "dispatch-test-map") };
    EXPECT_EQ(&model_ref, model.get());
    EXPECT_EQ(&map_ref, &map);
}

TEST(DataObjectDispatchTest, ExpectHelpersThrowOnUnsupportedTargetType)
{
    rg::AtomObject atom;
    EXPECT_THROW((void)rg::ExpectModelObject(atom, "dispatch-test-model"), std::runtime_error);
    EXPECT_THROW((void)rg::ExpectMapObject(atom, "dispatch-test-map"), std::runtime_error);
}

TEST(DataObjectDispatchTest, PainterTypedIngestionRequiresSupportedObjectTypes)
{
    auto model{ data_test::MakeModelWithBond() };
    rg::AtomObject atom;
    atom.SetSelectedFlag(true);
    atom.AddLocalPotentialEntry(std::make_unique<rg::LocalPotentialEntry>());

    rg::AtomPainter atom_painter;
    atom_painter.AddDataObject(&atom);
    atom_painter.AddReferenceDataObject(&atom, "ref");
    EXPECT_THROW(atom_painter.AddDataObject(model.get()), std::runtime_error);
    EXPECT_THROW(atom_painter.AddReferenceDataObject(model.get(), "ref"), std::runtime_error);

    rg::ModelPainter model_painter;
    model_painter.AddDataObject(model.get());
    model_painter.AddReferenceDataObject(model.get(), "ref");
    EXPECT_THROW(model_painter.AddDataObject(&atom), std::runtime_error);

    rg::GausPainter gaus_painter;
    gaus_painter.AddDataObject(model.get());
    gaus_painter.AddReferenceDataObject(model.get(), "ref");
    EXPECT_THROW(gaus_painter.AddDataObject(&atom), std::runtime_error);

    rg::ComparisonPainter comparison_painter;
    comparison_painter.AddDataObject(model.get());
    comparison_painter.AddReferenceDataObject(model.get(), "ref");
    EXPECT_THROW(comparison_painter.AddDataObject(&atom), std::runtime_error);

    rg::DemoPainter demo_painter;
    demo_painter.AddDataObject(model.get());
    demo_painter.AddReferenceDataObject(model.get(), "ref");
    EXPECT_THROW(demo_painter.AddDataObject(&atom), std::runtime_error);
}

TEST(DataObjectDispatchTest, PainterNullIngestThrows)
{
    rg::ModelPainter painter;
    EXPECT_THROW(painter.AddDataObject(nullptr), std::runtime_error);
    EXPECT_THROW(painter.AddReferenceDataObject(nullptr, "ref"), std::runtime_error);
}
