#include <gtest/gtest.h>

#include <rhbm_gem/data/object/ModelObject.hpp>

#include <memory>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "support/PublicHeaderSurfaceTestSupport.hpp"

namespace {

template <typename T, typename = void>
struct HasPublicRebuildSelectionIndex : std::false_type {};

template <typename T>
struct HasPublicRebuildSelectionIndex<
    T,
    std::void_t<decltype(std::declval<T &>().RebuildSelectionIndex())>> : std::true_type {};

template <typename T, typename = void>
struct HasPublicApplySymmetrySelection : std::false_type {};

template <typename T>
struct HasPublicApplySymmetrySelection<
    T,
    std::void_t<decltype(std::declval<T &>().ApplySymmetrySelection(false))>> : std::true_type {};

template <typename T, typename = void>
struct HasPublicSelectedAtomList : std::false_type {};

template <typename T>
struct HasPublicSelectedAtomList<
    T,
    std::void_t<decltype(std::declval<const T &>().GetSelectedAtomList())>> : std::true_type {};

template <typename T, typename = void>
struct HasPublicSetAtomList : std::false_type {};

template <typename T>
struct HasPublicSetAtomList<
    T,
    std::void_t<decltype(
        std::declval<T &>().SetAtomList(
            std::declval<std::vector<std::unique_ptr<rhbm_gem::AtomObject>>>()))>> : std::true_type {};

} // namespace

TEST(DataPublicSurfaceTest, DataPublicHeadersMatchApprovedSurface)
{
    const std::vector<std::string> expected{
        "data/io/DataRepository.hpp",
        "data/io/ModelMapFileIO.hpp",
        "data/object/AtomObject.hpp",
        "data/object/BondObject.hpp",
        "data/object/ChemicalComponentEntry.hpp",
        "data/object/MapObject.hpp",
        "data/object/ModelObject.hpp" };

    EXPECT_EQ(contract_test_support::CollectPublicHeadersForDomain("data"), expected);
}

TEST(DataPublicSurfaceTest, ModelObjectKeepsSelectionAndBuildWorkflowPrivate)
{
    static_assert(!HasPublicRebuildSelectionIndex<rhbm_gem::ModelObject>::value);
    static_assert(!HasPublicApplySymmetrySelection<rhbm_gem::ModelObject>::value);
    static_assert(!HasPublicSelectedAtomList<rhbm_gem::ModelObject>::value);
    static_assert(!HasPublicSetAtomList<rhbm_gem::ModelObject>::value);
    SUCCEED();
}
