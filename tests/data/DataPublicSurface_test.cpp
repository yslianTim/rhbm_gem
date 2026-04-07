#include <gtest/gtest.h>

#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/BondObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>

#include <filesystem>
#include <fstream>
#include <memory>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

#include "support/PublicHeaderSurfaceTestSupport.hpp"

namespace {

template <typename T, typename = void>
struct HasPublicRebuildSelection : std::false_type {};

template <typename T>
struct HasPublicRebuildSelection<
    T,
    std::void_t<decltype(std::declval<T &>().RebuildSelection())>> : std::true_type {};

template <typename T, typename = void>
struct HasPublicApplySymmetrySelection : std::false_type {};

template <typename T>
struct HasPublicApplySymmetrySelection<
    T,
    std::void_t<decltype(std::declval<T &>().ApplySymmetrySelection(false))>> : std::true_type {};

template <typename T, typename = void>
struct HasPublicSelectedAtoms : std::false_type {};

template <typename T>
struct HasPublicSelectedAtoms<
    T,
    std::void_t<decltype(std::declval<const T &>().GetSelectedAtoms())>> : std::true_type {};

template <typename T, typename = void>
struct HasPublicSetAtomList : std::false_type {};

template <typename T>
struct HasPublicSetAtomList<
    T,
    std::void_t<decltype(
        std::declval<T &>().SetAtomList(
            std::declval<std::vector<std::unique_ptr<rhbm_gem::AtomObject>>>()))>> : std::true_type {};

template <typename T, typename = void>
struct HasPublicFindChemicalComponentEntry : std::false_type {};

template <typename T>
struct HasPublicFindChemicalComponentEntry<
    T,
    std::void_t<decltype(std::declval<const T &>().FindChemicalComponentEntry(
        std::declval<ComponentKey>()))>> : std::true_type {};

template <typename T, typename = void>
struct HasPublicHasChemicalComponentEntry : std::false_type {};

template <typename T>
struct HasPublicHasChemicalComponentEntry<
    T,
    std::void_t<decltype(std::declval<const T &>().HasChemicalComponentEntry(
        std::declval<ComponentKey>()))>> : std::true_type {};

template <typename T, typename = void>
struct HasPublicChemicalComponentEntryMap : std::false_type {};

template <typename T>
struct HasPublicChemicalComponentEntryMap<
    T,
    std::void_t<decltype(std::declval<const T &>().GetChemicalComponentEntryMap())>>
    : std::true_type {};

template <typename T, typename = void>
struct HasPublicBuildKDTreeRoot : std::false_type {};

template <typename T>
struct HasPublicBuildKDTreeRoot<
    T,
    std::void_t<decltype(std::declval<T &>().BuildKDTreeRoot())>> : std::true_type {};

template <typename T, typename = void>
struct HasPublicPeekKDTreeRoot : std::false_type {};

template <typename T>
struct HasPublicPeekKDTreeRoot<
    T,
    std::void_t<decltype(std::declval<T &>().PeekKDTreeRoot())>> : std::true_type {};

template <typename T, typename = void>
struct HasPublicGetSelectedFlag : std::false_type {};

template <typename T>
struct HasPublicGetSelectedFlag<
    T,
    std::void_t<decltype(std::declval<const T &>().GetSelectedFlag())>> : std::true_type {};

template <typename T, typename = void>
struct HasStringResidueSetter : std::false_type {};

template <typename T>
struct HasStringResidueSetter<
    T,
    std::void_t<decltype(std::declval<T &>().SetResidue(std::declval<const std::string &>()))>>
    : std::true_type {};

template <typename T, typename = void>
struct HasStringElementSetter : std::false_type {};

template <typename T>
struct HasStringElementSetter<
    T,
    std::void_t<decltype(std::declval<T &>().SetElement(std::declval<const std::string &>()))>>
    : std::true_type {};

template <typename T, typename = void>
struct HasStringSpotSetter : std::false_type {};

template <typename T>
struct HasStringSpotSetter<
    T,
    std::void_t<decltype(std::declval<T &>().SetSpot(std::declval<const std::string &>()))>>
    : std::true_type {};

template <typename T, typename = void>
struct HasAdlMutableAnalysisData : std::false_type {};

template <typename T>
struct HasAdlMutableAnalysisData<
    T,
    std::void_t<decltype(MutableAnalysisData(std::declval<T &>()))>> : std::true_type {};

template <typename T, typename = void>
struct HasAdlReadAnalysisData : std::false_type {};

template <typename T>
struct HasAdlReadAnalysisData<
    T,
    std::void_t<decltype(ReadAnalysisData(std::declval<const T &>()))>> : std::true_type {};

template <typename T, typename = void>
struct HasAdlClearAnalysisFitStates : std::false_type {};

template <typename T>
struct HasAdlClearAnalysisFitStates<
    T,
    std::void_t<decltype(ClearAnalysisFitStates(std::declval<T &>()))>> : std::true_type {};

template <typename T, typename = void>
struct HasAdlFindAtomsInRange : std::false_type {};

template <typename T>
struct HasAdlFindAtomsInRange<
    T,
    std::void_t<decltype(FindAtomsInRange(
        std::declval<T &>(),
        std::declval<const rhbm_gem::AtomObject &>(),
        std::declval<double>()))>> : std::true_type {};

template <typename T, typename = void>
struct HasAdlOwnerModelOf : std::false_type {};

template <typename T>
struct HasAdlOwnerModelOf<
    T,
    std::void_t<decltype(OwnerModelOf(std::declval<const T &>()))>> : std::true_type {};

bool SourceTreeContains(
    const std::filesystem::path & root,
    std::string_view needle,
    const std::vector<std::filesystem::path> & ignored_paths = {})
{
    for (const auto & entry : std::filesystem::recursive_directory_iterator(root))
    {
        if (!entry.is_regular_file())
        {
            continue;
        }
        if (std::find(ignored_paths.begin(), ignored_paths.end(), entry.path()) != ignored_paths.end())
        {
            continue;
        }

        const auto extension{ entry.path().extension().string() };
        if (extension != ".cpp" && extension != ".hpp")
        {
            continue;
        }

        std::ifstream stream(entry.path());
        if (!stream.is_open())
        {
            continue;
        }

        const std::string content{
            std::istreambuf_iterator<char>(stream),
            std::istreambuf_iterator<char>() };
        if (content.find(needle) != std::string::npos)
        {
            return true;
        }
    }
    return false;
}

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

TEST(DataPublicSurfaceTest, ModelObjectExposesSelectionQueriesButKeepsBuildWorkflowPrivate)
{
    static_assert(HasPublicRebuildSelection<rhbm_gem::ModelObject>::value);
    static_assert(HasPublicApplySymmetrySelection<rhbm_gem::ModelObject>::value);
    static_assert(HasPublicSelectedAtoms<rhbm_gem::ModelObject>::value);
    static_assert(HasPublicFindChemicalComponentEntry<rhbm_gem::ModelObject>::value);
    static_assert(HasPublicHasChemicalComponentEntry<rhbm_gem::ModelObject>::value);
    static_assert(!HasPublicSetAtomList<rhbm_gem::ModelObject>::value);
    static_assert(!HasPublicChemicalComponentEntryMap<rhbm_gem::ModelObject>::value);
    static_assert(!HasPublicBuildKDTreeRoot<rhbm_gem::ModelObject>::value);
    static_assert(!HasPublicPeekKDTreeRoot<rhbm_gem::ModelObject>::value);
    static_assert(!HasPublicGetSelectedFlag<rhbm_gem::AtomObject>::value);
    static_assert(!HasPublicGetSelectedFlag<rhbm_gem::BondObject>::value);
    static_assert(!HasStringResidueSetter<rhbm_gem::AtomObject>::value);
    static_assert(!HasStringElementSetter<rhbm_gem::AtomObject>::value);
    static_assert(!HasStringSpotSetter<rhbm_gem::AtomObject>::value);
    static_assert(!HasAdlMutableAnalysisData<rhbm_gem::ModelObject>::value);
    static_assert(!HasAdlReadAnalysisData<rhbm_gem::ModelObject>::value);
    static_assert(!HasAdlClearAnalysisFitStates<rhbm_gem::ModelObject>::value);
    static_assert(!HasAdlFindAtomsInRange<rhbm_gem::ModelObject>::value);
    static_assert(!HasAdlOwnerModelOf<rhbm_gem::AtomObject>::value);
    static_assert(!HasAdlOwnerModelOf<rhbm_gem::BondObject>::value);
    SUCCEED();
}

TEST(DataPublicSurfaceTest, InternalTransitionWrappersAreNotAvailableOrIncluded)
{
    const auto project_root{ command_test::ProjectRootPath() };
    const std::vector<std::filesystem::path> ignored_test_files{
        project_root / "tests/data/DataPublicSurface_test.cpp" };
    EXPECT_FALSE(std::filesystem::exists(project_root / "src/core/detail/LocalPotentialAccess.hpp"));
    EXPECT_FALSE(std::filesystem::exists(project_root / "src/core/detail/GroupPotentialAccess.hpp"));
    EXPECT_FALSE(std::filesystem::exists(project_root / "src/data/detail/ModelAnalysisAccess.hpp"));
    EXPECT_FALSE(std::filesystem::exists(project_root / "src/data/detail/ModelObjectBuilder.hpp"));
    EXPECT_FALSE(std::filesystem::exists(project_root / "src/data/detail/ModelSelectionAccess.hpp"));
    EXPECT_FALSE(SourceTreeContains(project_root / "src", "LocalPotentialAccess.hpp"));
    EXPECT_FALSE(SourceTreeContains(project_root / "src", "GroupPotentialAccess.hpp"));
    EXPECT_FALSE(SourceTreeContains(project_root / "src", "ModelAnalysisAccess"));
    EXPECT_FALSE(SourceTreeContains(project_root / "src", "ModelObjectBuilder"));
    EXPECT_FALSE(SourceTreeContains(project_root / "src", "ModelSelectionAccess.hpp"));
    EXPECT_FALSE(SourceTreeContains(project_root / "tests", "ModelAnalysisAccess", ignored_test_files));
    EXPECT_FALSE(SourceTreeContains(project_root / "tests", "ModelObjectBuilder", ignored_test_files));
    EXPECT_FALSE(SourceTreeContains(project_root / "include", "ModelAnalysisAccess"));
    EXPECT_FALSE(SourceTreeContains(project_root / "include", "ModelObjectBuilder"));
    EXPECT_FALSE(SourceTreeContains(project_root / "include", "ModelSelectionAccess"));
}
