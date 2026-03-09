set(DATA_SOURCES
    dispatch/DataObjectDispatch.cpp
    io/DataObjectManager.cpp
    io/file/CCP4Format.cpp
    io/file/CifFormat.cpp
    io/file/FileFormatBackendFactory.cpp
    io/file/FileFormatRegistry.cpp
    io/file/FileProcessFactoryResolver.cpp
    io/file/MapAxisOrderHelper.cpp
    io/file/MapFileReader.cpp
    io/file/MapFileWriter.cpp
    io/file/MapObjectFactory.cpp
    io/file/ModelFileReader.cpp
    io/file/ModelFileWriter.cpp
    io/file/ModelObjectFactory.cpp
    io/file/MrcFormat.cpp
    io/file/PdbFormat.cpp
    io/sqlite/DataObjectDAOFactoryRegistry.cpp
    io/sqlite/DatabaseManager.cpp
    io/sqlite/ManagedStoreRegistry.cpp
    io/sqlite/MapObjectDAO.cpp
    io/sqlite/ModelAnalysisPersistence.cpp
    io/sqlite/ModelObjectDAO.cpp
    io/sqlite/ModelObjectDaoSqlite.cpp
    io/sqlite/ModelStructurePersistence.cpp
    migration/legacy_v1/LegacyModelObjectReader.cpp
    object/AtomClassifier.cpp
    object/AtomObject.cpp
    object/AtomicModelDataBlock.cpp
    object/BondClassifier.cpp
    object/BondObject.cpp
    object/ChemicalComponentEntry.cpp
    object/GroupPotentialEntry.cpp
    object/LocalPotentialEntry.cpp
    object/MapObject.cpp
    object/ModelObject.cpp
    object/PotentialEntryIterator.cpp
    schema/DatabaseSchemaManager.cpp
)

# Locate SQLite package
if(USE_SYSTEM_LIBS)
    find_package(SQLite3)
    if(SQLite3_FOUND)
        message(STATUS "Using system SQLite3")
        set(USE_SYSTEM_SQLITE3 TRUE)
    else()
        message(STATUS "System SQLite3 not found, will use third_party version")
        set(USE_SYSTEM_SQLITE3 FALSE)
    endif()
else()
    message(STATUS "Forced to use SQLite3 package in third_party libraries")
    set(USE_SYSTEM_SQLITE3 FALSE)
endif()

if(NOT USE_SYSTEM_SQLITE3)
    add_subdirectory(${PROJECT_SOURCE_DIR}/third_party/sqlite ${PROJECT_BINARY_DIR}/third_party/sqlite)
    if(NOT TARGET RHBM_GEM::sqlite3)
        add_library(RHBM_GEM::sqlite3 ALIAS sqlite3)
    endif()
endif()

add_library(data ${DATA_SOURCES})
target_include_directories(data PUBLIC
    "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)
target_include_directories(data PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/internal"
)
target_link_libraries(data PUBLIC
    utils
    $<$<BOOL:${USE_SYSTEM_SQLITE3}>:SQLite::SQLite3>
    $<$<NOT:$<BOOL:${USE_SYSTEM_SQLITE3}>>:RHBM_GEM::sqlite3>
)
target_link_libraries(data PRIVATE $<BUILD_INTERFACE:CompilerFlags>)

set_target_properties(data PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}"
    RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}"
    LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}"
)

install(
    TARGETS data
    EXPORT RHBM_GEMTargets
    RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
    LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
)

set(RHBM_GEM_USE_SYSTEM_SQLITE3 "${USE_SYSTEM_SQLITE3}" CACHE INTERNAL "RHBM_GEM uses system SQLite3" FORCE)
