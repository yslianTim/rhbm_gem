set(DATA_OBJECT_SOURCES
    dispatch/DataObjectDispatch.cpp
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
)

set(DATA_IO_FILE_SOURCES
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
)

set(DATA_IO_SQLITE_SOURCES
    io/sqlite/DataObjectDAOFactoryRegistry.cpp
    io/sqlite/DatabaseManager.cpp
    io/sqlite/ManagedStoreRegistry.cpp
    io/sqlite/MapObjectDAO.cpp
    io/sqlite/ModelAnalysisPersistence.cpp
    io/sqlite/ModelObjectDAO.cpp
    io/sqlite/ModelObjectDaoSqlite.cpp
    io/sqlite/ModelStructurePersistence.cpp
)

set(DATA_SCHEMA_MIGRATION_SOURCES
    schema/DatabaseSchemaManager.cpp
)

if(RHBM_GEM_LEGACY_V1_SUPPORT)
    list(APPEND DATA_SCHEMA_MIGRATION_SOURCES
        migration/legacy_v1/LegacyModelObjectReader.cpp
    )
endif()

set(DATA_FRONTDOOR_SOURCES
    io/DataObjectManager.cpp
)

# Locate SQLite package
if(USE_SYSTEM_LIBS)
    find_package(SQLite3 QUIET)
    if(SQLite3_FOUND)
        message(STATUS "Using system SQLite3")
        set(USE_SYSTEM_SQLITE3 TRUE)
    else()
        message(STATUS "System SQLite3 not found, will fetch via FetchContent")
        set(USE_SYSTEM_SQLITE3 FALSE)
    endif()
else()
    message(STATUS "Forced to use FetchContent for SQLite3")
    set(USE_SYSTEM_SQLITE3 FALSE)
endif()

if(NOT USE_SYSTEM_SQLITE3)
    rhbm_gem_fetch_sqlite3(RHBM_GEM_SQLITE3_SOURCE_DIR)

    if(NOT EXISTS "${RHBM_GEM_SQLITE3_SOURCE_DIR}/sqlite3.c"
        OR NOT EXISTS "${RHBM_GEM_SQLITE3_SOURCE_DIR}/sqlite3.h")
        message(FATAL_ERROR
            "Fetched SQLite3 source at '${RHBM_GEM_SQLITE3_SOURCE_DIR}' "
            "does not contain sqlite3.c/sqlite3.h.")
    endif()

    if(NOT TARGET sqlite3)
        add_library(sqlite3 STATIC "${RHBM_GEM_SQLITE3_SOURCE_DIR}/sqlite3.c")
        target_include_directories(sqlite3 PUBLIC
            "$<BUILD_INTERFACE:${RHBM_GEM_SQLITE3_SOURCE_DIR}>"
            "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
        )
        set_target_properties(sqlite3 PROPERTIES POSITION_INDEPENDENT_CODE ON)

        install(
            TARGETS sqlite3
            EXPORT RHBM_GEMTargets
            ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
            LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
            RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
        )
        install(
            FILES "${RHBM_GEM_SQLITE3_SOURCE_DIR}/sqlite3.h"
            DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
        )
    endif()

    if(NOT TARGET RHBM_GEM::sqlite3)
        add_library(RHBM_GEM::sqlite3 ALIAS sqlite3)
    endif()
endif()

function(rhbm_configure_data_object_target target_name)
    target_include_directories(${target_name} PRIVATE
        "${PROJECT_SOURCE_DIR}/include"
        "${CMAKE_CURRENT_SOURCE_DIR}"
        "${CMAKE_CURRENT_SOURCE_DIR}/internal"
    )
    target_link_libraries(${target_name} PRIVATE
        utils
        $<$<BOOL:${USE_SYSTEM_SQLITE3}>:SQLite::SQLite3>
        $<$<NOT:$<BOOL:${USE_SYSTEM_SQLITE3}>>:RHBM_GEM::sqlite3>
        $<BUILD_INTERFACE:CompilerFlags>
    )
    if(RHBM_GEM_LEGACY_V1_SUPPORT)
        target_compile_definitions(${target_name} PRIVATE RHBM_GEM_LEGACY_V1_SUPPORT)
    endif()
endfunction()

add_library(data_object OBJECT ${DATA_OBJECT_SOURCES})
rhbm_configure_data_object_target(data_object)

add_library(data_io_file OBJECT ${DATA_IO_FILE_SOURCES})
rhbm_configure_data_object_target(data_io_file)

add_library(data_io_sqlite OBJECT ${DATA_IO_SQLITE_SOURCES})
rhbm_configure_data_object_target(data_io_sqlite)

add_library(data_schema_migration OBJECT ${DATA_SCHEMA_MIGRATION_SOURCES})
rhbm_configure_data_object_target(data_schema_migration)

add_library(data
    ${DATA_FRONTDOOR_SOURCES}
    $<TARGET_OBJECTS:data_object>
    $<TARGET_OBJECTS:data_io_file>
    $<TARGET_OBJECTS:data_io_sqlite>
    $<TARGET_OBJECTS:data_schema_migration>
)
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
target_link_libraries(data PRIVATE
    data_object
    data_io_file
    data_io_sqlite
    data_schema_migration
    $<BUILD_INTERFACE:CompilerFlags>
)
if(RHBM_GEM_LEGACY_V1_SUPPORT)
    target_compile_definitions(data PRIVATE RHBM_GEM_LEGACY_V1_SUPPORT)
endif()

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
