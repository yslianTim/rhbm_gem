set(DATA_SOURCES
    AtomClassifier.cpp
    AtomObject.cpp
    AtomicModelDataBlock.cpp
    BondClassifier.cpp
    BondObject.cpp
    CCP4Format.cpp
    ChemicalComponentEntry.cpp
    CifFormat.cpp
    DataObjectDAOFactoryRegistry.cpp
    DatabaseSchemaManager.cpp
    DatabaseManager.cpp
    FileFormatBackendFactory.cpp
    FileFormatRegistry.cpp
    FileProcessFactoryRegistry.cpp
    GroupPotentialEntry.cpp
    LegacyModelObjectDAO.cpp
    LocalPotentialEntry.cpp
    MapFileReader.cpp
    MapFileWriter.cpp
    MapObject.cpp
    MapObjectDAO.cpp
    MapObjectFactory.cpp
    ModelFileReader.cpp
    ModelFileWriter.cpp
    ModelObject.cpp
    ModelObjectDAO.cpp
    ModelObjectDAOv2.cpp
    ModelObjectFactory.cpp
    model_io/ModelAnalysisPersistence.cpp
    model_io/ModelStructurePersistence.cpp
    MrcFormat.cpp
    PdbFormat.cpp
    PotentialEntryIterator.cpp
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
    "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/data>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/data>"
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
