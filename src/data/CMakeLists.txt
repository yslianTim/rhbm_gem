# File dependencies
file(GLOB SOURCE *.cpp)

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
    add_subdirectory(${CMAKE_SOURCE_DIR}/third_party/sqlite ${CMAKE_BINARY_DIR}/third_party/sqlite)
endif()

add_library(data SHARED ${SOURCE})
target_include_directories(data PUBLIC ${CMAKE_SOURCE_DIR}/include/data)
target_link_libraries(data PUBLIC
    utils
    $<$<BOOL:${USE_SYSTEM_SQLITE3}>:SQLite::SQLite3>
    $<$<NOT:$<BOOL:${USE_SYSTEM_SQLITE3}>>:sqlite3>
    CompilerFlags
)

set_target_properties(data PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
)
