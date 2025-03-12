# File dependencies
file(GLOB SOURCE *.cpp)

# Locate SQLite package
#find_package(SQLite3 REQUIRED)
find_package(SQLite3 QUIET)
if(SQLite3_FOUND)
    message(STATUS "Using system SQLite3")
    set(USE_SYSTEM_SQLITE3 TRUE)
else()
    message(STATUS "System SQLite3 not found, using third_party version")
    set(USE_SYSTEM_SQLITE3 FALSE)
    add_subdirectory(third_party/sqlite3)
endif()

add_library(data SHARED ${SOURCE})
target_include_directories(data PUBLIC ${CMAKE_SOURCE_DIR}/include/data)
target_link_libraries(data PUBLIC
    utils
    $<$<BOOL:${USE_SYSTEM_SQLITE3}>:SQLite::SQLite3>
    $<$<NOT:$<BOOL:${USE_SYSTEM_SQLITE3}>>:sqlite3>
)