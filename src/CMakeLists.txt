include(${CMAKE_CURRENT_SOURCE_DIR}/rhbm_gem_sources.cmake)

add_library(rhbm_gem ${RHBM_GEM_LIBRARY_SOURCES})
target_compile_features(rhbm_gem PUBLIC cxx_std_17)

target_include_directories(rhbm_gem PUBLIC
    "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)
target_include_directories(rhbm_gem PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/core"
    "${CMAKE_CURRENT_SOURCE_DIR}/core/internal"
    "${CMAKE_CURRENT_SOURCE_DIR}/data"
    "${CMAKE_CURRENT_SOURCE_DIR}/data/internal"
    "${CMAKE_CURRENT_SOURCE_DIR}/utils"
)

target_link_libraries(rhbm_gem PUBLIC rhbm_gem_dependencies)
target_link_libraries(rhbm_gem PRIVATE CompilerFlags)

if(WIN32 AND BUILD_SHARED_LIBS)
    set_target_properties(rhbm_gem PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS ON)
endif()
rhbm_gem_apply_library_install_rpath(rhbm_gem)

if(RHBM_GEM_ENABLE_EXPERIMENTAL_BOND_ANALYSIS)
    target_compile_definitions(rhbm_gem PRIVATE RHBM_GEM_ENABLE_EXPERIMENTAL_BOND_ANALYSIS)
endif()
if(RHBM_GEM_LEGACY_V1_SUPPORT)
    target_compile_definitions(rhbm_gem PRIVATE RHBM_GEM_LEGACY_V1_SUPPORT)
endif()

install(
    TARGETS rhbm_gem rhbm_gem_dependencies
    EXPORT RHBM_GEMTargets
    RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
    LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    INCLUDES DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
)

# Avoid target-directory collisions on case-insensitive filesystems
# (e.g. `rhbm_gem` library vs `RHBM_GEM` executable).
set(MAIN_EXEC_TARGET rhbm_gem_cli)
add_executable(${MAIN_EXEC_TARGET} main.cpp)
set_target_properties(${MAIN_EXEC_TARGET} PROPERTIES
    OUTPUT_NAME "${PROJECT_DISPLAY_NAME}"
)
target_link_libraries(${MAIN_EXEC_TARGET} PRIVATE
    rhbm_gem
    CompilerFlags
)
target_include_directories(${MAIN_EXEC_TARGET} PRIVATE "${PROJECT_BINARY_DIR}")
rhbm_gem_apply_executable_install_rpath(${MAIN_EXEC_TARGET})

install(
    TARGETS ${MAIN_EXEC_TARGET}
    RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
)

if(TARGET rhbm_gem_gui_dependencies)
    set(GUI_EXEC_TARGET ${PROJECT_NAME}_gui)
    add_executable(${GUI_EXEC_TARGET}
        gui/main.cpp
        gui/MainWindow.cpp
    )
    set_target_properties(${GUI_EXEC_TARGET} PROPERTIES
        OUTPUT_NAME "${PROJECT_DISPLAY_NAME}-GUI"
    )
    target_link_libraries(${GUI_EXEC_TARGET} PRIVATE
        rhbm_gem
        rhbm_gem_gui_dependencies
        CompilerFlags
    )
    target_include_directories(${GUI_EXEC_TARGET} PRIVATE
        "${PROJECT_SOURCE_DIR}/src/gui"
    )
    rhbm_gem_apply_executable_install_rpath(${GUI_EXEC_TARGET})

    install(
        TARGETS ${GUI_EXEC_TARGET}
        RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
    )
endif()
