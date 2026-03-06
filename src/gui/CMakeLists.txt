set(GUI_ADAPTER_SOURCES
    GuiCommandExecutor.cpp
)

add_library(gui_adapter ${GUI_ADAPTER_SOURCES})
target_include_directories(gui_adapter PUBLIC
    "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/gui>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/gui>"
)
target_link_libraries(gui_adapter PUBLIC core)
target_link_libraries(gui_adapter PRIVATE $<BUILD_INTERFACE:CompilerFlags>)
set_target_properties(gui_adapter PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}"
    RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}"
    LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}"
)

install(
    TARGETS gui_adapter
    EXPORT RHBM_GEMTargets
    RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
    LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
)

if(RHBM_GEM_BUILD_GUI)
    find_package(Qt6 REQUIRED COMPONENTS Core Widgets)

    set(GUI_SOURCES
        main.cpp
        MainWindow.cpp
    )

    set(GUI_EXEC_TARGET ${PROJECT_NAME}_gui)
    add_executable(${GUI_EXEC_TARGET} ${GUI_SOURCES})
    set_target_properties(${GUI_EXEC_TARGET} PROPERTIES
        OUTPUT_NAME "${PROJECT_DISPLAY_NAME}-GUI"
        RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}"
    )

    target_link_libraries(${GUI_EXEC_TARGET} PRIVATE
        gui_adapter
        Qt6::Core
        Qt6::Widgets
        CompilerFlags
    )
    target_include_directories(${GUI_EXEC_TARGET} PRIVATE
        "${PROJECT_SOURCE_DIR}/src/gui"
    )

    install(
        TARGETS ${GUI_EXEC_TARGET}
        RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
    )
endif()
