add_custom_target(lint_repo
    COMMAND ${CMAKE_COMMAND}
        -DPROJECT_SOURCE_DIR:PATH=${PROJECT_SOURCE_DIR}
        -DPROJECT_BINARY_DIR:PATH=${PROJECT_BINARY_DIR}
        -DRHBM_GEM_INCLUDE_INSTALL_SMOKE:BOOL=OFF
        -P ${PROJECT_SOURCE_DIR}/cmake/RepoGuards.cmake
    USES_TERMINAL
    COMMENT "Running repository lint checks (style/structure/hygiene)"
)

add_custom_target(lint_install_smoke
    COMMAND ${CMAKE_COMMAND}
        -DPROJECT_SOURCE_DIR:PATH=${PROJECT_SOURCE_DIR}
        -DPROJECT_BINARY_DIR:PATH=${PROJECT_BINARY_DIR}
        -P ${PROJECT_SOURCE_DIR}/cmake/InstallConsumerSmoke.cmake
    USES_TERMINAL
    COMMENT "Running install consumer smoke test"
)

function(rhbm_gem_add_lint_build_dependencies target_name)
    if(TARGET rhbm_gem)
        add_dependencies(${target_name} rhbm_gem)
    endif()
    if(TARGET rhbm_gem_cli)
        add_dependencies(${target_name} rhbm_gem_cli)
    endif()
    if(TARGET ${PROJECT_NAME}_gui)
        add_dependencies(${target_name} ${PROJECT_NAME}_gui)
    endif()
    if(TARGET rhbm_gem_module)
        add_dependencies(${target_name} rhbm_gem_module)
    endif()
endfunction()

rhbm_gem_add_lint_build_dependencies(lint_install_smoke)

add_custom_target(lint_all
    DEPENDS
        lint_repo
        lint_install_smoke
)
