function(rhbm_gem_install_fetched_headers)
    install(
        DIRECTORY "${RHBM_GEM_EIGEN_INCLUDE_DIR}/Eigen"
        DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
    )
    if(EXISTS "${RHBM_GEM_EIGEN_INCLUDE_DIR}/unsupported")
        install(
            DIRECTORY "${RHBM_GEM_EIGEN_INCLUDE_DIR}/unsupported"
            DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
        )
    endif()

    install(
        DIRECTORY "${RHBM_GEM_CLI11_INCLUDE_DIR}/CLI"
        DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
    )

    install(
        FILES "${RHBM_GEM_SQLITE3_SOURCE_DIR}/sqlite3.h"
        DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
    )

    install(
        DIRECTORY "${RHBM_GEM_BOOST_INCLUDE_DIR}/boost"
        DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
    )
endfunction()

# Install public headers and legal notices.
install(
    DIRECTORY "${PROJECT_SOURCE_DIR}/include/rhbm_gem/"
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/rhbm_gem"
    FILES_MATCHING
    PATTERN "*.hpp"
    PATTERN "*.def"
)
install(
    FILES
        "${PROJECT_SOURCE_DIR}/LICENSE"
        "${PROJECT_SOURCE_DIR}/THIRD_PARTY_NOTICES.md"
    DESTINATION "${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}"
)
install(
    PROGRAMS
        "${PROJECT_SOURCE_DIR}/resources/examples/cli/00_quickstart.sh"
        "${PROJECT_SOURCE_DIR}/resources/examples/cli/01_estimate_three_examples.sh"
        "${PROJECT_SOURCE_DIR}/resources/examples/cli/common.sh"
    DESTINATION "${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}/resources/examples/cli"
)
install(
    DIRECTORY "${PROJECT_SOURCE_DIR}/resources/examples/python/"
    DESTINATION "${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}/resources/examples/python"
    FILES_MATCHING
    PATTERN "*.py"
)
install(
    FILES "${PROJECT_SOURCE_DIR}/tests/fixtures/test_model.cif"
    DESTINATION "${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}/resources/examples/python/data"
)

# Temporary compatibility alias for the legacy installed example path.
install(
    DIRECTORY "${PROJECT_SOURCE_DIR}/resources/examples/python/"
    DESTINATION "${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}/examples/python"
    FILES_MATCHING
    PATTERN "*.py"
)
install(
    FILES "${PROJECT_SOURCE_DIR}/tests/fixtures/test_model.cif"
    DESTINATION "${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}/examples/python/data"
)

if(RHBM_GEM_DEP_PROVIDER STREQUAL "FETCH")
    rhbm_gem_install_fetched_headers()
endif()

set(RHBM_GEM_CMAKE_INSTALL_DIR "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}")

configure_package_config_file(
    "${PROJECT_SOURCE_DIR}/cmake/RHBM_GEMConfig.cmake.in"
    "${PROJECT_BINARY_DIR}/RHBM_GEMConfig.cmake"
    INSTALL_DESTINATION "${RHBM_GEM_CMAKE_INSTALL_DIR}"
)

write_basic_package_version_file(
    "${PROJECT_BINARY_DIR}/RHBM_GEMConfigVersion.cmake"
    VERSION "${PROJECT_VERSION}"
    COMPATIBILITY SameMajorVersion
)

install(
    EXPORT RHBM_GEMTargets
    FILE RHBM_GEMTargets.cmake
    NAMESPACE RHBM_GEM::
    DESTINATION "${RHBM_GEM_CMAKE_INSTALL_DIR}"
)

install(
    FILES
        "${PROJECT_BINARY_DIR}/RHBM_GEMConfig.cmake"
        "${PROJECT_BINARY_DIR}/RHBM_GEMConfigVersion.cmake"
    DESTINATION "${RHBM_GEM_CMAKE_INSTALL_DIR}"
)
