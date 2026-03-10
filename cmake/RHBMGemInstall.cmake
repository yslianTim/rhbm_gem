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
    DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/include/rhbm_gem/"
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/rhbm_gem"
    FILES_MATCHING
    PATTERN "*.hpp"
)
install(
    FILES
        "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE"
        "${CMAKE_CURRENT_SOURCE_DIR}/THIRD_PARTY_NOTICES.md"
    DESTINATION "${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}"
)
install(
    DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/examples/python/"
    DESTINATION "${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}/examples/python"
    FILES_MATCHING
    PATTERN "*.py"
)
install(
    FILES "${CMAKE_CURRENT_SOURCE_DIR}/tests/data/test_model.cif"
    DESTINATION "${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}/examples/python/data"
)

if(RHBM_GEM_DEP_PROVIDER STREQUAL "FETCH")
    rhbm_gem_install_fetched_headers()
endif()

set(RHBM_GEM_CMAKE_INSTALL_DIR "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}")

configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/RHBM_GEMConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/RHBM_GEMConfig.cmake"
    INSTALL_DESTINATION "${RHBM_GEM_CMAKE_INSTALL_DIR}"
)

write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/RHBM_GEMConfigVersion.cmake"
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
        "${CMAKE_CURRENT_BINARY_DIR}/RHBM_GEMConfig.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/RHBM_GEMConfigVersion.cmake"
    DESTINATION "${RHBM_GEM_CMAKE_INSTALL_DIR}"
)
