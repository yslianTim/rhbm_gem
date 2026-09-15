file(STRINGS "${INPUT_LIST}" _sources)
set(_source_hashes "")
foreach(_source IN LISTS _sources)
    file(RELATIVE_PATH _relative "${SOURCE_ROOT}" "${_source}")
    file(SHA256 "${_source}" _hash)
    string(APPEND _source_hashes "${_relative}\n${_hash}\n")
endforeach()
string(SHA256 _source_sha256 "${_source_hashes}")
file(SHA256 "${CONFIG_FILE}" _config_sha256)
string(SHA256 _build_sha256 "${_source_sha256}\n${_config_sha256}\n")
file(WRITE "${OUTPUT_FILE}.tmp" "#pragma once
#define RHBM_GEM_SIMULATION_VERSION \"${PROJECT_VERSION}\"
#define RHBM_GEM_SIMULATION_SOURCE_SHA256 \"${_source_sha256}\"
#define RHBM_GEM_SIMULATION_CONFIG_SHA256 \"${_config_sha256}\"
#define RHBM_GEM_SIMULATION_BUILD_SHA256 \"${_build_sha256}\"
")
configure_file("${OUTPUT_FILE}.tmp" "${OUTPUT_FILE}" COPYONLY)
file(REMOVE "${OUTPUT_FILE}.tmp")
