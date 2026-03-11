add_rhbm_gtest_target(rhbm_tests_data_io
    DOMAIN data
    INTENTS io contract
    SOURCES
        data/contract/PublicHeaderSurface_test.cpp
        data/io/DataObjectDispatchIterationArchitecture_test.cpp
        data/io/DataObjectManager_test.cpp
        data/io/MapFileIO_test.cpp
        data/schema/DataObjectIOContract_test.cpp
)

add_rhbm_gtest_target(rhbm_tests_data_schema
    DOMAIN data
    INTENTS schema
    SOURCES
        data/schema/DataObjectSchemaBootstrap_test.cpp
        data/schema/DataObjectSchemaValidation_test.cpp
)

if(RHBM_GEM_LEGACY_V1_SUPPORT)
    add_rhbm_gtest_target(rhbm_tests_data_migration
        DOMAIN data
        INTENTS migration schema
        SOURCES
            data/schema/DataObjectSchemaMigration_test.cpp
    )
endif()
