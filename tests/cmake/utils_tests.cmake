add_rhbm_gtest_target(rhbm_tests_utils_math
    DOMAIN utils
    INTENTS algorithm
    SOURCES
        utils/math/ArrayStats_test.cpp
        utils/math/EigenMatrixUtility_test.cpp
        utils/math/ElectricPotential_test.cpp
        utils/math/GausLinearTransformHelper_test.cpp
        utils/math/KDTreeAlgorithm_test.cpp
        utils/math/SphereSampler_test.cpp
)

add_rhbm_gtest_target(rhbm_tests_utils_domain
    DOMAIN utils
    INTENTS validation
    SOURCES
        utils/domain/AtomSelector_test.cpp
        utils/domain/ChemicalDataHelper_test.cpp
        utils/domain/ChimeraXHelper_test.cpp
        utils/domain/ComponentHelper_test.cpp
        utils/domain/FilePathHelper_test.cpp
        utils/domain/LoggerLevelBehavior_test.cpp
        utils/domain/LoggerOutput_test.cpp
        utils/domain/LoggerProgressTimer_test.cpp
        utils/domain/ROOTHelper_test.cpp
        utils/domain/StringHelperFormatting_test.cpp
        utils/domain/StringHelperParsing_test.cpp
        utils/domain/StringHelperTokenization_test.cpp
)

add_rhbm_gtest_target(rhbm_tests_utils_hrl
    DOMAIN utils
    INTENTS algorithm
    SOURCES
        utils/hrl/HRLAlphaTrainer_test.cpp
        utils/hrl/HRLDataTransform_test.cpp
        utils/hrl/HRLGroupEstimator_test.cpp
        utils/hrl/HRLModelAlgorithms_test.cpp
)
