#
# Just load the env and then run this with:
#
#   env CTEST_BUILD_CONFIGURATION_NAME=<build-config-name> \
#     ctest -V -S ctest_std_driver.cmake
#
# The default CTEST_BUILD_NAME will be set to
# $ENV{CTEST_BUILD_CONFIGURATION_NAME}.  Otherwise, it can be set with
# CTEST_BUILD_NAME=<build-name> in the env.
#

INCLUDE("${CTEST_SCRIPT_DIRECTORY}/TrilinosCTestDriverCore.sems.cmake")

SET(STD_CONFIG_FILE "${TRIBITS_PROJECT_ROOT}/cmake/std/$ENV{CTEST_BUILD_CONFIGURATION_NAME}.cmake")

MESSAGE("\nIncluding the file '${STD_CONFIG_FILE}' ...") 
INCLUDE("${STD_CONFIG_FILE}")

SET_DEFAULT_AND_FROM_ENV(CTEST_BUILD_NAME $ENV{CTEST_BUILD_CONFIGURATION_NAME})

SET(CTEST_BUILD_FLAGS "-j8 -i")
SET(CTEST_PARALLEL_LEVEL "8")
SET(Trilinos_CTEST_DO_ALL_AT_ONCE TRUE)
SET(Trilinos_ENABLE_CONFIGURE_TIMING ON)
SET(CTEST_EXPLICITLY_ENABLE_IMPLICITLY_ENABLED_PACKAGES FALSE)

SET(Trilinos_REPOSITORY_LOCATION "https://github.com/trilinos/Trilinos.git")
SET(Trilinos_BRANCH develop)

SET( EXTRA_CONFIGURE_OPTIONS
  "-C${STD_CONFIG_FILE}"
  "-DTrilinos_ENABLE_CONFIGURE_TIMING=ON"
  )

TRILINOS_SYSTEM_SPECIFIC_CTEST_DRIVER()