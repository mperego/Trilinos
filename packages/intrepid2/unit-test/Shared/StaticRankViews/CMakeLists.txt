# PROTOTYPE: static-rank Kokkos::View overloads vs. the DynRankView overloads that forward to them.
# The test uses Sacado Fad scalars alongside double.
IF(${PACKAGE_NAME}_ENABLE_Sacado)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    test_static_rank_views
    SOURCES test_static_rank_views.cpp
    NUM_MPI_PROCS 1
    PASS_REGULAR_EXPRESSION "TEST PASSED"
    ADD_DIR_TO_NAME
    )
ENDIF()
