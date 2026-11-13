# PROTOTYPE: DynRankView arguments vs. static-rank Kokkos::View arguments; see bench_static_rank_views.cpp.
IF(${PACKAGE_NAME}_ENABLE_Sacado)
  TRIBITS_ADD_EXECUTABLE(
    bench_static_rank_views
    SOURCES bench_static_rank_views.cpp
    ADD_DIR_TO_NAME
    )
ENDIF()
