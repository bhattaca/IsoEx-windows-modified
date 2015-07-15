cmake_minimum_required (VERSION 2.6)

include (ACGCommon)

acg_openmp ()

include_directories (
  ..
  ${OPENMESH_INCLUDE_DIRS}
  ${OPENMESH_INCLUDE_DIR}#old openmesh finder
  ${CMAKE_CURRENT_SOURCE_DIR}/../..
  ${CMAKE_CURRENT_BINARY_DIR}
)

link_directories (
  ${OPENMESH_LIBRARY_DIR}
)

# source code directories
set (directories 
  emc_poly
)

if (WIN32)
  add_definitions(
      -D_USE_MATH_DEFINES -DNOMINMAX
  )
endif ()

# collect all header,source and ui files
acg_append_files (headers "*.hh" ${directories})
acg_append_files (sources "*.cc" ${directories})

of_list_filter  ( headers )
of_list_filter  ( sources )


# remove template cc files from source file list
acg_drop_templates (sources)

acg_add_executable (EMCpoly ${sources} ${headers} )

target_link_libraries (EMCpoly
  OpenMeshCore
  OpenMeshTools
)

