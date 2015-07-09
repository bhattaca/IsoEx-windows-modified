cmake_minimum_required (VERSION 2.6)

project(IsoEx)

# add our macro directory to cmake search path
set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

include (ACGCommon)

acg_qt4 ()
acg_openmp ()

acg_get_version ()

include (ACGOutput)



find_package (OpenMesh)
if (NOT OPENMESH_FOUND )
  message (FATAL_ERROR "OPENMESH not found!")
endif ()


include_directories (
  ..
  ${OPENMESH_INCLUDE_DIRS}
  ${OPENMESH_INCLUDE_DIR}#old openmesh finder
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}
)

# generate dllexport macros on windows
if (WIN32)
  add_definitions(-DISOEXDLL)
endif ()


link_directories (
  ${OPENMESH_LIBRARY_DIR}
)

# source code directories
set (directories 
  .
  Extractors
  Grids
  Implicits
  Math
)

if (WIN32)
  add_definitions(
      -D_USE_MATH_DEFINES -DNOMINMAX
  )
endif ()

# collect all header,source and ui files
acg_append_files (headers "*.hh" ${directories})
acg_append_files (sources "*.cc" ${directories})
acg_append_files (ui "*.ui" ${directories})

macro (of_list_filter _list)
  if (WIN32)
   foreach (_element ${${_list}})
     if (_element MATCHES "gnuplot_i\\.(cc|hh)$")
       list (REMOVE_ITEM ${_list} ${_element})
     endif ()
   endforeach ()
  endif ()
endmacro ()

of_list_filter  ( headers )
of_list_filter  ( sources )


# remove template cc files from source file list
acg_drop_templates (sources)

# genereate uic and moc targets
acg_qt4_autouic (uic_targets ${ui})
acg_qt4_automoc (moc_targets ${headers})

acg_add_library (IsoEx SHARED ${uic_targets} ${sources} ${headers} ${moc_targets})

if (NOT APPLE)
	target_link_libraries (IsoEx 
	  OpenMeshCore
	  OpenMeshTools
	)
else(NOT APPLE)
	target_link_libraries (IsoEx 
	  OpenMeshCore
	  OpenMeshTools
	)
endif(NOT APPLE)

add_subdirectory(Apps)

# display results
acg_print_configure_header (ISOEX "IsoEx")

