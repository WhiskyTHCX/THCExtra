cmake_minimum_required (VERSION 2.6)
project (TovStars)
set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

set(Boost_NO_BOOST_CMAKE TRUE)
find_package(Boost COMPONENTS program_options filesystem system REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
set(LIBS ${LIBS} ${Boost_LIBRARIES})

include_directories ("${PROJECT_SOURCE_DIR}/tovlib")
add_subdirectory (tovlib)
set(LIBS ${LIBS} tov_static)

add_executable(tovstar tovstar.cpp)
target_link_libraries (tovstar ${LIBS})

add_executable(tovstars tovstars.cpp)
target_link_libraries (tovstars ${LIBS})

option (BUILD_PYTHON_BINDINGS "Build python bindings to tov libraries using SWIG" ON)
if (BUILD_PYTHON_BINDINGS)
  add_subdirectory (python_bindings)
endif (BUILD_PYTHON_BINDINGS)

install(TARGETS tovstar tovstars RUNTIME DESTINATION bin)
install(PROGRAMS tovstar_plot_seq  DESTINATION bin)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/examples DESTINATION share/doc/tovstar) 

option (BUILD_DOCUMENTATION "Create documentation using Doxygen" ON)
if(BUILD_DOCUMENTATION)
  find_package(Doxygen)
  if(DOXYGEN_FOUND)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
    add_custom_target(tovstar-doc ALL ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMENT "Generating API documentation with Doxygen" VERBATIM)
    install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc/ DESTINATION share/doc/tovstar) 
  endif(DOXYGEN_FOUND)
endif(BUILD_DOCUMENTATION)

