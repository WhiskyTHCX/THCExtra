find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})
set(TOVLIB_DEPS ${TOVLIB_DEPS} ${GSL_LIBRARIES})

find_package(Boost  REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

set(TOVLIB_SOURCES eos_barotropic.cpp lookuptab.cpp splines.cpp eos_barotropic_file.cpp modes.cpp sturm.cpp eos_polytropic.cpp odes.cpp tovseq.cpp eos_tabulated.cpp eos_piecewise_poly.cpp polynomial.cpp tovsol.cpp frobenius.cpp roots.cpp minima.cpp unitconv.cpp functors.cpp spharmonics.cpp)

add_library(tov_static STATIC ${TOVLIB_SOURCES})
target_link_libraries (tov_static ${TOVLIB_DEPS})

add_library(tov SHARED ${TOVLIB_SOURCES})
target_link_libraries (tov ${TOVLIB_DEPS})

install(TARGETS tov LIBRARY DESTINATION lib)
FILE(GLOB _tov_headers "${CMAKE_CURRENT_SOURCE_DIR}/*.h")

set (_tov_headers_real "")
foreach (_tov_header ${_tov_headers})
    get_filename_component(_tov_header_real "${_tov_header}" REALPATH)
    list (APPEND _tov_headers_real "${_tov_header_real}")
endforeach()
install(FILES ${_tov_headers_real} DESTINATION include)

find_package(Doxygen)
if(DOXYGEN_FOUND)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
  add_custom_target(tovlib-doc ALL ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMENT "Generating API documentation with Doxygen" VERBATIM)
  install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc/ DESTINATION share/doc/tovlib) 
endif(DOXYGEN_FOUND)

