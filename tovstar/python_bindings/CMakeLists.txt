FIND_PACKAGE(SWIG)
if(SWIG_FOUND)
  INCLUDE(${SWIG_USE_FILE})

  find_package(GSL REQUIRED)
  include_directories(${GSL_INCLUDE_DIRS})
  set(TOVLIB_DEPS ${TOVLIB_DEPS} ${GSL_LIBRARIES})

  find_package(Boost  REQUIRED)
  include_directories(${Boost_INCLUDE_DIRS})

  FIND_PACKAGE(PythonInterp)

  FIND_PACKAGE(PythonLibs REQUIRED)
  INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})

  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

  SET(CMAKE_SWIG_FLAGS "")

  SET_SOURCE_FILES_PROPERTIES(pytov.i PROPERTIES CPLUSPLUS ON)
  #SET_SOURCE_FILES_PROPERTIES(pytov.i PROPERTIES SWIG_FLAGS "-includeall")
  SWIG_ADD_MODULE(pytov python pytov.i eos_barotropic.cpp lookuptab.cpp splines.cpp eos_barotropic_file.cpp modes.cpp sturm.cpp eos_polytropic.cpp odes.cpp tovseq.cpp eos_tabulated.cpp eos_piecewise_poly.cpp polynomial.cpp tovsol.cpp frobenius.cpp roots.cpp minima.cpp unitconv.cpp functors.cpp spharmonics.cpp)
  SWIG_LINK_LIBRARIES(pytov ${PYTHON_LIBRARIES} ${TOVLIB_DEPS})


  INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/pytov.py ${CMAKE_CURRENT_BINARY_DIR}/_pytov.so DESTINATION lib/python${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}/site-packages)

endif(SWIG_FOUND)

