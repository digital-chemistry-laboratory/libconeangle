cmake_minimum_required(VERSION 3.17)
project(libconeangle LANGUAGES Fortran)

# Set flags for Windows compatibility
set(CMAKE_SHARED_LIBRARY_PREFIX lib)
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)

# Fetch stdlib
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_MAXIMUM_RANK 2)
set(BUILD_TESTING OFF)
set(WITH_QP OFF)
set(WITH_XDP OFF)
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/stdlib/CMakeLists.txt")
  add_subdirectory(stdlib)
else()
  set(stdlib-url "https://github.com/fortran-lang/stdlib")
  message(STATUS "Retrieving stdlib from ${stdlib-url}")
  include(FetchContent)
  FetchContent_Declare(
    fortran_stdlib
    GIT_REPOSITORY ${stdlib-url}
    GIT_TAG "HEAD"
  )
  FetchContent_GetProperties(fortran_stdlib)
  if(NOT fortran_stdlib_POPULATED)
    FetchContent_Populate(fortran_stdlib)
    add_subdirectory(${fortran_stdlib_SOURCE_DIR} ${fortran_stdlib_BINARY_DIR} EXCLUDE_FROM_ALL)
  endif()
endif()

# Set up shared library target
set(srcs)
add_subdirectory(src)
add_library(coneangle SHARED ${srcs})
target_link_libraries(coneangle fortran_stdlib)

# Install
if(SKBUILD)
  # Link statically against dependencies for Windows
  if(WIN32)
    if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
      set(CMAKE_Fortran_FLAGS "-static" ${CMAKE_Fortran_FLAGS})
    endif()
  endif()

  message(STATUS "The project is built using scikit-build")
  install(TARGETS coneangle LIBRARY DESTINATION libconeangle ARCHIVE DESTINATION libconeangle RUNTIME DESTINATION libconeangle) 
else()
  set_target_properties(coneangle PROPERTIES PUBLIC_HEADER "include/cone_angle.h")
  install(TARGETS coneangle LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} RUNTIME DESTINATION ${CMAKE_INSTALL_LIBDIR} PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/libconeangle)
endif()