cmake_minimum_required(VERSION 3.17)
project(libconeangle LANGUAGES Fortran)

# Fetch stdlib
set(CMAKE_MAXIMUM_RANK 2)
set(BUILD_TESTING OFF)
set(WITH_QP OFF)
set(WITH_XDP OFF)

set(stdlib-url "https://github.com/fortran-lang/stdlib")
message(STATUS "Retrieving stdlib from ${stdlib-url}")
include(FetchContent)
FetchContent_Declare(
  stdlib
  GIT_REPOSITORY ${stdlib-url}
  GIT_TAG "HEAD"
)
FetchContent_MakeAvailable(stdlib)

# Set compiler flags
set_target_properties(fortran_stdlib PROPERTIES POSITION_INDEPENDENT_CODE ON)
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU" AND CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 10)
  target_compile_options(fortran_stdlib PRIVATE "-fno-range-check")
endif()

# Set up shared library target
set(srcs)
add_subdirectory(src)
add_library(coneangle SHARED ${srcs})
target_link_libraries(coneangle fortran_stdlib)

# Install
if(SKBUILD)
  if(WIN32)
    # Set flags for Windows compatibility
    set_target_properties(coneangle PROPERTIES LIBRARY_PREFIX "lib")
    set_target_properties(coneangle WINDOWS_EXPORT_ALL_SYMBOLS ON)
    
    # Link statically against dependencies for Windows
    if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
      target_compile_options(coneangle "-static")
    endif()
  endif()

  message(STATUS "The project is built using scikit-build")
  install(TARGETS coneangle LIBRARY DESTINATION libconeangle ARCHIVE DESTINATION libconeangle RUNTIME DESTINATION libconeangle) 
else()
  set_target_properties(coneangle PROPERTIES PUBLIC_HEADER "include/cone_angle.h")
  install(TARGETS coneangle LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} RUNTIME DESTINATION ${CMAKE_INSTALL_LIBDIR} PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/libconeangle)
endif()