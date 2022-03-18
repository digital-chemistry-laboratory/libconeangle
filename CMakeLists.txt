cmake_minimum_required(VERSION 3.17)
project(libconeangle LANGUAGES C Fortran)

# Set flags for Windows compatibility
set(CMAKE_SHARED_LIBRARY_PREFIX lib)
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)

# Fetch stdlib
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/stdlib/CMakeLists.txt)
  add_subdirectory("stdlib")
else()
  set("stdlib-url" "https://github.com/fortran-lang/stdlib")
  message(STATUS "Retrieving stdlib from ${stdlib-url}")
  include(FetchContent)
  FetchContent_Declare(
    "fortran_stdlib"
    GIT_REPOSITORY "${stdlib-url}"
    GIT_TAG "HEAD"
  )
  FetchContent_GetProperties("fortran_stdlib")
  if(NOT fortran_stdlib_POPULATED)
    FetchContent_Populate("fortran_stdlib")
    add_subdirectory(${fortran_stdlib_SOURCE_DIR} ${fortran_stdlib_BINARY_DIR} EXCLUDE_FROM_ALL)
  endif()
endif()

# Set up shared library target
set(srcs)
add_subdirectory(src)
add_library("coneangle" SHARED ${srcs})
target_link_libraries("coneangle" "fortran_stdlib")
set_target_properties("coneangle" PROPERTIES PUBLIC_HEADER "include/cone_angle.h")

# Install
if(SKBUILD)
  message(STATUS "The project is built using scikit-build")
  install(TARGETS coneangle LIBRARY DESTINATION libconeangle ARCHIVE DESTINATION libconeangle RUNTIME DESTINATION libconeangle)
else()
  install(TARGETS coneangle LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} RUNTIME DESTINATION ${CMAKE_INSTALL_LIBDIR} PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/libconeangle)
endif()