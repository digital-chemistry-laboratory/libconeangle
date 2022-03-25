cmake_minimum_required(VERSION 3.17)
project(libconeangle LANGUAGES Fortran)

# Fetch stdlib
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
  # Set compiler flags
  set_target_properties(fortran_stdlib PROPERTIES POSITION_INDEPENDENT_CODE ON)
  if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU" AND CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 10)
    target_compile_options(fortran_stdlib PRIVATE "-fno-range-check")
  endif()
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
    set_target_properties(coneangle PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS ON)
    
    # Link statically against dependencies for Windows
    if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
      target_link_options(coneangle PRIVATE "-static")
    endif()
  endif()

  message(STATUS "The project is built using scikit-build")
  install(TARGETS coneangle LIBRARY DESTINATION libconeangle ARCHIVE DESTINATION libconeangle RUNTIME DESTINATION libconeangle) 
else()
  set_target_properties(coneangle PROPERTIES PUBLIC_HEADER "include/cone_angle.h")
  install(TARGETS coneangle LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} RUNTIME DESTINATION ${CMAKE_INSTALL_LIBDIR} PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/libconeangle)
endif()