cmake_minimum_required(VERSION 3.17)
project(libconeangle LANGUAGES C Fortran)

# Set flags for Windows compatibility
set(CMAKE_SHARED_LIBRARY_PREFIX lib)
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)

# Fetch stdlib and link
set("stdlib-url" "https://github.com/fortran-lang/stdlib")
message(STATUS "Retrieving fortran_stdlib from ${stdlib-url}")
include(FetchContent)
FetchContent_Declare(
  "fortran_stdlib"
  GIT_REPOSITORY "${stdlib-url}"
  GIT_TAG "HEAD"
)
FetchContent_MakeAvailable("fortran_stdlib")

set(srcs)
add_subdirectory(src)
add_library(coneangle SHARED ${srcs})
target_link_libraries(coneangle PUBLIC fortran_stdlib)

if(SKBUILD)
  message(STATUS "The project is built using scikit-build")
  install(TARGETS coneangle LIBRARY DESTINATION libconeangle ARCHIVE DESTINATION libconeangle RUNTIME DESTINATION libconeangle)
endif()
