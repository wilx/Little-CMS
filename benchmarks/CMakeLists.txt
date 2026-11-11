cmake_minimum_required(VERSION 3.16)
project(lcms_benchmarks LANGUAGES CXX)

add_library(benchmark_clock STATIC clock.cpp)
target_compile_features(benchmark_clock PUBLIC cxx_std_17)
target_include_directories(benchmark_clock PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")

find_package(OpenMP REQUIRED COMPONENTS CXX)
add_executable(lab2rgb lab2rgb.cpp)
target_include_directories(lab2rgb PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/../include")
target_link_libraries(lab2rgb PRIVATE benchmark_clock OpenMP::OpenMP_CXX ${CMAKE_DL_LIBS})
if(WIN32)
  # Match the calling convention of public LCMS DLL entry points under MSVC.
  target_compile_definitions(lab2rgb PRIVATE CMS_DLL)
endif()

find_package(Threads REQUIRED)
add_executable(test_clock test_clock.cpp)
target_link_libraries(test_clock PRIVATE benchmark_clock Threads::Threads)
enable_testing()
add_test(NAME benchmark_clocks COMMAND test_clock)
