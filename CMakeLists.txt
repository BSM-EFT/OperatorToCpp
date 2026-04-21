cmake_minimum_required(VERSION 3.15...4.0)
project(OpToCpp LANGUAGES CXX)


#  General setup
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Include paths
include_directories(${PROJECT_SOURCE_DIR}/include)


#  Detect platform + OpenMP setup
find_package(OpenMP)

if(OpenMP_CXX_FOUND)
    message(STATUS "OpenMP found, using ${OpenMP_CXX_LIB_NAMES}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(OPENMP_LIBS OpenMP::OpenMP_CXX)
else()
    message(WARNING "OpenMP not found — attempting manual setup")
    if(APPLE)
        message(STATUS "Using Apple Clang — adding Homebrew OpenMP flags manually")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Xpreprocessor -fopenmp -I/opt/homebrew/Cellar/libomp/21.1.3/include")
        link_directories(/opt/homebrew/Cellar/libomp/21.1.3/lib)
        set(OPENMP_LIBS omp)
    else()
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
        set(OPENMP_LIBS gomp)
    endif()
endif()


#  Precompiled headers
set(PCH_HEADER ${PROJECT_SOURCE_DIR}/include/pch.h)


#  Build the core library
file(GLOB LIB_SOURCES ${PROJECT_SOURCE_DIR}/lib/*.cpp ${PROJECT_SOURCE_DIR}/lib/*/*.cpp)
add_library(cpp_backend STATIC ${LIB_SOURCES})
target_include_directories(cpp_backend PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(cpp_backend PRIVATE ${OPENMP_LIBS})

if(EXISTS ${PCH_HEADER})
    target_precompile_headers(cpp_backend PRIVATE ${PCH_HEADER})
endif()


#  Build the C++ executable
add_executable(main.out ${PROJECT_SOURCE_DIR}/src/main.cpp)
target_link_libraries(main.out PRIVATE cpp_backend ${OPENMP_LIBS})

if(EXISTS ${PCH_HEADER})
    target_precompile_headers(main.out REUSE_FROM cpp_backend)
endif()


#  Pybind11 Integration
# set(PYBIND11_FINDPYTHON ON)
# find_package(pybind11 CONFIG REQUIRED)
include(FetchContent)
find_package(pybind11 CONFIG)

if(NOT pybind11_FOUND)
    message(STATUS "pybind11 not found locally. Fetching from GitHub...")
    FetchContent_Declare(
      pybind11
      GIT_REPOSITORY https://github.com/pybind/pybind11.git
      GIT_TAG        v3.0.4
    )
    FetchContent_MakeAvailable(pybind11)
endif()

pybind11_add_module(match_to_py ${PROJECT_SOURCE_DIR}/src/pyBindings.cpp)
target_link_libraries(match_to_py PRIVATE cpp_backend)

if(EXISTS ${PCH_HEADER})
    target_precompile_headers(match_to_py REUSE_FROM cpp_backend)
endif()

install(TARGETS match_to_py DESTINATION .)


#  Messages
message(STATUS "------------------------------------------")
message(STATUS "  Platform: ${CMAKE_SYSTEM_NAME}")
message(STATUS "  Compiler: ${CMAKE_CXX_COMPILER}")
message(STATUS "  C++ Flags: ${CMAKE_CXX_FLAGS}")
message(STATUS "  Build Type: ${CMAKE_BUILD_TYPE}")
message(STATUS "------------------------------------------")
