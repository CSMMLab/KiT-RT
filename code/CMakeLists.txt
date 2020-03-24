cmake_minimum_required( VERSION 3.5 )
project( RTSN LANGUAGES CXX VERSION 0.0.1 )
set( CMAKE_CXX_STANDARD 17 )

set( CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -march=native" )
set( CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -O3 -march=native -pg -no-pie" )
set( CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -Wall -Werror" )

set( CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin" )

file( GLOB_RECURSE SRCS RELATIVE ${CMAKE_SOURCE_DIR} "src/*.cpp" )
include_directories( ${CMAKE_SOURCE_DIR}/include )
add_executable( ${CMAKE_PROJECT_NAME} ${SRCS} )

find_package( OpenMP REQUIRED )
if( OPENMP_FOUND )
    set( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}" )
    set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}" )
endif()

find_package( MPI REQUIRED )
include_directories( ${MPI_INCLUDE_PATH} )

find_package( LAPACK REQUIRED )
include_directories( ${LAPACK_INCLUDE_DIR} )

include_directories( ${CMAKE_SOURCE_DIR}/../ext/cpptoml/include )
include_directories( ${CMAKE_SOURCE_DIR}/../ext/spdlog/include )

execute_process(
    COMMAND git rev-parse HEAD
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_HASH
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
add_definitions( '-DGIT_HASH="${GIT_HASH}"' )

set( CORE_LIBRARIES ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} -lstdc++fs )

target_link_libraries( ${CMAKE_PROJECT_NAME} ${CORE_LIBRARIES} )

if( BUILD_TESTING )
    add_subdirectory( ${CMAKE_SOURCE_DIR}/ext/Catch2 )
    include_directories( ${CMAKE_SOURCE_DIR}/ext/Catch2/single_include/catch2 )
    file( GLOB_RECURSE TEST_SRCS RELATIVE ${CMAKE_SOURCE_DIR} "tests/*.cpp" )
    list( REMOVE_ITEM SRCS "src/main.cpp" )
    add_executable( unit_tests ${SRCS} ${TEST_SRCS} )
    target_link_libraries( unit_tests ${CORE_LIBRARIES} )
endif()
