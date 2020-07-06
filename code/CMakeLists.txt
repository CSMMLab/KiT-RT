cmake_minimum_required( VERSION 3.12.4 )
project( RTSN VERSION 0.1.0 LANGUAGES CXX )
set( CMAKE_CXX_STANDARD 17 )
set( CMAKE_CXX_STANDARD_REQUIRED ON )

set( CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native -Wno-dev" )
set( CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -march=native -pg -no-pie" )
set( CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -Wall" )

set( CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin" )

file( GLOB_RECURSE SRCS RELATIVE ${CMAKE_SOURCE_DIR} "src/*.cpp" "include/*.h" )
include_directories( ${CMAKE_SOURCE_DIR}/include )
add_executable( ${CMAKE_PROJECT_NAME} ${SRCS} )

find_package( Python3 COMPONENTS Interpreter Development NumPy REQUIRED )
include_directories( ${Python3_INCLUDE_DIRS} ${Python3_NumPy_INCLUDE_DIRS} )
add_compile_definitions( RTSN_PYTHON_PATH="${CMAKE_SOURCE_DIR}/python" )

find_package( OpenMP REQUIRED )
if( OPENMP_FOUND )
    set( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}" )
    set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}" )
endif()

find_package( MPI REQUIRED )
include_directories( ${MPI_INCLUDE_PATH} )

find_package( LAPACK REQUIRED )
include_directories( ${LAPACK_INCLUDE_DIR} )

find_package( BLAS )
if( BLAS_FOUND )
    message( STATUS "Blaze: BLAS mode enabled" )
    add_compile_definitions( BLAZE_BLAS_MODE=1 )
    include_directories( ${BLAS_INCLUDE_DIR} )
else()
    message( STATUS "Blaze: BLAS mode disabled" )
    add_compile_definitions( BLAZE_BLAS_MODE=0 )
endif()

add_compile_definitions( BLAZE_USE_SHARED_MEMORY_PARALLELIZATION=0 )

set(flag 1)
if (WIN32)
    execute_process(COMMAND wmic cpu get L3CacheSize
        OUTPUT_VARIABLE tmp
        RESULT_VARIABLE flag
        ERROR_QUIET)
    if (flag)
        execute_process(COMMAND wmic cpu get L2CacheSize
            OUTPUT_VARIABLE tmp
            RESULT_VARIABLE flag
            ERROR_QUIET)
    endif (flag)
    if (flag)
        execute_process(COMMAND wmic cpu get L1CacheSize
            OUTPUT_VARIABLE tmp
            RESULT_VARIABLE flag
            ERROR_QUIET)
    endif (flag)
endif (WIN32)

if (UNIX)
    execute_process(COMMAND cat /sys/devices/system/cpu/cpu0/cache/index3/size
        OUTPUT_VARIABLE tmp
        RESULT_VARIABLE flag
        ERROR_QUIET)
    if (flag)
        execute_process(COMMAND cat /sys/devices/system/cpu/cpu0/cache/index2/size
            OUTPUT_VARIABLE tmp
            RESULT_VARIABLE flag
            ERROR_QUIET)
    endif (flag)
    if (flag)
        execute_process(COMMAND cat /sys/devices/system/cpu/cpu0/cache/index1/size
            OUTPUT_VARIABLE tmp
            RESULT_VARIABLE flag
            ERROR_QUIET)
    endif (flag)
endif (UNIX)

if (APPLE)
    execute_process(COMMAND sysctl -n hw.l3cachesize
        OUTPUT_VARIABLE tmp
        RESULT_VARIABLE flag
        ERROR_QUIET)
    if (flag)
        execute_process(COMMAND sysctl -n hw.l2cachesize
            OUTPUT_VARIABLE tmp
            RESULT_VARIABLE flag
            ERROR_QUIET)
    endif (flag)
    if (flag)
        execute_process(COMMAND sysctl -n hw.l1icachesize
            OUTPUT_VARIABLE tmp
            RESULT_VARIABLE flag
            ERROR_QUIET)
    endif (flag)

    if (flag EQUAL 0)
        math(EXPR tmp ${tmp}/1024)  # If successful convert to kibibytes to comply with rest
    endif (flag EQUAL 0)
endif (APPLE)

if (flag)
    message(WARNING "Cache size not found automatically. Using default value as cache size.")
    set(tmp "3072")
endif (flag)

string( REGEX MATCH "([0-9][0-9]+)" tmp ${tmp} )
math( EXPR BLAZE_CACHE_SIZE ${tmp}*1024 )
add_compile_definitions( BLAZE_CACHE_SIZE=${BLAZE_CACHE_SIZE}UL )
message( STATUS "Blaze: Automatic Cache Size Configuration = ${BLAZE_CACHE_SIZE} KiB" )
include_directories( ${CMAKE_SOURCE_DIR}/ext/blaze )

add_compile_definitions( METIS_EXPORT= )
set( DISABLE_PARMETIS_PROGRAMS ON )
set( ParMETIS_PATH ${CMAKE_CURRENT_SOURCE_DIR}/ext/parmetis )
include_directories( ${ParMETIS_PATH}/include )
include_directories( ${ParMETIS_PATH}/metis/include )
add_subdirectory( ${ParMETIS_PATH} )

find_package( VTK REQUIRED COMPONENTS vtkIOGeometry vtkFiltersCore HINTS /home/pia/vtk-inst/include/vtk-8.2)
include(${VTK_USE_FILE})

include_directories( ${CMAKE_SOURCE_DIR}/ext/cpptoml/include )
include_directories( ${CMAKE_SOURCE_DIR}/ext/spdlog/include )

execute_process(
    COMMAND git rev-parse HEAD
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_HASH
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
add_compile_definitions( GIT_HASH="${GIT_HASH}" )

set( CORE_LIBRARIES ${Python3_LIBRARIES} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} ${MPI_LIBRARIES} ${VTK_LIBRARIES} OpenMP::OpenMP_CXX parmetis -lstdc++fs )

target_link_libraries( ${CMAKE_PROJECT_NAME} ${CORE_LIBRARIES} )

if( BUILD_TESTING )
    add_compile_definitions( BUILD_TESTING=1 )
    add_subdirectory( ${CMAKE_SOURCE_DIR}/ext/Catch2 )
    include_directories( ${CMAKE_SOURCE_DIR}/ext/Catch2/single_include/catch2 )
    file( GLOB_RECURSE TEST_SRCS RELATIVE ${CMAKE_SOURCE_DIR} "tests/*.cpp" )
    list( REMOVE_ITEM SRCS "src/main.cpp" )
    add_executable( unit_tests ${SRCS} ${TEST_SRCS} )
    target_link_libraries( unit_tests ${CORE_LIBRARIES} )
endif()
