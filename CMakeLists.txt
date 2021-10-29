cmake_minimum_required( VERSION 3.16 )
project( KiT-RT VERSION 0.1.0 LANGUAGES CXX )

### OPTIONS #####################################
option( BUILD_TESTING "builds all available unit tests" OFF )
option( BUILD_GUI "additionally builds a user interface" OFF )
option( BUILD_DOC "builds Doxygen and Sphinx documentation" OFF )
option( BUILD_UNITY "enables unity build for faster compile times" ON )
option( BUILD_CODE_COV "enables compiler option required for code coverage analysis" OFF )
#################################################


### COMPILER ####################################
set( CMAKE_CXX_STANDARD 17 )
set( CMAKE_CXX_STANDARD_REQUIRED ON )
set( KITRT_RELEASE_OPTIONS -march=native -w )
set( KITRT_RELWITHDEBINFO_OPTIONS -march=native -pg -no-pie )
set( KITRT_DEBUG_OPTIONS -Wall -Wextra -Wpedantic )
if( BUILD_UNITY AND NOT BUILD_CODE_COV )
    message( STATUS "Unity build enabled" )
    set( CMAKE_UNITY_BUILD ON )
    set( CMAKE_UNITY_BUILD_BATCH_SIZE 0 )
endif()
#################################################


### LIBRARIES ###################################
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})

find_package( OpenMP REQUIRED )

find_package( MPI REQUIRED )
include_directories( ${MPI_INCLUDE_PATH} )

find_package( LAPACK REQUIRED )
include_directories( ${LAPACK_INCLUDE_DIR} )

find_package( BLAS )
# TODO change Hints for VTK and Python before merging
find_package( VTK COMPONENTS vtkIOGeometry vtkFiltersCore REQUIRED HINTS ~/VTK-install)

find_package( Python3 COMPONENTS Interpreter Development NumPy REQUIRED HINTS ~/anaconda3/lib/python3.8)
include_directories( ${Python3_INCLUDE_DIRS} ${Python3_NumPy_INCLUDE_DIRS} )
add_compile_definitions( KITRT_PYTHON_PATH="${CMAKE_SOURCE_DIR}/python" )

add_compile_definitions( BLAZE_USE_SHARED_MEMORY_PARALLELIZATION=0 )
message( STATUS "Blaze: Shared memory parallelization disabled" )
if( BLAS_FOUND )
    message( STATUS "Blaze: BLAS mode enabled" )
    add_compile_definitions( BLAZE_BLAS_MODE=1 )
    include_directories( ${BLAS_INCLUDE_DIR} )
else()
    message( STATUS "Blaze: BLAS mode disabled" )
    add_compile_definitions( BLAZE_BLAS_MODE=0 )
endif()
include( blaze-cache-config )
include_directories( ${CMAKE_SOURCE_DIR}/ext/blaze )

include_directories( ${CMAKE_SOURCE_DIR}/ext/cpptoml/include )
include_directories( ${CMAKE_SOURCE_DIR}/ext/spdlog/include )

set( CORE_LIBRARIES ${Python3_LIBRARIES} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} ${MPI_LIBRARIES} ${VTK_LIBRARIES} OpenMP::OpenMP_CXX -lstdc++fs )
#################################################


### MISC ########################################
execute_process(
    COMMAND git rev-parse HEAD
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_HASH
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
add_compile_definitions( GIT_HASH="${GIT_HASH}" )
#################################################


### BUILD KIT-RT ################################
file( GLOB_RECURSE SRCS RELATIVE ${CMAKE_SOURCE_DIR} "src/*.cpp" "include/*.h" )
set( EXCLUDE_DIR "/gui/" )
foreach( TMP_PATH ${SRCS} )
    string( FIND ${TMP_PATH} ${EXCLUDE_DIR} EXCLUDE_DIR_FOUND )
    if( NOT ${EXCLUDE_DIR_FOUND} EQUAL -1 )
        list( REMOVE_ITEM SRCS ${TMP_PATH} )
    endif()
endforeach( TMP_PATH )
include_directories( ${CMAKE_SOURCE_DIR}/include )
add_executable( ${CMAKE_PROJECT_NAME} ${SRCS} ${EXT_SRCS} )
target_link_libraries( ${CMAKE_PROJECT_NAME} ${CORE_LIBRARIES} )
target_compile_options( ${CMAKE_PROJECT_NAME} PUBLIC "$<$<CONFIG:DEBUG>:${KITRT_DEBUG_OPTIONS}>" )
target_compile_options( ${CMAKE_PROJECT_NAME} PUBLIC "$<$<CONFIG:RELWITHDEBINFO>:${KITRT_RELWITHDEBINFO_OPTIONS}>" )
target_compile_options( ${CMAKE_PROJECT_NAME} PUBLIC "$<$<CONFIG:RELEASE>:${KITRT_RELEASE_OPTIONS}>" )
#################################################


### BUILD UNIT TESTS ############################
if( BUILD_TESTING )
    include( CTest )
    set( CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/ext/Catch2/contrib ${CMAKE_MODULE_PATH} )
    include( Catch )
    set( CATCH_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/ext/Catch2/single_include/catch2 )
    add_compile_definitions( BUILD_TESTING )
    add_compile_definitions( TESTS_PATH="${CMAKE_SOURCE_DIR}/tests/" )
    add_library( Catch INTERFACE )
    target_include_directories( Catch INTERFACE ${CATCH_INCLUDE_DIR} )
    file( GLOB_RECURSE TEST_SRCS RELATIVE ${CMAKE_SOURCE_DIR} "tests/*.cpp" )
    list( REMOVE_ITEM SRCS "src/main.cpp" )
    add_executable( unit_tests ${TEST_SRCS} ${SRCS} ${EXT_SRCS} )
    target_compile_definitions( unit_tests PUBLIC BUILD_TESTING )
    target_compile_definitions( unit_tests PUBLIC TESTS_PATH="${CMAKE_SOURCE_DIR}/tests/" )
    target_link_libraries( unit_tests Catch ${CORE_LIBRARIES} )
    target_compile_options( unit_tests PUBLIC "$<$<CONFIG:DEBUG>:${KITRT_DEBUG_OPTIONS}>" )
    if( BUILD_CODE_COV )
        if( CMAKE_COMPILER_IS_GNUCXX )
            set( CODE_COVERAGE_OPTIONS --coverage -g -O0 -w )
            target_compile_options( unit_tests PUBLIC "$<$<CONFIG:DEBUG>:${CODE_COVERAGE_OPTIONS}>" )
            target_link_libraries( unit_tests Catch gcov )
        else()
            message( FATAL_ERROR "Code coverage is currently only supported for gcc!" )
        endif()
    else()
        target_link_libraries( unit_tests Catch )
    endif()
    target_compile_options( unit_tests PUBLIC "$<$<CONFIG:RELWITHDEBINFO>:${KITRT_RELWITHDEBINFO_OPTIONS}>" )
    target_compile_options( unit_tests PUBLIC "$<$<CONFIG:RELEASE>:${KITRT_RELEASE_OPTIONS}>" )
    catch_discover_tests( unit_tests )
endif()
#################################################


### BUILD GUI ###################################
if( BUILD_GUI )
    find_package( Qt5 COMPONENTS Core Widgets Gui OpenGL REQUIRED )
    find_package( VTK COMPONENTS vtkGUISupportQt vtkInteractionStyle vtkInteractionWidgets vtkRenderingAnnotation vtkRenderingCore vtkRenderingFreeType vtkRenderingOpenGL2 REQUIRED )
    include_directories( ${CMAKE_SOURCE_DIR}/include/gui )
    file( GLOB_RECURSE GUI_SRCS RELATIVE ${CMAKE_SOURCE_DIR} "src/gui/*" )
    file( GLOB_RECURSE GUI_INCL RELATIVE ${CMAKE_SOURCE_DIR} "include/gui/*" )
    add_executable( ${CMAKE_PROJECT_NAME}_gui ${SRCS} ${EXT_SRCS} ${GUI_SRCS} ${GUI_INCL})
    target_compile_definitions( ${CMAKE_PROJECT_NAME}_gui PUBLIC BUILD_GUI )
    set_target_properties( ${CMAKE_PROJECT_NAME}_gui PROPERTIES AUTOMOC ON )
    set_target_properties( ${CMAKE_PROJECT_NAME}_gui PROPERTIES AUTOUIC ON )
    set_target_properties( ${CMAKE_PROJECT_NAME}_gui PROPERTIES AUTORCC ON )
    set( GUI_LIBRARIES Qt5::Core Qt5::Widgets Qt5::Gui Qt5::OpenGL ${VTK_LIBRARIES} )
    target_link_libraries( ${CMAKE_PROJECT_NAME}_gui ${CORE_LIBRARIES} ${GUI_LIBRARIES} )
    set_target_properties( ${CMAKE_PROJECT_NAME}_gui PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin" )
    target_compile_options( ${CMAKE_PROJECT_NAME}_gui PUBLIC "$<$<CONFIG:DEBUG>:${KITRT_DEBUG_OPTIONS}>" )
    target_compile_options( ${CMAKE_PROJECT_NAME}_gui PUBLIC "$<$<CONFIG:RELWITHDEBINFO>:${KITRT_RELWITHDEBINFO_OPTIONS}>" )
    target_compile_options( ${CMAKE_PROJECT_NAME}_gui PUBLIC "$<$<CONFIG:RELEASE>:${KITRT_RELEASE_OPTIONS}>" )
endif()
#################################################


### BUILD DOCUMENTATION #########################
if( BUILD_DOC )
    find_package( Doxygen REQUIRED )
    find_package( Sphinx REQUIRED )

    file( GLOB_RECURSE DOC_HEADER RELATIVE ${CMAKE_SOURCE_DIR} "include/*.h" )

    set( DOXYGEN_INPUT_DIR ${PROJECT_SOURCE_DIR} )
    set( DOXYGEN_OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/docs/doxygen )
    set( DOXYGEN_INDEX_FILE ${DOXYGEN_OUTPUT_DIR}/xml/index.xml )
    set( DOXYFILE_IN ${CMAKE_CURRENT_SOURCE_DIR}/doc/Doxyfile.in )
    set( DOXYFILE_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile )
    configure_file( ${DOXYFILE_IN} ${DOXYFILE_OUT} @ONLY )
    file( MAKE_DIRECTORY ${DOXYGEN_OUTPUT_DIR} )
    add_custom_command( OUTPUT ${DOXYGEN_INDEX_FILE}
                        DEPENDS ${DOC_HEADER}
                        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYFILE_OUT}
                        MAIN_DEPENDENCY ${DOXYFILE_OUT} ${DOXYFILE_IN}
                        COMMENT "Generating documentaion with Doxygen"
                        VERBATIM )
    add_custom_target( Doxygen ALL DEPENDS ${DOXYGEN_INDEX_FILE} )

    set( SPHINX_SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/doc )
    set( SPHINX_BUILD ${CMAKE_CURRENT_BINARY_DIR}/docs/sphinx )
    set( SPHINX_INDEX_FILE ${SPHINX_BUILD}/index.html )
    add_custom_command( OUTPUT ${SPHINX_INDEX_FILE}
                        COMMAND ${SPHINX_EXECUTABLE} ${SPHINX_SOURCE} ${SPHINX_BUILD}
                        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                        DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/doc/index.rst ${DOXYGEN_INDEX_FILE}
                        MAIN_DEPENDENCY ${SPHINX_SOURCE}/conf.py
                        COMMENT "Converting documentation from Doxygen to Sphinx" )

    add_custom_target( Sphinx ALL DEPENDS ${SPHINX_INDEX_FILE} )
endif()
#################################################
