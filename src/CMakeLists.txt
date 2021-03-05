cmake_minimum_required(VERSION 3.16)

include(gen_headers)
include(gen_test_target)

option(BUILD_TESTS "Whether to build the test programs associated with the RPN_Comm library")

file(GLOB PROJECT_INCLUDE_FILES RPN_COMM*.h RPN_COMM*.inc RPN_COMM*.hf)
file(GLOB PROJECT_C_FILES RPN_COMM*.c)
file(GLOB PROJECT_F_FILES RPN_COMM*.F90  RPN_COMM*.f90)

message(STATUS "RPN_COMM_VERSION_STRING = ${RPN_COMM_VERSION_STRING}")
message(STATUS "RPN_COMM_VERSION_MARKER = ${RPN_COMM_VERSION_MARKER}")

add_library(${PROJECT_NAME} STATIC ${PROJECT_INCLUDE_FILES} ${PROJECT_C_FILES} ${PROJECT_F_FILES})

configure_file(RPN_COMM_version.F90.in ${CMAKE_CURRENT_SOURCE_DIR}/RPN_COMM_version.F90 @ONLY)
target_compile_definitions(${PROJECT_NAME} PUBLIC MARKER=RPN_COMM_marker_${RPN_COMM_VERSION_MARKER})
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

gen_headers(SCRIPT ${CMAKE_CURRENT_SOURCE_DIR}/gen_interfaces.sh INPUT_FILES ${PROJECT_F_FILES} OUTPUT_SUFFIX inc OUTPUT_FILENAME RPN_COMM_interfaces_int.inc)

add_custom_command(
    OUTPUT RPN_COMM_interfaces.inc
    COMMAND cat RPN_COMM_interfaces_int.inc | grep -v '^\#' > RPN_COMM_interfaces.inc
    DEPENDS RPN_COMM_interfaces_int.inc
)
add_custom_target(RPN_COMM_interfaces DEPENDS RPN_COMM_interfaces.inc)
add_dependencies(${PROJECT_NAME} RPN_COMM_interfaces)

if (BUILD_TESTS)
    gen_test_target(GENERATE_F_PROGRAM rpn_comm INPUT_FILES TEST_000.F90)
    gen_test_target(GENERATE_F_PROGRAM rpn_comm INPUT_FILES TEST_001.F90 LINK_LIBRARIES ${PROJECT_NAME})
    gen_test_target(GENERATE_F_PROGRAM rpn_comm INPUT_FILES TEST_002.F90 LINK_LIBRARIES ${PROJECT_NAME})
    gen_test_target(GENERATE_F_PROGRAM rpn_comm INPUT_FILES TEST_003.F90 LINK_LIBRARIES ${PROJECT_NAME})
    gen_test_target(GENERATE_F_PROGRAM rpn_comm INPUT_FILES TEST_004.F90 LINK_LIBRARIES ${PROJECT_NAME})
    gen_test_target(GENERATE_F_PROGRAM rpn_comm INPUT_FILES TEST_005.F90)
    gen_test_target(GENERATE_F_PROGRAM rpn_comm INPUT_FILES TEST_006.F90 LINK_LIBRARIES ${PROJECT_NAME})
    gen_test_target(GENERATE_F_PROGRAM rpn_comm INPUT_FILES TEST_007.F90 LINK_LIBRARIES ${PROJECT_NAME})
    gen_test_target(GENERATE_F_PROGRAM rpn_comm INPUT_FILES TEST_008.F90 LINK_LIBRARIES ${PROJECT_NAME})
    gen_test_target(GENERATE_F_PROGRAM rpn_comm INPUT_FILES TEST_009.F90 LINK_LIBRARIES ${PROJECT_NAME})
    gen_test_target(GENERATE_F_PROGRAM rpn_comm INPUT_FILES TEST_010.F90 LINK_LIBRARIES ${PROJECT_NAME})
    gen_test_target(GENERATE_F_PROGRAM rpn_comm INPUT_FILES TEST_011.F90 LINK_LIBRARIES ${PROJECT_NAME})
    gen_test_target(GENERATE_F_PROGRAM rpn_comm INPUT_FILES TEST_012.F90 LINK_LIBRARIES ${PROJECT_NAME})
    # gen_test_target(GENERATE_F_PROGRAM rpn_comm INPUT_FILES TEST_013.F90)
endif()