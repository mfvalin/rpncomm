cmake_minimum_required(VERSION 3.16)

include(gen_headers)
include(gen_test_target)

option(RPN_BUILD_TESTS "Whether to build the test programs associated with the RPN_Comm library")
option(RPN_BUILD_STUBS "Whether to build the stubs and, if building the serial version of RPN_Comm, include them in the library")

# Files to compile
file(GLOB PROJECT_INCLUDE_FILES RPN_COMM*.h RPN_COMM*.inc RPN_COMM*.hf)
file(GLOB PROJECT_C_FILES RPN_COMM*.c)
file(GLOB PROJECT_F_FILES RPN_COMM*.F90  RPN_COMM*.f90)

#message(STATUS "RPN_COMM_VERSION_STRING = ${RPN_COMM_VERSION_STRING}")
#message(STATUS "RPN_COMM_VERSION_MARKER = ${RPN_COMM_VERSION_MARKER}")

if (RPN_BUILD_STUBS)
   add_subdirectory(stubs)
   set(PROJECT_STUB_FILES "${RPN_STUB_SOURCE_FILES}")
endif()

# Main target and some custom compilation configuration
add_library(${PROJECT_NAME} STATIC
   ${PROJECT_INCLUDE_FILES}
   ${PROJECT_C_FILES}
   ${PROJECT_F_FILES}
   ${PROJECT_STUB_FILES}
)

configure_file(RPN_COMM_version.F90.in ${CMAKE_CURRENT_SOURCE_DIR}/RPN_COMM_version.F90 @ONLY)
target_compile_definitions(${PROJECT_NAME} PUBLIC MARKER=RPN_COMM_marker_${RPN_COMM_VERSION_MARKER})
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_BINARY_DIR})

# Targets for generating header files
add_custom_command(
   OUTPUT RPN_COMM_interfaces_int.inc
   COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/gen_interfaces.sh ${PROJECT_F_FILES} > RPN_COMM_interfaces_int.inc
   DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/gen_interfaces.sh ${PROJECT_F_FILES}
   )
add_custom_command(
   OUTPUT RPN_COMM_interfaces.inc
   COMMAND cat RPN_COMM_interfaces_int.inc | grep -v '^\#' > RPN_COMM_interfaces.inc
   DEPENDS RPN_COMM_interfaces_int.inc
)
add_custom_target(RPN_COMM_interfaces DEPENDS RPN_COMM_interfaces.inc)
add_dependencies(${PROJECT_NAME} RPN_COMM_interfaces)

# Adjust library file name and setup for a symbolic link
set(FULL_LIB_NAME ${PROJECT_NAME}_${PROJECT_VERSION}${MPI_LIBRARY_VERSION})
set(SHORT_LIB_FILE_NAME lib${PROJECT_NAME}${MPI_LIBRARY_VERSION}.a)
set_target_properties(${PROJECT_NAME} PROPERTIES OUTPUT_NAME ${FULL_LIB_NAME})
add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
   COMMAND ln -sfv $<TARGET_FILE_NAME:${PROJECT_NAME}> ${SHORT_LIB_FILE_NAME}
)

# Installation
set(DISTINCLUDES
   RPN_COMM.inc
   rpn_comm.inc 
   RPN_COMM_types.inc
   RPN_COMM_constants.inc
   RPN_COMM_ftoc.inc
   RPN_COMM_is_null.inc
   ${CMAKE_CURRENT_BINARY_DIR}/RPN_COMM_interfaces.inc
)

install(TARGETS ${PROJECT_NAME})
install(FILES ${DISTINCLUDES} DESTINATION include)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${SHORT_LIB_FILE_NAME} DESTINATION lib)

# The tests
if (RPN_BUILD_TESTS)
    gen_test_target(GENERATE_F_PROGRAM rpn_comm INPUT_FILES TEST_000.F90)
    gen_test_target(GENERATE_F_PROGRAM rpn_comm INPUT_FILES TEST_001.F90 LINK_LIBRARIES ${PROJECT_NAME})
    gen_test_target(GENERATE_F_PROGRAM rpn_comm INPUT_FILES TEST_002.F90 LINK_LIBRARIES ${PROJECT_NAME})
    # gen_test_target(GENERATE_F_PROGRAM rpn_comm INPUT_FILES TEST_003.F90 LINK_LIBRARIES ${PROJECT_NAME})
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
