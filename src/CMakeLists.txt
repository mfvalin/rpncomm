cmake_minimum_required(VERSION 3.16)

option(RPN_BUILD_TESTS "Whether to build the test programs associated with the RPN_Comm library")
option(RPN_BUILD_STUBS "Whether to build the stubs and, if building the serial version of RPN_Comm, include them in the library")

# Files to compile
set(RPN_VERSION_FILE ${CMAKE_CURRENT_BINARY_DIR}/RPN_COMM_version.F90)
file(GLOB PROJECT_INCLUDE_FILES RPN_*.h RPN_COMM*.inc RPN_*.hf c_binding_extras.hf iso_c_binding_extras.hf)
file(GLOB PROJECT_C_FILES RPN_*.c)
file(GLOB PROJECT_F_FILES RPN_*.F90  RPN_COMM*.f90)
#file(GLOB PROJECT_MPI_F_FILES RPN_MPI*.F90)
list(APPEND PROJECT_F_FILES ${RPN_VERSION_FILE})

#message(STATUS "(EC) RPN_COMM_VERSION_STRING = ${RPN_COMM_VERSION_STRING}")
#message(STATUS "(EC) RPN_COMM_VERSION_MARKER = ${RPN_COMM_VERSION_MARKER}")

if (RPN_BUILD_STUBS)
   add_subdirectory(stubs)
   set(PROJECT_STUB_FILES "${RPN_STUB_SOURCE_FILES}")
endif()

# Main target and some custom compilation configuration
add_library(${PROJECT_NAME} STATIC
   ${PROJECT_INCLUDE_FILES}
   ${PROJECT_C_FILES}
   ${PROJECT_F_FILES}
   ${PROJECT_MPI_F_FILES}
   ${PROJECT_STUB_FILES}
)

configure_file(RPN_COMM_version.F90.in ${RPN_VERSION_FILE} @ONLY)
target_compile_definitions(${PROJECT_NAME} PUBLIC MARKER=RPN_COMM_marker_${RPN_COMM_VERSION_MARKER})
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
target_compile_options(
  ${PROJECT_NAME}
  PRIVATE
  $<$<COMPILE_LANG_AND_ID:Fortran,Intel>:-diag-disable 5268>
)
# Targets for generating header files
add_executable(c_binding_extras c_binding_extras.c)
add_custom_command(
  OUTPUT c_binding_extras.hf
  COMMAND ${CMAKE_CURRENT_BINARY_DIR}/c_binding_extras > ${CMAKE_CURRENT_SOURCE_DIR}/c_binding_extras.hf 
  COMMENT "Generating c_binding_extras.hf"
  DEPENDS c_binding_extras
  BYPRODUCTS c_binding_extras.hf
  USES_TERMINAL
)
add_custom_target(c_binding_extras_include DEPENDS c_binding_extras.hf)
add_dependencies(${PROJECT_NAME} c_binding_extras_include)

set(GEN_INTERFACES_SCRIPT ${CMAKE_CURRENT_SOURCE_DIR}/../scripts/gen_interfaces.sh)
add_custom_command(
   OUTPUT RPN_COMM_interfaces_int.inc
   COMMAND ${GEN_INTERFACES_SCRIPT} ${PROJECT_F_FILES} > RPN_COMM_interfaces_int.inc
   DEPENDS ${GEN_INTERFACES_SCRIPT} ${PROJECT_F_FILES}
)
add_custom_command(
   OUTPUT RPN_COMM_interfaces.inc
   COMMAND cat RPN_COMM_interfaces_int.inc | grep -v '^\#' > RPN_COMM_interfaces.inc
   DEPENDS RPN_COMM_interfaces_int.inc
)
add_custom_target(RPN_COMM_interfaces DEPENDS RPN_COMM_interfaces.inc)
add_dependencies(${PROJECT_NAME} RPN_COMM_interfaces)

set(FULL_LIB_NAME ${PROJECT_NAME}_${PROJECT_VERSION}${MPI_LIBRARY_VERSION})
set(SHORT_LIB_FILE_NAME lib${PROJECT_NAME}${MPI_LIBRARY_VERSION}.a)
# Adjust library name to include MPI version
set_target_properties(${PROJECT_NAME} PROPERTIES OUTPUT_NAME ${FULL_LIB_NAME})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
   # Sanity check for library content
   COMMAND ar t $<TARGET_FILE_NAME:${PROJECT_NAME}> | sort -u > objects.lst
   COMMAND sort -u ${CMAKE_CURRENT_SOURCE_DIR}/REFERENCE.lst > sorted_reference.lst
   COMMAND diff sorted_reference.lst objects.lst

   # Symbolic link to latest build
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
   set(TEST_SOURCES
      TEST_000.F90
      TEST_001.F90
      TEST_002.F90
      # TEST_003.F90
      TEST_004.F90
      TEST_005.F90
      TEST_006.F90
      TEST_007.F90
      TEST_008.F90
      TEST_009.F90
      TEST_010.F90
      TEST_011.F90
      TEST_012.F90
      # TEST_013.F90
   )

   function(get_test_name SOURCE_FILE TEST_NAME)
      get_filename_component(BASE_NAME ${SOURCE_FILE} NAME_WLE)
      set(${TEST_NAME} ${PROJECT_NAME}_${BASE_NAME} PARENT_SCOPE)
   endfunction(get_test_name)

   # Build the tests
   foreach(TEST_SOURCE ${TEST_SOURCES})
      get_test_name(${TEST_SOURCE} TEST_TARGET)


      set(CMAKE_RPN_TEST_FUNCTION_NAME ${TEST_TARGET})
      set(TEST_PROGRAM_FILE ${CMAKE_CURRENT_BINARY_DIR}/${TEST_SOURCE}_main.f90)
      file(WRITE ${TEST_PROGRAM_FILE} "program MAIN_TEST_PROGRAM\n  call ${CMAKE_RPN_TEST_FUNCTION_NAME}\n  stop\nend program MAIN_TEST_PROGRAM\n")

      add_executable(${TEST_TARGET} ${TEST_SOURCE} ${TEST_PROGRAM_FILE})
      target_link_libraries(${TEST_TARGET} ${PROJECT_NAME})
   endforeach(TEST_SOURCE)
endif()
