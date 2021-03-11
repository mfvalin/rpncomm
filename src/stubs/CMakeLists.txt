cmake_minimum_required(VERSION 3.16)

# Setup variables
set(C_STUB_SOURCE rpn_comm_c_stubs.c)
set(FORTRAN_STUB_SOURCE rpn_comm_fortran_stubs.F90)

set(MPI_STUB_INCLUDE ${CMAKE_CURRENT_BINARY_DIR}/mpi_stub.h)
set(MAKE_MPI_STUB_INCLUDE make_mpi_stub_include)
set(MAKE_MPI_STUB_INCLUDE_SOURCE ${MAKE_MPI_STUB_INCLUDE}.f90)

set(RPN_STUB_SOURCE_FILES ${C_STUB_SOURCE} ${FORTRAN_STUB_SOURCE})

# Library targets and compile configuration
if (NOT RPN_SERIAL)
   set(STUB_LIBRARY ${PROJECT_NAME}_stubs)
   set(FULL_STUB_LIB_NAME ${STUB_LIBRARY}_${PROJECT_VERSION}${MPI_LIBRARY_VERSION})
   set(SHORT_STUB_LIB_FILE_NAME lib${STUB_LIBRARY}${MPI_LIBRARY_VERSION}.a)

   add_library(${STUB_LIBRARY} STATIC ${RPN_STUB_SOURCE_FILES})
   target_include_directories(${STUB_LIBRARY} PRIVATE ${CMAKE_CURRENT_BINARY_DIR})

   set_target_properties(${STUB_LIBRARY} PROPERTIES OUTPUT_NAME ${FULL_STUB_LIB_NAME})
   add_custom_command(TARGET ${STUB_LIBRARY} POST_BUILD
      COMMAND ln -sfv $<TARGET_FILE_NAME:${STUB_LIBRARY}> ${SHORT_STUB_LIB_FILE_NAME}
   )

   install(TARGETS ${STUB_LIBRARY})
   install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${SHORT_STUB_LIB_FILE_NAME} DESTINATION lib)
endif()

add_executable(${MAKE_MPI_STUB_INCLUDE} ${MAKE_MPI_STUB_INCLUDE_SOURCE})

# Generate the program that generates the header used by stubs
add_custom_command(
   OUTPUT ${MAKE_MPI_STUB_INCLUDE_SOURCE}
   COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/gen_make_mpi_stub_include.sh > ${MAKE_MPI_STUB_INCLUDE_SOURCE}
   DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/gen_make_mpi_stub_include.sh
)

# Generate the header used by stubs
add_custom_command(
   OUTPUT ${MPI_STUB_INCLUDE}
   COMMAND ${MAKE_MPI_STUB_INCLUDE} > ${MPI_STUB_INCLUDE}
   DEPENDS ${MAKE_MPI_STUB_INCLUDE}
)

# Setup dependencies on the generated files
set_property(SOURCE ${FORTRAN_STUB_SOURCE} APPEND PROPERTY OBJECT_DEPENDS ${MPI_STUB_INCLUDE})
set_property(SOURCE ${C_STUB_SOURCE} APPEND PROPERTY OBJECT_DEPENDS ${MPI_STUB_INCLUDE})


