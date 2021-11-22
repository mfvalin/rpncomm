cmake_minimum_required(VERSION 3.16)

message(STATUS "(EC) Generating rpncomm Makefile")

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake_rpn)

include(ec_init)
include(ec_parse_manifest)

ec_parse_manifest()

project(${NAME} VERSION ${VERSION} DESCRIPTION "${DESCRIPTION}") 

set(RPN_COMM_VERSION_STRING "${VERSION}")
string(REPLACE "." "" RPN_COMM_VERSION_MARKER ${VERSION})

option(RPN_SERIAL "Whether the library will avoid MPI (not fully implemented, DO NOT USE)")

if (RPN_SERIAL)
   message(WARNING "(EC) The RPN_SERIAL option is not fully implemented. It will be deactivated.")
   set(RPN_SERIAL FALSE)
endif()

enable_language(Fortran)
enable_language(C)
include(ec_compiler_presets)

find_package(MPI)
if (MPI_FOUND)
   set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
   set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
endif()

if (NOT RPN_SERIAL)
   find_program(RPN_MPIRUN mpirun)
   if (${RPN_MPIRUN} STREQUAL "RPN_MPIRUN-NOTFOUND")
      message(FATAL_ERROR "(EC) mpirun not found. Won't be able to configure the version.")
   else()
      execute_process(
         COMMAND ${RPN_MPIRUN} --version 2>&1
         COMMAND grep -w MPI
         COMMAND sed -e "s/.* //" -e "s/)//"
         OUTPUT_VARIABLE MPI_LIBRARY_VERSION
         OUTPUT_STRIP_TRAILING_WHITESPACE)

     set(MPI_LIBRARY_VERSION "_${MPI_LIBRARY_VERSION}")
   endif()
endif()

add_subdirectory(src)

