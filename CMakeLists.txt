# SWAT+ cmake, req min v3.20 (ifx, etc. support)
###############################################################
if("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}" OR
   "${CMAKE_SOURCE_DIR}/src" STREQUAL "${CMAKE_BINARY_DIR}")
   message(FATAL_ERROR "'in-source' builds not permitted, e.g. use 'cmake -B build' instead.")
endif()

cmake_minimum_required(VERSION 3.20)

project(swatplus VERSION 0.2 LANGUAGES Fortran)

include(CTest)

if(NOT TAG)
    find_package(Git REQUIRED)
    execute_process(
      COMMAND ${GIT_EXECUTABLE} describe --tags
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
      OUTPUT_VARIABLE TAG
      OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    if(TAG STREQUAL "")
    # message(FATAL_ERROR "No tag info found.")
    # this happens when git is not installed or not in the path
      set(TAG "unknown")
    endif()
endif()


if (UNIX)
    if(CMAKE_Fortran_COMPILER_ID STREQUAL Intel)
        set(fdialect        "-free -fpe0 -traceback -diag-disable=10448")
        set(fdebug          "-warn all")
        set(frelease        "-O")
        set(FFC             "ifo")
        link_libraries("-static")
    elseif(CMAKE_Fortran_COMPILER_ID STREQUAL IntelLLVM)
        set(fdialect        "-free -fpe0 -traceback")
        # set(fdebug          "-warn all -O0")
        set(fdebug          "-O0")
        set(frelease        "-O")
        set(FFC             "ifx")
        link_libraries("-static")
    elseif(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
        set(fdialect        "-fcheck=all -ffpe-trap=invalid,zero,overflow,underflow -fimplicit-none -ffree-line-length-none -fbacktrace -finit-local-zero -fno-unsafe-math-optimizations -frounding-math -fsignaling-nans")
        set(fdebug          "-Wall")
        set(frelease        "-O")
        set(FFC             "gnu")
        if(NOT APPLE)
            link_libraries("-static")
        endif()
    endif()
elseif(WIN32)
    if(CMAKE_Fortran_COMPILER_ID STREQUAL Intel)
        set(fdialect        "/free /fpe0 /traceback /Qdiag-disable:all /Qdiag-disable:remarks")
        set(fdebug          "/warn:all")
        set(frelease        "/O")
        set(FFC             "ifo")
        link_libraries("-static")
    elseif(CMAKE_Fortran_COMPILER_ID STREQUAL IntelLLVM)
        set(fdialect        "/free /fpe0 /traceback")
        set(fdebug          "/warn:all")
        set(frelease        "/O")
        set(FFC             "ifx")
        link_libraries("-static")
    elseif(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
        set(fdialect        "-fcheck=all -ffpe-trap=invalid,zero,overflow,underflow -fimplicit-none -ffree-line-length-none -fbacktrace -finit-local-zero -fno-unsafe-math-optimizations -frounding-math -fsignaling-nans")
        set(fdebug          "-Wall ")
        set(FFC             "gnu")
        set(frelease        "-O")
    endif()
endif()

string(TOLOWER ${CMAKE_HOST_SYSTEM_PROCESSOR} ARCH)
if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    set(AR "lin_${ARCH}")
elseif(CMAKE_SYSTEM_NAME STREQUAL "Windows")
    set(AR "win_${ARCH}")
elseif(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
    set(AR "mac_${ARCH}")
else()
    set(AR "unknown")
endif()

if(CMAKE_BUILD_TYPE STREQUAL "Release")
    set(TY "-Rel")
elseif(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(TY "-Dbg")
else()
    set(TY "")
endif()

# SWAT Version number
set(SWAT_VERSION ${TAG})
set(SWATPLUS_EXE "swatplus-${SWAT_VERSION}-${FFC}-${AR}${TY}")

# Set this to 'TRUE' to see the fortran command on compile
set(CMAKE_VERBOSE_MAKEFILE FALSE)

set(CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} ${fdialect}")
set(CMAKE_Fortran_FLAGS_DEBUG   "${CMAKE_Fortran_FLAGS_DEBUG} ${fdebug}")
# set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} ${frelease}")
set(CMAKE_Fortran_FLAGS_RELEASE "${frelease}")

#############################################################################
# Build
# list (SORT _variableNames)
# foreach (_variableName ${_variableNames})
#     message(STATUS "${_variableName}=${${_variableName}}")
# endforeach()

if(EXISTS "${PROJECT_SOURCE_DIR}/src/main.f90.in")
    string(TIMESTAMP ISO   "%Y-%m-%d %H:%M:%S")
    string(TIMESTAMP TODAY "%Y-%m-%d")  # e.g. produces 2023-12-07
    string(TIMESTAMP YEAR  "%Y")

    configure_file(
        "${PROJECT_SOURCE_DIR}/src/main.f90.in"
        "${PROJECT_SOURCE_DIR}/src/main.f90"
    )
endif()

if(EXISTS "${PROJECT_SOURCE_DIR}/ford.md.in")
    configure_file(
        "${PROJECT_SOURCE_DIR}/ford.md.in"
        "${PROJECT_SOURCE_DIR}/ford.md"
    )   
endif()


file(GLOB sources src/*.f90)
add_executable(${SWATPLUS_EXE} ${sources})

#############################################################################
# Install
install(TARGETS ${SWATPLUS_EXE} DESTINATION RUNTIME DESTINATION ${PROJECT_SOURCE_DIR})


#############################################################################
# Testing
find_package(Python3 REQUIRED)

set(spcheck "${PROJECT_SOURCE_DIR}/test/spcheck.py")
set(exe_path "${PROJECT_BINARY_DIR}/${SWATPLUS_EXE}")
set(test_dir "${PROJECT_BINARY_DIR}/data")
set(ref_dir  "${PROJECT_SOURCE_DIR}/data")

# error tolerances
set(rerr "0.01")
set(aerr "1e-8")

add_test(Ames_sub1   ${Python3_EXECUTABLE} ${spcheck} ctest ${exe_path} ${ref_dir}/Ames_sub1 ${test_dir} --abserr ${aerr} --relerr ${rerr})
# add_test(Ithaca_sub6 python3 ${check_py} ${exe_path} ${ref_dir}/Ithaca_sub6 ${test_dir} ${aerr} ${rerr})


#############################################################################
# Summary
message( "-- Verbose Makefile   ${CMAKE_VERBOSE_MAKEFILE}")
message( "-- Generator          ${CMAKE_GENERATOR_TOOLSET}")
message( "-- VS Platform        ${CMAKE_VS_PLATFORM_TOOLSET}")
message( "-- Build Type         ${CMAKE_BUILD_TYPE}")
message( "-- Compiler           ${CMAKE_Fortran_COMPILER} (${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION})")
message( "-- FFLAGS            ${CMAKE_Fortran_FLAGS}")
message( "--  RELEASE           ${CMAKE_Fortran_FLAGS_RELEASE}")
message( "--  DEBUG             ${CMAKE_Fortran_FLAGS_DEBUG}")
message( "-- PROJECT_DIR        ${PROJECT_SOURCE_DIR}")
message( "-- BINARY_DIR         ${PROJECT_BINARY_DIR}")
message( "-- Install Prefix     ${CMAKE_INSTALL_PREFIX}")
message( "-- Executable         ${SWATPLUS_EXE}, rev/tag: ${SWAT_VERSION}")
message( "-- FORD project       ${PROJECT_SOURCE_DIR}/ford.md")
message( "-- Test Data Dir      ${ref_dir}")
message( "--  Abs Error         ${aerr}")
message( "--  Rel Error         ${rerr}")
