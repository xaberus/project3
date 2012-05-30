cmake_minimum_required (VERSION 2.8.6)

project(projekt)

include(FindPkgConfig)

#set(CMAKE_VERBOSE_MAKEFILE true)

#include_directories()

#pkg_check_modules(PKG REQUIRED "cairo" "fftw3")
pkg_check_modules(PKG REQUIRED "fftw3")

set(CFLAGS_LIST
  "-std=gnu99"
  "-Wall"
  "-Wextra"
  #"-Werror"
  "-pedantic"
  "-ggdb"
  "-O0"
  #"-O1"
  #"-O3"
  "-msse4"
  #"-fno-inline"
  "-fPIC"
  #"-Wno-braced-groups"
  #"-fprofile-arcs"
  #"-ftest-coverage"
  #"-fcatch-undefined-behavior"
  #"-ftrapv"
  ${PKG_CFLAGS}
)

add_definitions(-DTEST)
add_definitions(-D_XOPEN_SOURCE=700)
add_definitions(-D_GNU_SOURCE)
#add_definitions(-DTESTPROF="cpuprof")
add_definitions(-DBROOT="${CMAKE_CURRENT_SOURCE_DIR}")

set(LFLAGS_LIST
  "-Wl,--no-undefined"
  "-Wl,--warn-unresolved-symbols"
  "-pthread"
  "-rdynamic"
  ${PKG_LDFLAGS}
)

include_directories(${libnih_SOURCE_DIR})

string (REPLACE ";" " " CFLAGS "${CFLAGS_LIST}")
string (REPLACE ";" " " LFLAGS "${LFLAGS_LIST}")

add_executable(projekt3
  ${projekt_SOURCE_DIR}/hashmap.c
  ${projekt_SOURCE_DIR}/projekt.c
  ${projekt_SOURCE_DIR}/splitop.c
  ${projekt_SOURCE_DIR}/carray.c
  ${projekt_SOURCE_DIR}/expr_parser.c
  ${projekt_SOURCE_DIR}/lexer.c
  ${projekt_SOURCE_DIR}/expr.c
)
set_target_properties(projekt3 PROPERTIES COMPILE_FLAGS ${CFLAGS})
set_target_properties(projekt3 PROPERTIES LINK_FLAGS ${LFLAGS})
target_link_libraries(projekt3 ${PKG_LIBRARIES} fftw3_threads m)

#add_custom_target(chop
#  $<TARGET_FILE:butcher> $<TARGET_FILE:nih> --bexec $<TARGET_FILE:bexec> $(BFLAGS)
#  DEPENDS nih butcher bexec
#)
