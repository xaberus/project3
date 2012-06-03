cmake_minimum_required (VERSION 2.8.6)
project(projekt)
include(FindPkgConfig)
find_package(Doxygen REQUIRED)
find_program(FFMPEG ffmpeg)
find_program(GNUPLOT gnuplot)
find_program(PDFLATEX pdflatex)

set(PRESETS
  "simple"
  "square"
  "harmosca"
)

#set(CMAKE_VERBOSE_MAKEFILE true)

#include_directories()

pkg_check_modules(PKG REQUIRED "fftw3" "lua")
pkg_check_modules(CAIRO REQUIRED "cairo")

set(CFLAGS_LIST
  "-std=gnu99"
  "-Wall"
  "-Wextra"
  "-Werror"
  "-pedantic"
  "-ggdb"
  #"-O0"
  #"-O1"
  "-O3"
  "-msse4"
  #"-fno-inline"
  #"-fPIC"
  #"-Wno-braced-groups"
  #"-fprofile-arcs"
  #"-ftest-coverage"
  #"-fcatch-undefined-behavior"
  #"-ftrapv"
  ${PKG_CFLAGS}
  ${CAIRO_CFLAGS}
)

add_definitions(-DTEST)
add_definitions(-D_XOPEN_SOURCE=700)
add_definitions(-D_GNU_SOURCE)

if(${CAIRO_FOUND})
  add_definitions(-DUSE_CAIRO)
endif()

#add_definitions(-DTESTPROF="cpuprof")
add_definitions(-DBROOT="${CMAKE_CURRENT_SOURCE_DIR}")

set(LFLAGS_LIST
  "-Wl,--no-undefined"
  "-Wl,--warn-unresolved-symbols"
  "-pthread"
  "-rdynamic"
  ${PKG_LDFLAGS}
  ${CAIRO_LDFLAGS}
)

include_directories(${libnih_SOURCE_DIR})

string (REPLACE ";" " " CFLAGS "${CFLAGS_LIST}")
string (REPLACE ";" " " LFLAGS "${LFLAGS_LIST}")

add_executable(projekt3
  ${projekt_SOURCE_DIR}/simulation.c
  ${projekt_SOURCE_DIR}/projekt.c
  ${projekt_SOURCE_DIR}/splitop.c
  ${projekt_SOURCE_DIR}/carray.c
  ${projekt_SOURCE_DIR}/array.c
)
set_target_properties(projekt3 PROPERTIES COMPILE_FLAGS ${CFLAGS})
set_target_properties(projekt3 PROPERTIES LINK_FLAGS ${LFLAGS})
target_link_libraries(projekt3 ${PKG_LIBRARIES} fftw3_threads m)

foreach(PRESET ${PRESETS})
  message("adding ${PRESET}")
  if(${CAIRO_FOUND})
    add_custom_command(
      OUTPUT
        "${projekt_SOURCE_DIR}/${PRESET}/apsi.dat"
        "${projekt_SOURCE_DIR}/${PRESET}/corr.dat"
        "${projekt_SOURCE_DIR}/${PRESET}/dftcorr.dat"
        "${projekt_SOURCE_DIR}/${PRESET}/pot.dat"
        "${projekt_SOURCE_DIR}/${PRESET}/theoenrg.dat"
        "${projekt_SOURCE_DIR}/${PRESET}/video/image00000.png"
      COMMAND
        mkdir -p "${projekt_SOURCE_DIR}/${PRESET}"
      COMMAND
        $<TARGET_FILE:projekt3> "${projekt_SOURCE_DIR}/${PRESET}.lua"
      MAIN_DEPENDENCY
        "${projekt_SOURCE_DIR}/${PRESET}.lua")
  else()
    add_custom_command(
      OUTPUT
        "${projekt_SOURCE_DIR}/${PRESET}/apsi.dat"
        "${projekt_SOURCE_DIR}/${PRESET}/corr.dat"
        "${projekt_SOURCE_DIR}/${PRESET}/dftcorr.dat"
        "${projekt_SOURCE_DIR}/${PRESET}/pot.dat"
        "${projekt_SOURCE_DIR}/${PRESET}/theoenrg.dat"
      COMMAND
        mkdir -p "${projekt_SOURCE_DIR}/${PRESET}"
      COMMAND
        $<TARGET_FILE:projekt3> "${projekt_SOURCE_DIR}/${PRESET}.lua"
      MAIN_DEPENDENCY
        "${projekt_SOURCE_DIR}/${PRESET}.lua")
  endif()

  if(${CAIRO_FOUND})
    add_custom_command(
      OUTPUT
        "${projekt_SOURCE_DIR}/${PRESET}/movie.webm"
      COMMAND
        ${FFMPEG} -y -r 20 -b 9600 -i
          "${projekt_SOURCE_DIR}/${PRESET}/video/image%05d.png"
          "${projekt_SOURCE_DIR}/${PRESET}/movie.webm"
      MAIN_DEPENDENCY
        "${projekt_SOURCE_DIR}/${PRESET}/video/image00000.png")
  endif()

###########################################################################

  #"${projekt_SOURCE_DIR}/${PRESET}/init.png"
  add_custom_command(
    OUTPUT
      "${projekt_SOURCE_DIR}/${PRESET}/init.pdf"
    COMMAND
      ${GNUPLOT} "${projekt_SOURCE_DIR}/init.gp"
    WORKING_DIRECTORY
      "${projekt_SOURCE_DIR}/${PRESET}/"
    DEPENDS
      "${projekt_SOURCE_DIR}/init.gp"
      "${projekt_SOURCE_DIR}/${PRESET}/apsi.dat"
      "${projekt_SOURCE_DIR}/${PRESET}/pot.dat")

  add_custom_command(
    OUTPUT
      "${projekt_SOURCE_DIR}/${PRESET}/corr.pdf"
    COMMAND
      ${GNUPLOT} "${projekt_SOURCE_DIR}/corr.gp"
    WORKING_DIRECTORY
      "${projekt_SOURCE_DIR}/${PRESET}/"
    DEPENDS
      "${projekt_SOURCE_DIR}/corr.gp"
      "${projekt_SOURCE_DIR}/${PRESET}/corr.dat")

  add_custom_command(
    OUTPUT
      "${projekt_SOURCE_DIR}/${PRESET}/spektrum.pdf"
    COMMAND
      ${GNUPLOT} "${projekt_SOURCE_DIR}/spektrum.gp"
    WORKING_DIRECTORY
      "${projekt_SOURCE_DIR}/${PRESET}/"
    DEPENDS
      "${projekt_SOURCE_DIR}/spektrum.gp"
      "${projekt_SOURCE_DIR}/${PRESET}/dftcorr.dat"
      "${projekt_SOURCE_DIR}/${PRESET}/theoenrg.dat")

###########################################################################

  if(${CAIRO_FOUND})
    add_custom_target(${PRESET}
      DEPENDS
        "${projekt_SOURCE_DIR}/${PRESET}/apsi.dat"
        "${projekt_SOURCE_DIR}/${PRESET}/corr.dat"
        "${projekt_SOURCE_DIR}/${PRESET}/dftcorr.dat"
        "${projekt_SOURCE_DIR}/${PRESET}/pot.dat"
        "${projekt_SOURCE_DIR}/${PRESET}/theoenrg.dat"
        "${projekt_SOURCE_DIR}/${PRESET}/movie.webm"

        "${projekt_SOURCE_DIR}/${PRESET}/init.pdf"
        "${projekt_SOURCE_DIR}/${PRESET}/corr.pdf"
        "${projekt_SOURCE_DIR}/${PRESET}/spektrum.pdf")
    else()
    add_custom_target(${PRESET}
      DEPENDS
        "${projekt_SOURCE_DIR}/${PRESET}.lua"
        "${projekt_SOURCE_DIR}/${PRESET}/apsi.dat"
        "${projekt_SOURCE_DIR}/${PRESET}/corr.dat"
        "${projekt_SOURCE_DIR}/${PRESET}/dftcorr.dat"
        "${projekt_SOURCE_DIR}/${PRESET}/pot.dat"
        "${projekt_SOURCE_DIR}/${PRESET}/theoenrg.dat"

        "${projekt_SOURCE_DIR}/${PRESET}/init.pdf"
        "${projekt_SOURCE_DIR}/${PRESET}/corr.pdf"
        "${projekt_SOURCE_DIR}/${PRESET}/spektrum.pdf")
    endif()
endforeach()

add_custom_command(
  OUTPUT
    "${projekt_SOURCE_DIR}/adoc"
  COMMAND
    ${DOXYGEN_EXECUTABLE} "${projekt_SOURCE_DIR}/Doxyfile"
  WORKING_DIRECTORY
    "${projekt_SOURCE_DIR}"
  COMMENT
    "Generating API documentation with Doxygen"
  VERBATIM)

add_custom_command(
  OUTPUT
    "${projekt_SOURCE_DIR}/doc/projekt3.pdf"
  COMMAND
    ${PDFLATEX} "${projekt_SOURCE_DIR}/doc/projekt3.tex"
  WORKING_DIRECTORY
    "${projekt_SOURCE_DIR}/doc"
  DEPENDS
    "${projekt_SOURCE_DIR}/doc/projekt3.tex"
  COMMENT
    "Generating LaTeX doc"
  VERBATIM)

add_custom_target(doc
  DEPENDS
    "${projekt_SOURCE_DIR}/adoc"
    "${projekt_SOURCE_DIR}/doc/projekt3.pdf"
  COMMENT
    "Generating documentation")