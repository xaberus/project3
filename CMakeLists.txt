cmake_minimum_required (VERSION 2.8.6)
project(projekt)

include(FindPkgConfig)
find_package(Doxygen REQUIRED)
find_program(FFMPEG ffmpeg)
find_program(GNUPLOT gnuplot REQUIRED)
find_program(PDFLATEX pdflatex REQUIRED)
find_program(LUA lua REQUIRED)

set(BADNESS
  "badness"
)

set(PRESETS
  "simple"
  "square"
  "square1"
  "harmosca"
  ${BADNESS}
)

#set(CMAKE_VERBOSE_MAKEFILE true)

#include_directories()

pkg_check_modules(PKG REQUIRED "fftw3" "lua")
#pkg_check_modules(CAIRO REQUIRED "cairo")

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

add_library(projekt STATIC
  ${projekt_SOURCE_DIR}/simulation.c
  ${projekt_SOURCE_DIR}/splitop.c
  ${projekt_SOURCE_DIR}/numerov.c
  ${projekt_SOURCE_DIR}/carray.c
  ${projekt_SOURCE_DIR}/array.c
  ${projekt_SOURCE_DIR}/peaks.c
)
set_target_properties(projekt PROPERTIES COMPILE_FLAGS ${CFLAGS})
set_target_properties(projekt PROPERTIES LINK_FLAGS ${LFLAGS})

add_executable(simulate ${projekt_SOURCE_DIR}/simulate.c)
set_target_properties(simulate PROPERTIES COMPILE_FLAGS ${CFLAGS})
set_target_properties(simulate PROPERTIES LINK_FLAGS ${LFLAGS})
target_link_libraries(simulate ${PKG_LIBRARIES} fftw3_threads m projekt)

add_executable(evaluate ${projekt_SOURCE_DIR}/evaluate.c)
set_target_properties(evaluate PROPERTIES COMPILE_FLAGS ${CFLAGS})
set_target_properties(evaluate PROPERTIES LINK_FLAGS ${LFLAGS})
target_link_libraries(evaluate ${PKG_LIBRARIES} fftw3_threads m projekt)

add_executable(testpeaks ${projekt_SOURCE_DIR}/testpeaks.c)
set_target_properties(testpeaks PROPERTIES COMPILE_FLAGS ${CFLAGS})
target_link_libraries(testpeaks ${PKG_LIBRARIES} m projekt)

#add_executable(tabgen
#  ${projekt_SOURCE_DIR}/tabgen.c
#)
#set_target_properties(tabgen PROPERTIES COMPILE_FLAGS ${CFLAGS})
#target_link_libraries(tabgen ${PKG_LIBRARIES} m)

foreach(PRESET ${PRESETS})
  message("adding ${PRESET}")

  add_custom_command(
    OUTPUT
      "${projekt_SOURCE_DIR}/${PRESET}/result.dat"
    COMMAND
      $<TARGET_FILE:simulate>
        "${projekt_SOURCE_DIR}/${PRESET}.lua"
        "${projekt_SOURCE_DIR}/${PRESET}/result.dat"
    MAIN_DEPENDENCY
      "${projekt_SOURCE_DIR}/${PRESET}.lua")

  if(${CAIRO_FOUND})
    add_custom_command(
      OUTPUT
        "${projekt_SOURCE_DIR}/${PRESET}/apsi.dat"
        "${projekt_SOURCE_DIR}/${PRESET}/corr.dat"
        "${projekt_SOURCE_DIR}/${PRESET}/dftcorr.dat"
        "${projekt_SOURCE_DIR}/${PRESET}/pot.dat"
        "${projekt_SOURCE_DIR}/${PRESET}/spectrum.dat"
        "${projekt_SOURCE_DIR}/${PRESET}/numen.dat"
        "${projekt_SOURCE_DIR}/${PRESET}/stats"
        "${projekt_SOURCE_DIR}/${PRESET}/video/image00000.png"
      COMMAND
        $<TARGET_FILE:evaluate>
          "${projekt_SOURCE_DIR}/${PRESET}.lua"
          "${projekt_SOURCE_DIR}/${PRESET}/result.dat"
      MAIN_DEPENDENCY
        "${projekt_SOURCE_DIR}/${PRESET}.lua"
        "${projekt_SOURCE_DIR}/${PRESET}/result.dat")
  else()
    add_custom_command(
      OUTPUT
        "${projekt_SOURCE_DIR}/${PRESET}/apsi.dat"
        "${projekt_SOURCE_DIR}/${PRESET}/corr.dat"
        "${projekt_SOURCE_DIR}/${PRESET}/dftcorr.dat"
        "${projekt_SOURCE_DIR}/${PRESET}/pot.dat"
        "${projekt_SOURCE_DIR}/${PRESET}/spectrum.dat"
        "${projekt_SOURCE_DIR}/${PRESET}/numen.dat"
        "${projekt_SOURCE_DIR}/${PRESET}/stats"
      COMMAND
        $<TARGET_FILE:evaluate>
          "${projekt_SOURCE_DIR}/${PRESET}.lua"
          "${projekt_SOURCE_DIR}/${PRESET}/result.dat"
      MAIN_DEPENDENCY
        "${projekt_SOURCE_DIR}/${PRESET}.lua"
        "${projekt_SOURCE_DIR}/${PRESET}/result.dat")
  endif()

  set(CONFS ${CONFS} "${projekt_SOURCE_DIR}/${PRESET}.lua")

  add_custom_command(
    OUTPUT
      "${projekt_SOURCE_DIR}/${PRESET}/spectrum.tex.inc"
    COMMAND
      ${LUA} "${projekt_SOURCE_DIR}/tabgen.lua"
        "${projekt_SOURCE_DIR}/${PRESET}.lua"
        "${projekt_SOURCE_DIR}/${PRESET}/spectrum.dat"
        "${projekt_SOURCE_DIR}/${PRESET}/numen.dat"
        >"${projekt_SOURCE_DIR}/${PRESET}/spectrum.tex.inc"
    MAIN_DEPENDENCY
      "${projekt_SOURCE_DIR}/tabgen.lua"
      "${projekt_SOURCE_DIR}/${PRESET}/spectrum.dat"
      "${projekt_SOURCE_DIR}/${PRESET}/numen.dat")

  set(SPECTRA ${SPECTRA} "${projekt_SOURCE_DIR}/${PRESET}/spectrum.tex.inc")

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
      "${projekt_SOURCE_DIR}/${PRESET}/dftcorr.dat")

###########################################################################

  if(${CAIRO_FOUND})
    add_custom_target(${PRESET}
      DEPENDS
        "${projekt_SOURCE_DIR}/${PRESET}.lua"
        "${projekt_SOURCE_DIR}/${PRESET}/apsi.dat"
        "${projekt_SOURCE_DIR}/${PRESET}/corr.dat"
        "${projekt_SOURCE_DIR}/${PRESET}/dftcorr.dat"
        "${projekt_SOURCE_DIR}/${PRESET}/pot.dat"
        "${projekt_SOURCE_DIR}/${PRESET}/spectrum.dat"
        "${projekt_SOURCE_DIR}/${PRESET}/movie.webm"
        "${projekt_SOURCE_DIR}/${PRESET}/numen.dat"
        "${projekt_SOURCE_DIR}/${PRESET}/spectrum.tex.inc"

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
        "${projekt_SOURCE_DIR}/${PRESET}/spectrum.dat"
        "${projekt_SOURCE_DIR}/${PRESET}/numen.dat"
        "${projekt_SOURCE_DIR}/${PRESET}/spectrum.tex.inc"

        "${projekt_SOURCE_DIR}/${PRESET}/init.pdf"
        "${projekt_SOURCE_DIR}/${PRESET}/corr.pdf"
        "${projekt_SOURCE_DIR}/${PRESET}/spektrum.pdf")
  endif()
  set(PLOTS ${PLOTS}
    "${projekt_SOURCE_DIR}/${PRESET}/init.pdf"
    "${projekt_SOURCE_DIR}/${PRESET}/corr.pdf"
    "${projekt_SOURCE_DIR}/${PRESET}/spektrum.pdf")
endforeach()

foreach(BAD ${BADNESS})
  message("adding ${BAD}")

  add_custom_command(
    OUTPUT
      "${projekt_SOURCE_DIR}/${BAD}/bspektrum.pdf"
    COMMAND
      ${GNUPLOT} "${projekt_SOURCE_DIR}/bspektrum.gp"
    WORKING_DIRECTORY
      "${projekt_SOURCE_DIR}/${BAD}/"
    DEPENDS
      "${projekt_SOURCE_DIR}/bspektrum.gp"
      "${projekt_SOURCE_DIR}/${BAD}/dftcorr.dat"
      "${projekt_SOURCE_DIR}/${BAD}/spectrum.dat"
      "${projekt_SOURCE_DIR}/${BAD}/numen.dat")
  set(PLOTS ${PLOTS} "${projekt_SOURCE_DIR}/${BAD}/bspektrum.pdf")
  set(BDEP
    ${BDEP} "${projekt_SOURCE_DIR}/${BAD}/bspektrum.pdf"
    "${projekt_SOURCE_DIR}/${BAD}/spectrum.tex.inc")
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
    echo "remove the comments for release!!!!"
  COMMAND
    sleep 1
  COMMAND
    ${PDFLATEX} "${projekt_SOURCE_DIR}/doc/projekt3.tex"
  WORKING_DIRECTORY
    "${projekt_SOURCE_DIR}/doc"
  DEPENDS
    "${projekt_SOURCE_DIR}/doc/projekt3.tex"
    #${PLOTS}
    ${SPECTRA}
    ${CONFS}
  COMMENT
    "Generating LaTeX doc"
  VERBATIM)

add_custom_target(bad
  DEPENDS ${BDEP})

add_custom_target(doc
  DEPENDS
    "${projekt_SOURCE_DIR}/adoc"
    "${projekt_SOURCE_DIR}/doc/projekt3.pdf"
  COMMENT
    "Generating documentation")
