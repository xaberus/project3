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

set(CMP "cmp" "cmp1")

set(PRESETS
  "test"

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
  ${projekt_SOURCE_DIR}/source/simulation.c
  ${projekt_SOURCE_DIR}/source/splitop.c
  ${projekt_SOURCE_DIR}/source/numerov.c
  ${projekt_SOURCE_DIR}/source/carray.c
  ${projekt_SOURCE_DIR}/source/array.c
  ${projekt_SOURCE_DIR}/source/peaks.c
  ${projekt_SOURCE_DIR}/source/akima.c
  ${projekt_SOURCE_DIR}/source/squares.c
)
set_target_properties(projekt PROPERTIES COMPILE_FLAGS ${CFLAGS})
set_target_properties(projekt PROPERTIES LINK_FLAGS ${LFLAGS})

add_executable(simulate ${projekt_SOURCE_DIR}/source/simulate.c)
set_target_properties(simulate PROPERTIES COMPILE_FLAGS ${CFLAGS})
set_target_properties(simulate PROPERTIES LINK_FLAGS ${LFLAGS})
target_link_libraries(simulate ${PKG_LIBRARIES} fftw3_threads m projekt)

add_executable(evaluate ${projekt_SOURCE_DIR}/source/evaluate.c)
set_target_properties(evaluate PROPERTIES COMPILE_FLAGS ${CFLAGS})
set_target_properties(evaluate PROPERTIES LINK_FLAGS ${LFLAGS})
target_link_libraries(evaluate ${PKG_LIBRARIES} fftw3_threads m projekt)

#add_executable(tabgen
#  ${projekt_SOURCE_DIR}/tabgen.c
#)
#set_target_properties(tabgen PROPERTIES COMPILE_FLAGS ${CFLAGS})
#target_link_libraries(tabgen ${PKG_LIBRARIES} m)

set(LUA_SCRIPTS "${projekt_SOURCE_DIR}/scripts/?.lua;${projekt_SOURCE_DIR}/presets/?.inc")
set(SCRIPTS "${projekt_SOURCE_DIR}/scripts")

macro(get_lua_var cfg var cvar)
  execute_process(
    COMMAND
      env "LUA_PATH=${LUA_SCRIPTS}" ${LUA} -e "loadfile(\"${cfg}\")(); print(assert(${var}))"
    OUTPUT_VARIABLE _get_var
    RESULT_VARIABLE _get_result
    OUTPUT_STRIP_TRAILING_WHITESPACE)
  if (_get_result)
    message(SEND_ERROR "Could not get ${LUA} variable '${var}' from ${cfg}")
  else()
    #message("${cfg}:${var} = '${_get_var}'")
    set(${cvar} ${_get_var})
  endif()
endmacro()

foreach(PRESET ${PRESETS})
  message("adding ${PRESET}")

  set(PRESET_LUA "${projekt_SOURCE_DIR}/presets/${PRESET}.lua")

  get_lua_var(${PRESET_LUA} "config.output.dir" CONFIG_DIR)
  get_lua_var(${PRESET_LUA} "config.output.apsi" CONFIG_APSI)
  get_lua_var(${PRESET_LUA} "config.output.pot" CONFIG_POT)
  get_lua_var(${PRESET_LUA} "config.output.corr" CONFIG_CORR)
  get_lua_var(${PRESET_LUA} "config.output.dftcorr" CONFIG_DFTCORR)
  get_lua_var(${PRESET_LUA} "config.output.spectrum" CONFIG_SPECTRUM)
  get_lua_var(${PRESET_LUA} "config.output.numen" CONFIG_NUMEN)
  get_lua_var(${PRESET_LUA} "config.output.splen" CONFIG_SPLEN)
  get_lua_var(${PRESET_LUA} "config.output.aken" CONFIG_AKEN)

  set(PRESET_DIR "${projekt_SOURCE_DIR}/${CONFIG_DIR}")

  add_custom_command(
    OUTPUT
      "${PRESET_DIR}/result.dat"
    COMMAND
      mkdir -p "${PRESET_DIR}"
    COMMAND
      env
        "'LUA_PATH=${LUA_SCRIPTS}'"
        $<TARGET_FILE:simulate>
        ${PRESET_LUA}
        "${PRESET_DIR}/result.dat"
    MAIN_DEPENDENCY
      ${PRESET_LUA})

  if(${CAIRO_FOUND})
    add_custom_command(
      OUTPUT
        "${PRESET_DIR}/${CONFIG_APSI}"
        "${PRESET_DIR}/${CONFIG_CORR}"
        "${PRESET_DIR}/${CONFIG_DFTCORR}"
        "${PRESET_DIR}/${CONFIG_POT}"
        "${PRESET_DIR}/${CONFIG_SPECTRUM}"
        "${PRESET_DIR}/${CONFIG_NUMEN}"
        "${PRESET_DIR}/stats"
        "${PRESET_DIR}/video/image00000.png"
      COMMAND
        env "'LUA_PATH=${LUA_SCRIPTS}'"
          $<TARGET_FILE:evaluate>
          ${PRESET_LUA}
          "${PRESET_DIR}/result.dat"
      MAIN_DEPENDENCY
        ${PRESET_LUA}
        "${PRESET_DIR}/result.dat")
  else()
    add_custom_command(
      OUTPUT
        "${PRESET_DIR}/${CONFIG_APSI}"
        "${PRESET_DIR}/${CONFIG_CORR}"
        "${PRESET_DIR}/${CONFIG_DFTCORR}"
        "${PRESET_DIR}/${CONFIG_POT}"
        "${PRESET_DIR}/${CONFIG_SPECTRUM}"
        "${PRESET_DIR}/${CONFIG_NUMEN}"
        "${PRESET_DIR}/stats"
      COMMAND
        env "'LUA_PATH=${LUA_SCRIPTS}'" $<TARGET_FILE:evaluate>
          ${PRESET_LUA}
          "${PRESET_DIR}/result.dat"
      MAIN_DEPENDENCY
        ${PRESET_LUA}
        "${PRESET_DIR}/result.dat")
  endif()

  set(CONFS ${CONFS} ${PRESET_LUA})

  add_custom_command(
    OUTPUT
      "${PRESET_DIR}/spectrum.tex.inc"
    COMMAND
      env "'LUA_PATH=${LUA_SCRIPTS}'" ${LUA} "${SCRIPTS}/tabgen.lua"
        ${PRESET_LUA}
        "${PRESET_DIR}/${CONFIG_SPECTRUM}"
        "${PRESET_DIR}/${CONFIG_NUMEN}"
        >"${PRESET_DIR}/spectrum.tex.inc"
    MAIN_DEPENDENCY
      "${projekt_SOURCE_DIR}/tabgen.lua"
      "${PRESET_DIR}/${CONFIG_SPECTRUM}"
      "${PRESET_DIR}/${CONFIG_NUMEN}")

  set(SPECTRA ${SPECTRA} "${PRESET_DIR}/spectrum.tex.inc")

  if(${CAIRO_FOUND})
    add_custom_command(
      OUTPUT
        "${PRESET_DIR}/movie.webm"
      COMMAND
        ${FFMPEG} -y -r 20 -b 9600 -i
          "${PRESET_DIR}/video/image%05d.png"
          "${PRESET_DIR}/movie.webm"
      MAIN_DEPENDENCY
        "${PRESET_DIR}/video/image00000.png")
  endif()

###########################################################################

  #"${PRESET_DIR}/init.png"
  add_custom_command(
    OUTPUT
      "${PRESET_DIR}/init.pdf"
    COMMAND
      ${GNUPLOT} "${SCRIPTS}/init.gp"
    WORKING_DIRECTORY
      "${PRESET_DIR}/"
    DEPENDS
      "${SCRIPTS}/init.gp"
      "${PRESET_DIR}/${CONFIG_APSI}"
      "${PRESET_DIR}/${CONFIG_POT}")

  add_custom_command(
    OUTPUT
      "${PRESET_DIR}/corr.pdf"
    COMMAND
      ${GNUPLOT} "${SCRIPTS}/corr.gp"
    WORKING_DIRECTORY
      "${PRESET_DIR}/"
    DEPENDS
      "${SCRIPTS}/corr.gp"
      "${PRESET_DIR}/${CONFIG_CORR}")

  add_custom_command(
    OUTPUT
      "${PRESET_DIR}/spektrum.pdf"
    COMMAND
      ${GNUPLOT} "${SCRIPTS}/spektrum.gp"
    WORKING_DIRECTORY
      "${PRESET_DIR}/"
    DEPENDS
      "${SCRIPTS}/spektrum.gp"
      "${PRESET_DIR}/${CONFIG_DFTCORR}")

###########################################################################

  add_custom_target(${PRESET}_compile
    DEPENDS
      "${PRESET_DIR}/result.dat")

  if(${CAIRO_FOUND})
    add_custom_target(${PRESET}
      DEPENDS
        ${PRESET_LUA}
        "${PRESET_DIR}/${CONFIG_APSI}"
        "${PRESET_DIR}/${CONFIG_CORR}"
        "${PRESET_DIR}/${CONFIG_DFTCORR}"
        "${PRESET_DIR}/${CONFIG_POT}"
        "${PRESET_DIR}/${CONFIG_SPECTRUM}"
        "${PRESET_DIR}/movie.webm"
        "${PRESET_DIR}/${CONFIG_NUMEN}"
        "${PRESET_DIR}/spectrum.tex.inc"

        "${PRESET_DIR}/init.pdf"
        "${PRESET_DIR}/corr.pdf"
        "${PRESET_DIR}/spektrum.pdf")
  else()
    add_custom_target(${PRESET}
      DEPENDS
        ${PRESET_LUA}
        "${PRESET_DIR}/${CONFIG_APSI}"
        "${PRESET_DIR}/${CONFIG_CORR}"
        "${PRESET_DIR}/${CONFIG_DFTCORR}"
        "${PRESET_DIR}/${CONFIG_POT}"
        "${PRESET_DIR}/${CONFIG_SPECTRUM}"
        "${PRESET_DIR}/${CONFIG_NUMEN}"
        "${PRESET_DIR}/spectrum.tex.inc"

        "${PRESET_DIR}/init.pdf"
        "${PRESET_DIR}/corr.pdf"
        "${PRESET_DIR}/spektrum.pdf")
  endif()
  set(PLOTS ${PLOTS}
    "${PRESET_DIR}/init.pdf"
    "${PRESET_DIR}/corr.pdf"
    "${PRESET_DIR}/spektrum.pdf")
endforeach()

foreach(CM ${CMP})
  set(PRESET_LUA "${projekt_SOURCE_DIR}/presets/${CM}.lua")

  get_lua_var(${PRESET_LUA} "config.output.dir" CONFIG_DIR)

  set(PRESET_DIR "${projekt_SOURCE_DIR}/${CONFIG_DIR}")

  add_custom_command(
    OUTPUT
      "${PRESET_DIR}/result.dat"
    COMMAND
      mkdir -p "${PRESET_DIR}"
    COMMAND
      env
        "'LUA_PATH=${LUA_SCRIPTS}'"
        $<TARGET_FILE:simulate>
        ${PRESET_LUA}
        "${PRESET_DIR}/result.dat"
    MAIN_DEPENDENCY
      ${PRESET_LUA})

   add_custom_target(${CM}
    COMMAND
      env "'LUA_PATH=${LUA_SCRIPTS}'"
        $<TARGET_FILE:evaluate>
        ${PRESET_LUA}
        "${PRESET_DIR}/result.dat"
    DEPENDS "${PRESET_DIR}/result.dat")
endforeach()

foreach(BAD ${BADNESS})
  message("adding badness ${BAD}")

  get_lua_var(${PRESET_LUA} "config.output.dir" CONFIG_DIR)

  set(BAD_DIR "${projekt_SOURCE_DIR}/${CONFIG_DIR}")

  add_custom_command(
    OUTPUT
      "${BAD_DIR}/bspektrum.pdf"
    COMMAND
      ${GNUPLOT} "${SCRIPTS}/bspektrum.gp"
    WORKING_DIRECTORY
      "${BAD_DIR}/"
    DEPENDS
      "${SCRIPTS}/bspektrum.gp"
      "${BAD_DIR}/${CONFIG_DFTCORR}"
      "${BAD_DIR}/${CONFIG_SPECTRUM}"
      "${BAD_DIR}/${CONFIG_NUMEN}")
  set(PLOTS ${PLOTS} "${projekt_SOURCE_DIR}/${BAD}/bspektrum.pdf")
  set(BDEP
    ${BDEP} "${BAD_DIR}/bspektrum.pdf"
    "${BAD_DIR}/spectrum.tex.inc")
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
    ${PLOTS}
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
