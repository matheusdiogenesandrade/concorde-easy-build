# CPLEX_ROOT_DIR = IBM CPLEX Studio install root (directory containing cplex/ and concert/).
# If CPLEX_ROOT_DIR is set in the environment, it overrides the CMake cache (fixes stale /path/to/... entries).
# Otherwise: cache / -D, then other env vars, then derive from CPLEX_DIR, then common Linux paths.

if(DEFINED ENV{CPLEX_ROOT_DIR} AND NOT "$ENV{CPLEX_ROOT_DIR}" STREQUAL "")
    set(CPLEX_ROOT_DIR "$ENV{CPLEX_ROOT_DIR}" CACHE PATH "IBM CPLEX Studio root (contains cplex/ and concert/)" FORCE)
endif()

if(NOT CPLEX_ROOT_DIR)
    if(DEFINED ENV{CPLEX_STUDIO_DIR} AND NOT "$ENV{CPLEX_STUDIO_DIR}" STREQUAL "")
        set(CPLEX_ROOT_DIR "$ENV{CPLEX_STUDIO_DIR}")
    elseif(DEFINED ENV{CPLEX_DIR} AND NOT "$ENV{CPLEX_DIR}" STREQUAL "")
        get_filename_component(CPLEX_ROOT_DIR "$ENV{CPLEX_DIR}" DIRECTORY)
    endif()
endif()
if(NOT CPLEX_ROOT_DIR AND DEFINED CPLEX_DIR AND NOT "${CPLEX_DIR}" STREQUAL "")
    get_filename_component(CPLEX_ROOT_DIR "${CPLEX_DIR}" DIRECTORY)
endif()
if(NOT CPLEX_ROOT_DIR AND UNIX AND NOT APPLE)
    foreach(_ilog /opt/ibm/ILOG /opt/IBM/ILOG)
        file(GLOB _studios "${_ilog}/CPLEX_Studio*")
        list(LENGTH _studios _n)
        if(_n GREATER 0)
            list(SORT _studios)
            list(REVERSE _studios)
            list(GET _studios 0 CPLEX_ROOT_DIR)
            break()
        endif()
    endforeach()
endif()

if(CPLEX_ROOT_DIR)
    set(CPLEX_ROOT_DIR "${CPLEX_ROOT_DIR}" CACHE PATH "IBM CPLEX Studio root (contains cplex/ and concert/)")
endif()

if(NOT CPLEX_ROOT_DIR)
    message(FATAL_ERROR
        "FindCPLEX (concorde-easy-build): CPLEX_ROOT_DIR is not set and could not be inferred.\n"
        "Set your CPLEX Studio install root (the folder that contains cplex/ and concert/), for example:\n"
        "  export CPLEX_ROOT_DIR=/opt/ibm/ILOG/CPLEX_Studio2211\n"
        "  cmake -B build -S . ...\n"
        "or pass -DCPLEX_ROOT_DIR=/opt/ibm/ILOG/CPLEX_Studio2211\n"
        "If you only have CPLEX_DIR pointing at .../cplex, set CPLEX_ROOT_DIR to its parent directory.")
endif()

get_filename_component(_cplex_studio "${CPLEX_ROOT_DIR}" ABSOLUTE)

string(TOLOWER "${_cplex_studio}" _cplex_lc)
if(_cplex_lc MATCHES "/path/to")
    message(FATAL_ERROR
        "CPLEX_ROOT_DIR='${_cplex_studio}' looks like a documentation placeholder (e.g. /path/to/...).\n"
        "Replace it with your real IBM CPLEX Studio install directory (the one that contains cplex/ and concert/).")
endif()

if(NOT EXISTS "${_cplex_studio}")
    message(FATAL_ERROR
        "CPLEX_ROOT_DIR='${_cplex_studio}' does not exist on disk.\n"
        "Fix the path, export it in the same shell as cmake, or pass -DCPLEX_ROOT_DIR=...")
endif()

if(NOT EXISTS "${_cplex_studio}/cplex" OR NOT EXISTS "${_cplex_studio}/concert")
    message(FATAL_ERROR
        "CPLEX_ROOT_DIR='${_cplex_studio}' is not a CPLEX Studio root: expected subdirectories 'cplex/' and 'concert/'.\n"
        "Often the mistake is pointing at .../cplex — use its parent directory instead (the CPLEX_Studio... folder).")
endif()

find_path(CPLEX_INCLUDE_DIR
    NAMES ilcplex/cplex.h
    PATHS ${_cplex_studio}
    PATH_SUFFIXES include cplex/include
)

find_path(CONCERT_INCLUDE_DIR
    NAMES ilconcert/iloenv.h
    PATHS ${_cplex_studio}
    PATH_SUFFIXES include concert/include
)

file(GLOB CPLEX_LIB_PATHS "${_cplex_studio}/cplex/lib/*/static_pic")

find_library(CPLEX_LIBRARY
    NAMES cplex
    PATHS ${CPLEX_LIB_PATHS}
)

find_library(ILOCPLEX_LIBRARY
    NAMES ilocplex
    PATHS ${CPLEX_LIB_PATHS}
)

file(GLOB CONCERT_LIB_PATHS "${_cplex_studio}/concert/lib/*/static_pic")

find_library(CONCERT_LIBRARY
    NAMES concert
    PATHS ${CONCERT_LIB_PATHS}
)

mark_as_advanced(CPLEX_FOUND CPLEX_INCLUDE_DIR CONCERT_INCLUDE_DIR CPLEX_LIBRARY ILOCPLEX_LIBRARY CONCERT_LIBRARY)

if(NOT CPLEX_INCLUDE_DIR
    OR NOT CONCERT_INCLUDE_DIR
    OR NOT CPLEX_LIBRARY
    OR NOT ILOCPLEX_LIBRARY
    OR NOT CONCERT_LIBRARY)
    file(GLOB _cw_cplex_pic_dirs "${_cplex_studio}/cplex/lib/*/static_pic")
    file(GLOB _cw_concert_pic_dirs "${_cplex_studio}/concert/lib/*/static_pic")
    message(FATAL_ERROR
        "FindCPLEX (concorde-easy-build): could not find CPLEX headers or *static_pic* libraries.\n"
        "  CPLEX_ROOT_DIR (absolute): ${_cplex_studio}\n"
        "  CPLEX_INCLUDE_DIR: ${CPLEX_INCLUDE_DIR}\n"
        "  CONCERT_INCLUDE_DIR: ${CONCERT_INCLUDE_DIR}\n"
        "  CPLEX_LIBRARY: ${CPLEX_LIBRARY}\n"
        "  ILOCPLEX_LIBRARY: ${ILOCPLEX_LIBRARY}\n"
        "  CONCERT_LIBRARY: ${CONCERT_LIBRARY}\n"
        "  Globs: cplex .../static_pic -> ${_cw_cplex_pic_dirs}\n"
        "         concert .../static_pic -> ${_cw_concert_pic_dirs}\n"
        "Concorde links the static PIC archives. Install CPLEX with static libraries, or ensure the "
        "architecture folder under cplex/lib/ matches this machine (e.g. x86-64_linux/static_pic).")
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CPLEX
    REQUIRED_VARS CPLEX_INCLUDE_DIR CONCERT_INCLUDE_DIR CPLEX_LIBRARY ILOCPLEX_LIBRARY CONCERT_LIBRARY
)

if(CPLEX_FOUND)
    set(CPLEX_INCLUDE_DIRS ${CPLEX_INCLUDE_DIR} ${CONCERT_INCLUDE_DIR})
    set(CPLEX_LIBRARIES ${CONCERT_LIBRARY} ${ILOCPLEX_LIBRARY} ${CPLEX_LIBRARY})
endif()

if(CPLEX_FOUND AND NOT TARGET "Cplex::cplex")
    # GLOBAL so the parent concorde_wrapper directory (sibling of nested project(concorde)) can link it.
    add_library(Cplex::cplex INTERFACE IMPORTED GLOBAL)
    set_target_properties(Cplex::cplex PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${CPLEX_INCLUDE_DIRS}"
        INTERFACE_LINK_LIBRARIES "${CPLEX_LIBRARIES}"
    )
endif()
