#
# Try to find GLFW library and include path.
# Once done this will define
#
# GLFW_FOUND
# GLFW_INCLUDE_DIR
# GLFW_LIBRARIES
#

FIND_PATH(GLFW_INCLUDE_DIR GLFW/glfw3.h
  PATHS
    /home/taodu/research/libigl/external/glfw/include
    /usr/local/include
    /usr/X11/include
    /usr/include
    /opt/local/include
    NO_DEFAULT_PATH
    )

FIND_LIBRARY( GLFW_LIBRARIES NAMES glfw3
  PATHS
    /home/taodu/research/libigl/external/glfw/src
    /home/taodu/research/libigl/external/glfw/lib/x64
    /usr/lib
    /usr/local
    /usr/X11
    /usr
    PATH_SUFFIXES
    a
    lib64
    lib
    NO_DEFAULT_PATH
)

SET(GLFW_FOUND "NO")
IF (GLFW_INCLUDE_DIR AND GLFW_LIBRARIES)
	SET(GLFW_FOUND "YES")
ENDIF (GLFW_INCLUDE_DIR AND GLFW_LIBRARIES)

if(GLFW_FOUND)
  message(STATUS "Found GLFW: ${GLFW_INCLUDE_DIR}")
else(GLFW_FOUND)
  message(FATAL_ERROR "could NOT find GLFW")
endif(GLFW_FOUND)
