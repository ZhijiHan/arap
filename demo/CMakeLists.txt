cmake_minimum_required(VERSION 2.6)
project(demo)

include("../CMakeLists.shared")

set(SOURCES
${PROJECT_SOURCE_DIR}/main.cc
${PROJECT_SOURCE_DIR}/adaptadmmfixedsolver.cc
${PROJECT_SOURCE_DIR}/adaptadmmfreesolver.cc
${PROJECT_SOURCE_DIR}/admmfixedsolver.cc
${PROJECT_SOURCE_DIR}/admmfreesolver.cc
${PROJECT_SOURCE_DIR}/arapsolver.cc
${PROJECT_SOURCE_DIR}/arapbenchmarksolver.cc
${PROJECT_SOURCE_DIR}/energy.cc
${PROJECT_SOURCE_DIR}/solver.cc
)

add_executable(${PROJECT_NAME}_bin ${SOURCES} ${SHARED_SOURCES})
target_link_libraries(${PROJECT_NAME}_bin ${SHARED_LIBRARIES})
