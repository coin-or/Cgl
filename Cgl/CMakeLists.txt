set(NAME "CGL")

# PTHREAD
if(PTHREADS_FOUND)
  set(${NAME}_PTHREADS "1" CACHE INTERNAL "Use pthread")
endif()

set(COIN_${NAME}_CHECKLEVEL "0" CACHE INTERNAL
  "${NAME} check level")
set(COIN_${NAME}_VERBOSITY "0" CACHE INTERNAL
  "${NAME} verbosity level")
configure_file(config.h.cmake.in config.h)
configure_file(config_cgl.h.cmake.in config_cgl.h)

set(_SRCS
  src/CglMixedIntegerRounding/CglMixedIntegerRoundingTest.cpp
  src/CglMixedIntegerRounding/CglMixedIntegerRounding.cpp
  src/CglDuplicateRow/CglDuplicateRow.cpp
  src/CglCutGenerator.cpp
  src/CglPreProcess/CglPreProcess.cpp
  src/CglProbing/CglProbing.cpp
  src/CglProbing/CglProbingTest.cpp
  src/CglMixedIntegerRounding2/CglMixedIntegerRounding2.cpp
  src/CglMixedIntegerRounding2/CglMixedIntegerRounding2Test.cpp
  src/CglLandP/CglLandPUtils.cpp
  src/CglLandP/CglLandPTest.cpp
  src/CglLandP/CglLandPTabRow.cpp
  src/CglLandP/CglLandPSimplex.cpp
  src/CglLandP/CglLandPValidator.cpp
  src/CglLandP/CglLandPMessages.cpp
  src/CglLandP/CglLandP.cpp
  src/CglOddHole/CglOddHoleTest.cpp
  src/CglOddHole/CglOddHole.cpp
  src/CglRedSplit/CglRedSplitParam.cpp
  src/CglRedSplit/CglRedSplit.cpp
  src/CglRedSplit/CglRedSplitTest.cpp
  src/CglAllDifferent/CglAllDifferent.cpp
  src/CglKnapsackCover/CglKnapsackCoverTest.cpp
  src/CglKnapsackCover/CglKnapsackCover.cpp
  src/CglSimpleRounding/CglSimpleRoundingTest.cpp
  src/CglSimpleRounding/CglSimpleRounding.cpp
  src/CglZeroHalf/CglZeroHalfTest.cpp
  src/CglZeroHalf/Cgl012cut.cpp
  src/CglZeroHalf/CglZeroHalf.cpp
  src/CglRedSplit2/CglRedSplit2Test.cpp
  src/CglRedSplit2/CglRedSplit2Param.cpp
  src/CglRedSplit2/CglRedSplit2.cpp
  src/CglMessage.cpp
  src/CglStored.cpp
  src/CglParam.cpp
  src/CglResidualCapacity/CglResidualCapacity.cpp
  src/CglResidualCapacity/CglResidualCapacityTest.cpp
  src/CglTwomir/CglTwomirTest.cpp
  src/CglTwomir/CglTwomir.cpp
  src/CglFlowCover/CglFlowCoverTest.cpp
  src/CglFlowCover/CglFlowCover.cpp
  src/CglClique/CglCliqueHelper.cpp
  src/CglClique/CglCliqueTest.cpp
  src/CglClique/CglClique.cpp
  src/CglTreeInfo.cpp
  src/CglLiftAndProject/CglLiftAndProject.cpp
  src/CglGomory/CglGomoryTest.cpp
  src/CglGomory/CglGomory.cpp
  src/CglGMI/CglGMI.cpp
  src/CglGMI/CglGMIParam.cpp)

set(_HDRS
  src/CglConfig.h
  src/CglMixedIntegerRounding/CglMixedIntegerRounding.hpp
  src/CglDuplicateRow/CglDuplicateRow.hpp
  src/CglStored.hpp
  src/CglPreProcess/CglPreProcess.hpp
  src/CglProbing/CglProbing.hpp
  src/CglMixedIntegerRounding2/CglMixedIntegerRounding2.hpp
  src/CglLandP/CglLandP.hpp
  src/CglLandP/CglLandPUtils.hpp
  src/CglLandP/CglLandPValidator.hpp
  src/CglLandP/CglLandPTabRow.hpp
  src/CglLandP/CglLandPMessages.hpp
  src/CglLandP/CglLandPSimplex.hpp
  src/CglOddHole/CglOddHole.hpp
  src/CglRedSplit/CglRedSplitParam.hpp
  src/CglRedSplit/CglRedSplit.hpp
  src/CglAllDifferent/CglAllDifferent.hpp
  src/CglKnapsackCover/CglKnapsackCover.hpp
  src/CglSimpleRounding/CglSimpleRounding.hpp
  src/CglMessage.hpp
  src/CglZeroHalf/CglZeroHalf.hpp
  src/CglZeroHalf/Cgl012cut.hpp
  src/CglRedSplit2/CglRedSplit2Param.hpp
  src/CglRedSplit2/CglRedSplit2.hpp
  src/CglCutGenerator.hpp
  src/CglResidualCapacity/CglResidualCapacity.hpp
  src/CglTwomir/CglTwomir.hpp
  src/CglFlowCover/CglFlowCover.hpp
  src/CglClique/CglClique.hpp
  src/CglLiftAndProject/CglLiftAndProject.hpp
  src/CglGomory/CglGomory.hpp
  src/CglGMI/CglGMIParam.hpp
  src/CglGMI/CglGMI.hpp
  src/CglParam.hpp
  src/CglTreeInfo.hpp)

include(GNUInstallDirs)

add_library(Cgl ${_SRCS} ${_HDRS})
target_include_directories(Cgl PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/CglClique>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/CglDuplicateRow>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/CglFlowCover>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/CglGMI>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/CglGomory>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/CglKnapsackCover>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/CglLandP>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/CglMixedIntegerRounding>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/CglMixedIntegerRounding2>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/CglPreProcess>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/CglProbing>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/CglOddHole>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/CglRedSplit>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/CglRedSplit2>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/CglResidualCapacity>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/CglSimpleRounding>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/CglTwomir>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/CglZeroHalf>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
  $<INSTALL_INTERFACE:include/coin>)
target_compile_definitions(Cgl
  PUBLIC HAVE_CONFIG_H
  PRIVATE CGL_BUILD)
if(CMAKE_VERSION VERSION_LESS "3.8.2")
  set_property(TARGET Cgl PROPERTY CXX_STANDARD 11)
  set_property(TARGET Cgl PROPERTY CXX_STANDARD_REQUIRED ON)
else()
  target_compile_features(Cgl PUBLIC cxx_std_11)
endif()
if(APPLE)
  set_target_properties(Cgl PROPERTIES INSTALL_RPATH "@loader_path")
elseif(UNIX)
  set_target_properties(Cgl PROPERTIES INSTALL_RPATH "$ORIGIN")
endif()
target_link_libraries(Cgl PUBLIC
  Coin::CoinUtils
  Coin::Osi
  Coin::OsiClp)
set_target_properties(Cgl PROPERTIES
  PUBLIC_HEADER "${_HDRS};${CMAKE_CURRENT_BINARY_DIR}/config_cgl.h"
  VERSION ${PROJECT_VERSION}
  SOVERSION ${PROJECT_VERSION_MAJOR})
add_library(Coin::Cgl ALIAS Cgl)

# Install
install(TARGETS Cgl
  EXPORT ${PROJECT_NAME}Targets
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/coin
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  )
