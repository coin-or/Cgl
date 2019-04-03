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
  CglAllDifferent/CglAllDifferent.cpp
  CglBKClique/CglBKClique.cpp
  CglClique/CglClique.cpp
  CglClique/CglCliqueHelper.cpp
  CglClique/CglCliqueTest.cpp
  CglCliqueMerging/CglCliqueMerging.cpp
  CglCommon/CglCutGenerator.cpp
  CglCommon/CglMessage.cpp
  CglCommon/CglParam.cpp
  CglCommon/CglStored.cpp
  CglCommon/CglTreeInfo.cpp
  CglDuplicateRow/CglDuplicateRow.cpp
  CglFlowCover/CglFlowCover.cpp
  CglFlowCover/CglFlowCoverTest.cpp
  CglGMI/CglGMI.cpp
  CglGMI/CglGMIParam.cpp
  CglGomory/CglGomory.cpp
  CglGomory/CglGomoryTest.cpp
  CglKnapsackCover/CglKnapsackCover.cpp
  CglKnapsackCover/CglKnapsackCoverTest.cpp
  CglLandP/CglLandP.cpp
  CglLandP/CglLandPMessages.cpp
  CglLandP/CglLandPSimplex.cpp
  CglLandP/CglLandPTabRow.cpp
  CglLandP/CglLandPTest.cpp
  CglLandP/CglLandPUtils.cpp
  CglLandP/CglLandPValidator.cpp
  CglLiftAndProject/CglLiftAndProject.cpp
  CglMixedIntegerRounding/CglMixedIntegerRounding.cpp
  CglMixedIntegerRounding/CglMixedIntegerRoundingTest.cpp
  CglMixedIntegerRounding2/CglMixedIntegerRounding2.cpp
  CglMixedIntegerRounding2/CglMixedIntegerRounding2Test.cpp
  CglOddHole/CglOddHole.cpp
  CglOddHole/CglOddHoleTest.cpp
  CglOddHoleWC/CglOddHoleWC.cpp
  CglOddHoleWC/oddhs.cpp
  CglOddHoleWC/spaths.cpp
  CglPreProcess/CglPreProcess.cpp
  CglProbing/CglProbing.cpp
  CglProbing/CglProbingTest.cpp
  CglRedSplit/CglRedSplit.cpp
  CglRedSplit/CglRedSplitParam.cpp
  CglRedSplit/CglRedSplitTest.cpp
  CglRedSplit2/CglRedSplit2.cpp
  CglRedSplit2/CglRedSplit2Param.cpp
  CglRedSplit2/CglRedSplit2Test.cpp
  CglResidualCapacity/CglResidualCapacity.cpp
  CglResidualCapacity/CglResidualCapacityTest.cpp
  CglSimpleRounding/CglSimpleRounding.cpp
  CglSimpleRounding/CglSimpleRoundingTest.cpp
  CglTwomir/CglTwomir.cpp
  CglTwomir/CglTwomirTest.cpp
  CglZeroHalf/Cgl012cut.cpp
  CglZeroHalf/CglZeroHalf.cpp
  CglZeroHalf/CglZeroHalfTest.cpp
  )

set(_HDRS
  CglAllDifferent/CglAllDifferent.hpp
  CglBKClique/CglBKClique.hpp
  CglClique/CglClique.hpp
  CglCliqueMerging/CglCliqueMerging.hpp
  CglCommon/CglConfig.h
  CglCommon/CglCutGenerator.hpp
  CglCommon/CglMessage.hpp
  CglCommon/CglParam.hpp
  CglCommon/CglStored.hpp
  CglCommon/CglTreeInfo.hpp
  CglDuplicateRow/CglDuplicateRow.hpp
  CglFlowCover/CglFlowCover.hpp
  CglGMI/CglGMI.hpp
  CglGMI/CglGMIParam.hpp
  CglGomory/CglGomory.hpp
  CglKnapsackCover/CglKnapsackCover.hpp
  CglLandP/CglLandP.hpp
  CglLandP/CglLandPMessages.hpp
  CglLandP/CglLandPSimplex.hpp
  CglLandP/CglLandPTabRow.hpp
  CglLandP/CglLandPUtils.hpp
  CglLandP/CglLandPValidator.hpp
  CglLiftAndProject/CglLiftAndProject.hpp
  CglMixedIntegerRounding/CglMixedIntegerRounding.hpp
  CglMixedIntegerRounding2/CglMixedIntegerRounding2.hpp
  CglOddHole/CglOddHole.hpp
  CglOddHoleWC/CglOddHoleWC.hpp
  CglOddHoleWC/oddhs.h
  CglOddHoleWC/spaths.h
  CglPreProcess/CglPreProcess.hpp
  CglProbing/CglProbing.hpp
  CglRedSplit/CglRedSplit.hpp
  CglRedSplit/CglRedSplitParam.hpp
  CglRedSplit2/CglRedSplit2.hpp
  CglRedSplit2/CglRedSplit2Param.hpp
  CglResidualCapacity/CglResidualCapacity.hpp
  CglSimpleRounding/CglSimpleRounding.hpp
  CglTwomir/CglTwomir.hpp
  CglZeroHalf/Cgl012cut.hpp
  CglZeroHalf/CglZeroHalf.hpp
  )

add_library(Cgl ${_SRCS} ${_HDRS})
target_include_directories(Cgl PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/CglCommon>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/CglClique>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/CglDuplicateRow>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/CglFlowCover>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/CglGMI>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/CglGomory>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/CglKnapsackCover>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/CglLandP>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/CglMixedIntegerRounding>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/CglMixedIntegerRounding2>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/CglPreProcess>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/CglProbing>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/CglOddHole>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/CglRedSplit>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/CglRedSplit2>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/CglResidualCapacity>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/CglSimpleRounding>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/CglTwomir>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/CglZeroHalf>
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
include(GNUInstallDirs)
install(TARGETS Cgl
  EXPORT ${PROJECT_NAME}Targets
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/coin
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  )
