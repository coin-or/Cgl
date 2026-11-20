
add_executable(CglTest "unitTest.cpp")
target_include_directories(CglTest PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
  $<INSTALL_INTERFACE:include/coin>)
target_compile_definitions(CglTest
  PUBLIC HAVE_CONFIG_H
  PRIVATE CGL_BUILD COIN_HAS_CGL)
if(CMAKE_VERSION VERSION_LESS "3.8.2")
  set_property(TARGET CglTest PROPERTY CXX_STANDARD 11)
  set_property(TARGET CglTest PROPERTY CXX_STANDARD_REQUIRED ON)
else()
  target_compile_features(CglTest PUBLIC cxx_std_11)
endif()
if(APPLE)
  set_target_properties(CglTest PROPERTIES
    INSTALL_RPATH "@loader_path/../${CMAKE_INSTALL_LIBDIR};@loader_path")
elseif(UNIX)
  cmake_path(RELATIVE_PATH CMAKE_INSTALL_FULL_LIBDIR
    BASE_DIRECTORY ${CMAKE_INSTALL_FULL_BINDIR}
    OUTPUT_VARIABLE libdir_relative_path)
  set_target_properties(CglTest PROPERTIES
    INSTALL_RPATH "$ORIGIN/${libdir_relative_path}:$ORIGIN")
endif()
target_link_libraries(CglTest PRIVATE
  Coin::CoinUtils
  Coin::Osi
  Coin::Cgl)
set_target_properties(CglTest PROPERTIES
  VERSION ${PROJECT_VERSION}
  SOVERSION ${PROJECT_VERSION_MAJOR})
add_executable(Coin::CglTest ALIAS CglTest)

add_test(NAME CglTest COMMAND CglTest)
