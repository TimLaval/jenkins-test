
cmake_minimum_required(VERSION 2.8)

if(POLICY CMP0042)
  cmake_policy(SET CMP0042 NEW)
endif()

PROJECT(DermScan)


include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
include_directories( Source )

add_subdirectory( DS_ABCD )
#add_subdirectory( DS_FocusAnalysis )
#add_subdirectory( DS_FusedRGB )
#add_subdirectory( DS_LesionSegmentation )
