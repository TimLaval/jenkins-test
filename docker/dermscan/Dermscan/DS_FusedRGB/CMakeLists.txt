cmake_minimum_required(VERSION 2.8)

cmake_policy(SET CMP0043 NEW)

set( APP_NAME DS_FusedRGB )

PROJECT(${APP_NAME})

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

option(DermScan_USE_VTK_WITH_QT "Turn this On in case you have VTK built with Qt support")


if(DermScan_USE_VTK_WITH_QT)
    set(CMAKE_AUTOMOC ON)
    find_package(Qt5Widgets REQUIRED QUIET)
    qt5_wrap_ui(UISrcs ${UI_FILES} )
endif()


include_directories( ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_BINARY_DIR}
    ${PROJECT_SOURCE_DIR}/../Source
)

set( ${APP_NAME}_SOURCES
    ${APP_NAME}.cpp
    ${PROJECT_SOURCE_DIR}/../Source/bdBIP.cpp
    ${PROJECT_SOURCE_DIR}/../Source/bdGeometry.cpp
    ${PROJECT_SOURCE_DIR}/../Source/bdGIP.cpp
    ${PROJECT_SOURCE_DIR}/../Source/bdImage.cpp
    ${PROJECT_SOURCE_DIR}/../Source/bdImageSliceTag.cpp
    ${PROJECT_SOURCE_DIR}/../Source/bdMelanoma.cpp
    ${PROJECT_SOURCE_DIR}/../Source/bdObject.cpp
    ${PROJECT_SOURCE_DIR}/../Source/bdRegion2D.cpp
    ${PROJECT_SOURCE_DIR}/../Source/bdRegion3D.cpp
    ${PROJECT_SOURCE_DIR}/../Source/bdRegularGrid.cpp
    ${PROJECT_SOURCE_DIR}/../Source/bdString.cpp
    ${PROJECT_SOURCE_DIR}/../Source/bdVoxel.cpp
    #${PROJECT_SOURCE_DIR}/../Source/vbdBasicImage.cpp
    ${PROJECT_SOURCE_DIR}/../Source/vbdMelanoma.cpp
)

set( ${APP_NAME}_HEADERS
    ${PROJECT_SOURCE_DIR}/../Source/bdArray.h
    ${PROJECT_SOURCE_DIR}/../Source/bdBIP.h
    ${PROJECT_SOURCE_DIR}/../Source/bdGeometry.h
    ${PROJECT_SOURCE_DIR}/../Source/bdGIP.h
    ${PROJECT_SOURCE_DIR}/../Source/bdImage.h
    ${PROJECT_SOURCE_DIR}/../Source/bdImageSliceTag.h
    ${PROJECT_SOURCE_DIR}/../Source/bdList.h
    ${PROJECT_SOURCE_DIR}/../Source/bdMatrix.h
    ${PROJECT_SOURCE_DIR}/../Source/bdMelanoma.h
    ${PROJECT_SOURCE_DIR}/../Source/bdObject.h
    ${PROJECT_SOURCE_DIR}/../Source/bdRegion2D.h
    ${PROJECT_SOURCE_DIR}/../Source/bdRegion3D.h
    ${PROJECT_SOURCE_DIR}/../Source/bdRegularGrid.h
    ${PROJECT_SOURCE_DIR}/../Source/bdStandard.h
    ${PROJECT_SOURCE_DIR}/../Source/bdString.h
    ${PROJECT_SOURCE_DIR}/../Source/bdVoxel.h
    #${PROJECT_SOURCE_DIR}/../Source/vbdBasicImage.h
    ${PROJECT_SOURCE_DIR}/../Source/vbdMelanoma.h
)

set( ${APP_NAME}_LIBS
)

#set( MACOSX_RESOURCE_FILES ${CMAKE_SOURCE_DIR}/Files/Icons/dermscan.png )
set( MACOSX_BUNDLE_INFO_STRING "@EXE_NAME@ - Version ${VERSION}" )
set( MACOSX_BUNDLE_BUNDLE_VERSION ${VERSION} )
#set( MACOSX_BUNDLE_ICON_FILE "dermscan.png" )
set( MACOSX_BUNDLE_GUI_IDENTIFIER "com.sri.aic" )
set( MACOSX_BUNDLE_BUNDLE_NAME "@EXE_NAME@" )

add_executable(${APP_NAME} MACOSX_BUNDLE
    ${${APP_NAME}_SOURCES}
    ${${APP_NAME}_HEADERS}
    #${MOC_SOURCE_FILES}
    #${UI_SOURCE_FILES}
    #${RESOURCE_SOURCE_FILES}
    #${HEADER_FILES}
    ${MACOSX_RESOURCE_FILES}
)

if(DermScan_USE_VTK_WITH_QT)
    qt5_use_modules(${APP_NAME} Core Gui)
endif()

target_link_libraries(${APP_NAME}
    ${VTK_LIBRARIES}
    #${APP_NAME}_LIBS
)

set_source_files_properties(
      ${MACOSX_RESOURCE_FILES}
      PROPERTIES
      MACOSX_PACKAGE_LOCATION Resources
)

install(TARGETS ${APP_NAME} DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)

install(CODE "
	include(BundleUtilities)
	fixup_bundle( \"${CMAKE_INSTALL_PREFIX}/bin/${APP_NAME}.app\" \"\" \"${CMAKE_INSTALL_PREFIX}/bin\" )
	" COMPONENT Runtime)
