project (abcd)

find_package(Doxygen)

if(DOXYGEN_FOUND)
    add_subdirectory(doxygen)
endif()
