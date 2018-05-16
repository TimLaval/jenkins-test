In order to build the code for your platform download CMake:

https://cmake.org/download/

and download VTK:

http://www.vtk.org/download/

Build the VTK using CMake and your IDE, some guidelines are here:

http://www.vtk.org/Wiki/VTK/Configure_and_Build

After you have built VTK, do the same with this project (project folder DermScam_CMAKE contains the CMakeLists.txt file). Subprojects will build executables which you can use from command line and hence write batch commands to process multiple images. The project "install" will install the built app to the folder you specified in the CMake configuration under CMAKE_INSTALL_PREFIX.