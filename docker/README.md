### Docker Image for the ABCD
Contains VTK, GoogleTest, lcov, cppcheck, doxygen and python3
Is able to launch ABCD's code for unit testing, coverage and technical documentation, python tool for the segmentation annotation.


Contains .sh files for each execution.
    build_abcd.sh to use cmake on ABCD code (first thing to do)
    testing_and_coverage.sh to start unit testing and see the coverage on it
    analysis.sh to start status code analysis with cppcheck
    segmentation_tool.sh to start python application on segmentation
    documentation.sh to start technical documentation with Doxygen.