#!/bin/bash

rm -rf abcd_build
mkdir abcd_build
cd abcd_build
cmake ../gtest-jenkins
make test
