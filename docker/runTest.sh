#!/bin/bash

rm -rf abcd_build
mkdir abcd_build
cd abcd_build
cmake /work/abcd/gtest-jenkins
make
