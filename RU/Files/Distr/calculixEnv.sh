#!/bin/bash

#
# Change this path to your 
#
export CalculiXHome=/home/matukha/Документы/UniCourses/MeshingCourse/trunk/data/Files/Distr

export LaMathHome=$CalculiXHome/LaMath
export MPI_COMPILE_FLAGS="-I$MPI_ARCH_PATH/include"
export MPI_LINK_FLAGS="-L$MPI_ARCH_PATH/lib -lmpi"

export PATH="$CalculiXHome/CalculiXCrunchiX/:$CalculiXHome/CalculiXGraphiX/:$PATH"
