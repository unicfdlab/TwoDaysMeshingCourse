#!/bin/bash

source calculixEnv.sh

#
# Clean and build ARPACK
#

cd $LaMathHome/ARPACK
make clean
make all

#
# Clean and build SPOOLES
#

cd $LaMathHome/SPOOLES
make clean
make global

#
# Clean and build CalculiX CrunchiX
#
cd $CalculiXHome/CalculiXCrunchiX/src
rm -rf *.o
make


#
# Clean and build CalculiX GraphiX
#
cd $CalculiXHome/CalculiXGraphiX/cgx_2.6/src
rm -rf *.o
make

cd $CalculixHome

#
#END-OF-FILE
#



