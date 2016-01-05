#!/bin/bash

m4 beam-linear.m4 > beam-linear.fbd
cgx_2.6.1 -b beam-linear.fbd 

#
#END-OF-FILE
#

