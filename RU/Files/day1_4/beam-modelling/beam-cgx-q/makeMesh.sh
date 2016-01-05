#!/bin/bash

m4 beam-quad.m4 > beam-quad.fbd
cgx_2.6.1 -b beam-quad.fbd 

#
#END-OF-FILE
#

