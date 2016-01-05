##!/usr/bin/python
# -*- coding: utf-8 -*-

# author: Michael Kalugin <shaman@ispras.ru>
# date: 28.09.2012
# This file is licensed under GPL

import re
import sys
from abc import ABCMeta, abstractmethod


def map_point(line):
	values = line.split()
	
	if len(values) == 4:
		return "(%s %s %s)\n" % (values[1], values[2], values[3])
	else:
		print "unknown line format: %s \n" % line

def map_edge(line):
	values = line.split()
	if len(values) == 4:
		if values[1] == "102":
			try:
				int(values[2])
				int(values[3])
				return "(%s %s)\n" % (int(values[2]) -1, int(values[3])-1)
			except Exception as e:
				print str(e) + " at line: " + line
		else:
			print "wrong line, should contain only integers: %s \n" % line
	else:
		print "unknown line format: %s \n" % line	

if __name__ == "__main__":
        try:
                if len(sys.argv) != 3:
                        print "script @1 @2 \n@1 input file in dat \n@2 ouput file emesh"
                        sys.exit(1)

                lines = open(sys.argv[1]).read().split("\n")
		
		quantity = lines[0].split()
		NPoints = int(quantity[0])
		NEdges = int(quantity[1])

		result_points = ""

		for line in lines[1:NPoints + 1]:
			line = line.strip()
			if len(line) > 0:
				result_points += map_point(line)

		result_edges = ""
		for line in lines[NPoints + 1:]:
			line = line.strip()
			if len(line) > 0:
				result_edges += map_edge(line) 

                output = open(sys.argv[2],"w+")
                #
                # OpenFOAM header
                #

		output.write("/*--------------------------------*- C++ -*----------------------------------*\\\n")
		output.write("| =========                 |                                                 |\n")
		output.write("| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n")
		output.write("|  \\\\    /   O peration     | Version:  2.1.0                                 |\n")
		output.write("|   \\\\  /    A nd           | Web:      www.OpenFOAM.org                      |\n")
		output.write("|    \\\\/     M anipulation  |                                                 |\n")
		output.write("\\*---------------------------------------------------------------------------*/\n")
		output.write("FoamFile\n")
		output.write("{\n")
		output.write("    version     2.0;\n")
		output.write("    format      ascii;\n")
		output.write("    class       featureEdgeMesh;\n")
		output.write("    location    \"constant/triSurface\";\n")
		output.write("    object      %s;\n" % sys.argv[2])
		output.write("}\n")
		output.write("// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n\n")


		output.write("// points:\n\n")

		output.write(quantity[0] + "\n(\n")
                output.write(result_points)
		output.write(")\n\n\n// edges:\n\n")

		output.write(quantity[1] + "\n(\n")
		output.write(result_edges)
		output.write(")\n\n")

		output.write("// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n\n")


		output.write("// ************************************************************************* //\n");

                output.close()

        except Exception as e:
                print str(e)


