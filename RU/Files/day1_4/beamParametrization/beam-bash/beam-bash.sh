#!/bin/bash

#
# initial parameters
#
toMeters="1.0;"
cx=0.0
cy=0.5
cz=0.5
l=1.0
h=10.0
xdivs=10
ydivs=2
zdivs=2
#
# block vertices
#
minx=$cx
maxx=`echo "$cx+$h" | bc`
miny=`echo "$cy-$l*0.5" | bc`
maxy=`echo "$cy+$l*0.5" | bc`
minz=`echo "$cz-$l*0.5" | bc`
maxz=`echo "$cz+$l*0.5" | bc`

#
# block 1
#
export j=0
export n=8
export block1=""
while [ `echo "$j<$n" | bc` -eq 1 ]
do
    block1="$block1 $j"
    j=`echo "$j+1" | bc`
done
#
# write blockMesh header
#
cat << BLOCKMESHHEADER
/*--------------------------------*- C++ -*----------------------------------*\ 
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  2.1.1                                 |
|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
 FoamFile
 {
    version     2.0;
    format      ascii;
    class       dictionary;
    object      blockMeshDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

convertToMeters $toMeters

vertices
(
    ($minx $miny $minz)
    ($maxx $miny $minz)
    ($maxx $maxy $minz)
    ($minx $maxy $minz)

    ($minx $miny $maxz)
    ($maxx $miny $maxz)
    ($maxx $maxy $maxz)
    ($minx $maxy $maxz)
);

blocks
(
    hex ($block1) ($xdivs $ydivs $zdivs) simpleGrading (1 1 1)
);

edges
(
);

boundary
(
    zfaces
    {
	type	patch;
	faces
	(
	    (0 1 2 3)
	    (4 5 6 7)
	);
    }
    
    yfaces
    {
	type	patch;
	faces
	(
	    (0 1 5 4)
	    (2 3 7 6)
	);
    }
    
    traction
    {
	type	patch;
	faces
	(
	    (1 2 6 5)
	);
    }
    
    fixed
    {
	type	patch;
	faces
	(
	    (0 3 7 4)
	);
    }
);

mergePatchPairs
(
);

BLOCKMESHHEADER



echo "// ************************************************************************* //"


#
#END-OF-FILE
#

