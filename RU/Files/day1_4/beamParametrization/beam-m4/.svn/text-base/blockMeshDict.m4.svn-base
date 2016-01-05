dnl>
dnl> Initial parameters
dnl>
define(toMeters, 1.0)dnl>
define(cx, 0.0)dnl>
define(cy, 0.5)dnl>
define(cz, 0.5)dnl>
define(l, 1.0)dnl>
define(h,10.0)dnl>
define(xdivs, 10)dnl>
define(ydivs, 2)dnl>
define(zdivs, 2)dnl>
dnl>
dnl> Functions definitions
dnl>
changecom(//) dnl>
changequote([,]) dnl>
define(calc, [esyscmd(perl -e 'print ($1)')]) dnl>
define(VCOUNT, 0)  dnl>
define(vlabel, [[// ]point VCOUNT ($1) define($1, VCOUNT)define([VCOUNT], incr(VCOUNT))])  dnl>
define(vert,
    [($1 $2 $3)] dnl>
    [vlabel($4)] dnl>
    ) dnl> 
define(hex2D, hex (b$1 b$2 b$3 b$4 f$1 f$2 f$3 f$4)) dnl>
define(quad, ($1 $2 $3 $4))dnl>
define(defpatch,
    $1
    [{]
	[type]	[patch];
	faces
	(
	    $2
	);
    [}]
    ) dnl>
dnl>
dnl> Some characteristic lengths
dnl>
define(minx,cx)dnl>
define(maxx,calc(minx+h))dnl>
define(miny,calc(cy-l*0.5))dnl>
define(maxy,calc(cy+l*0.5))dnl>
define(minz,calc(cz-l*0.5))dnl>
define(maxz,calc(cz+l*0.5))dnl>


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
    [format]      ascii;
    class       dictionary;
    object      blockMeshDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

convertToMeters toMeters;

vertices
(
    vert(minx, miny, minz, bv1)
    vert(maxx, miny, minz, bv2)
    vert(maxx, maxy, minz, bv3)
    vert(minx, maxy, minz, bv4)
    
    vert(minx, miny, maxz, fv1)
    vert(maxx, miny, maxz, fv2)
    vert(maxx, maxy, maxz, fv3)
    vert(minx, maxy, maxz, fv4)
);

blocks
(
    hex2D(v1, v2, v3, v4) ( xdivs ydivs zdivs ) simpleGrading (1 1 1)
);

edges
(
);

boundary
(
    defpatch(zfaces, [quad(bv1,bv2,bv3,bv4) quad(fv1,fv2,fv3,fv4)])
    defpatch(yfaces, [quad(bv1,bv2,fv2,fv1) quad(bv3,bv4,fv4,fv3)])
    defpatch(fixed,[quad(bv1,bv4,fv4,fv1)])
    defpatch(traction,[quad(bv2,bv3,fv3,fv2)])
);

mergePatchPairs
(
);

// ************************************************************************* //
