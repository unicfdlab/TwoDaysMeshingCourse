dnl>
dnl> Initial parameters
dnl>
define(toMeters, 1.0)dnl>
define(l,2.0)dnl>
define(r,0.5)dnl>
define(zmax,0.5)dnl>
define(sdivs, 10)dnl>
define(bdivs, 20)dnl>
define(zdivs, 1)dnl>
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
define(hex3D, hex ($1 $2 $3 $4 $5 $6 $7 $8)) dnl>
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
define(x1,calc(r*sin(3.14159265359/4)))dnl>
define(x2,calc(2*r*sin(3.14159265359/4)))dnl>
define(x3,calc(2*r))dnl>
define(y1,calc(r*sin(3.14159265359/4)))dnl>
define(y2,calc(2*r*cos(3.14159265359/4)))dnl>
define(y3,calc(2*r))dnl>


define(createArc, arc [$1 $2] )dnl>
define(arcx1,calc(r*cos(3.14159265359/8)))dnl>
define(arcy1,calc(r*sin(3.14159265359/8)))dnl>
define(arcx2,calc(r*sin(3.14159265359/8)))dnl>
define(arcy2,calc(r*cos(3.14159265359/8)))dnl>
define(arcx3,calc(2*r*cos(3.14159265359/8)))dnl>
define(arcy3,calc(2*r*sin(3.14159265359/8)))dnl>
define(arcx4,calc(2*r*sin(3.14159265359/8)))dnl>
define(arcy4,calc(2*r*cos(3.14159265359/8)))dnl>


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
    vert(r, 0, 0, bv1)
    vert(x3, 0, 0, bv2)
    vert(l, 0, 0, bv3)
    vert(l, y2, 0, bv4)
    vert(x2, y2, 0, bv5)
    vert(x1, y1, 0, bv6)
    vert(l, l, 0, bv7)
    vert(x2, l, 0, bv8)
    vert(0, l, 0, bv9)
    vert(0, y3, 0, bv10)
    vert(0, r, 0, bv11)

    vert(r, 0, zmax, fv1)
    vert(x3, 0, zmax, fv2)
    vert(l, 0, zmax, fv3)
    vert(l, y2, zmax, fv4)
    vert(x2, y2, zmax, fv5)
    vert(x1, y1, zmax, fv6)
    vert(l, l, zmax, fv7)
    vert(x2, l, zmax, fv8)
    vert(0, l, zmax, fv9)
    vert(0, y3, zmax, fv10)
    vert(0, r, zmax, fv11)
    
);

blocks
(
    hex2D( v6, v5, v10, v11) ( sdivs sdivs zdivs ) simpleGrading (1 1 1)
    hex2D( v1, v2,  v5,  v6) ( sdivs sdivs zdivs ) simpleGrading (1 1 1)
    hex2D( v2, v3,  v4,  v5) ( bdivs sdivs zdivs ) simpleGrading (1 1 1)
    hex2D( v5, v4,  v7,  v8) ( bdivs bdivs zdivs ) simpleGrading (1 1 1)
    hex2D(v10, v5,  v8,  v9) ( sdivs bdivs zdivs ) simpleGrading (1 1 1)
);

edges
(
    createArc(bv1,bv6)   (arcx1 arcy1 0)
    createArc(bv6,bv11)   (arcx2 arcy2 0)
    createArc(bv2,bv5)   (arcx3 arcy3 0)
    createArc(bv5,bv10)   (arcx4 arcy4 0)

    createArc(fv1,fv6)   (arcx1 arcy1 zmax)
    createArc(fv6,fv11)   (arcx2 arcy2 zmax)
    createArc(fv2,fv5)   (arcx3 arcy3 zmax)
    createArc(fv5,fv10)   (arcx4 arcy4 zmax)


);

boundary
(
    defpatch(left, [quad(bv9,bv10,fv10,fv9) quad(bv10,bv11,fv11,fv10)])
    defpatch(right, [quad(bv3,bv4,fv4,fv3) quad(bv4,bv7,fv7,fv4)])
    defpatch(down,[quad(bv1,bv2,fv2,fv1) quad(bv2,bv3,fv3,fv2)])
    defpatch(up,[quad(bv8,bv9,fv9,fv8) quad(bv7,bv8,fv8,fv7)])
    defpatch(hole,[quad(bv11,bv6,fv6,fv11) quad(bv6,bv1,fv1,fv6)])
);

mergePatchPairs
(
);

// ************************************************************************* //
