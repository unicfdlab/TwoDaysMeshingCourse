dnl>
dnl> MESH DESCRIPTION FOR THE BEAM LOAD
dnl>
dnl>
dnl> Initial parameters
dnl>
define(cx, 0.0)dnl>
define(cy, 0.0)dnl>
define(cz, 0.0)dnl>
define(l, 1.0)dnl>
define(h,10.0)dnl>
define(xdivs, 10)dnl>
define(ydivs, 4)dnl>
define(zdivs, 4)dnl>
changequote([,]) dnl>
define(calc, [esyscmd(perl -e 'print ($1)')]) dnl>
dnl>
dnl> Intermediate variables
dnl>
define(minx,cx)dnl>
define(maxx,calc(minx+h))dnl>
define(miny,calc(cy-l*0.5))dnl>
define(maxy,calc(cy+l*0.5))dnl>
define(minz,calc(cz-l*0.5))dnl>
define(maxz,calc(cz+l*0.5))dnl>
dnl>
**
**
** POINTS
**
**
PNT Pbv1 minx miny minz
PNT Pbv2 maxx miny minz
PNT Pbv3 maxx maxy minz
PNT Pbv4 minx maxy minz
**
PNT Pfv1 minx miny maxz
PNT Pfv2 maxx miny maxz
PNT Pfv3 maxx maxy maxz
PNT Pfv4 minx maxy maxz
**
**
** LINES
**
**
LINE Lb1 Pbv1 Pbv2 xdivs
LINE Lb2 Pbv2 Pbv3 ydivs
LINE Lb3 Pbv3 Pbv4 xdivs
Line Lb4 Pbv4 Pbv1 ydivs
**
LINE Lf1 Pfv1 Pfv2 xdivs
LINE Lf2 Pfv2 Pfv3 ydivs
LINE Lf3 Pfv3 Pfv4 xdivs
Line Lf4 Pfv4 Pfv1 ydivs
**
LINE Lbf1 Pbv1 Pfv1 zdivs
LINE Lbf2 Pbv2 Pfv2 zdivs
LINE Lbf3 Pbv3 Pfv3 zdivs
LINE Lbf4 Pbv4 Pfv4 zdivs
**
**
** SURFACES
**
**
SURF SZ1 Lb1 Lb2 Lb3 Lb4
SURF SZ2 Lf1 Lf2 Lf3 Lf4
SURF SY1 Lb1 Lbf2 Lf1 Lbf1
SURF SY2 Lb3 Lbf4 Lf3 Lbf3
SURF SX1 Lb4 Lbf4 Lf4 Lbf1
SURF SX2 Lb2 Lbf3 Lf2 Lbf2
**
**
** BODIES
**
**
BODY B1 SZ1 SZ2 SY1 SY2 SX1 SX2
**
**
** MAKE SETS FOR NODES
**
**
** FIXED SURFACE
seta fixed_ae SX1
seta fixed_ae Lb4
seta fixed_ae Lf4
seta fixed_ae Lbf1
seta fixed_ae Lbf4
seta fixed_ae Pfv1
seta fixed_ae Pbv1
seta fixed_ae Pfv4
seta fixed_ae Pbv4
** TRACTION SURFACE
**seta traction_ae SX2
**seta traction_ae Lb2
**seta traction_ae Lf2
**seta traction_ae Lbf2
**seta traction_ae Lbf3
seta traction_ae Pfv2
seta traction_ae Pfv3
seta traction_ae Pbv2
seta traction_ae Pbv3
** YFACES SURFACES
seta yfaces_ae SY1
seta yfaces_ae SY2
comp yfaces_ae do
** ZFACES SURFACES
seta zfaces_ae SZ1
seta zfaces_ae SZ2
comp zfaces_ae do
**
**
** MESHING ALL WITH HE8
**
**
elty all he8
mesh all
seta fixed n fixed_ae
seta traction n traction_ae
seta yfaces n yfaces_ae
seta zfaces n zfaces_ae
**
**
** SEND DATA TO FILES
**
**
send all abq
send fixed abq
send traction abq
send yfaces abq
send zfaces abq
**
**END-OF-FILE
**

