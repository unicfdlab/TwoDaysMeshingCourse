#!/usr/bin/perl
#
# initial parameters
#
$toMeters=1.0;
$cx=0.0;
$cy=0.5;
$cz=0.5;
$l=1.0;
$h=10.0;
$xdivs=10;
$ydivs=2;
$zdivs=2;
#
# block vertices
#
$minx=$cx;
$maxx=$cx+$h;
$miny=$cy-$l*0.5;
$maxy=$cy+$l*0.5;
$minz=$cz-$l*0.5;
$maxz=$cz+$l*0.5;

%points = ();

$blockMeshHeader = <<BLOCKMESHHEADER;
/*--------------------------------*- C++ -*----------------------------------*\ 
| =========                |                                                 |
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
BLOCKMESHHEADER

print $blockMeshHeader;

&convertToMeters($toMeters);

&newLine();
&startPoints();
&definePoint($minx,$miny,$minz,'bv1');
&definePoint($maxx,$miny,$minz,'bv2');
&definePoint($maxx,$maxy,$minz,'bv3');
&definePoint($minx,$maxy,$minz,'bv4');

&definePoint($minx,$miny,$minz,'fv1');
&definePoint($maxx,$miny,$minz,'fv2');
&definePoint($maxx,$maxy,$minz,'fv3');
&definePoint($minx,$maxy,$minz,'fv4');
&endPoints();

&newLine();
&startBlocks();
&defineHexaBlock('bv1', 'bv2', 'bv3', 'bv4', 'fv1', 'fv2', 'fv3', 'fv4');
&hexaDivs($xdivs, $ydivs, $zdivs);
&endBlocks();

&newLine();
&startEdges();
&endEdges();

&newLine();
&startBoundary();
&definePatch('zfaces',(quad2D('bv1', 'bv2', 'bv3', 'bv4'), quad2D('fv1', 'fv2', 'fv3', 'fv4')));
&definePatch('yfaces',(quad2D('bv1', 'bv2', 'fv2', 'fv1'), quad2D('bv3', 'bv4', 'fv4', 'fv3')));
&definePatch('fixed',quad2D('bv1','bv4','fv4','fv1'));
&definePatch('traction',quad2D('bv2','bv3','fv3','fv2'));
&endBoundary();

&newLine();
&startMergePatchPairs();
&endMergePatchPairs();

$i=0;
while ($i++ < 10)
{
    print "$i\n";
}

print "// ************************************************************************* //\n";
sub newLine()
{
    print "\n";
}
sub definePoint()
{
    $x = $_[0];
    $y = $_[1];
    $z = $_[2];
    $pn= $_[3];
    
    @ks = keys %points;
    $pidx = @ks;
    
    if (exists($points{$pn}))
    {
	print "Found point $pn with duplicating name at index $pidx\n";
	exit ("Found point $pn with duplicating name at index $pidx");
    }
    
    $points{$pn} = $pidx;
    
    print "\t\( $x $y $z \) // point $pidx ($pn)\n";
}
sub startPoints()
{
    print "vertices\n(\n";
}
sub endPoints()
{
    print ");\n";
}

sub hexaDivs()
{
    print " ($_[0] $_[1] $_[2]) simpleGrading (1 1 1)\n";
}
sub defineHexaBlock()
{
    print "\thex ($points{$_[0]} $points{$_[1]} $points{$_[2]} $points{$_[3]} $points{$_[4]} $points{$_[5]} $points{$_[6]} $points{$_[7]}) ";
}
sub startBlocks()
{
    print "blocks\n(\n";
}
sub endBlocks()
{
    print ");\n";
}

sub startEdges()
{
    print "edges\n(\n";
}
sub endEdges()
{
    print ");\n";
}

sub startBoundary()
{
    print "boundary\n(\n";
}
sub endBoundary()
{
    print ");\n";
}
sub quad2D()
{
    return "($points{$_[0]} $points{$_[1]} $points{$_[2]} $points{$_[3]})";
}
sub definePatch()
{
    print "\t$_[0]\n";
    print "\t{\n";
    print "\t\t"."type\t"."patch;\n";
    print "\t\t"."faces\n";
    print "\t\t(\n";
    my @quads = @_;
    shift(@quads);
    for $quad (@quads)
    {
	print "\t\t\t"."$quad\n";
    }
    print "\t\t);\n";
    print "\t}\n";
}

sub startMergePatchPairs()
{
    print "mergePatchPairs\n(\n";
}
sub endMergePatchPairs()
{
    print ");\n";
}

sub convertToMeters()
{
    print "convertToMeters $_[0];\n";
}
#
#END-OF-FILE
#

