Include "Parameters.geo";

RECOMBINE = 0;
If (MeshType == TRI)
	Printf(Str("In TRI"));
ElseIf (MeshType == QUAD)
	Printf(Str("In QUAD"));
	RECOMBINE = 1;
EndIf



lc = 1;

L = 1;
H = 1;

/*
test = "Test";
Printf(test);
Printf("%f",TRI);
Printf("%g",QUAD);
Printf("%e",MeshType);
Printf(Str("test"));
*/


// Geometry Specification

Point(1) = {-L,-H,-0,lc};
Point(2) = {+L,-H,-0,lc};
Point(3) = {-L,+H,-0,lc};
Point(4) = {+L,+H,-0,lc};

Line(1001) = {1,2};
Line(1002) = {3,4};
Line(2001) = {1,3};
Line(2002) = {2,4};

Transfinite Line{1001:1002} = 2 Using Progression 1;
Transfinite Line{2001:2002} = 2 Using Progression 1;

Line Loop (4001) = {1001,2002,-1002,-2001};

Plane Surface(4001) = {4001};

Transfinite Surface {4001};
If (RECOMBINE)
	Recombine Surface{4001};
EndIf



// Physical parameters for '.msh' file

Physical Line(10001) = {1001:1002,2001:2002};

Physical Surface(9401) = 4001;



// Visualization in gmsh

Color Black{ Surface{4001}; }
Geometry.Color.Points = Black;
