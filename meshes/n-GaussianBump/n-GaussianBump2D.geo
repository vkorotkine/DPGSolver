Include "../Parameters.geo";
//Geom_2BEXP = GEOM_2BEXP_1; MeshLevel = 1; MeshType = QUAD; MeshCurving = CURVED; PDEName = EULER;

// Modifiable Parameters



// Geometry Specification
N = 1e2;

a = 0.0625;
c0 = 0.2;
c = c0/2^(Geom_2BEXP/2.0);
l = Sqrt(-2.0*c0^2*Log(EPS/a));
h = 3*a;


xL = -l; x = xL; yL = a*Exp(-x^2/(2*c^2));
xR =  l; x = xR; yR = a*Exp(-x^2/(2*c^2));
xM =  0; x = xM; yM = a*Exp(-x^2/(2*c^2));

Point(1) = {xL,yL,0,lc};
Point(2) = {xR,yR,0,lc};
Point(3) = {xL,h,0,lc};
Point(4) = {xR,h,0,lc};
Point(5) = {xM,yM,0,lc};
Point(6) = {xM,h,0,lc};

pListBL[0] = 1;
For i In {1:N-1}
    x = -l + i/N*l;
    y = a*Exp(-x^2/(2*c^2));

    pListBL[i] = newp;
    Point(pListBL[i]) = {x,y,0,lc};
EndFor
pListBL[N] = 5;

pListBR[0] = 5;
For i In {1:N-1}
    x = 0 + i/N*l;
    y = a*Exp(-x^2/(2*c^2));

    pListBR[i] = newp;
    Point(pListBR[i]) = {x,y,0,lc};
EndFor
pListBR[N] = 2;

Spline(1001) = pListBL[];
Spline(1002) = pListBR[];
Line(1003) = {3,6};
Line(1004) = {6,4};
Line(1005) = {1,3};
Line(1006) = {2,4};
Line(1007) = {5,6};

Transfinite Line {1001:1004} = 3*2^(MeshLevel)+1 Using Progression 1;
Transfinite Line {1005:1007} = 1*2^(MeshLevel)+1 Using Progression 1;

Line Loop (4001) = {1001,1007,-1003,-1005};
Line Loop (4002) = {1002,1006,-1004,-1007};

Plane Surface(4001) = {4001};
Plane Surface(4002) = {4002};

Transfinite Surface{4001};
Transfinite Surface{4002} Right;


If (MeshType == QUAD)
	Recombine Surface{4001,4002};
ElseIf (MeshType == MIXED2D)
	Recombine Surface{4001};
EndIf



// Physical Parameters for '.msh' file
If (PDEName == EULER) // Euler
    Physical Line(BC_STRAIGHT+BC_TOTAL_TP)     = {1005};
    Physical Line(BC_STRAIGHT+BC_BACKPRESSURE) = {1006};
    Physical Line(BC_STRAIGHT+BC_SLIPWALL)     = {1003,1004};
    Physical Line(BC_CURVED  +BC_SLIPWALL)     = {1001,1002};
EndIf

Physical Surface(9401) = {4001:4002};



// Visualization in gmsh

Color Black{ Surface{4001:4002}; }
Geometry.Color.Points = Black;
Coherence;
