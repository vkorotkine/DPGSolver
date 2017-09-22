// Modifiable Parameters

Refine = 0;
Angle = 30;

lc = 0.6/2.0^Refine;

l = 1.0;
h = l;
w = 0.5*l;

r = 0.75*l;
t = 2*Pi-Angle/180*Pi;

xR = r*Cos(t);
yR = r*Sin(t);


// Geometry Specification

Point(0) = {-w,0,0,lc};
Point(1) = {0,0,0,lc};
Point(2) = {xR,yR,0,lc};
Point(3) = {-w,h,0,lc};
Point(4) = {0,h,0,lc};
Point(5) = {xR,h,0,lc};


Line(1001) = {0,3};
Line(1002) = {1,4};
Line(1003) = {2,5};

Line(1004) = {0,1};
Line(1005) = {1,2};
Line(1006) = {3,4};
Line(1007) = {4,5};



Transfinite Line {1001:1003} = 1*2^(Refine)+1 Using Progression 1;
Transfinite Line {1004:1007} = 1*2^(Refine)+1 Using Progression 1;


Line Loop (4001) = {1004,1002,-1006,-1001};
Line Loop (4002) = {1005,1003,-1007,-1002};


Plane Surface(4001) = {4001};
Plane Surface(4002) = {4002};

Transfinite Surface{4001:4002};
//Recombine Surface{4001:40025};




// Physical Parameters for '.msh' file

Physical Line(10005) = {1001};      // Straight Supersonic Inflow
Physical Line(10006) = {1003};      // Straight Supersonic Outflow
Physical Line(10002) = {1004:1007}; // Straight SlipWall


Physical Surface(9401) = {4001};
Physical Surface(9402) = {4002};



// Visualization in gmsh

Color Black{ Surface{4001:4002}; }
Geometry.Color.Points = Black;