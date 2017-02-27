// Modifiable Parameters

Refine = 0;

lc = 0.6/2.0^Refine;

l = 1.0;
h = 3*l;
w = l;

a = 2*l;
b = l;



t = Pi/4.0;
r = Sqrt(1.0/((Cos(t)/a)^2.0+(Sin(t)/b)^2.0));


// Geometry Specification

Point(0) = {0,0,0,lc};
Point(1) = {-l,b,0,lc};
Point(2) = {0,b,0,lc};
Point(3) = {r*Cos(t),r*Sin(t),0,lc};
Point(4) = {a,0,0,lc};
Point(5) = {a,-l,0,lc};

Point(6) = {-l,b+h,0,lc};
Point(7) = {0,b+h,0,lc};
Point(8) = {l,b+h,0,lc};
Point(9) = {a+w,b+h,0,lc};
Point(10) = {a+w,b+l,0,lc};
Point(11) = {a+w,0,0,lc};
Point(12) = {a+w,-l,0,lc};



Line(1001)    = {1,6};
Line(1002)    = {2,7};
Line(1003)    = {3,8};
Line(1004)    = {3,10};
Line(1005)    = {4,11};
Line(1006)    = {5,12};

Line(1007)    = {1,2};
Ellipse(1008) = {2,0,4,3}; // start, centre, point on major-axis, end
Ellipse(1009) = {3,0,4,4};
Line(1010)    = {4,5};

Line(1011)    = {6,7};
Line(1012)    = {7,8};
Line(1013)    = {8,9};
Line(1014)    = {9,10};
Line(1015)    = {10,11};
Line(1016)    = {11,12};


Transfinite Line {1001:1006} = 1*2^(Refine)+1 Using Progression 1;
Transfinite Line {1007:1016} = 1*2^(Refine)+1 Using Progression 1;


Line Loop (4001) = {1007,1002,-1011,-1001};
Line Loop (4002) = {1008,1003,-1012,-1002};
Line Loop (4003) = {1004,-1014,-1013,-1003};
Line Loop (4004) = {1009,1005,-1015,-1004};
Line Loop (4005) = {1010,1006,-1016,-1005};


Plane Surface(4001) = {4001};
Plane Surface(4002) = {4002};
Plane Surface(4003) = {4003};
Plane Surface(4004) = {4004};
Plane Surface(4005) = {4005};

Transfinite Surface{4001:4005};
//Recombine Surface{4001:4005};




// Physical Parameters for '.msh' file

Physical Line(10005) = {1001};                // Straight Supersonic Inflow
Physical Line(10006) = {1006,1014:1016};      // Straight Supersonic Outflow
Physical Line(10002) = {1007,1010,1011:1013}; // Straight SlipWall
Physical Line(20002) = {1008:1009};           // Curved   SlipWall


Physical Surface(9401) = {4001};
Physical Surface(9402) = {4002};
Physical Surface(9403) = {4003};
Physical Surface(9404) = {4004};
Physical Surface(9405) = {4005};




// Visualization in gmsh

Color Black{ Surface{4001:4005}; }
Geometry.Color.Points = Black;