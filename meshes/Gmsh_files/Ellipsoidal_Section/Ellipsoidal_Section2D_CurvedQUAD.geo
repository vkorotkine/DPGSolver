// Modifiable Parameters

Refine = 0;

lc = 3/2.0^Refine;


aIn  = 0.5;
aOut = 1.0;

bIn  = 0.5;
bOut = 3.0;

t = Pi/4.0;

r_e = Sqrt(1.0/((Cos(t)/aOut)^2.0+(Sin(t)/bOut)^2.0));


// Geometry Specification

Point(1) = {aIn,0,0,lc};
Point(2) = {aOut,0,0,lc};
Point(3) = {0,bIn,0,lc};
Point(4) = {aIn*Cos(t),aIn*Sin(t),0,lc};
Point(5) = {0,bOut,0,lc};
Point(6) = {r_e*Cos(t),r_e*Sin(t),0,lc};
Point(7) = {0,0,0,lc};


Line(1001)    = {1,2};
Line(1002)    = {3,5};
Ellipse(1003) = {1,7,3,4}; // start, centre, point on major-axis, end
Ellipse(1004) = {4,7,3,3};
Ellipse(1005) = {6,7,5,5};
Ellipse(1006) = {6,7,5,2};
Line(1007)    = {4,6};



Transfinite Line {1001:1007} = 2 Using Progression 1;

Line Loop (4001) = {1001,-1006,-1007,-1003};
Line Loop (4002) = {-1004,1007,1005,-1002};

Plane Surface(4001) = {4001};
Plane Surface(4002) = {4002};

Recombine Surface{4001:4002};




// Physical Parameters for '.msh' file

Physical Line(10011) = {1002}; // Straight Dirichlet
Physical Line(10012) = {1001}; // Straight Neumann
Physical Line(20011) = {1003:1004}; // Curved Dirichlet
Physical Line(20012) = {1005:1006}; // Curved Neumann

Physical Surface(9401) = {4001};
Physical Surface(9402) = {4002};



// Visualization in gmsh

Color Black{ Surface{4001:4002}; }
Geometry.Color.Points = Black;
