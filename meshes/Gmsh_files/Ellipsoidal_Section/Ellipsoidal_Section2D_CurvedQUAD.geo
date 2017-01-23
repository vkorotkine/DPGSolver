// Modifiable Parameters

Refine = 0;

lc = 3/2.0^Refine;

aIn  = 0.25; aOut = 0.75;
bIn  = 0.75; bOut = 2.0;

//aIn  = 0.25;  aOut = 1.5;
//bIn  = 0.375; bOut = 2.0;



// Geometry Specification

Point(1) = {aIn,0,0,lc};
Point(2) = {aOut,0,0,lc};
Point(3) = {0,bIn,0,lc};
//Point(4) = {aIn*Cos(Pi/4.0),bIn*Sin(Pi/4.0),0,lc};
Point(5) = {0,bOut,0,lc};
//Point(6) = {aOut*Cos(Pi/5.0),bOut*Sin(Pi/5.0),0,lc};
Point(7) = {0,0,0,lc};

Line(1001)    = {1,2};
Line(1002)    = {3,5};
Ellipse(1003) = {1,7,3,3}; // start, centre, point on major-axis, end
//Ellipse(1004) = {4,7,3,1};
Ellipse(1005) = {2,7,5,5};
//Ellipse(1006) = {6,7,5,2};



Transfinite Line {1001:1003,1005} = 2 Using Progression 1;

Line Loop (4001) = {1001,1005,-1002,-1003};

Plane Surface(4001) = {4001};

Recombine Surface(4001);




// Physical Parameters for '.msh' file

Physical Line(10011) = {1002}; // Straight Dirichlet
Physical Line(10012) = {1001}; // Straight Neumann
Physical Line(20011) = {1003}; // Curved Dirichlet
Physical Line(20012) = {1005}; // Curved Neumann

Physical Surface(9401) = {4001};




// Visualization in gmsh

Color Black{ Surface{4001}; }
Geometry.Color.Points = Black;
