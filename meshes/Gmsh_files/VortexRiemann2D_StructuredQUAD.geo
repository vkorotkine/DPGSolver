/*********************************************************************
 *
 *  Gmsh
 *
 *  Riemann/Periodic Vortex 2D
 *
 *********************************************************************/





TO BE EDITED FOR NEW CONVENTION.





/*
Use ‘Transfinite Surface’ for structured tri mesh; for an unstructured mesh, comment it.
Uncomment ‘Recombine Surface{1}’ to convert to a quad mesh.
See page 103 of the user manual for explanation of numbering in the .msh file.
*/

Refine = 3;

lc = 1;

L = 8;
H = 8;

Point(1) = {-L,-H,-0,lc};
Point(2) = {+L,-H,-0,lc};
Point(3) = {-L,+H,-0,lc};
Point(4) = {+L,+H,-0,lc};

Line(101) = {1,2};
Line(102) = {3,4};
Line(103) = {1,3};
Line(104) = {2,4};

Transfinite Line{101:102} = 2^(Refine)+1 Using Progression 1;
Transfinite Line{103:104} = 2^(Refine)+1 Using Progression 1;

Line Loop (1001) = {101,104,-102,-103};

Plane Surface(201) = {1001};

Transfinite Surface {201};
Recombine Surface{201};

// Apply an elliptic smoother to the grid
// Mesh.Smoothing = 100;

Color Black{ Surface{201}; }
Geometry.Color.Points = Black;

// Physical parameters for '.msh' file
Physical Surface(100) = 201;

//Riemann
Physical Point(1003) = {1:4};     // Riemann

Physical Line(1003)  = {101:104}; // Riemann


// Periodic
/*
Physical Point(2051) = {1,3}; // Periodic (-x)
Physical Point(2052) = {2,4}; // Periodic (+x)
Physical Point(2053) = {1,2}; // Periodic (-y)
Physical Point(2054) = {3,4}; // Periodic (+y)

Physical Line(2051) = {103}; // Periodic -x
Physical Line(2052) = {104}; // Periodic +x
Physical Line(2053) = {101}; // Periodic -y
Physical Line(2054) = {102}; // Periodic +y


// Periodic Indicator
Periodic Line {101} = {-102};
Periodic Line {103} = {-104};
*/
