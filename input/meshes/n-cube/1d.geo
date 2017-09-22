Include "../Parameters.geo";
//MeshLevel = 2; PDEName = POISSON;

// Geometry Specification
L = 1;

Point(1) = {-L,-0,-0,lc};
Point(2) = {+L,-0,-0,lc};

Line(1001) = {1,2};

Transfinite Line{1001} = 2^(MeshLevel)+1 Using Progression 1;




// Physical parameters for '.msh' file
BC_Straight = BC_STEP_SC;
BC_Base = BC_Straight;

If (PDEName == POISSON)
	Physical Point(BC_Base+BC_DIRICHLET) = {1,2};
EndIf

Physical Line(9101) = 1001;



// Visualization in gmsh

Color Black{ Line{1001}; }
Geometry.Color.Points = Black;
