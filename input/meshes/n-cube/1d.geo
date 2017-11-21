Include "../parameters.geo";
//MESH_LEVEL = 0; PDE_NAME = ADVECTION; GEOM_ADV = GEOM_ADV_XL;

// Geometry Specification
L = 1;

Point(1) = {-L,0,0,lc};
Point(2) = {+L,0,0,lc};

Line(1001) = {1,2};

Transfinite Line{1001} = 2^(MESH_LEVEL)+1 Using Progression 1;



// Physical parameters for '.msh' file
BC_Straight = BC_STEP_SC;
BC_Base = BC_Straight;

If (PDE_NAME == ADVECTION)
	If (GEOM_ADV == GEOM_ADV_XL)
		Physical Point(BC_Base+BC_INFLOW)  = {1};
		Physical Point(BC_Base+BC_OUTFLOW) = {2};
	Else
		Error("Unsupported GEOM_ADV: %d",GEOM_ADV); Exit;
	EndIf
ElseIf (PDE_NAME == POISSON)
	Physical Point(BC_Base+BC_DIRICHLET) = {1,2};
EndIf

Physical Line(9101) = 1001;



// Visualization in gmsh

Color Black{ Line{1001}; }
Geometry.Color.Points = Black;
