Include "../parameters.geo";
//mesh_level = 0; pde_name = ADVECTION; geom_adv = GEOM_ADV_XL;

// Geometry Specification
l = 1;

Point(1) = {-l,0,0,lc};
Point(2) = {+l,0,0,lc};

Line(1001) = {1,2};

Transfinite Line{1001} = 2^(mesh_level)+1 Using Progression 1;



// Physical parameters for '.msh' file
bc_straight = BC_STEP_SC;
bc_base = bc_straight;

If (pde_name == ADVECTION)
	If (geom_adv == GEOM_ADV_XL)
		Physical Point(bc_base+BC_INFLOW)  = {1};
		Physical Point(bc_base+BC_OUTFLOW) = {2};
	Else
		Error("Unsupported geom_adv: %d",geom_adv); Exit;
	EndIf
ElseIf (pde_name == POISSON)
	Physical Point(bc_base+BC_DIRICHLET) = {1,2};
Else
	Error("Unsupported pde_name: %d",pde_name); Exit;
EndIf

Physical Line(9101) = 1001;



// Visualization in gmsh

Color Black{ Line{1001}; }
Geometry.Color.Points = Black;
