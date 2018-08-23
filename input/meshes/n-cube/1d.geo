Include "../parameters.geo";
//mesh_level = 2; pde_name = ADVECTION; geom_adv = GEOM_ADV_PERIODIC;
//mesh_level = 2; pde_name = BURGERS_INVISCID; geom_adv = GEOM_ADV_PERIODIC;

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
		Physical Point(bc_base+BC_UPWIND)  = {1};
		Physical Point(bc_base+BC_OUTFLOW) = {2};
	ElseIf (geom_adv == GEOM_ADV_PERIODIC)
		Physical Point(bc_base+PERIODIC_XL) = {1};
		Physical Point(bc_base+PERIODIC_XR) = {2};

		// Periodic Indicator (Slave = Master). Note: There are no "Periodic Point"s.
		Line(1) = {1,1};
		Line(2) = {2,2};

		Periodic Line {2} = {1}; // Periodic (x)
	Else
		Error("Unsupported geom_adv: %d",geom_adv); Exit;
	EndIf
ElseIf (pde_name == DIFFUSION)
	Physical Point(bc_base+BC_DIRICHLET) = {1};
	Physical Point(bc_base+BC_NEUMANN)   = {2};
ElseIf (pde_name == BURGERS_INVISCID)
	If (geom_adv == GEOM_ADV_PERIODIC)
		Physical Point(bc_base+PERIODIC_XL) = {1};
		Physical Point(bc_base+PERIODIC_XR) = {2};

		// Periodic Indicator (Slave = Master). Note: There are no "Periodic Point"s.
		Line(1) = {1,1};
		Line(2) = {2,2};

		Periodic Line {2} = {1}; // Periodic (x)
	Else
		Error("Unsupported geom_adv: %d",geom_adv); Exit;
	EndIf
Else
	Error("Unsupported pde_name: %d",pde_name); Exit;
EndIf

Physical Line(9101) = 1001;



// Visualization in gmsh

Color Black{ Line{1001}; }
Geometry.Color.Points = Black;
