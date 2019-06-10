Include "../parameters.geo";
//mesh_level = 2; mesh_type = MIXED; mesh_domain = STRAIGHT; pde_name = ADVECTION; geom_adv = GEOM_ADV_XL; geom_unaligned = 1;
//mesh_level = 2; mesh_type = QUAD; mesh_domain = PARAMETRIC; pde_name = EULER; geom_adv = GEOM_ADV_INTERNAL; geom_unaligned = 0; geom_ar = 0.5;


// ********************************************************************
// Gmsh file used to create an unstructured grid. Works currently for the 
// airfoil and internal channel flow case. 
//
// This file will use higher refinement levels at the points (-1,-1), 
// (0,-1) and (1,-1) as these correspond to the leading and trailing
// edge of the airfoil for the o-grid. 
//
// Set local refinement by using characteristic length. The base 
// characteristic length and local refinement will depend on the mesh level.
// ********************************************************************


// ===================================
//        Geometry Specification
// ===================================

// Geometry Specification
l = 1;
h = 1;
lc = 0.125;

Point(1) = {-l,    -h,-0,lc};  // lower left
Point(2) = {-0.5*l,-h,-0,lc};
Point(3) = {0.0,   -h,-0,lc};
Point(4) = {0.5*l, -h,-0,lc};
Point(5) = {+l,    -h,-0,lc};  // lower right

Point(6) = {-l,    +h,-0,lc};  // upper left
Point(7) = {-0.5*l,+h,-0,lc};
Point(8) = {0.0,   +h,-0,lc};
Point(9) = {0.5*l, +h,-0,lc};
Point(10) ={l,     +h,-0,lc};  // upper right


// All lines are from left to right or bottom to top
//
//         300x
//        ------
//  200x |      | 400x
//       |      | 
//        ------
//         100x

Line(1001) = {1,2};  
Line(1002) = {2,3};
Line(1003) = {3,4};
Line(1004) = {4,5};

Line(2001) = {1,6};

Line(3001) = {6,7};
Line(3002) = {7,8};
Line(3003) = {8,9};
Line(3004) = {9,10};

Line(4001) = {5,10};


// Create the line loop and plane surface
Line Loop (5001) = {1001, 1002, 1003, 1004, 4001, -3004, -3003, -3002, -3001, -2001};
Plane Surface (5001) = {5001};


// ===================================
//         Mesh Refinement
// ===================================

// Perform the refinement

If (pde_name == EULER)
	If (geom_adv == GEOM_NURBS_AIRFOIL_O_GRID)
		If (mesh_level == 0)
			
			// TODO: Implement
			Error("Unsupported mesh_type: %d",mesh_type); Exit;

		ElseIf (mesh_level == 1)
			
			// TODO: Implement
			Error("Unsupported mesh_type: %d",mesh_type); Exit;

		ElseIf (mesh_level == 2)

			// Create transifinite lines so that periodic BCs can be done
			//Transfinite Line {2001, 4001} = 2*2^(mesh_level+1)+1 Using Progression 1.2;

			// Tests
			Characteristic Length {1,3,5} = 0.0125;
			Characteristic Length {2,4} = 0.05;

		ElseIf (mesh_level == 3)

			// TODO: Implement
			Error("Unsupported mesh_type: %d",mesh_type); Exit;

		ElseIf (mesh_level == 4)

			// TODO: Implement
			Error("Unsupported mesh_type: %d",mesh_type); Exit;

		Else

			// TODO: Implement
			Error("Unsupported mesh_type: %d",mesh_type); Exit;

		EndIf
	EndIf
EndIf

// ===================================
//       Boundary Conditions
// ===================================

// TODO: Be able to handle the Channel case too

// Physical parameters for '.msh' file
bc_straight =   BC_STEP_SC;
bc_curved   = 2*BC_STEP_SC;
If (mesh_domain == STRAIGHT)
	bc_base = bc_straight;
Else
	bc_base = bc_curved;
EndIf

If (pde_name == EULER)
	If (geom_adv == GEOM_NURBS_AIRFOIL_O_GRID)
		
		// Periodic left and right x face
		Physical Line(bc_base+PERIODIC_XL) = {2001};
	 	Physical Line(bc_base+PERIODIC_XR) = {4001};
		Periodic Line {4001} = {2001}; // Periodic (x)

		// Slipwall for bottom face and reimann for top
		Physical Line (bc_base+BC_SLIPWALL)= {1001, 1002, 1003, 1004};
		Physical Line (bc_base+BC_RIEMANN) = {3001, 3002, 3003, 3004};
	
	ElseIf (geom_adv == GEOM_ADV_INTERNAL)

		Physical Line(bc_straight+BC_BACKPRESSURE) = {4001};
		Physical Line(bc_straight+BC_TOTAL_TP)     = {2001};
		Physical Line(bc_base+BC_SLIPWALL)         = {1001,1002,1003,1004};
		Physical Line(bc_straight+BC_SLIPWALL)     = {3001,3002,3003,3004};

	Else

		Error("Unsupported geom_adv: %d",geom_adv); Exit;

	EndIf

Else
	Error("Unsupported pde_name: %d",pde_name); Exit;
EndIf


Physical Surface(9401) = 5001;


// Can create either a TRI or QUAD unstructured mesh
If (mesh_type == TRI)
	// Do nothing.
ElseIf (mesh_type == QUAD)
	Recombine Surface{5001};
ElseIf (mesh_type == MIXED)
	Error("Unsupported mesh_type: %d",mesh_type); Exit;
Else
	Error("Unsupported mesh_type: %d",mesh_type); Exit;
EndIf


// Visualization in gmsh
Color Black{ Surface{5001}; }
Geometry.Color.Points = Black;
