Include "../parameters.geo";
//mesh_domain = BLENDED; mesh_level = 0; mesh_type = MIXED; pde_name = EULER; pde_spec = STEADY_SUPERSONIC_VORTEX;

// Geometry Specification
If (pde_spec == STEADY_SUPERSONIC_VORTEX)
	Include "../../input_files/euler/steady/supersonic_vortex/geometry_parameters.geo";
EndIf

If (mesh_domain == PARAMETRIC)
	Point(1) = {r_i,0,0,lc};
	Point(2) = {r_o,0,0,lc};
	Point(3) = {0,r_i,0,lc};
	Point(4) = {r_i,r_i,0,lc};
	Point(5) = {0,r_o,0,lc};
	Point(6) = {r_o,r_o,0,lc};

	Line(1001) = {1,2};
	Line(1002) = {3,5};
	Line(1003) = {4,3};
	Line(1004) = {4,1};
	Line(1005) = {6,5};
	Line(1006) = {6,2};
	Line(1007) = {4,6};
ElseIf (mesh_domain == BLENDED)
	Point(1) = {r_i,0,0,lc};
	Point(2) = {r_o,0,0,lc};
	Point(3) = {0,r_i,0,lc};
	Point(4) = {Sqrt(0.5)*r_i,Sqrt(0.5)*r_i,0,lc};
	Point(5) = {0,r_o,0,lc};
	Point(6) = {Sqrt(0.5)*r_o,Sqrt(0.5)*r_o,0,lc};
	Point(7) = {0,0,0,lc};

	Line(1001)   = {1,2};
	Line(1002)   = {3,5};
	Circle(1003) = {4,7,3};
	Circle(1004) = {4,7,1};
	Circle(1005) = {6,7,5};
	Circle(1006) = {6,7,2};
	Line(1007)   = {4,6};
Else
	Error("Unsupported mesh_domain: %d",mesh_domain); Exit;
EndIf


// Include something for aspect ratio: 1.0, 2.5, 5.0, 20.0
Transfinite Line {1003:1006}      = 2^(mesh_level)+1 Using Progression 1;
Transfinite Line {1001,1002,1007} = 2^(mesh_level)+1 Using Progression 1;


Line Loop (4001) = {1007,1005,-1002,-1003};
Line Loop (4002) = {-1007,1004,1001,-1006};

Plane Surface(4001) = {4001};
Plane Surface(4002) = {4002};

Transfinite Surface{4001} Left;
Transfinite Surface{4002} Right;

If (mesh_type == MIXED)
	Recombine Surface{4002};
ElseIf (mesh_type == QUAD)
	Recombine Surface{4001,4002};
EndIf



// Physical parameters for '.msh' file
BC_Straight =   BC_STEP_SC;
BC_Curved   = 2*BC_STEP_SC;

If (pde_name == EULER)
	Physical Line (1*BC_STEP_SC+BC_RIEMANN)  = {1001,1002};
	Physical Line (2*BC_STEP_SC+BC_SLIPWALL) = {1003:1004};
	Physical Line (3*BC_STEP_SC+BC_SLIPWALL) = {1005:1006};
EndIf

Physical Surface(9401) = {4001};
Physical Surface(9402) = {4002};



// Visualization in gmsh

Color Black{ Surface{4001:4002}; }
Geometry.Color.Points = Black;
