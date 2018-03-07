Include "../parameters.geo";
//mesh_domain = PARAMETRIC; mesh_level = 0; mesh_type = MIXED; pde_name = DIFFUSION; pde_spec = STEADY_DEFAULT; geom_ar = 5;

// Geometry Specification
If (pde_spec == STEADY_SUPERSONIC_VORTEX)
	Include "../../input_files/euler/steady/supersonic_vortex/geometry_parameters.geo";
ElseIf (pde_spec == STEADY_DEFAULT)
	Include "../../input_files/diffusion/steady/default/geometry_parameters.geo";
Else
	Error("Unsupported pde_spec: %d",pde_spec); Exit;
EndIf
Printf("r_i, r_o: %g %g",r_i,r_o);

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

aspect_ratio = geom_ar;

Printf("aspect_ratio ~= %g.",aspect_ratio);
If (aspect_ratio == 1.0)
	Transfinite Line {1003:1004}      = 5*2^(mesh_level)+1  Using Progression 1;
	Transfinite Line {1005:1006}      = 5*2^(mesh_level)+1  Using Progression 1;
	Transfinite Line {1001,1002,1007} = 2*2^(mesh_level)+1  Using Progression 1;
ElseIf (aspect_ratio == 2.5)
	Transfinite Line {1003:1004}      = 1*2^(mesh_level)+1  Using Progression 1;
	Transfinite Line {1005:1006}      = 1*2^(mesh_level)+1  Using Progression 1;
	Transfinite Line {1001,1002,1007} = 1*2^(mesh_level)+1  Using Progression 1;
ElseIf (aspect_ratio == 5.0)
	Transfinite Line {1003:1004}      = 2*2^(mesh_level)+1  Using Progression 1;
	Transfinite Line {1005:1006}      = 2*2^(mesh_level)+1  Using Progression 1;
	Transfinite Line {1001,1002,1007} = 5*2^(mesh_level)+1  Using Progression 1;
ElseIf (aspect_ratio == 20.0)
	Transfinite Line {1003:1004}      = 1*2^(mesh_level)+1  Using Progression 1;
	Transfinite Line {1005:1006}      = 1*2^(mesh_level)+1  Using Progression 1;
	Transfinite Line {1001,1002,1007} = 10*2^(mesh_level)+1 Using Progression 1;
Else
    Error("Unsupported aspect_ratio: %d",aspect_ratio); Exit;
EndIf


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

If (pde_name == DIFFUSION)
	Physical Line (1*BC_STEP_SC+BC_DIRICHLET)    = {1001,1002};
	Physical Line (2*BC_STEP_SC+BC_NEUMANN)      = {1003:1004};
	Physical Line (3*BC_STEP_SC+BC_NEUMANN_ALT1) = {1005:1006};
ElseIf (pde_name == EULER)
	Physical Line (1*BC_STEP_SC+BC_RIEMANN)  = {1001,1002};
	Physical Line (2*BC_STEP_SC+BC_SLIPWALL) = {1003:1004};
	Physical Line (3*BC_STEP_SC+BC_SLIPWALL) = {1005:1006};
Else
	Error("Unsupported pde_name: %d",pde_name); Exit;
EndIf

Physical Surface(9401) = {4001};
Physical Surface(9402) = {4002};



// Visualization in gmsh

Color Black{ Surface{4001:4002}; }
Geometry.Color.Points = Black;
