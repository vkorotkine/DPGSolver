Include "../parameters.geo";
//mesh_level = 0; mesh_type = MIXED; mesh_domain = PARAMETRIC; pde_name = EULER; geom_conformal = GEOM_CONFORMAL_FULL; pde_spec = STEADY_JOUKOWSKI; geom_ar = GEOM_AR_4;

// Geometry Specification
If (pde_name == EULER && pde_spec == STEADY_JOUKOWSKI)
	Include "../../input_files/euler/steady/joukowski/geometry_parameters.geo";
ElseIf (pde_name == NAVIER_STOKES && pde_spec == STEADY_JOUKOWSKI)
	Include "../../input_files/navier_stokes/steady/joukowski/geometry_parameters.geo";
Else
	Error("Unsupported pde_name, pde_spec: %d, %d",pde_name,pde_spec); Exit;
EndIf
Printf("s_offset = %g.",s_offset);

b = s_offset;
l = 1;
h = 1;

all_transfinite  = 1;
bl_recombined    = 1;
all_recombined   = 1;
use_bump         = 0;
t_progression_x  = 3.0;
t_bump_x         = 0.15;
t_progression_j  = 1.2;

// Allows different (A)spect (R)atio elements
aspect_ratio = geom_ar;
If (aspect_ratio == 1)
	t_progression_y = 4.0;
	n_y = 2^(mesh_level+2)+1;
ElseIf (aspect_ratio == 2)
	t_progression_y = 5.5;
	n_y = 2^(mesh_level+2)+1;
ElseIf (aspect_ratio == 4)
	t_progression_y = 4.0;
	n_y = 2^(mesh_level+2)+2;
Else
    Error("Unsupported aspect_ratio: %d",aspect_ratio); Exit;
EndIf
Printf("aspect_ratio ~= %g.",aspect_ratio);

Point(0) = {-l,b,-0,lc};
Point(1) = {0, b,-0,lc};
Point(2) = {+l,b,-0,lc};
Point(3) = {-l,b+h,-0,lc};
Point(4) = {0, b+h,-0,lc};
Point(5) = {+l,b+h,-0,lc};

Line(1000) = {0,1};
Line(1001) = {1,2};
Line(1002) = {3,4};
Line(1003) = {4,5};

Line(2000) = {0,3};
Line(2001) = {1,4};
Line(2002) = {2,5};


Line Loop(4000) = {1000,2001,-1002,-2000};
Line Loop(4001) = {1001,2002,-1003,-2001};

Plane Surface(4000) = {4000};
Plane Surface(4001) = {4001};

Transfinite Line{2000:2002} = n_y Using Progression t_progression_y;
If (use_bump)
	Transfinite Line{1000,1002} = 2^(mesh_level+2)+1 Using Bump t_bump_x;
Else
	Transfinite Line{1000,1002} = 2^(mesh_level+2)+1 Using Progression t_progression_j;
EndIf
Transfinite Line{1001,1003} = 2^(mesh_level+2)+1 Using Progression t_progression_x;

If (all_transfinite)
	Transfinite Surface {4000:4001};
Else
	Transfinite Surface {4000};
EndIf

If (bl_recombined)
	Recombine Surface{4000};
EndIf
If (all_recombined)
	Recombine Surface{4001};
EndIf


// Physical parameters for '.msh' file
bc_s = 1*BC_STEP_SC;
bc_c = 2*BC_STEP_SC;

If (geom_conformal == GEOM_CONFORMAL_HALF)
	Physical Line(bc_c+BC_RIEMANN)  = {1002};
	Physical Line(bc_s+BC_RIEMANN)  = {1003,2002};
	Physical Line(bc_s+BC_SLIPWALL) = {2000,1001};
	If (pde_name == EULER)
		Physical Line(bc_c+BC_SLIPWALL) = {1000};
	ElseIf (pde_name == NAVIER_STOKES)
		Physical Line(bc_c+BC_NOSLIP_ADIABATIC) = {1000};
	Else
		Error("Unsupported pde_name: %d",pde_name); Exit;
	EndIf
EndIf

Physical Surface(9400) = 4000;
Physical Surface(9401) = 4001;



// Visualization in gmsh

Color Black{ Surface{4000:4001}; }
Geometry.Color.Points = Black;

If (geom_conformal == GEOM_CONFORMAL_FULL)
	Physical Point(9998) = {0,0,0,lc}; // Required to fix the numbering of the symmetric entities below.

	Symmetry{ 0.0,-1.0,0.0,0.0 }{Duplicata{Surface{4000:4001};}}

	Transfinite Line{-10003,10001,10006} = n_y Using Progression t_progression_y;
	If (use_bump)
		Transfinite Line{10000,-10002} = 2^(mesh_level+2)+1 Using Bump t_bump_x;
	Else
		Transfinite Line{10000,-10002} = 2^(mesh_level+2)+1 Using Progression t_progression_j;
	EndIf

	Transfinite Line{10005,-10007}       = 2^(mesh_level+2)+1 Using Progression t_progression_x;
	If (all_transfinite)
		Transfinite Surface {9999,10004};
	Else
		Transfinite Surface {9999};
	EndIf
	If (bl_recombined)
		Recombine Surface{9999};
	EndIf
	If (all_recombined)
		Recombine Surface{10004};
	EndIf

	Physical Line(bc_c+BC_RIEMANN) = {1002,10002};
	Physical Line(bc_s+BC_RIEMANN) = {1003,2002,10007,10006};

	Physical Line(bc_s+PERIODIC_XL_REFLECTED_Y) = {2000};
	Physical Line(bc_s+PERIODIC_XR_REFLECTED_Y) = {10003};
	Physical Line(bc_s+PERIODIC_YL)             = {1001};
	Physical Line(bc_s+PERIODIC_YR)             = {10005};

	// Periodic Indicator (Slave = Master)
	Periodic Line {10003} = {-2000}; // Periodic (x)
	Periodic Line {10005} = {1001};  // Periodic (y)

	If (pde_name == EULER)
		Physical Line(bc_c+BC_SLIPWALL) = {1000,10000};
	ElseIf (pde_name == NAVIER_STOKES)
		Physical Line(bc_c+BC_NOSLIP_ADIABATIC) = {1000,10000};
	Else
		Error("Unsupported pde_name: %d",pde_name); Exit;
	EndIf

	Physical Surface(9402) = -9999;
	Physical Surface(9403) = -10004;


	Color Black{ Surface{9999,10004}; }
	Geometry.Color.Points = Black;
EndIf
