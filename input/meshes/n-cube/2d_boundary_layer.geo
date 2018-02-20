Include "../parameters.geo";
//mesh_level = 1; mesh_type = MIXED; mesh_domain = PARAMETRIC; pde_name = EULER; geom_conformal = GEOM_CONFORMAL_FULL; pde_spec = STEADY_JOUKOWSKI;

// Geometry Specification
If (pde_name == EULER && pde_spec == STEADY_JOUKOWSKI)
	Include "../../input_files/euler/steady/joukowski/geometry_parameters.geo";
Else
	Error("Unsupported pde_name, pde_spec: %d, %d",pde_name,pde_spec); Exit;
EndIf
Printf("s_offset = %g.",s_offset);

b = s_offset;
l = 1;
h = 1;
bl = 0.1*h;

t_progression_bl = 1.2;
t_progression_x  = 1.4;
t_progression_y  = 1.4;
t_bump_x         = 0.2;

Point(0) = {-l,b,-0,lc};
Point(1) = {0, b,-0,lc};
Point(2) = {+l,b,-0,lc};
Point(3) = {-l,b+bl,-0,lc};
Point(4) = {0, b+bl,-0,lc};
Point(5) = {+l,b+bl,-0,lc};
Point(6) = {-l,b+h,-0,lc};
Point(7) = {0, b+h,-0,lc};
Point(8) = {+l,b+h,-0,lc};

Line(1000) = {0,1};
Line(1001) = {1,2};
Line(1002) = {3,4};
Line(1003) = {4,5};
Line(1004) = {6,7};
Line(1005) = {7,8};

Line(2000) = {0,3};
Line(2001) = {3,6};
Line(2002) = {1,4};
Line(2003) = {4,7};
Line(2004) = {2,5};
Line(2005) = {5,8};


Line Loop(4000) = {1000,2002,-1002,-2000};
Line Loop(4001) = {1002,2003,-1004,-2001};
Line Loop(4002) = {1001,2004,-1003,-2002};
Line Loop(4003) = {1003,2005,-1005,-2003};

Plane Surface(4000) = {4000};
Plane Surface(4001) = {4001};
Plane Surface(4002) = {4002};
Plane Surface(4003) = {4003};

all_transfinite = 1;
If (all_transfinite)
	Transfinite Line{2000,2002,2004} = 2^(mesh_level+1)+1 Using Progression t_progression_bl;
	Transfinite Line{2001,2003,2005} = 2^(mesh_level+1)+0 Using Progression t_progression_y;
	Transfinite Line{1000,1002,1004} = 2^(mesh_level+1)+1 Using Bump t_bump_x;
	Transfinite Line{1001,1003,1005} = 2^(mesh_level+1)+1 Using Progression t_progression_x;
	Transfinite Surface {4000:4003};
Else
	Transfinite Line{2000,2002} = 2^(mesh_level+1)+1 Using Progression t_progression_bl;
	Transfinite Line{2001,2003} = 2^(mesh_level+1)+0 Using Progression t_progression_y;
	Transfinite Line{1000,1002,1004} = 2^(mesh_level+1)+1 Using Bump t_bump_x;
	Transfinite Line{1001,1003} = 2^(mesh_level+1)+0 Using Progression t_progression_x;
	Transfinite Surface {4000:4001};
EndIf

all_recombined = 0;
Recombine Surface{4000};
If (all_recombined)
	Recombine Surface{4001,4002,4003};
EndIf


// Physical parameters for '.msh' file
bc_s = 1*BC_STEP_SC;
bc_c = 2*BC_STEP_SC;

If (geom_conformal == GEOM_CONFORMAL_HALF)
	Physical Line(bc_c+BC_RIEMANN)  = {1004};
	Physical Line(bc_s+BC_RIEMANN)  = {1005,2004,2005};
	Physical Line(bc_s+BC_SLIPWALL) = {2000,2001,1001};
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
Physical Surface(9402) = 4002;
Physical Surface(9403) = 4003;



// Visualization in gmsh

Color Black{ Surface{4000:4003}; }
Geometry.Color.Points = Black;

If (geom_conformal == GEOM_CONFORMAL_FULL)
	Physical Point(20102) = {0,0,0,lc}; // Required to fix the numbering of the symmetric entities below.

	Symmetry{ 0.0,-1.0,0.0,0.0 }{Duplicata{Surface{4000:4003};}}
	If (all_transfinite)
		Transfinite Line{20107,20105,20115} = 2^(mesh_level+1)+1 Using Progression t_progression_bl;
		Transfinite Line{20112,20110,20120} = 2^(mesh_level+1)+0 Using Progression t_progression_y;
		Transfinite Line{20104,20106,20111} = 2^(mesh_level+1)+1 Using Bump t_bump_x;
		Transfinite Line{20114,-20116,-20121} = 2^(mesh_level+1)+1 Using Progression t_progression_x;
		Transfinite Surface {20103,20108,20113,20118};
	Else
		Transfinite Line{20107,20105} = 2^(mesh_level+1)+1 Using Progression t_progression_bl;
		Transfinite Line{20112,20110} = 2^(mesh_level+1)+0 Using Progression t_progression_y;
		Transfinite Line{20104,20106,20111} = 2^(mesh_level+1)+1 Using Bump t_bump_x;
		Transfinite Line{20114,-20116} = 2^(mesh_level+1)+0 Using Progression t_progression_x;
		Transfinite Surface {20103,20108};
	EndIf
	Recombine Surface{20103};
	If (all_recombined)
		Recombine Surface{20108,20113,20118};
	EndIf

	Physical Line(bc_c+BC_RIEMANN) = {1004,20111};
	Physical Line(bc_s+BC_RIEMANN) = {1005,2004,2005,20121,20115,20120};

	Physical Line(bc_s+PERIODIC_XL_REFLECTED_Y) = {2000:2001};
	Physical Line(bc_s+PERIODIC_XR_REFLECTED_Y) = {20107,20112};
	Physical Line(bc_s+PERIODIC_YL)             = {1001};
	Physical Line(bc_s+PERIODIC_YR)             = {20114};

	// Periodic Indicator (Slave = Master)
	Periodic Line {20107} = {-2000}; // Periodic (x)
	Periodic Line {20112} = {-2001};
	Periodic Line {20114} = {1001};  // Periodic (y)

	If (pde_name == EULER)
		Physical Line(bc_c+BC_SLIPWALL) = {1000,20104};
	ElseIf (pde_name == NAVIER_STOKES)
		Physical Line(bc_c+BC_NOSLIP_ADIABATIC) = {1000,20104};
	Else
		Error("Unsupported pde_name: %d",pde_name); Exit;
	EndIf

	Physical Surface(9404) = -20103;
	Physical Surface(9405) = -20108;
	Physical Surface(9406) = -20113;
	Physical Surface(9407) = -20118;


	Color Black{ Surface{20103,20108,20113,20118}; }
	Geometry.Color.Points = Black;
EndIf
