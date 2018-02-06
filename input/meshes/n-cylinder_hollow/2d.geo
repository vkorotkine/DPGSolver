Include "../parameters.geo";
//mesh_domain = PARAMETRIC; mesh_type = MIXED; pde_name = NAVIERSTOKES; mesh_level = 0; geom_bc = GEOM_BC_ADIABATIC_O; pde_spec = STEADY_TAYLOR_COUETTE;


// Geometry Specification
If (pde_spec == STEADY_TAYLOR_COUETTE)
	Include "../../input_files/navier_stokes/steady/taylor_couette/geometry_parameters.geo";
Else
	Error("Unsupported pde_spec: %d",pde_spec); Exit;
EndIf
Printf("r_i, r_o = %g, %g",r_i,r_o);

If (mesh_domain == PARAMETRIC)
	ind_p = 1; r = r_i;
	Point(ind_p) = {-r,-r,0,lc}; ind_p++;
	Point(ind_p) = { r,-r,0,lc}; ind_p++;
	Point(ind_p) = {-r, r,0,lc}; ind_p++;
	Point(ind_p) = { r, r,0,lc}; ind_p++; r = r_o;
	Point(ind_p) = {-r,-r,0,lc}; ind_p++;
	Point(ind_p) = { r,-r,0,lc}; ind_p++;
	Point(ind_p) = {-r, r,0,lc}; ind_p++;
	Point(ind_p) = { r, r,0,lc}; ind_p++;


	ind_l = 1001; ind_p = 0;
	Line(ind_l) = {1+ind_p,2+ind_p}; ind_l++;
	Line(ind_l) = {3+ind_p,4+ind_p}; ind_l++;
	Line(ind_l) = {1+ind_p,3+ind_p}; ind_l++;
	Line(ind_l) = {2+ind_p,4+ind_p}; ind_l++; ind_p = 4;
	Line(ind_l) = {1+ind_p,2+ind_p}; ind_l++;
	Line(ind_l) = {3+ind_p,4+ind_p}; ind_l++;
	Line(ind_l) = {1+ind_p,3+ind_p}; ind_l++;
	Line(ind_l) = {2+ind_p,4+ind_p}; ind_l++; ind_p = 1;
	Line(ind_l) = {ind_p,ind_p+4};   ind_l++; ind_p++;
	Line(ind_l) = {ind_p,ind_p+4};   ind_l++; ind_p++;
	Line(ind_l) = {ind_p,ind_p+4};   ind_l++; ind_p++;
	Line(ind_l) = {ind_p,ind_p+4};   ind_l++; ind_p++;
ElseIf (mesh_domain == BLENDED)
	ind_p = 0;
	Point(ind_p) = {0,0,0,lc}; ind_p++;
	r = r_i;
	t = 5.0/4.0*Pi; Point(ind_p) = {r*Cos(t),r*Sin(t),0,lc}; ind_p++;
	t = 7.0/4.0*Pi; Point(ind_p) = {r*Cos(t),r*Sin(t),0,lc}; ind_p++;
	t = 3.0/4.0*Pi; Point(ind_p) = {r*Cos(t),r*Sin(t),0,lc}; ind_p++;
	t = 1.0/4.0*Pi; Point(ind_p) = {r*Cos(t),r*Sin(t),0,lc}; ind_p++;
	r = r_o;
	t = 5.0/4.0*Pi; Point(ind_p) = {r*Cos(t),r*Sin(t),0,lc}; ind_p++;
	t = 7.0/4.0*Pi; Point(ind_p) = {r*Cos(t),r*Sin(t),0,lc}; ind_p++;
	t = 3.0/4.0*Pi; Point(ind_p) = {r*Cos(t),r*Sin(t),0,lc}; ind_p++;
	t = 1.0/4.0*Pi; Point(ind_p) = {r*Cos(t),r*Sin(t),0,lc}; ind_p++;


	ind_l = 1001; ind_p = 0;
	Circle(ind_l) = {1+ind_p,0,2+ind_p}; ind_l++;
	Circle(ind_l) = {3+ind_p,0,4+ind_p}; ind_l++;
	Circle(ind_l) = {1+ind_p,0,3+ind_p}; ind_l++;
	Circle(ind_l) = {2+ind_p,0,4+ind_p}; ind_l++; ind_p = 4;
	Circle(ind_l) = {1+ind_p,0,2+ind_p}; ind_l++;
	Circle(ind_l) = {3+ind_p,0,4+ind_p}; ind_l++;
	Circle(ind_l) = {1+ind_p,0,3+ind_p}; ind_l++;
	Circle(ind_l) = {2+ind_p,0,4+ind_p}; ind_l++; ind_p = 1;
	Line(ind_l) = {ind_p,ind_p+4}; ind_l++; ind_p++;
	Line(ind_l) = {ind_p,ind_p+4}; ind_l++; ind_p++;
	Line(ind_l) = {ind_p,ind_p+4}; ind_l++; ind_p++;
	Line(ind_l) = {ind_p,ind_p+4}; ind_l++; ind_p++;
Else
	Error("Unsupported mesh_domain: %d",mesh_domain); Exit;
EndIf


// Allows different (A)spect (R)atio elements
aspect_ratio = 1.0;
If (aspect_ratio == 1.0)
	Transfinite Line {1001:1008} = 2^(mesh_level+1)+1 Using Progression 1;
	Transfinite Line {1009:1012} = 2^(mesh_level)+1 Using Progression 1;
ElseIf (aspect_ratio == 2.0)
	Transfinite Line {1001:1008} = 2^(mesh_level)+1 Using Progression 1;
	Transfinite Line {1009:1012} = 2^(mesh_level)+1 Using Progression 1;
ElseIf (aspect_ratio == 4.0)
	Transfinite Line {1001:1008} = 2^(mesh_level)+1 Using Progression 1;
	Transfinite Line {1009:1012} = 2^(mesh_level+1)+1 Using Progression 1;
ElseIf (aspect_ratio == 8.0)
	Transfinite Line {1001:1008} = 2^(mesh_level)+1 Using Progression 1;
	Transfinite Line {1009:1012} = 2^(mesh_level+2)+1 Using Progression 1;
ElseIf (aspect_ratio == 16.0)
	Transfinite Line {1001:1008} = 2^(mesh_level)+1 Using Progression 1;
	Transfinite Line {1009:1012} = 2^(mesh_level+3)+1 Using Progression 1;
EndIf

ind_l = 4001;
Line Loop (ind_l) = {1005,-1010,-1001,1009}; ind_l++;
Line Loop (ind_l) = {1002,1012,-1006,-1011}; ind_l++;
Line Loop (ind_l) = {1003,1011,-1007,-1009}; ind_l++;
Line Loop (ind_l) = {1008,-1012,-1004,1010}; ind_l++;

ind_l = 4001;
For (1:4)
//	Plane Surface(ind_l) = {ind_l}; ind_l++; // Makes a more equilateral mesh
	Plane Surface(ind_l) = {ind_l}; Transfinite Surface{ind_l}; ind_l++;
EndFor


If (mesh_type == MIXED)
	Recombine Surface{4001,4003};
ElseIf (mesh_type == QUAD)
	Recombine Surface{4001:4004};
EndIf



// Physical parameters for '.msh' file
bc_curved   = 2*BC_STEP_SC;

If (pde_name == NAVIER_STOKES)
	Physical Line (bc_curved+BC_NOSLIP_ALL_ROTATING) = {1001:1004};
	If (geom_bc == GEOM_BC_ADIABATIC_O)
		Physical Line (bc_curved+BC_NOSLIP_ADIABATIC) = {1005:1008};
	ElseIf (geom_bc == GEOM_BC_DIABATIC_O)
		Physical Line (bc_curved+BC_NOSLIP_DIABATIC) = {1005:1008};
	Else
		Error("Unsupported geom_bc: %d",geom_bc); Exit;
Else
	Error("Unsupported pde_name: %d",pde_name); Exit;
EndIf

ind_v = 9401; ind_s = 4001;
For (1:4)
	Physical Surface(ind_v) = {ind_s}; ind_v++; ind_s++;
EndFor



// Visualization in gmsh

Color Black{ Surface{4001:4004}; }
Geometry.Color.Points = Black;
