// Modifiable Parameters

Refine = 0;

lc = 0.5/2.0^Refine;

// RefType Options: 0 (Std. transfinite), 1 (Special transfinite), 2 (Std.)
RefType = 1;

LC_ratio = 2;
Nc = 1e3;

// radii determined here assumed l = 1.0;
If (LC_ratio == 1)
	r = 2.54161;
ElseIf (LC_ratio == 2)
	r = 1.52613;
ElseIf (LC_ratio == 10)
	r = 0.363683;
ElseIf (LC_ratio == 100)
	r = 0.0380061;
ElseIf (LC_ratio == 1000)
	r = 0.0038178;
EndIf

// Geometry Specification
Nl = Nc*LC_ratio;

h = 0.5;
l = 1.0;


xR = l*Cos(15*Pi/180);
yR = l*Sin(15*Pi/180);

Xc = -r/Tan(82.5*Pi/180);
Yc = r;


Point(1) = {-l,0,0,lc};
Point(2) = {xR,yR,0,lc};
Point(3) = {-l,h,0,lc};
Point(4) = {xR,h,0,lc};

pListB[0] = 1;
For i In {1:Nl}
	x = -l + i/Nl*(l+Xc);
	y = 0.0;

	pListB[i] = newp;
	Point(pListB[i]) = {x,y,0,lc};
EndFor

For i In {1:Nc}
	t = 3/2*Pi+Pi*15/180*i/Nc;
	x = Xc+r*Cos(t);
	y = Yc+r*Sin(t);

	pListB[Nl+i] = newp;
	Point(pListB[Nl+i]) = {x,y,0,lc};
EndFor

xL = Xc+r*Cos(3/2*Pi+Pi*15/180);
yL = Yc+r*Sin(3/2*Pi+Pi*15/180);
For i In {1:Nl-1}
	x = xL+(xR-xL)*i/Nl;
	y = yL+(yR-yL)*i/Nl;

	pListB[Nl+Nc+i] = newp;
	Point(pListB[Nl+Nc+i]) = {x,y,0,lc};
EndFor
pListB[2*Nl+Nc-1] = 2;

Spline(1001) = pListB[];
Line(1002) = {3,4};
Line(1003) = {1,3};
Line(1004) = {2,4};

If (RefType == 0)
	Transfinite Line {1001,1002} = 3*2^(Refine)+1 Using Progression 1;
	Transfinite Line {1003,1004} = 1*2^(Refine)+1 Using Progression 1;
ElseIf (RefType == 1)
	Transfinite Line {1001,1002} = 4*2^(Refine)   Using Progression 1;
	Transfinite Line {1003,1004} = 1*2^(Refine)+1 Using Progression 1;
EndIf

Line Loop (4001) = {1001,1004,-1002,-1003};

Plane Surface(4001) = {4001};

If (RefType == 0 || RefType == 1)
	Transfinite Surface{4001};
EndIf

//Recombine Surface{4001};



// Physical Parameters for '.msh' file

Physical Line(10011) = {1002};      // Straight Dirichlet
Physical Line(10012) = {1003,1004}; // Straight Neumann
Physical Line(20011) = {1001};      // Curved Dirichlet
//Physical Line(20012) = {};          // Curved Neumann

Physical Surface(9401) = {4001};



// Visualization in gmsh

Color Black{ Surface{4001:4002}; }
Geometry.Color.Points = Black;