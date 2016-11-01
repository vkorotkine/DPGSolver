Mesh.Algorithm = 6;

Refine = 0;

lc = .55/2.0^Refine;

Point(1) = {0.0,0.0,0.0,lc};
Point(2) = {1,0.0,0.0,lc};
Point(3) = {0,1,0.0,lc};
Point(4) = {0,0,1,lc};

Circle(1001) = {2,1,3};
Circle(1002) = {4,1,3};
Circle(1003) = {2,1,4};

Line(1004) = {1,2};
Line(1005) = {1,3};
Line(1006) = {1,4};

Line Loop(4008) = {1004,1003,-1006};
Line Loop(4010) = {1005,-1002,-1006};
Line Loop(4012) = {1004,1001,-1005};
Line Loop(4017) = {-1002,-1003,1001};

Plane Surface(4009) = {4008};
Plane Surface(4011) = {4010};
Plane Surface(4013) = {4012};
Ruled Surface(4018) = {4017};

Surface Loop(7029) = {4009,4011,4013,4018};
Volume(7030) = {7029};

Physical Surface(10011) = {4009,4011};
Physical Surface(10012) = {4013};
Physical Surface(20012) = {4018};

Physical Volume(9701) = {7030};
