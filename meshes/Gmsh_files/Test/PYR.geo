// Modifiable Parameters
lc = 1; // Not used.

L = 1;
H = 1;
W = 1;



// Geometry Specification

Point(1) = {-L,-H,-W,lc};
Point(2) = { 0,-H,-W,lc};
Point(3) = {-L, 0,-W,lc};
Point(4) = { 0, 0,-W,lc};

Point(5) = {-L,-H, 0,lc};
Point(6) = { 0,-H, 0,lc};
Point(7) = {-L, 0, 0,lc};
Point(8) = { 0, 0, 0,lc};

Line(1001) = {1,2};
Line(1002) = {3,4};
Line(1003) = {1+4,2+4};
Line(1004) = {3+4,4+4};

Line(2001) = {1,3};
Line(2002) = {2,4};
Line(2003) = {1+4,3+4};
Line(2004) = {2+4,4+4};

Line(3001) = {1,1+4};
Line(3002) = {2,2+4};
Line(3003) = {3,3+4};
Line(3004) = {4,4+4};

Transfinite Line{1001:1004} = 2 Using Progression 1;
Transfinite Line{2001:2004} = 2 Using Progression 1;
Transfinite Line{3001:3004} = 2 Using Progression 1;

Line Loop (4001) = {1001,2002,-1002,-2001};
Line Loop (4002) = {1003,2004,-1004,-2003};
Line Loop (5001) = {1001,3002,-1003,-3001};
Line Loop (5002) = {1002,3004,-1004,-3003};
Line Loop (6001) = {2001,3003,-2003,-3001};
Line Loop (6002) = {2002,3004,-2004,-3002};

Plane Surface(4001) = {4001};
Plane Surface(4002) = {4002};
Plane Surface(5001) = {5001};
Plane Surface(5002) = {5002};
Plane Surface(6001) = {6001};
Plane Surface(6002) = {6002};

Transfinite Surface {4001:4002,5001:5002,6001:6002};

Surface Loop (7001) = {4001:4002,5001:5002,6001:6002};

Volume(7001) = {7001};


// Form 8 octants for the initial mesh
Symmetry{ -1.0,0.0,0.0,0.0 }{Duplicata{Volume{7001};}}

Transfinite Line{7004,7006,7011,7009,7007,7022,7012,7017} = 2 Using Progression 1;
Transfinite Surface {7003,7008,7013,7018,7023};


Symmetry{ 0.0,-1.0,0.0,0.0 }{Duplicata{Volume{7001,7002};}}

Transfinite Line{7034,7029,7032,7027,7065,7060} = 2 Using Progression 1;
Transfinite Line{7031,7026,7039,7037,7062,7057,7070} = 2 Using Progression 1;
Transfinite Surface {7025,7030,7045,7050,7035,7076,7061,7056,7066};


Symmetry{ 0.0,0.0,-1.0,0.0 }{Duplicata{Volume{7001:7002,7024,7055};}}
Transfinite Line{7154,7152,7185,7097,7095,7128,7092,7090,7123,7141,7172,7144,7142,7175,7081,7112,7082,7080,7113,7079,7110} = 2 Using Progression 1;
Transfinite Surface {7150,7181,7160,7165,7191,7093,7124,7098,7103,
                     7129,7088,7119,7140,7171,7078,7109};


Recombine   Surface {4001,7003,7025,7056,5001,7013,7088,7119,7023,7076,7129,7191};


// Physical parameters for '.msh' file

Physical Surface(10001) = {4001,5001,6001,7003,7013,7023,7025,7056,7045,7035,7066,7076,
                           7140,7171,7078,7109,7150,7181,7160,7191,7098,7129,7088,7119};

Physical Volume(9701) = {7001:7002,7024,7055,7077,7108,7139,7170};



// Visualization in gmsh

Color Black{ Surface{4001:4002,5001:5002,6001:6002}; }
Color Black{ Volume{7001:7002,7024,7055,7077,7108,7139,7170}; }
Geometry.Color.Points = Black;