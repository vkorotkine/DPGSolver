lc = 0.75;
lc_edges = 0.05;

Point(1) = {0.5, 0.0, 0, lc_edges};
Point(2) = {4.11462305875860e-01, -6.98572814764988e-03, 0, lc};
Point(3) = {2.60549469047750e-01, -4.28144270425459e-02, 0, lc};
Point(4) = {-1.46065162300945e-02, -6.12445446570306e-02, 0, lc};
Point(5) = {-2.22070939017145e-01, -8.77768858230900e-02, 0, lc_edges};
Point(6) = {-5.79487813867564e-01, 2.68678493546853e-17, 0, lc_edges};
Point(7) = {-2.22070939017144e-01, 8.77768858230900e-02, 0, lc_edges};
Point(8) = {-1.46065162300946e-02, 6.12445446570307e-02, 0, lc};
Point(9) = {2.60549469047750e-01, 4.28144270425459e-02, 0, lc};
Point(10) = {4.11462305875859e-01, 6.98572814764991e-03, 0, lc};

// Create the NURBS Airfoil representation
Nurbs(1) = {1,2,3,4,5,6,7,8,9,10,1} Knots {-1, -1, -1, -1, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5, 0.75, 1.0, 1.0, 1.0, 1.0} Order 3;

Point(11) = {0, -20, -0, 1.0};
Point(12) = {0, 0, -0, 1.0};
Point(13) = {0, 20, -0, 1.0};

Circle(2) = {13, 12, 11};
Circle(3) = {11, 12, 13};

Line Loop(1) = {2, 3};
Line Loop(2) = {1};
Plane Surface(1) = {1, 2};

Point(14) = {20, 0, 0, 1.0};
Line(4) = {1, 14};
Line{4} In Surface{1};

// Airfoil Surface
Physical Line(1001) = {1};

// Connecting Line (periodic boundary)
Physical Line(2001) = {4};

// Outer Boundary
Physical Line(3001) = {2, 3};

// Inner domain
Physical Surface(4001) = {1};
