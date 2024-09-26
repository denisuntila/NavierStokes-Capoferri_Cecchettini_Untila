// Domain size
Lx = 2.2;		// Domain width
Ly = 0.41;		// Domain height

// Coordinates of the center of the circle
Dx = 0.2;
Dy = 0.2;

// Radius of the circle
r = 0.05;

h = 0.01;

Point(1) = {0, 0, 0, h};
Point(2) = {Lx, 0, 0, h};
Point(3) = {Lx, Ly, 0, h};
Point(4) = {0, Ly, 0, h};

Point(5) = {Dx, Dy, 0, h};
Point(6) = {Dx - r, Dy, 0, h};
Point(7) = {Dx + r, Dy, 0, h};


Line(8) = {1, 2};
Line(9) = {2, 3};
Line(10) = {3, 4};
Line(11) = {4, 1};


Line Loop(1) = {8, 9, 10, 11};

Circle(12) = {6, 5, 7};
Circle(13) = {7, 5, 6};

Curve Loop(2) = {12, 13};

Plane Surface(1) = {1, 2};

Physical Line(0) = {8};
Physical Line(1) = {9};
Physical Line(2) = {10};
Physical Line(3) = {11};
Physical Curve(4) = {12, 13};
Physical Surface(10) = {1};

Mesh 2;




