// Domain sizes
Lx = 2.50;
Ly = 0.41;
Lz = 0.41;

Dx = 0.45;
Dy = 0.20;

r = 0.05;


h = 0.05;



// Points
Point(1) = {0, 0, 0, h};
Point(2) = {Lx, 0, 0, h};
Point(3) = {Lx, Ly, 0, h};
Point(4) = {0, Ly, 0, h};

Point(5) = {Dx, Dy, 0, h};
Point(6) = {Dx - r, Dy, 0, h};
Point(7) = {Dx + r, Dy, 0, h};

Point(8) = {0, 0, Lz, h};
Point(9) = {Lx, 0, Lz, h};
Point(10) = {Lx, Ly, Lz, h};
Point(11) = {0, Ly, Lz, h};

Point(12) = {Dx, Dy, Lz, h};
Point(13) = {Dx - r, Dy, Lz, h};
Point(14) = {Dx + r, Dy, Lz, h};


// Lines
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};


Circle(5) = {6, 5, 7};
Circle(6) = {7, 5, 6};


Line(7) = {8, 9};
Line(8) = {9, 10};
Line(9) = {10, 11};
Line(10) = {11, 8};


Circle(11) = {13, 12, 14};
Circle(12) = {14, 12, 13};


// Links between the two surfaces
Line(13) = {6, 13};
Line(14) = {7, 14};

Line(15) = {1, 8};
Line(16) = {2, 9};
Line(17) = {3, 10};
Line(18) = {4, 11};


Line Loop(1) = {1, 2, 3, 4};
Line Loop(2) = {5, 6};
Line Loop(3) = {7, 8, 9, 10};
Line Loop(4) = {11, 12};

Line Loop(5) = {5, 14, -11, -13};
Line Loop(6) = {6, 13, -12, -14};

Line Loop(7) = {15, -10, -18, 4};
Line Loop(8) = {16, 8, -17, -2};
Line Loop(9) = {1, 16, -7, -15};
Line Loop(10) = {3, 18, -9, -17};

// Surfaces
Surface(1) = {5};	// Half cylinder
Surface(2) = {6};

Plane Surface(3) = {1, 2};	// Lateral faces
Plane Surface(4) = {3, 4};

Plane Surface(5) = {7}; 	// Inlet Surface
Plane Surface(6) = {8}; 	// Outlet Surface
Plane Surface(7) = {9};		// Down surface
Plane Surface(8) = {10};	// Up Surface


// Volume
Surface Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8};
Volume(1) = {1};


// Physical entities
Physical Surface(0) = {3, 4};
Physical Surface(1) = {6};
Physical Surface(2) = {7, 8};
Physical Surface(3) = {5};
Physical Surface(4) = {1, 2};
Physical Volume(10) = {1};


Mesh 3;


