// Domain sizes
Lx = 2.50;
Ly = 0.41;
Lz = 0.41;

Dx = 0.45;
Dy = 0.15;

S = 0.1;


h = 0.2;


Point(1) = {0, 0, 0, h};
Point(2) = {Lx, 0, 0, h};
Point(3) = {Lx, Ly, 0, h};
Point(4) = {0, Ly, 0, h};

Point(5) = {Dx, Dy, 0, h};
Point(6) = {Dx + S, Dy, 0, h};
Point(7) = {Dx + S, Dy + S, 0, h};
Point(8) = {Dx, Dy + S, 0, h};


Point(21) = {0, 0, Lz, h};
Point(22) = {Lx, 0, Lz, h};
Point(23) = {Lx, Ly, Lz, h};
Point(24) = {0, Ly, Lz, h};

Point(25) = {Dx, Dy, Lz, h};
Point(26) = {Dx + S, Dy, Lz, h};
Point(27) = {Dx + S, Dy + S, Lz, h};
Point(28) = {Dx, Dy + S, Lz, h};


Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};


Line(21) = {21, 22};
Line(22) = {22, 23};
Line(23) = {23, 24};
Line(24) = {24, 21};

Line(25) = {25, 26};
Line(26) = {26, 27};
Line(27) = {27, 28};
Line(28) = {28, 25};


Line(31) = {1, 21};
Line(32) = {2, 22};
Line(33) = {3, 23};
Line(34) = {4, 24};

Line(35) = {5, 25};
Line(36) = {6, 26};
Line(37) = {7, 27};
Line(38) = {8, 28};


Curve Loop(1) = {1, 2, 3, 4};
Curve Loop(2) = {5, 6, 7, 8};


Curve Loop(3) = {21, 22, 23, 24};
Curve Loop(4) = {25, 26, 27, 28};


Curve Loop(11) = {4, 31, -24, -34};
Curve Loop(12) = {8, 35, -28, -38};
Curve Loop(13) = {1, 32, -21, -31};
Curve Loop(14) = {5, 36, -25, -35};
Curve Loop(15) = {2, 33, -22, -32};
Curve Loop(16) = {6, 37, -26, -36};
Curve Loop(17) = {3, 34, -23, -33};
Curve Loop(18) = {7, 38, -27, -37};


Plane Surface(1) = {11};
Plane Surface(2) = {12};
Plane Surface(3) = {13};
Plane Surface(4) = {14};
Plane Surface(5) = {15};
Plane Surface(6) = {16};
Plane Surface(7) = {17};
Plane Surface(8) = {18};


Plane Surface(9) = {1, 2};
Plane Surface(10) = {3, 4};

Surface Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
Volume(1) = {1};

Physical Surface(0) = {9, 10};
Physical Surface(1) = {5};
Physical Surface(2) = {3, 7};
Physical Surface(3) = {1};
Physical Surface(4) = {2, 4, 6, 8};
Physical Volume(10) = {1};



// e() = Extrude{0, 0, Lz}{ Surface{1}; };

// Physical Volume(1) = {e(1)};






