// Domain size
Lx =     2.2000;
Ly =     0.4100;

// Coordinates of the center of the airfoil
Dx =     0.4000;
Dy =     0.2000;

h =     0.0100;

// Points
Point(0) = {0, 0, 0, h};
Point(1) = {Lx, 0, 0, h};
Point(2) = {Lx, Ly, 0, h};
Point(3) = {0, Ly, 0, h};

Point(4) = {Dx +0.24369, Dy -0.05441, 0, h};
Point(5) = {Dx +0.22036, Dy -0.04540, 0, h};
Point(6) = {Dx +0.19691, Dy -0.03662, 0, h};
Point(7) = {Dx +0.14975, Dy -0.01976, 0, h};
Point(8) = {Dx +0.10227, Dy -0.00380, 0, h};
Point(9) = {Dx +0.05447, Dy +0.01123, 0, h};
Point(10) = {Dx +0.00635, Dy +0.02524, 0, h};
Point(11) = {Dx -0.04212, Dy +0.03807, 0, h};
Point(12) = {Dx -0.09132, Dy +0.04916, 0, h};
Point(13) = {Dx -0.11613, Dy +0.05378, 0, h};
Point(14) = {Dx -0.14110, Dy +0.05766, 0, h};
Point(15) = {Dx -0.16622, Dy +0.06066, 0, h};
Point(16) = {Dx -0.19152, Dy +0.06251, 0, h};
Point(17) = {Dx -0.20426, Dy +0.06286, 0, h};
Point(18) = {Dx -0.21705, Dy +0.06261, 0, h};
Point(19) = {Dx -0.23002, Dy +0.06121, 0, h};
Point(20) = {Dx -0.23654, Dy +0.05993, 0, h};
Point(21) = {Dx -0.24359, Dy +0.05480, 0, h};
Point(22) = {Dx -0.23818, Dy +0.04806, 0, h};
Point(23) = {Dx -0.23230, Dy +0.04498, 0, h};
Point(24) = {Dx -0.22036, Dy +0.04036, 0, h};
Point(25) = {Dx -0.20832, Dy +0.03658, 0, h};
Point(26) = {Dx -0.19626, Dy +0.03325, 0, h};
Point(27) = {Dx -0.17206, Dy +0.02732, 0, h};
Point(28) = {Dx -0.14783, Dy +0.02195, 0, h};
Point(29) = {Dx -0.12359, Dy +0.01689, 0, h};
Point(30) = {Dx -0.09934, Dy +0.01199, 0, h};
Point(31) = {Dx -0.05082, Dy +0.00233, 0, h};
Point(32) = {Dx -0.00197, Dy -0.00728, 0, h};
Point(33) = {Dx +0.04697, Dy -0.01672, 0, h};
Point(34) = {Dx +0.09598, Dy -0.02618, 0, h};
Point(35) = {Dx +0.14506, Dy -0.03573, 0, h};
Point(36) = {Dx +0.19421, Dy -0.04541, 0, h};
Point(37) = {Dx +0.21883, Dy -0.05028, 0, h};
Point(38) = {Dx +0.24350, Dy -0.05518, 0, h};


// Lines
Line(0) = {0, 1};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 0};

Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 9};
Line(9) = {9, 10};
Line(10) = {10, 11};
Line(11) = {11, 12};
Line(12) = {12, 13};
Line(13) = {13, 14};
Line(14) = {14, 15};
Line(15) = {15, 16};
Line(16) = {16, 17};
Line(17) = {17, 18};
Line(18) = {18, 19};
Line(19) = {19, 20};
Line(20) = {20, 21};
Line(21) = {21, 22};
Line(22) = {22, 23};
Line(23) = {23, 24};
Line(24) = {24, 25};
Line(25) = {25, 26};
Line(26) = {26, 27};
Line(27) = {27, 28};
Line(28) = {28, 29};
Line(29) = {29, 30};
Line(30) = {30, 31};
Line(31) = {31, 32};
Line(32) = {32, 33};
Line(33) = {33, 34};
Line(34) = {34, 35};
Line(35) = {35, 36};
Line(36) = {36, 37};
Line(37) = {37, 38};
Line(38) = {38, 4};


// Loops
Line Loop(1) = {0, 1, 2, 3};

Line Loop(2) = {4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38};


// Surfaces
Plane Surface(0) = {1, 2};


// Physical entities
Physical Line(0) = {0};
Physical Line(1) = {1};
Physical Line(2) = {2};
Physical Line(3) = {3};
Physical Line(4) = {4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38};
Physical Surface(10) = {0};


Mesh 2;
