// Gmsh project created on Tue Mar 14 16:46:07 2023
SetFactory("OpenCASCADE");
R = 1.0;
R1 = 0.10;
R2 = 0.05;
lc = R2/10;
Lc = 0.1;

Point(1) = { 0.0, 0.0, 0.0, lc};
Point(2) = { R1 * Cos(0), R2 * Sin(0), 0.0, lc };
Point(3) = { R1 * Cos(Pi/2), R2 * Sin(Pi/2), 0.0, lc };
Point(4) = { R1 * Cos(Pi), R2 * Sin(Pi), 0.0, lc };
Point(6) = { R1 * Cos(Pi/6), R2 * Sin(Pi/6), 0.0, lc };
Point(7) = { R1 * Cos(5*Pi/6), R2 * Sin(5*Pi/6), 0.0, lc };
Point(5) = { R1 * Cos(-Pi/2), R2 * Sin(-Pi/2), 0.0, lc };
Point(8) = { R1 * Cos(-5*Pi/6), R2 * Sin(-5*Pi/6), 0.0, lc };
Point(9) = { R1 * Cos(-Pi/6), R2 * Sin(-Pi/6), 0.0, lc };

Ellipse(1) = { 2, 1, 2, 6 };
Ellipse(2) = { 6, 1, 2, 3 };
Ellipse(3) = { 3, 1, 2, 7 };
Ellipse(4) = { 4, 1, 2, 7 };
Ellipse(5) = { 4, 1, 2, 8 };
Ellipse(6) = { 8, 1, 2, 5 };
Ellipse(7) = { 5, 1, 2, 9 };
Ellipse(8) = { 9, 1, 2, 2 };

Curve Loop(20) = {1,2,3,4, 5, 6, 7, 8}; 


//+
Point(10) = {R, 0, 0, Lc};
//+
Point(11) = {0, R, 0, Lc};
//+
Point(12) = {0, -R, 0, Lc};
//+
Point(13) = {-R, 0, 0, Lc};

Circle(9) = {10, 1, 11};
Circle(10) = {11, 1, 12};
Circle(11) = {12, 1, 13};
Circle(12) = {13, 1, 10};


Curve Loop(21) = {9, 10, 11, 12};

Plane Surface(22) = {21, 20}; 
