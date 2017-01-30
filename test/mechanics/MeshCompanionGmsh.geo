Point(1) = {-1.3, -0.3, -1, 2.0};
Point(2) = {0.5, -0.9, 0, 2.0};
Point(3) = {0.4, 0, 1, 2.0};
Point(4) = {-1, -1.1, -0, 2.0};
Line(1) = {1, 4};
Line(2) = {4, 2};
Line(3) = {2, 3};
Line(4) = {3, 1};
Line(5) = {1, 2};
Line(6) = {3, 4};
//+
Line Loop(7) = {1, -6, 4};
//+
Plane Surface(8) = {7};
//+
Line Loop(9) = {5, 3, 4};
//+
Plane Surface(10) = {9};
//+
Line Loop(11) = {3, 6, 2};
//+
Plane Surface(12) = {11};
//+
Line Loop(13) = {1, 2, -5};
//+
Plane Surface(14) = {13};
//+
Surface Loop(15) = {10, 14, 8, 12};
//+
Volume(16) = {15};
//+
Physical Volume(17) = {16};
