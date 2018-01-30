lx = 4;
ly = 4;
r = 1;

Mesh.SecondOrderIncomplete = 1;

NTrainsfinite = 16;

cos45 = 0.70710678118654757;

Point(1) = {0,0,0,0};
Point(2) = {r,0,0,0};
Point(3) = {lx,0,0,0};
Point(4) = {lx,ly,0,0};
Point(5) = {0, ly, 0, 0};
Point(6) = {0, r, 0,0};
Point(7) = {cos45*r,cos45*r,0,0};

Line(1) = {2,3};
Line(2) = {3,4};
Line(3) = {4,5};
Line(4) = {5,6};
Circle(5) = {6,1,7};
Circle(6) = {7,1,2};
Line(7) = {7, 4};


//+
Line Loop(8) = {6, 1, 2, -7};
//+
Plane Surface(9) = {8};
//+
Line Loop(10) = {5, 7, 3, 4};
//+
Plane Surface(11) = {10};
//+
Transfinite Line {1, 2, 3, 4, 5, 6, 7} = NTrainsfinite;
//+
Transfinite Surface {11};
//+
Transfinite Surface {9};
//+
Recombine Surface {11, 9};
//+
Physical Surface("Plate") = {11, 9};
//+
Physical Line("Top") = {3};
//+
Physical Line("Right") = {2};
//+
Physical Line("Left") = {4};
//+
Physical Line("Bottom") = {1};
