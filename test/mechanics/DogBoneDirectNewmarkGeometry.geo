//mesh using gmsh -2 -order 2 CoarseScale.geo
Mesh.ChacoSeed = 1;
Mesh.SubdivisionAlgorithm=3;
//D for the dog bone specimen
D = 0.05;
//mesh size
lc = 0.01;
//angle for control points
phi=0.01;

//define points
Point(1)  = { 0. , 1.5*D , 0.0 , lc};
Point(2)  = { 0. , 1.25*D, 0.0 , lc};
Point(3)  = { 0. -0.525*D+Cos(phi)*0.725*D, 0.75*D+Sin(phi)*0.725*D, 0.0 , lc};
Point(4)  = { 0. -0.525*D+Cos(phi)*0.725*D, 0.75*D-Sin(phi)*0.725*D, 0.0 , lc};
Point(5)  = { 0. , 0.25*D, 0.0 , lc};
Point(6)  = { 0. , 0, 0.0 , lc};
Point(7)  = { 0. -0.525*D , 0.75*D, 0.0 , lc};

Point(8)  = { D , 1.5*D , 0.0 , lc};
Point(9)  = { D , 1.25*D, 0.0 , lc};
Point(10)  = { D+0.525*D-Cos(phi)*0.725*D, 0.75*D+Sin(phi)*0.725*D, 0.0 , lc};
Point(11)  = { D+0.525*D-Cos(phi)*0.725*D, 0.75*D-Sin(phi)*0.725*D, 0.0 , lc};
Point(12)  = { D , 0.25*D, 0.0 , lc};
Point(13)  = { D , 0, 0.0 , lc};
Point(14)  = { D +0.525*D , 0.75*D, 0.0 , lc};


Line(1)  = {1,2};
Circle(2)  = {3,7,2};
Circle(3)  = {4,7,3};
Circle(4)  = {5,7,4};
Line(5)  = {5,6};
Line(6)  = {6,13};
Line(7)  = {13,12};
Circle(8)  = {12,14,11};
Circle(9)  = {11,14,10};
Circle(10) = {10,14,9};
Line(11) = {9,8};
Line(12) = {8,1}; 

Line Loop(1) = {1,-2,-3,-4,5,6,7,8,9,10,11,12};
Plane Surface(1) = {1} ;

Physical Surface(100) = {1};
