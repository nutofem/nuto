//mesh using gmsh -2 -order 2 ImportGmsh.geo

//width of the bounding box
w = 100;
//height of the bounding box
h = 100;
//mesh size
lc = 1;

//define points
Point(1) = { 0.0, 0.0 , 0.0 , lc};
Point(2) = { w, 0.0 , 0.0 , lc};
Point(3) = { w, h , 0.0 , lc};
Point(4) = { 0.0, h , 0.0 , lc};

//define ellipsoides
//neight and width
he = 50;
we = 50;

Point(5) = { w/2, h/2 , 0.0 , lc};
Point(6) = { w/2, h/2+he/2 , 0.0 , lc};
Point(7) = { w/2+we/2, h/2 , 0.0 , lc};
Point(8) = { w/2, h/2-he/2 , 0.0 , lc};
Point(9) = { w/2-we/2, h/2 , 0.0 , lc};
Point(10) = { w/2, h/2 , 10.0 , lc};

Line(1)  = {1 ,2};
Line(2)  = {2 ,3};
Line(3)  = {3 ,4};
Line(4)  = {4 ,1};

to_rad = 3.141592654/180;
Rotate { { 0, 0., 1.0}, { w/2, h/2, 0.0}, 0*to_rad}
{
Ellipse(5) = {6,5,7,7};
Ellipse(6) = {7,5,7,8};
Ellipse(7) = {8,5,7,9};
Ellipse(8) = {9,5,7,6};
}

//internal surface
Line Loop(6) = {5,6,7,8};
Plane Surface(1) = {6};

//external surface with holes
Line Loop(7) = {1,2,3,4};
Plane Surface(2) = {7,6};

//physical surface to distinguish the groups
Physical Surface(101) = {1};
Physical Surface(102) = {2};
