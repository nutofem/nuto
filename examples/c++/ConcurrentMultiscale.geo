//mesh using gmsh -2 -order 1 ConcurrentMultiscale.geo
Mesh.ChacoSeed = 2;
Mesh.SubdivisionAlgorithm=2;
//width of the bounding box
w = 100;
//height of the bounding box
h = 100;
//mesh size
lc = 100;

//define points
Point(1) = { 0.0, 0.0 , 0.0 , lc};
Point(2) = { w, 0.0 , 0.0 , lc};
Point(3) = { w, h , 0.0 , lc};
Point(4) = { 0.0, h , 0.0 , lc};

Line(1)  = {1 ,2};
Line(2)  = {2 ,3};
Line(3)  = {3 ,4};
Line(4)  = {4 ,1};
Transfinite Line{1,3} = w/lc+1;
Transfinite Line{2,4} = h/lc+1;

Line Loop(1) = {1,2,3,4};
Plane Surface(1) = {1} ;
Transfinite Surface{1} = {1,2,3,4};
Recombine Surface {1};

//Extrude {0,0,0.15} {Surface{1}; Layers{1}; Recombine; }
Physical Surface(100) = {1};
