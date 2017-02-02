cl__1 = 5;
Point(1) = {0, 0, 0, 5};
Point(2) = {10, 0, 0, 5};
Point(3) = {10, 1, 0, 5};
Point(4) = {0, 1, 0, 5};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Transfinite Line {1:4} = 2;
Line Loop(6) = {1, 2, 3, 4};
Plane Surface(6) = {6};
Transfinite Surface {6};
//Recombine Surface {6};



Extrude {0,0,2} 
{
Surface{6}; 
Layers {1};
}

Transfinite Volume {1};
Physical Volume(999) = {1};