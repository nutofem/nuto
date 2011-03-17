//mesh using gmsh -2 -order 1 ConcurrentMultiscale.geo
Mesh.ChacoSeed = 1;
Mesh.SubdivisionAlgorithm=3;
//width of the bounding box
w = 100;
//height of the bounding box
h = 100;
//mesh size
lc = 5;

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


to_rad = 3.141592654/180;

Function myHole 

   p1 = newp; Point(p1) = { x, y , 0.0 , lc};
   p2 = newp; Point(p2) = { x, y+he/2 , 0.0 , lc};
   p3 = newp; Point(p3) = { x+we/2, y , 0.0 , lc};
   p4 = newp; Point(p4) = { x, y-he/2 , 0.0 , lc};
   p5 = newp; Point(p5) = { x-we/2, y , 0.0 , lc};

   e1 = newreg; 
   e2 = e1+1; 
   e3 = e2+1; 
   e4 = e3+1; 
   Rotate { { 0, 0., 1.0}, { x, y, 0.0}, angle*to_rad}
   {
   Ellipse(e1) = {p2,p1,p3,p3};
   Ellipse(e2) = {p3,p1,p3,p4};
   Ellipse(e3) = {p4,p1,p3,p5};
   Ellipse(e4) = {p5,p1,p3,p2};
   }

   l1 = newreg; Line Loop(l1) = {e1,e2,e3,e4};

   theloops[t] = l1 ; 
   t = t+1;

Return

//generate ellipse
t = 1; 
x = w/2;
y = h/2;
angle = Rand(180.);
he = 8;
we = 8;
Call myHole;

x = w*0.15;
y = h*(0.5-0.1);
angle = Rand(180.);
he = 8;
we = 8;
Call myHole;

x = w*0.85;
y = h*(0.5+0.1);
angle = Rand(180.);
he = 8;
we = 8;
Call myHole;

x = w*0.2;
y = h*(0.2);
angle = Rand(180.);
he = 8;
we = 8;
Call myHole;

x = w*0.65;
y = h*0.7;
angle = Rand(180.);
he = 8;
we = 8;
Call myHole;

x = w*0.7;
y = h*(0.15);
angle = Rand(180.);
he = 8;
we = 8;
Call myHole;

x = w*0.3;
y = h*0.83;
angle = Rand(180.);
he = 8;
we = 8;
Call myHole;


theloops[0] = newreg;
Line Loop(theloops[0]) = {1,2,3,4};

pall = news;
Plane Surface(pall) = {theloops[]} ;

Physical Surface(100) = {pall};

//Line Loop(1) = {1,2,3,4};
//Plane Surface(1) = {1} ;

//Transfinite Surface{1} = {1,2,3,4};
//Recombine Surface {1};

//Physical Surface(100) = {1};
