//mesh using gmsh -2 -order 2 ConcurrentMultiscale.geo
Mesh.ChacoSeed = 2;
//width of the bounding box
w = 100;
//height of the bounding box
h = 100;
//mesh size
lc = 10;

//define points
Point(1) = { 0.0, 0.0 , 0.0 , lc};
Point(2) = { w, 0.0 , 0.0 , lc};
Point(3) = { w, h , 0.0 , lc};
Point(4) = { 0.0, h , 0.0 , lc};

//define ellipsoides
//height and width

numEllipse1D = 1;

he = 12;
we = 12;

variationDiameter = 0.2;
variationPosition = 1;

deltaPos = ((w-numEllipse1D*we)/numEllipse1D)+we;

Line(1)  = {1 ,2};
Line(2)  = {2 ,3};
Line(3)  = {3 ,4};
Line(4)  = {4 ,1};

to_rad = 3.141592654/180;


Function myHole 

   p1 = newp; Point(p1) = { x_rand, y_rand , 0.0 , lc};
   p2 = newp; Point(p2) = { x_rand, y_rand+he_rand/2 , 0.0 , lc};
   p3 = newp; Point(p3) = { x_rand+we_rand/2, y_rand , 0.0 , lc};
   p4 = newp; Point(p4) = { x_rand, y_rand-he_rand/2 , 0.0 , lc};
   p5 = newp; Point(p5) = { x_rand-we_rand/2, y_rand , 0.0 , lc};

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

   //p1 = news; Plane Surface(p1) = {l1};

   theloops[t] = l1 ; 
   Printf("Hole %g (center = {%g,%g}, width = %g height = %g angle = %g)",
	 t, x_rand, y_rand, we_rand, he_rand, angle) ;

Return


// We can use a `For' loop to generate five holes in the cube:

x = deltaPos/2 ; 

t = 0;
For tx In {1:numEllipse1D}
y = deltaPos/2 ;
For ty In {1:numEllipse1D}

  t += 1; 

  he_rand = he*(1.-variationDiameter*0.5+Rand(variationDiameter));
  we_rand = we*(1.-variationDiameter*0.5+Rand(variationDiameter));
  
  x_rand = x-variationPosition*(0.5+Rand(1.));
  y_rand = y-variationPosition*(0.5+Rand(1.));
   
  angle = Rand(180.);

  Call myHole ;

  y += deltaPos ; 
EndFor
  x += deltaPos ; 
EndFor

x = deltaPos ; 
For tx In {1:numEllipse1D-1}
y = deltaPos ;
For ty In {1:numEllipse1D-1}

  t += 1; 

  he_rand = he*(1.-variationDiameter*0.5+Rand(variationDiameter));
  we_rand = we*(1.-variationDiameter*0.5+Rand(variationDiameter));
  
  x_rand = x-variationPosition*(0.5+Rand(1.));
  y_rand = y-variationPosition*(0.5+Rand(1.));
   
  angle = Rand(180.);

  Call myHole ;

  y += deltaPos ; 
EndFor
  x += deltaPos ; 
EndFor


// We can then define the surface loop for the exterior surface of the
// cube:

//external surface with holes
theloops[0] = newreg;
Line Loop(theloops[0]) = {1,2,3,4};

pall = news;
Plane Surface(pall) = {theloops[]} ;

// We finally define a physical volume for the elements discretizing
// the cube, without the holes (whose elements were already tagged
// with numbers 1 to 5 in the `For' loop):

//Physical Surface(101) = {};
Physical Surface(102) = {pall};
