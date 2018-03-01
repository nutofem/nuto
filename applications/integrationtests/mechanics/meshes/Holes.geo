Mesh.Algorithm = 1;
Mesh.Optimize = 1;
Mesh.Smoothing = 2;
Function MySphere 
  meshCircleR = meshCircle; 
  p1 = newp; Point(p1) = {xC,  yC,  0,  meshCircleR} ; 
  p2 = newp; Point(p2) = {xC+R,yC,  0,  meshCircleR} ; 
  p3 = newp; Point(p3) = {xC-R,yC,  0,  meshCircleR} ; 
 
  c1  = newreg; Circle(c1)  = {p2,p1,p3}; c2  = newreg; Circle(c2)  = {p3,p1,p2}; 
 
  l1 = newreg; Line Loop(l1) = {c1,c2};
  s1 = newreg; Plane Surface(s1) = {l1}; 
   
  theLoops[t] = l1; 
  theAggregates[t] = s1; 
 
Return 
 
xS = 0; xE = 16;
yS = 0; yE = 16;
meshSpecimen = 1.5;
// defines a box-shaped specimen 
// by start coordinates <xyz>S 
// and end coordinates  <xyz>E 

// points: 
p0 = newp; Point(p0) = {xS, yS, 0, meshSpecimen}; 
p1 = newp; Point(p1) = {xE, yS, 0, meshSpecimen}; 
p2 = newp; Point(p2) = {xE, yE, 0, meshSpecimen}; 
p3 = newp; Point(p3) = {xS, yE, 0, meshSpecimen}; 

// lines 
l0 = newreg; Line(l0) = {p0, p1}; 
l1 = newreg; Line(l1) = {p1, p2}; 
l2 = newreg; Line(l2) = {p2, p3}; 
l3 = newreg; Line(l3) = {p3, p0}; 

// lineloops and surfaces 
// the index says which coordinate is constant 
box = newreg; Line Loop(box) = { l0, l1, l2, l3}; 


meshCircle = 1.5; 
t = 0;
xC = 5.00986; yC = 8.1089;
R = 2.94444; 
Call MySphere; 
 
 
t = 1;
xC = 12.3876; yC = 8.379;
R = 2.16157; 
Call MySphere; 
 
t = 2;
xC = 8.3876; yC = 4.379;
R = 1.16157; 
Call MySphere; 
 
 
volNr = newreg; 
Plane Surface(volNr) = {box, theLoops[]};
Physical Surface("matrix") = volNr;
Physical Surface("aggregates") = {theAggregates[]};
Physical Line("bottom") = {l0};
Physical Line("right") = {l1};
Physical Line("top") = {l2};
Physical Line("left") = {l3};
