Mesh.Algorithm = 6;
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
 
xS = 0; xE = 32;
yS = 0; yE = 16;
meshSpecimen = 0.75;
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


meshCircle = 0.75; 
t = 0;
xC = 24.1899; yC = 11.1245;
R = 2.3458; 
Call MySphere; 
 
 
t = 1;
xC = 12.0583; yC = 7.68818;
R = 1.30966; 
Call MySphere; 
 
 
t = 2;
xC = 3.63097; yC = 10.2585;
R = 2.46499; 
Call MySphere; 
 
 
t = 3;
xC = 9.00742; yC = 12.872;
R = 2.62253; 
Call MySphere; 
 
 
t = 4;
xC = 23.7289; yC = 3.0223;
R = 2.34118; 
Call MySphere; 
 
 
t = 5;
xC = 29.2447; yC = 2.76744;
R = 2.11593; 
Call MySphere; 
 
 
t = 6;
xC = 18.6219; yC = 6.56865;
R = 2.03546; 
Call MySphere; 
 
 
t = 7;
xC = 29.2319; yC = 8.70573;
R = 1.81797; 
Call MySphere; 
 
 
t = 8;
xC = 13.7849; yC = 2.73947;
R = 2.0221; 
Call MySphere; 
 
 
t = 9;
xC = 17.9597; yC = 13.1144;
R = 1.78413; 
Call MySphere; 
 
 
t = 10;
xC = 4.62145; yC = 4.14454;
R = 1.24145; 
Call MySphere; 
 
 
t = 11;
xC = 29.792; yC = 14.0883;
R = 1.09586; 
Call MySphere; 
 
 
t = 12;
xC = 2.25181; yC = 1.81313;
R = 0.860428; 
Call MySphere; 
 
 
t = 13;
xC = 1.58892; yC = 5.53815;
R = 0.920209; 
Call MySphere; 
 
 
t = 14;
xC = 7.09248; yC = 6.69391;
R = 0.991408; 
Call MySphere; 
 
 
t = 15;
xC = 9.76648; yC = 3.95881;
R = 1.02543; 
Call MySphere; 
 
 
volNr = newreg; 
Plane Surface(volNr) = {box, theLoops[]};
Physical Surface(newreg) = volNr;
Physical Surface(newreg) = {theAggregates[]}; 
