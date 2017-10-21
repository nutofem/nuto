Mesh.Algorithm = 6;
Mesh.HighOrderOptimize = 1;
Mesh.Optimize = 2;
Mesh.Smoothing = 2;

Function MySphere 
  meshCircleR = meshCircle; 
  p1 = newp; Point(p1) = {xC,  yC,  zC,  meshCircleR} ; 
  p2 = newp; Point(p2) = {xC+R,yC,  zC,  meshCircleR} ; 
  p3 = newp; Point(p3) = {xC,  yC+R,zC,  meshCircleR} ; 
  p4 = newp; Point(p4) = {xC,  yC,  zC+R,meshCircleR} ; 
  p5 = newp; Point(p5) = {xC-R,yC,  zC,  meshCircleR} ; 
  p6 = newp; Point(p6) = {xC,  yC-R,zC,  meshCircleR} ; 
  p7 = newp; Point(p7) = {xC,  yC,  zC-R,meshCircleR} ; 
 
  c1  = newreg; Circle(c1)  = {p2,p1,p7}; c2  = newreg; Circle(c2)  = {p7,p1,p5}; 
  c3  = newreg; Circle(c3)  = {p5,p1,p4}; c4  = newreg; Circle(c4)  = {p4,p1,p2}; 
  c5  = newreg; Circle(c5)  = {p2,p1,p3}; c6  = newreg; Circle(c6)  = {p3,p1,p5}; 
  c7  = newreg; Circle(c7)  = {p5,p1,p6}; c8  = newreg; Circle(c8)  = {p6,p1,p2}; 
  c9  = newreg; Circle(c9)  = {p7,p1,p3}; c10 = newreg; Circle(c10) = {p3,p1,p4}; 
  c11 = newreg; Circle(c11) = {p4,p1,p6}; c12 = newreg; Circle(c12) = {p6,p1,p7}; 
 
  l1 = newreg; Line Loop(l1) = {c5,c10,c4};   Ruled Surface(newreg) = {l1}; 
  l2 = newreg; Line Loop(l2) = {c9,-c5,c1};   Ruled Surface(newreg) = {l2}; 
  l3 = newreg; Line Loop(l3) = {c12,-c8,-c1}; Ruled Surface(newreg) = {l3}; 
  l4 = newreg; Line Loop(l4) = {c8,-c4,c11};  Ruled Surface(newreg) = {l4}; 
  l5 = newreg; Line Loop(l5) = {-c10,c6,c3};  Ruled Surface(newreg) = {l5}; 
  l6 = newreg; Line Loop(l6) = {-c11,-c3,c7}; Ruled Surface(newreg) = {l6}; 
  l7 = newreg; Line Loop(l7) = {-c2,-c7,-c12};Ruled Surface(newreg) = {l7}; 
  l8 = newreg; Line Loop(l8) = {-c6,-c9,c2};  Ruled Surface(newreg) = {l8}; 
   
  theLoops[t] = newreg; 
 
  Surface Loop(theLoops[t]) = {l8+1,l5+1,l1+1,l2+1,l3+1,l7+1,l6+1,l4+1}; 
 
  thehole = newreg; 
  Volume(thehole) = theLoops[t]; 
  theAggregates[t] = thehole; 
 
Return 
 
xS = 0; xE = 24;
yS = 0; yE = 24;
zS = 0; zE = 24;
meshSpecimen = 4; 
// defines a box-shaped specimen 
// by start coordinates <xyz>S 
// and end coordinates  <xyz>E 

// top points: z = zE 
p0 = newp; Point(p0) = {xS, yS, zE, meshSpecimen}; 
p1 = newp; Point(p1) = {xS, yE, zE, meshSpecimen}; 
p2 = newp; Point(p2) = {xE, yE, zE, meshSpecimen}; 
p3 = newp; Point(p3) = {xE, yS, zE, meshSpecimen}; 
// bottom points z = zS 
p4 = newp; Point(p4) = {xS, yS, zS, meshSpecimen}; 
p5 = newp; Point(p5) = {xS, yE, zS, meshSpecimen}; 
p6 = newp; Point(p6) = {xE, yE, zS, meshSpecimen}; 
p7 = newp; Point(p7) = {xE, yS, zS, meshSpecimen}; 

// top lines z = zS 
lT0 = newreg; Line(lT0) = {p0, p1}; 
lT1 = newreg; Line(lT1) = {p1, p2}; 
lT2 = newreg; Line(lT2) = {p2, p3}; 
lT3 = newreg; Line(lT3) = {p3, p0}; 
// bottom lines z = zE 
lB0 = newreg; Line(lB0) = {p4, p5}; 
lB1 = newreg; Line(lB1) = {p5, p6}; 
lB2 = newreg; Line(lB2) = {p6, p7}; 
lB3 = newreg; Line(lB3) = {p7, p4}; 
// connection zS --> zE 
lC0 = newreg; Line(lC0) = {p0, p4}; 
lC1 = newreg; Line(lC1) = {p1, p5}; 
lC2 = newreg; Line(lC2) = {p2, p6}; 
lC3 = newreg; Line(lC3) = {p3, p7}; 

// lineloops and surfaces 
// the index says which coordinate is constant 
lxS = newreg; Line Loop(lxS) = {-lT3, lC3, lB3,-lC0}; Plane Surface(newreg) = {lxS}; 
lxE = newreg; Line Loop(lxE) = {-lT1, lC1, lB1,-lC2}; Plane Surface(newreg) = {lxE}; 

lyS = newreg; Line Loop(lyS) = {-lT0, lC0, lB0,-lC1}; Plane Surface(newreg) = {lyS}; 
lyE = newreg; Line Loop(lyE) = {-lT2, lC2, lB2,-lC3}; Plane Surface(newreg) = {lyE}; 

lzS = newreg; Line Loop(lzS) = { lT0, lT1, lT2, lT3}; Plane Surface(newreg) = {lzS}; 
lzE = newreg; Line Loop(lzE) = {-lB3,-lB2,-lB1,-lB0}; Plane Surface(newreg) = {lzE}; 

theLoops[0] = newreg; 
Surface Loop(theLoops[0]) = {lxS+1, lyS+1, lxE+1, lyE+1, lzS+1, lzE+1}; 
theBox = newreg; 
Volume(theBox) = theLoops[0]; 


meshCircle = 4; 
t = 1;
xC = 15.9936; yC = 8.9257; zC = 15.9935;
R = 5.87829; 
Call MySphere; 
 
 
t = 2;
xC = 10.5554; yC = 17.7632; zC = 6.23681;
R = 4.10868; 
Call MySphere; 
 
 
t = 3;
xC = 6.16492; yC = 6.16492; zC = 6.16496;
R = 4.03682; 
Call MySphere; 
 
 
volNr = newreg; 
Volume(volNr) = {theLoops[]};
Physical Volume(newreg) = volNr;
Physical Volume(newreg) = {theAggregates[]};
