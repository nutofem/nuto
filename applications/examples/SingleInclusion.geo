Mesh.Algorithm = 6;
Mesh.Optimize = 1;
Mesh.Smoothing = 2;
Mesh.SecondOrderIncomplete = 1;

Function MySphere 
  meshCircleR = meshCircle; 
  p1 = newp; Point(p1) = {xC,  yC,  0,  meshCircleR} ; 
  p2 = newp; Point(p2) = {xC+R,yC,  0,  meshCircleR} ; 
  p3 = newp; Point(p3) = {xC-R,yC,  0,  meshCircleR} ; 
  
  p2t = newp; Point(p2t) = {xC+R-thickness,yC,  0,  meshCircleR} ; 
  p3t = newp; Point(p3t) = {xC-R+thickness,yC,  0,  meshCircleR} ; 
 
  c1  = newreg; Circle(c1)  = {p2, p1,p3};  c2  = newreg; Circle(c2)  = {p3,p1,p2}; 
  c1t = newreg; Circle(c1t) = {p2t,p1,p3t}; c2t = newreg; Circle(c2t) = {p3t,p1,p2t}; 
  l2 = newreg; Line(l2) = {p2, p2t}; l3 = newreg; Line(l3) = {p3t, p3}; 
 
  llOuter  = newreg; Line Loop(llOuter) = {c1,c2};
  llInner  = newreg; Line Loop(llInner) = {c1t,c2t};

  ll2 = newreg; Line Loop(ll2) = {c1, -l3, - c1t, -l2};
  ll3 = newreg; Line Loop(ll3) = {c2, l2, -c2t, l3};

  sInner  = newreg; Plane Surface(sInner) = {llInner}; 
  
  sInterface2  = newreg; Plane Surface(sInterface2) = {ll2};
  sInterface3  = newreg; Plane Surface(sInterface3) = {ll3};

  numTrans = 3.1415 * R / meshCircleR;
  If (numTrans < 4)
    numTrans = 4;
  EndIf
  Printf("%g", numTrans);
  Transfinite Line {c1, c1t, c2, c2t} = numTrans Using Progression 1;
  Transfinite Line {l2, l3} = 1 Using Progression 1;

  Transfinite Surface {sInterface2};
  Transfinite Surface {sInterface3};
  Recombine Surface {sInterface2};
  Recombine Surface {sInterface3};

  theInterfaces[2*t] = sInterface2;
  theInterfaces[2*t+1] = sInterface3;
  theLoops[t] = llOuter; 
  theAggregates[t] = sInner; 
 
Return 
 
xS = 0; xE = 60;
yS = 0; yE = 60;
meshSpecimen = 2;
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

thickness = 0.1;

meshCircle = 2; 

t = 0;
xC = 30; yC = 20;
R = 10; 
Call MySphere; 
 
volNr = newreg; 
Plane Surface(volNr) = {box, theLoops[]};
Physical Surface("Matrix") = volNr;
Physical Surface("Aggregates") = {theAggregates[]}; 
Physical Surface("Interfaces") = {theInterfaces[]}; 

