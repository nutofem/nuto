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
 
xS = 0; xE = 100;
yS = 0; yE = 20;
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


t = 0;
xC = 39.1614; yC = 12.4495;
R = 3.54062; 
meshCircle = 0.75; 
Call MySphere; 
 
 
t = 1;
xC = 67.4044; yC = 8.89679;
R = 4.50113; 
meshCircle = 0.75; 
Call MySphere; 
 
 
t = 2;
xC = 86.6112; yC = 6.51903;
R = 4.38525; 
meshCircle = 0.75; 
Call MySphere; 
 
 
t = 3;
xC = 73.8829; yC = 15.1759;
R = 3.15768; 
meshCircle = 0.75; 
Call MySphere; 
 
 
t = 4;
xC = 23.2516; yC = 11.8574;
R = 4.24237; 
meshCircle = 0.75; 
Call MySphere; 
 
 
t = 5;
xC = 7.81364; yC = 11.7894;
R = 4.0292; 
meshCircle = 0.75; 
Call MySphere; 
 
 
t = 6;
xC = 55.4486; yC = 4.2353;
R = 3.50636; 
meshCircle = 0.75; 
Call MySphere; 
 
 
t = 7;
xC = 4.31764; yC = 4.23527;
R = 2.82626; 
meshCircle = 0.75; 
Call MySphere; 
 
 
t = 8;
xC = 78.25; yC = 4.23578;
R = 2.48265; 
meshCircle = 0.75; 
Call MySphere; 
 
 
t = 9;
xC = 53.8921; yC = 13.2986;
R = 1.79448; 
meshCircle = 0.75; 
Call MySphere; 
 
 
t = 10;
xC = 96.1336; yC = 3.98436;
R = 2.98526; 
meshCircle = 0.75; 
Call MySphere; 
 
 
t = 11;
xC = 29.8907; yC = 6.97704;
R = 2.78934; 
meshCircle = 0.75; 
Call MySphere; 
 
 
t = 12;
xC = 59.1212; yC = 13.4204;
R = 2.02667; 
meshCircle = 0.75; 
Call MySphere; 
 
 
t = 13;
xC = 47.6001; yC = 16.4655;
R = 2.83093; 
meshCircle = 0.75; 
Call MySphere; 
 
 
t = 14;
xC = 81.3938; yC = 16.6574;
R = 2.62123; 
meshCircle = 0.75; 
Call MySphere; 
 
 
t = 15;
xC = 73.6194; yC = 2.95672;
R = 0.600524; 
meshCircle = 0.75; 
Call MySphere; 
 
 
t = 16;
xC = 16.1691; yC = 13.6074;
R = 1.85086; 
meshCircle = 0.75; 
Call MySphere; 
 
 
t = 17;
xC = 65.1391; yC = 16.7073;
R = 1.79685; 
meshCircle = 0.75; 
Call MySphere; 
 
 
t = 18;
xC = 90.1312; yC = 17.3604;
R = 2.12893; 
meshCircle = 0.75; 
Call MySphere; 
 
 
t = 19;
xC = 11.6531; yC = 17.5443;
R = 1.11057; 
meshCircle = 0.75; 
Call MySphere; 
 
 
t = 20;
xC = 94.9391; yC = 17.4841;
R = 1.9041; 
meshCircle = 0.75; 
Call MySphere; 
 
 
t = 21;
xC = 22.0484; yC = 2.40833;
R = 0.431616; 
meshCircle = 0.75; 
Call MySphere; 
 
 
t = 22;
xC = 19.402; yC = 17.6549;
R = 1.8427; 
meshCircle = 0.75; 
Call MySphere; 
 
 
t = 23;
xC = 35.4389; yC = 6.38053;
R = 1.82835; 
meshCircle = 0.75; 
Call MySphere; 
 
 
t = 24;
xC = 75.7422; yC = 9.92858;
R = 1.18538; 
meshCircle = 0.75; 
Call MySphere; 
 
 
t = 25;
xC = 88.6594; yC = 13.0737;
R = 1.56315; 
meshCircle = 0.75; 
Call MySphere; 
 
 
t = 26;
xC = 15.3406; yC = 17.926;
R = 0.692373; 
meshCircle = 0.75; 
Call MySphere; 
 
 
t = 27;
xC = 29.5654; yC = 17.9509;
R = 1.45124; 
meshCircle = 0.75; 
Call MySphere; 
 
 
t = 28;
xC = 69.3191; yC = 1.87099;
R = 1.37151; 
meshCircle = 0.75; 
Call MySphere; 
 
 
t = 29;
xC = 48.3148; yC = 11.448;
R = 0.276267; 
meshCircle = 0.75; 
Call MySphere; 
 
 
t = 30;
xC = 56.5345; yC = 17.961;
R = 0.739889; 
meshCircle = 0.75; 
Call MySphere; 
 
 
t = 31;
xC = 50.0321; yC = 7.03844;
R = 0.968781; 
meshCircle = 0.75; 
Call MySphere; 
 
 
t = 32;
xC = 52.3766; yC = 18.1175;
R = 1.25176; 
meshCircle = 0.75; 
Call MySphere; 
 
 
t = 33;
xC = 97.6537; yC = 14.4607;
R = 1.21248; 
meshCircle = 0.75; 
Call MySphere; 
 
 
t = 34;
xC = 44.7213; yC = 10.6616;
R = 0.805416; 
meshCircle = 0.75; 
Call MySphere; 
 
 
t = 35;
xC = 60.4234; yC = 18.2315;
R = 1.13282; 
meshCircle = 0.75; 
Call MySphere; 
 
 
t = 36;
xC = 74.0167; yC = 7.32936;
R = 0.826813; 
meshCircle = 0.75; 
Call MySphere; 
 
 
t = 37;
xC = 42.6387; yC = 5.47365;
R = 0.908775; 
meshCircle = 0.75; 
Call MySphere; 
 
 
t = 38;
xC = 35.2186; yC = 2.07236;
R = 1.28951; 
meshCircle = 0.75; 
Call MySphere; 
 
 
t = 39;
xC = 25.3321; yC = 18.3112;
R = 1.22056; 
meshCircle = 0.75; 
Call MySphere; 
 
 
t = 40;
xC = 42.1716; yC = 1.64117;
R = 0.677341; 
meshCircle = 0.75; 
Call MySphere; 
 
 
t = 41;
xC = 95.2957; yC = 8.87054;
R = 0.759614; 
meshCircle = 0.75; 
Call MySphere; 
 
 
t = 42;
xC = 11.8291; yC = 6.53579;
R = 1.18948; 
meshCircle = 0.75; 
Call MySphere; 
 
 
t = 43;
xC = 59.7514; yC = 7.97864;
R = 1.18998; 
meshCircle = 0.75; 
Call MySphere; 
 
 
t = 44;
xC = 1.93613; yC = 9.84664;
R = 1.10774; 
meshCircle = 0.75; 
Call MySphere; 
 
 
t = 45;
xC = 78.5744; yC = 12.2125;
R = 1.1245; 
meshCircle = 0.75; 
Call MySphere; 
 
 
t = 46;
xC = 33.1179; yC = 10.9544;
R = 1.13991; 
meshCircle = 0.75; 
Call MySphere; 
 
 
t = 47;
xC = 92.2191; yC = 13.7729;
R = 1.15701; 
meshCircle = 0.75; 
Call MySphere; 
 
 
t = 48;
xC = 33.0557; yC = 18.4578;
R = 1.16058; 
meshCircle = 0.75; 
Call MySphere; 
 
 
t = 49;
xC = 68.6567; yC = 18.2109;
R = 0.947757; 
meshCircle = 0.75; 
Call MySphere; 
 
 
t = 50;
xC = 18.0202; yC = 7.8458;
R = 1.1369; 
meshCircle = 0.75; 
Call MySphere; 
 
 
t = 51;
xC = 29.497; yC = 12.0092;
R = 1.14807; 
meshCircle = 0.75; 
Call MySphere; 
 
 
t = 52;
xC = 82.806; yC = 1.5144;
R = 0.865991; 
meshCircle = 0.75; 
Call MySphere; 
 
 
t = 53;
xC = 25.5805; yC = 4.4079;
R = 0.398269; 
meshCircle = 0.75; 
Call MySphere; 
 
 
t = 54;
xC = 38.8197; yC = 18.245;
R = 0.77283; 
meshCircle = 0.75; 
Call MySphere; 
 
 
t = 55;
xC = 55.6767; yC = 9.95038;
R = 0.766822; 
meshCircle = 0.75; 
Call MySphere; 
 
 
t = 56;
xC = 14.3628; yC = 8.46428;
R = 0.968248; 
meshCircle = 0.75; 
Call MySphere; 
 
 
t = 57;
xC = 31.0859; yC = 2.00462;
R = 0.945563; 
meshCircle = 0.75; 
Call MySphere; 
 
 
t = 58;
xC = 11.9288; yC = 1.88199;
R = 1.07205; 
meshCircle = 0.75; 
Call MySphere; 
 
 
t = 59;
xC = 90.8952; yC = 1.96565;
R = 0.941423; 
meshCircle = 0.75; 
Call MySphere; 
 
 
t = 60;
xC = 9.81539; yC = 3.96909;
R = 1.02851; 
meshCircle = 0.75; 
Call MySphere; 
 
 
t = 61;
xC = 63.2173; yC = 1.52892;
R = 0.981403; 
meshCircle = 0.75; 
Call MySphere; 
 
 
t = 62;
xC = 61.4572; yC = 4.32861;
R = 1.03193; 
meshCircle = 0.75; 
Call MySphere; 
 
 
t = 63;
xC = 50.6028; yC = 1.79541;
R = 0.999344; 
meshCircle = 0.75; 
Call MySphere; 
 
 
volNr = newreg; 
Plane Surface(volNr) = {box, theLoops[]};
Physical Surface(newreg) = volNr;
Physical Surface(newreg) = {theAggregates[]}; 
