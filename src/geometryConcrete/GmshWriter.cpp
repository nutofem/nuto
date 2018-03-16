#include "geometryConcrete/GmshWriter.h"
#include <fstream>

const std::string sAggregate = R"(
Function Aggregate
  p0 = newp; Point(p0) = {xC,  yC,  zC,  meshCircle} ; 
  p1 = newp; Point(p1) = {xC+R,yC,  zC,  meshCircle} ; 
  p2 = newp; Point(p2) = {xC,  yC+R,zC,  meshCircle} ; 
  p3 = newp; Point(p3) = {xC,  yC,  zC+R,meshCircle} ; 
  p4 = newp; Point(p4) = {xC-R,yC,  zC,  meshCircle} ; 
  p5 = newp; Point(p5) = {xC,  yC-R,zC,  meshCircle} ; 
  p6 = newp; Point(p6) = {xC,  yC,  zC-R,meshCircle} ; 
 
  c0  = newl; Circle(c0)  = {p1,p0,p6}; c1  = newl; Circle(c1)  = {p6,p0,p4}; 
  c2  = newl; Circle(c2)  = {p4,p0,p3}; c3  = newl; Circle(c3)  = {p3,p0,p1}; 
  c4  = newl; Circle(c4)  = {p1,p0,p2}; c5  = newl; Circle(c5)  = {p2,p0,p4}; 
  c6  = newl; Circle(c6)  = {p4,p0,p5}; c7  = newl; Circle(c7)  = {p5,p0,p1}; 
  c8  = newl; Circle(c8)  = {p6,p0,p2}; c9  = newl; Circle(c9) =  {p2,p0,p3}; 
  c10 = newl; Circle(c10) = {p3,p0,p5}; c11 = newl; Circle(c11) = {p5,p0,p6}; 
 
  l0 = newll; Line Loop(l0) = { c4,  c9, c3};   
  l1 = newll; Line Loop(l1) = { c8, -c4, c0};  
  l2 = newll; Line Loop(l2) = { c11,-c7,-c0};
  l3 = newll; Line Loop(l3) = { c7, -c3, c10};
  l4 = newll; Line Loop(l4) = {-c9,  c5, c2};
  l5 = newll; Line Loop(l5) = {-c10,-c2, c6};
  l6 = newll; Line Loop(l6) = {-c1, -c6,-c11};
  l7 = newll; Line Loop(l7) = {-c5, -c8, c1};
  
  s0 = news; Ruled Surface(s0) = {l0};
  s1 = news; Ruled Surface(s1) = {l1};
  s2 = news; Ruled Surface(s2) = {l2};
  s3 = news; Ruled Surface(s3) = {l3};
  s4 = news; Ruled Surface(s4) = {l4};
  s5 = news; Ruled Surface(s5) = {l5};
  s6 = news; Ruled Surface(s6) = {l6};
  s7 = news; Ruled Surface(s7) = {l7};
  
  theInnerInterface = newsl;
  theAggregate = newv; 
  Surface Loop(theInnerInterface) = {s0, s1, s2, s3, s4, s5, s6, s7};  
  Volume(theAggregate) = theInnerInterface; 
  theAggregates[i] = theAggregate; 
  
  theOuterInterface[i] = theInnerInterface; 
Return 
)";

const std::string sAggregateWithInterface = R"(
Function Aggregate
  R = R - 0.5 * h;
  p0 = newp; Point(p0) = {xC,  yC,  zC,  meshCircle} ; 
  p1 = newp; Point(p1) = {xC+R,yC,  zC,  meshCircle} ; 
  p2 = newp; Point(p2) = {xC,  yC+R,zC,  meshCircle} ; 
  p3 = newp; Point(p3) = {xC,  yC,  zC+R,meshCircle} ; 
  p4 = newp; Point(p4) = {xC-R,yC,  zC,  meshCircle} ; 
  p5 = newp; Point(p5) = {xC,  yC-R,zC,  meshCircle} ; 
  p6 = newp; Point(p6) = {xC,  yC,  zC-R,meshCircle} ; 
 
  c0  = newl; Circle(c0)  = {p1,p0,p6}; c1  = newl; Circle(c1)  = {p6,p0,p4}; 
  c2  = newl; Circle(c2)  = {p4,p0,p3}; c3  = newl; Circle(c3)  = {p3,p0,p1}; 
  c4  = newl; Circle(c4)  = {p1,p0,p2}; c5  = newl; Circle(c5)  = {p2,p0,p4}; 
  c6  = newl; Circle(c6)  = {p4,p0,p5}; c7  = newl; Circle(c7)  = {p5,p0,p1}; 
  c8  = newl; Circle(c8)  = {p6,p0,p2}; c9  = newl; Circle(c9) =  {p2,p0,p3}; 
  c10 = newl; Circle(c10) = {p3,p0,p5}; c11 = newl; Circle(c11) = {p5,p0,p6}; 
 
  l0 = newll; Line Loop(l0) = { c4,  c9, c3};   
  l1 = newll; Line Loop(l1) = { c8, -c4, c0};  
  l2 = newll; Line Loop(l2) = { c11,-c7,-c0};
  l3 = newll; Line Loop(l3) = { c7, -c3, c10};
  l4 = newll; Line Loop(l4) = {-c9,  c5, c2};
  l5 = newll; Line Loop(l5) = {-c10,-c2, c6};
  l6 = newll; Line Loop(l6) = {-c1, -c6,-c11};
  l7 = newll; Line Loop(l7) = {-c5, -c8, c1};
  
  s0 = news; Ruled Surface(s0) = {l0};
  s1 = news; Ruled Surface(s1) = {l1};
  s2 = news; Ruled Surface(s2) = {l2};
  s3 = news; Ruled Surface(s3) = {l3};
  s4 = news; Ruled Surface(s4) = {l4};
  s5 = news; Ruled Surface(s5) = {l5};
  s6 = news; Ruled Surface(s6) = {l6};
  s7 = news; Ruled Surface(s7) = {l7};
  
  theInnerInterface = newsl;
  theAggregate = newv; 
  Surface Loop(theInnerInterface) = {s0, s1, s2, s3, s4, s5, s6, s7};  
  Volume(theAggregate) = theInnerInterface; 
  theAggregates[i] = theAggregate; 
  
  Ex0[] = Extrude { Surface {s0};  Layers {{1}, {h}}; Recombine;};
  Ex1[] = Extrude { Surface {s1};  Layers {{1}, {h}}; Recombine;};
  Ex2[] = Extrude { Surface {s2};  Layers {{1}, {h}}; Recombine;};
  Ex3[] = Extrude { Surface {s3};  Layers {{1}, {h}}; Recombine;};
  Ex4[] = Extrude { Surface {s4};  Layers {{1}, {h}}; Recombine;};
  Ex5[] = Extrude { Surface {s5};  Layers {{1}, {h}}; Recombine;};
  Ex6[] = Extrude { Surface {s6};  Layers {{1}, {h}}; Recombine;};
  Ex7[] = Extrude { Surface {s7};  Layers {{1}, {h}}; Recombine;};

  theInterfaces[8*i + 0] = Ex0[1];
  theInterfaces[8*i + 1] = Ex1[1];
  theInterfaces[8*i + 2] = Ex2[1];
  theInterfaces[8*i + 3] = Ex3[1];
  theInterfaces[8*i + 4] = Ex4[1];
  theInterfaces[8*i + 5] = Ex5[1];
  theInterfaces[8*i + 6] = Ex6[1];
  theInterfaces[8*i + 7] = Ex7[1];

  theOuterInterface[i] = newsl; 
  Surface Loop(theOuterInterface[i]) = {Ex0[0],Ex1[0],Ex2[0],Ex3[0],Ex4[0],Ex5[0],Ex6[0],Ex7[0]};
Return 
)";


const std::string sBox3D = R"(
Function Box3D
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

  theBox = newreg; 
  Surface Loop(theBox) = {lxS+1, lyS+1, lxE+1, lyE+1, lzS+1, lzE+1}; 
Return
)";

const std::string sCylinder = R"(
Function Cylinder3D
    // base cirle
    p0 = newp; Point(p0) = {0, 0, 0, meshSpecimen};
    p1 = newp; Point(p1) = {0, radius, 0, meshSpecimen};
    p2 = newp; Point(p2) = {0, -radius, 0, meshSpecimen};
    c0 = newl; Circle(c0) = {p1, p0, p2};
    c1 = newl; Circle(c1) = {p2, p0, p1};
    baseLoop = newll; Line Loop(baseLoop) = {c0, c1};
    baseSurface = news; Plane Surface(baseSurface) = {baseLoop};
    // top cirle
    p3 = newp; Point(p3) = {0, 0, height, meshSpecimen};
    p4 = newp; Point(p4) = {0, radius, height, meshSpecimen};
    p5 = newp; Point(p5) = {0, -radius, height, meshSpecimen};
    c2 = newl; Circle(c2) = {p4, p3, p5};
    c3 = newl; Circle(c3) = {p5, p3, p4};
    topLoop = newll; Line Loop(topLoop) = {c2, c3};
    topSurface = news; Plane Surface(topSurface) = {topLoop};
    // barrell surfaces
    verticalOne = newl; Line(verticalOne) = {p2, p5};
    verticalTwo = newl; Line(verticalTwo) = {p1, p4};
    barrelLoopOne = newll; Line Loop(barrelLoopOne) = {c0, verticalOne, -c2, -verticalTwo};
    barrelLoopTwo = newll; Line Loop(barrelLoopTwo) = {c1, verticalTwo, -c3, -verticalOne};
    barrelSurfaceOne = news; Ruled Surface(barrelSurfaceOne) = {barrelLoopOne};
    barrelSurfaceTwo = news; Ruled Surface(barrelSurfaceTwo) = {barrelLoopTwo};
    // surface loop
    theBox = newreg;
    Surface Loop(theBox) = {barrelSurfaceOne, barrelSurfaceTwo, baseSurface, topSurface};
Return
)";

const std::string sBox2D = R"(
Function Box2D
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

  theBox = newreg; Line Loop(theBox) = { l0, l1, l2, l3};
Return
)";

const std::string sAggregate2D = R"(
Function Aggregate2D
  p1 = newp; Point(p1) = {xC,  yC,  0,  meshCircle};
  p2 = newp; Point(p2) = {xC+R,yC,  0,  meshCircle};
  p3 = newp; Point(p3) = {xC-R,yC,  0,  meshCircle};
  c1 = newreg; Circle(c1) = {p2,p1,p3}; 
  c2 = newreg; Circle(c2) = {p3,p1,p2};
  l1 = newreg; Line Loop(l1) = {c1,c2};
  s1 = newreg; Plane Surface(s1) = {l1};
  theOuterInterface[i] = l1;
  theAggregates[i] = s1;
Return
)";

const std::string sAggregate2DWithInterface = R"(
Function Aggregate2D
  p1 = newp; Point(p1) = {xC,  yC,  0,  meshCircle};
  p2 = newp; Point(p2) = {xC+R,yC,  0,  meshCircle};
  p3 = newp; Point(p3) = {xC,yC+R,  0,  meshCircle};
  p4 = newp; Point(p4) = {xC-R,yC,  0,  meshCircle};
  p5 = newp; Point(p5) = {xC,yC-R,  0,  meshCircle};

  c1 = newreg; Circle(c1) = {p2,p1,p3}; 
  c2 = newreg; Circle(c2) = {p3,p1,p4};
  c3 = newreg; Circle(c3) = {p4,p1,p5};
  c4 = newreg; Circle(c4) = {p5,p1,p2};
  l1 = newreg; Line Loop(l1) = {c1,c2,c3,c4};
  theInnerInterface = newreg; Plane Surface(theInnerInterface) = {l1};
  theAggregates[i] = theInnerInterface;

  Ex0[] = Extrude { Line {c1}; Layers {{1}, {h}}; Recombine;};
  Ex1[] = Extrude { Line {c2}; Layers {{1}, {h}}; Recombine;};
  Ex2[] = Extrude { Line {c3}; Layers {{1}, {h}}; Recombine;};
  Ex3[] = Extrude { Line {c4}; Layers {{1}, {h}}; Recombine;};

  theInterfaces[4*i + 0] = Ex0[1];
  theInterfaces[4*i + 1] = Ex1[1];
  theInterfaces[4*i + 2] = Ex2[1];
  theInterfaces[4*i + 3] = Ex3[1];

  theOuterInterface[i] = newsl;
  Line Loop(theOuterInterface[i]) = {Ex0[0],Ex1[0],Ex2[0],Ex3[0]};
Return
)";

void Header(std::ostream& out, NuTo::GmshWriter::Options opt)
{
    out << "Mesh.RecombinationAlgorithm = " << opt.meshRecombinationAlgorithm << ";\n";
    out << "Mesh.Algorithm = " << opt.meshAlg << ";\n";
    out << "Mesh.Algorithm3D = " << opt.meshAlg3D << ";\n";
    out << "Mesh.Optimize = " << opt.meshOptimize << ";\n";
    out << "Mesh.Smoothing = " << opt.meshSmoothing << ";\n";
    out << opt.additionalOptions << "\n";
    out << "meshSpecimen=" << opt.meshSizeMatrix << ";\n";
    out << "meshCircle=" << opt.meshSizeAggregates << ";\n";
}

void Footer3D(std::ostream& out, NuTo::GmshWriter::Options opt)
{
    out << "theMatrix = newv;";
    out << "Volume(theMatrix) = {theBox, theOuterInterface[]};\n";
    out << "Physical Volume(\"" << opt.physicalGroupMatrix << "\") = {theMatrix};\n";
    out << "Physical Volume(\"" << opt.physicalGroupAggregates << "\") = {theAggregates[]};\n";
    if (opt.interfaceThickness != 0)
        out << "Physical Volume(\"" << opt.physicalGroupInterfaces << "\") = {theInterfaces[]};\n";
}

void Footer2D(std::ostream& out, NuTo::GmshWriter::Options opt)
{
    out << "theMatrix = news;";
    out << "Plane Surface(theMatrix) = {theBox, theOuterInterface[]};\n";
    out << "Physical Surface(\"" << opt.physicalGroupMatrix << "\") = {theMatrix};\n";
    out << "Physical Surface(\"" << opt.physicalGroupAggregates << "\") = {theAggregates[]};\n";
    if (opt.interfaceThickness != 0)
        out << "Physical Surface(\"" << opt.physicalGroupInterfaces << "\") = {theInterfaces[]};\n";
}

void Aggregates3D(std::ostream& out, const Eigen::MatrixX4d& a, double thickness)
{
    // write functions to file
    if (thickness == 0)
    {
        out << sAggregate;
    }
    else
    {
        out << sAggregateWithInterface;
        out << "h = " << thickness << ";\n";
    }

    // write function call code
    for (int i = 0; i < a.rows(); ++i)
    {
        out << " i=" << i << ";";
        out << "xC=" << a(i, 0) << "; ";
        out << "yC=" << a(i, 1) << "; ";
        out << "zC=" << a(i, 2) << "; ";
        out << " R=" << a(i, 3) << ";\n";
        out << "Call Aggregate;\n";
    }
}

void NuTo::GmshWriter::Write(std::ostream& out, Box2D box, const Eigen::MatrixX3d& aggregates, Options opt)
{
    Header(out, opt);
    out << sBox2D;
    if (opt.interfaceThickness == 0)
    {
        out << sAggregate2D;
    }
    else
    {
        out << sAggregate2DWithInterface;
        out << "h = " << opt.interfaceThickness << ";\n";
    }

    // write function call code
    for (int i = 0; i < aggregates.rows(); ++i)
    {
        out << "i=" << i << ";";
        out << "xC=" << aggregates(i, 0) << "; ";
        out << "yC=" << aggregates(i, 1) << "; ";
        out << "R=" << aggregates(i, 2) << ";\n";
        out << "Call Aggregate2D;\n";
    }

    out << "xS=" << box.mStart.x() << ";";
    out << "yS=" << box.mStart.y() << ";\n";
    out << "xE=" << box.mEnd.x() << ";";
    out << "yE=" << box.mEnd.y() << ";\n";
    out << "Call Box2D;\n";

    Footer2D(out, opt);
}


void NuTo::GmshWriter::Write(std::ostream& out, Box3D box, const Eigen::MatrixX4d& aggregates, Options opt)
{
    Header(out, opt);
    out << sBox3D;
    Aggregates3D(out, aggregates, opt.interfaceThickness);

    out << "xS=" << box.mStart.x() << ";";
    out << "yS=" << box.mStart.y() << ";";
    out << "zS=" << box.mStart.z() << ";\n";
    out << "xE=" << box.mEnd.x() << ";";
    out << "yE=" << box.mEnd.y() << ";";
    out << "zE=" << box.mEnd.z() << ";\n";
    out << "Call Box3D;\n";

    Footer3D(out, opt);
}

void NuTo::GmshWriter::Write(std::ostream& out, Cylinder cylinder, const Eigen::MatrixX4d& aggregates, Options opt)
{
    Header(out, opt);
    out << sCylinder;
    Aggregates3D(out, aggregates, opt.interfaceThickness);

    out << "radius=" << cylinder.mRadius << ";\n";
    out << "height=" << cylinder.mHeight << ";\n";
    out << "Call Cylinder3D;\n";

    Footer3D(out, opt);
}

void NuTo::GmshWriter::Write(std::string filename, Box2D box, const Eigen::MatrixX3d& aggregates, Options opt)
{
    std::ofstream out(filename);
    Write(out, box, aggregates, opt);
}
void NuTo::GmshWriter::Write(std::string filename, Box3D box, const Eigen::MatrixX4d& aggregates, Options opt)
{
    std::ofstream out(filename);
    Write(out, box, aggregates, opt);
}
void NuTo::GmshWriter::Write(std::string filename, Cylinder cylinder, const Eigen::MatrixX4d& aggregates, Options opt)
{
    std::ofstream out(filename);
    Write(out, cylinder, aggregates, opt);
}
