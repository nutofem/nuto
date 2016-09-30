#include <iostream>
#include "nuto/mechanics/structures/StructureBase.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"

#include <eigen3/Eigen/Core>

#include "nuto/math/FullMatrix.h"
#include "nuto/math/SparseMatrixCSRGeneral.h"
#include "nuto/math/SparseMatrixCSRSymmetric.h"
#include "nuto/math/SparseDirectSolverMUMPS.h"
#include "nuto/math/SparseDirectSolverPardiso.h"

#include "nuto/mechanics/integrationtypes/IntegrationType1D2NLobatto3Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType1D2NLobatto4Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType1D2NLobatto5Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType2D4NLobatto9Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType2D4NLobatto16Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType2D4NLobatto25Ip.h"

#include "nuto/mechanics/IGA/BSplineCurve.h"
#include "nuto/mechanics/IGA/BSplineSurface.h"

#include <boost/filesystem.hpp>

#define PRINTRESULT true

/*
 |>*----*----*----*----*  -->F
*/
NuTo::Structure* buildStructure1D(double& DisplacementCorrect, int refinements)
{
    /** paramaters **/
    double  YoungsModulus = 20000.;
    double  Area = 100.0*0.1;
    double  Length = 1000.0;
    double  Force = 1.;

    DisplacementCorrect = (Force*Length)/(Area*YoungsModulus);

    /** Structure 1D **/
    NuTo::Structure* myStructure = new NuTo::Structure(1);

#ifdef _OPENMP
    int numThreads = 1;
    myStructure->SetNumProcessors(numThreads);
#endif

    /** create section **/
    int Section = myStructure->SectionCreate("Truss");
    myStructure->SectionSetArea(Section, Area);

    /** create material law **/
    int Material = myStructure->ConstitutiveLawCreate("LINEAR_ELASTIC_ENGINEERING_STRESS");
    myStructure->ConstitutiveLawSetParameterDouble(Material, NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, YoungsModulus);

    /** create IGA structure **/

    // create IGA curve
    int degree = 3;
    int numPoints = 10;
    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> points(numPoints,1);
    for(int i = 0; i < numPoints; i++)
    {
        points(i,0) = i*Length/(numPoints-1);
    }

    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> AInv;

    NuTo::BSplineCurve curve(degree, points, AInv);

    for(int i = 0; i < refinements; i++)
    {
        curve.DuplicateKnots();
    }

    if (PRINTRESULT) std::cout << "Control points: \n" << curve.GetControlPoints() << std::endl;
    if (PRINTRESULT) std::cout << "Knots: \n" << curve.GetKnotVector() << std::endl;

    if (PRINTRESULT) std::cout << "Control points: \n" << curve.GetControlPoints() << std::endl;
    if (PRINTRESULT) std::cout << "Knots: \n" << curve.GetKnotVector() << std::endl;

    int num = 100;
    NuTo::FullVector<double, Eigen::Dynamic> parameters(num);
    double inc = (1./num);
    for(int i = 1; i < num-1; i++) parameters[i] = i*inc;
    parameters[num-1] = 1.;

    std::set<NuTo::Node::eDof> setOfDOFS;
    setOfDOFS.insert(NuTo::Node::COORDINATES);
    setOfDOFS.insert(NuTo::Node::DISPLACEMENTS);

    /** create nodes **/
    for(int i = 0; i < curve.GetNumControlPoints(); i++)
    {
        myStructure->NodeCreateDOFs(i, setOfDOFS, curve.GetControlPoint(i));
    }

    int interpolationType = myStructure->InterpolationTypeCreate("IGA1D");

    Eigen::VectorXi vecDegree(1);
    vecDegree(0) = degree;

    std::vector<Eigen::VectorXd> vecKnots;
    vecKnots.push_back(curve.GetKnotVector());
    myStructure->InterpolationTypeAdd(interpolationType, NuTo::Node::COORDINATES, NuTo::Interpolation::eTypeOrder::SPLINE, vecDegree, vecKnots, curve.GetWeights());
    myStructure->InterpolationTypeAdd(interpolationType, NuTo::Node::DISPLACEMENTS, NuTo::Interpolation::eTypeOrder::SPLINE, vecDegree, vecKnots, curve.GetWeights());

    /** create elements **/
    NuTo::FullVector<int, Eigen::Dynamic> elementIncidence(degree+1);
    for(int element = 0; element < curve.GetNumIGAElements(); element++)
    {
        elementIncidence = curve.GetElementControlPointIDs(element);
        myStructure->ElementCreate(interpolationType, elementIncidence, curve.GetElementKnots(element), curve.GetElementKnotIDs(element));
        if (PRINTRESULT) std::cout << "Incidence: \n" << elementIncidence << std::endl;
    }

    myStructure->ElementTotalSetSection(Section);
    myStructure->ElementTotalSetConstitutiveLaw(Material);

    /** set boundary conditions and loads **/
    NuTo::FullVector<double,Eigen::Dynamic> direction(1);
    direction(0) = 1;
    // first node is fixed
    myStructure->ConstraintLinearSetDisplacementNode(0, direction, 0.0);
    myStructure->SetNumLoadCases(1);
    myStructure->LoadCreateNodeForce(0, curve.GetNumControlPoints()-1, direction, Force);

    /** start analysis **/

    myStructure->CalculateMaximumIndependentSets();
    myStructure->NodeInfo(10);
    myStructure->NodeBuildGlobalDofs();

    return myStructure;
}

Eigen::Vector2d rightPressure(Eigen::Vector2d pos)
{
    double stress = -10.;
    double Hight = 5.;
    Eigen::Vector2d stressRet(0, stress * pos(1) * (pos(1) - Hight) );
    std::cout << "Load_Upper: " << stressRet << std::cout;
    return stressRet;
}

NuTo::BSplineSurface buildRect2D(double x0, double y0, double Height, double Length)
{
    /** Knots and control points **/
    int numElementsX = 1;
    int numElementsY = 1;

    Eigen::Vector2i degree(2,2);

    int numKnotsX = 2*(degree(0)+1) + numElementsX - 1;
    int numKnotsY = 2*(degree(1)+1) + numElementsY - 1;

    Eigen::VectorXd knotsX(numKnotsX);
    Eigen::VectorXd knotsY(numKnotsY);

    for (int i = 0; i <= degree(0); i++) knotsX(i) = 0.;
    for (int i = degree(0) + 1; i <= degree(0) + numElementsX - 1 ; i++) knotsX(i) = knotsX(i-1) + 1./numElementsX;
    for (int i = degree(0) + numElementsX ; i < numKnotsX; i++) knotsX(i) = 1;

    for (int i = 0; i <= degree(1); i++) knotsY(i) = 0.;
    for (int i = degree(1) + 1; i <= degree(1) + numElementsY - 1 ; i++) knotsY(i) = knotsY(i-1) + 1./numElementsY;
    for (int i = degree(1) + numElementsY ; i < numKnotsY; i++) knotsY(i) = 1;

    int numControlPointsX = (numKnotsX -1) - degree(0);
    int numControlPointsY = (numKnotsY -1) - degree(1);

    Eigen::Matrix<Eigen::VectorXd, Eigen::Dynamic, Eigen::Dynamic> controlPoints(numControlPointsY, numControlPointsX);

    double incrx = (Length/(numControlPointsX-1));
    double incry = (Height/(numControlPointsY-1));
    for(int i = 0; i < numControlPointsY; i++)
    {
        for(int j = 0; j < numControlPointsX; j++)
        {
            controlPoints(i,j) = Eigen::Vector2d(x0 + incrx*j, y0 + incry*i);
            if (PRINTRESULT) std::cout << "(" << controlPoints(i, j)(0) << "," << controlPoints(i, j)(1) << ")" << "  ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

    Eigen::MatrixXd weights(numControlPointsY, numControlPointsX);
    weights.setOnes(numControlPointsY, numControlPointsX);

    return NuTo::BSplineSurface(degree, knotsX, knotsY, controlPoints, weights);
}

void contactPatchTest()
{
    NuTo::Structure* myStructure = new NuTo::Structure(2);

    /** parameters **/
    double  YoungsModulus = 20000.;
    double  PoissonRatio = 0.3;
    double  Height = 5.;
    // there should be no dependency, because the pressure at the boundary is predefined
    double  Thickness = 2.123548;
    double  Length = 10;

    NuTo::BSplineSurface slave  = buildRect2D(0,      0,  Height, Length);
    NuTo::BSplineSurface master = buildRect2D(0, -Height, Height, Length);

    for(int i = 0; i < 3; i++)
    {
        slave.DuplicateKnots(0);
        slave.DuplicateKnots(1);
    }

    for(int i = 0; i < 1; i++)
    {
        master.DuplicateKnots(0);
        master.DuplicateKnots(1);
    }

    std::set<NuTo::Node::eDof> setOfDOFS;
    setOfDOFS.insert(NuTo::Node::COORDINATES);
    setOfDOFS.insert(NuTo::Node::DISPLACEMENTS);

    int groupNodesSlave  = myStructure->GroupCreate("Nodes");
    int groupNodesMaster = myStructure->GroupCreate("Nodes");

    int groupElementsSlave  = myStructure->GroupCreate("Elements");
    int groupElementsMaster = myStructure->GroupCreate("Elements");

    slave.buildIGAStructure (*myStructure, setOfDOFS, groupElementsSlave, groupNodesSlave);
    master.buildIGAStructure(*myStructure, setOfDOFS, groupElementsMaster, groupNodesMaster);

    /** create section **/
    int mySection = myStructure->SectionCreate("Plane_Stress");
    myStructure->SectionSetThickness(mySection,Thickness);

    /** create constitutive law **/
    int myMatLin = myStructure->ConstitutiveLawCreate("LINEAR_ELASTIC_ENGINEERING_STRESS");
    myStructure->ConstitutiveLawSetParameterDouble(myMatLin, NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, YoungsModulus);
    myStructure->ConstitutiveLawSetParameterDouble(myMatLin, NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, PoissonRatio);

    int groupNodesSlaveContact  = myStructure->GroupCreate("Nodes");
    int groupNodesMasterContact  = myStructure->GroupCreate("Nodes");

    auto LambdaGetBoundaryNodeslower = [](NuTo::NodeBase* rNodePtr) -> bool
                                {
                                    double Tol = 1.e-6;
                                    if (rNodePtr->GetNum(NuTo::Node::COORDINATES)>0)
                                    {
                                        double y = rNodePtr->Get(NuTo::Node::COORDINATES)[1];
                                        if ((y >= 0 - Tol && y <= 0 + Tol))
                                        {
                                            return true;
                                        }
                                    }
                                    return false;
                                };  // LambdaGetBoundaryNodeslower

    myStructure->GroupAddNodeFunction(groupNodesSlaveContact, groupNodesSlave,  LambdaGetBoundaryNodeslower);
    myStructure->GroupAddNodeFunction(groupNodesMasterContact, groupNodesMaster,  LambdaGetBoundaryNodeslower);


    NuTo::FullVector<int,Eigen::Dynamic> rMembers;
    myStructure->NodeGroupGetMembers(groupNodesSlaveContact,  rMembers);
    std::cout << "groupNodesSlaveContact: \n" << rMembers.transpose() << std::endl;

    myStructure->NodeGroupGetMembers(groupNodesMasterContact,  rMembers);
    std::cout << "groupNodesMasterContact: \n" << rMembers.transpose() << std::endl;

    int groupElementsSlaveContact   = myStructure->GroupCreate("Elements");
    int groupElementsMasterContact  = myStructure->GroupCreate("Elements");

    myStructure->GroupAddElementsFromNodes(groupElementsMasterContact, groupNodesMasterContact, false);
    myStructure->GroupAddElementsFromNodes(groupElementsSlaveContact, groupNodesSlaveContact, false);

    myStructure->ElementGroupGetMembers(groupElementsSlaveContact,  rMembers);
    std::cout << "groupElementsSlaveContact: \n" << rMembers.transpose() << std::endl;

    myStructure->ElementGroupGetMembers(groupElementsMasterContact,  rMembers);
    std::cout << "groupNodesMasterContact: \n" << rMembers.transpose() << std::endl;

    std::vector<int> groups;
    groups.push_back(groupElementsSlaveContact);
    groups.push_back(groupElementsMasterContact);
    myStructure->Contact(groups);
}

/*
   ||>*----*----*----*----*
      |    |    |    |    | -->
      |    |    |    |    |
   ||>*----*----*----*----*     Sigma
      |    |    |    |    | -->
      |    |    |    |    |
   ||>*----*----*----*----*
      ^
*/
NuTo::Structure* constantStress(double& DisplacementCorrect, int refinements, const std::string &resultDir)
{
    /** parameters **/
    double  YoungsModulus = 20000.;
    double  PoissonRatio = 0.3;
    double  Height = 5.;
    // there should be no dependency, because the pressure at the boundary is predefined
    double  Thickness = 2.123548;
    double  Length = 10;
    double  Stress = 10.;

    DisplacementCorrect = (Stress*Length)/YoungsModulus;

    NuTo::BSplineSurface surface = buildRect2D(0, 0, Height, Length);

    /** Structure 2D **/
    NuTo::Structure* myStructure = new NuTo::Structure(2);

#ifdef _OPENMP
    int numThreads = 1;
    myStructure->SetNumProcessors(numThreads);
#endif

    /** create section **/
    int mySection = myStructure->SectionCreate("Plane_Stress");
    myStructure->SectionSetThickness(mySection,Thickness);

    /** create constitutive law **/
    int myMatLin = myStructure->ConstitutiveLawCreate("LINEAR_ELASTIC_ENGINEERING_STRESS");
    myStructure->ConstitutiveLawSetParameterDouble(myMatLin,NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS,YoungsModulus);
    myStructure->ConstitutiveLawSetParameterDouble(myMatLin,NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO,PoissonRatio);

    int numX = 100;
    int numY = 100;
    Eigen::MatrixXd parameters(numX*numY, 2);
    int count = 0;
    double incrX = 1./(numX-1);
    double incrY = 1./(numY-1);
    for(int i = 0; i < numX; i++)
    {
        for(int j = 0; j < numY; j++)
        {
            parameters(count, 0) = incrX*i;
            parameters(count, 1) = incrY*j;
            count++;
        }
    }

    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> plateWithHole = surface.SurfacePoints(parameters);
    boost::filesystem::path resultFileName(resultDir);
    std::string ident("Rectangle");
    resultFileName /= ident+".dat";
    plateWithHole.WriteToFile(resultFileName.string(), "  ");

    for(int i = 0; i < surface.GetNumControlPoints(1); i++)
    {
        for(int j = 0; j < surface.GetNumControlPoints(0); j++)
        {
            if (PRINTRESULT) std::cout << "(" << surface.GetControlPoint(i, j)(0) << "," <<  surface.GetControlPoint(i, j)(1) << ")" << "  ";
        }
    }


    for(int i = 0; i < surface.GetNumControlPoints(1); i++)
    {
        for(int j = 0; j < surface.GetNumControlPoints(0); j++)
        {
            if (PRINTRESULT) std::cout << "(" << surface.GetControlPoint(i, j)(0) << "," <<  surface.GetControlPoint(i, j)(1) << ")" << "  ";
        }
        std::cout << std::endl ;
    }
    std::cout << std::endl ;

    for(int i = 0; i < refinements; i++)
    {
        surface.DuplicateKnots(0);
        surface.DuplicateKnots(1);
    }

    if (PRINTRESULT) std::cout << "KnotsX: \n" << surface.GetKnotVector(0) << std::endl;
    if (PRINTRESULT) std::cout << "KnotsY: \n" << surface.GetKnotVector(1)  << std::endl << std::endl;

    for(int i = 0; i < surface.GetNumControlPoints(1); i++)
    {
        for(int j = 0; j < surface.GetNumControlPoints(0); j++)
        {
            if (PRINTRESULT) std::cout << "(" << surface.GetControlPoint(i, j)(0) << "," <<  surface.GetControlPoint(i, j)(1) << ")" << "  ";
        }
        std::cout << std::endl ;
    }

    std::set<NuTo::Node::eDof> setOfDOFS;
    setOfDOFS.insert(NuTo::Node::COORDINATES);
    setOfDOFS.insert(NuTo::Node::DISPLACEMENTS);

    int groupNodes  = myStructure->GroupCreate("Nodes");
    int groupElements = myStructure->GroupCreate("Elements");

    surface.buildIGAStructure(*myStructure, setOfDOFS, groupElements, groupNodes);

    myStructure->Info();

    /** assign constitutive law **/
    myStructure->ElementTotalSetConstitutiveLaw(myMatLin);
    myStructure->ElementTotalSetSection(mySection);

    /** Boundary condition **/

    // Dirichlet
    NuTo::FullVector<double,Eigen::Dynamic> direction(2);


    if (PRINTRESULT) std::cout <<  "Nodes fixed left: ";

    direction << 1, 0;
    for(int i = 0; i < surface.GetNumControlPoints(1); i++)
    {
        myStructure->ConstraintLinearSetDisplacementNode(i*surface.GetNumControlPoints(0), direction, 0.0);
        if (PRINTRESULT) std::cout << i*surface.GetNumControlPoints(0) << ", ";
    }

    if (PRINTRESULT) std::cout << std::endl;

    if (PRINTRESULT) std::cout <<  "Nodes fixed bottom: ";

    direction << 0, 1;
    for(int i = 0; i <  surface.GetNumControlPoints(0); i++)
    {
        myStructure->ConstraintLinearSetDisplacementNode(i, direction, 0.0);
        if (PRINTRESULT) std::cout << i << ", ";
    }
    if (PRINTRESULT) std::cout << std::endl;



    // right boundary
    int groupNumberNodesLeft = myStructure->GroupCreate("NODES");

    if (PRINTRESULT) std::cout <<  "Pressure on nodes: ";

    for(int i = 1; i <=  surface.GetNumControlPoints(1); i++)
    {
        myStructure->GroupAddNode(groupNumberNodesLeft, i*surface.GetNumControlPoints(0) - 1);
        if (PRINTRESULT) std::cout << i* surface.GetNumControlPoints(0) - 1 << ", ";
    }
    if (PRINTRESULT) std::cout << std::endl;

    if (PRINTRESULT) std::cout <<  "Pressure on elements: ";

    int groupNumberElementsLeft = myStructure->GroupCreate("ELEMENTS");
    for(int i = 1; i <= surface.GetNumIGAElements(1); i++)
    {
        myStructure->GroupAddElement(groupNumberElementsLeft, i*surface.GetNumIGAElements(0) - 1);
        if (PRINTRESULT) std::cout << i*surface.GetNumIGAElements(0) - 1 << ", ";
    }
    if (PRINTRESULT) std::cout << std::endl;

    //NuTo::FullVector<double, 2> ForceVecRight({Force,0.});
    myStructure->SetNumLoadCases(1);
    //myStructure->LoadSurfaceConstDirectionCreate2D(0, groupNumberElementsLeft, groupNumberNodesLeft, ForceVecRight);

    myStructure->LoadSurfacePressureCreate2D(0, groupNumberElementsLeft, groupNumberNodesLeft, -Stress);

//    ///////////////// stress function on the right ////////////////////////

//    std::function<Eigen::Vector2d(Eigen::Vector2d)> stress_right  = rightPressure;

//    myStructure->LoadSurfacePressureFunctionCreate2D(0, groupNumberElementsLeft,  groupNumberNodesLeft, stress_right);

//    ///////////////////////////////////////////////////////////////////////



    myStructure->CalculateMaximumIndependentSets();
    myStructure->NodeBuildGlobalDofs();

    myStructure->NodeInfo(10);

    return myStructure;
}

Eigen::VectorXd exact_plate_hole_disp(Eigen::VectorXd pt)
{
    double a = 1.0;
    double T = -1.;
    double E = 1.e5;
    double nu = 0.3;
    Eigen::VectorXd disp(2);
    double r     = std::sqrt(pt(0)*pt(0) + pt(1)*pt(1));
    double theta = std::atan(pt(1)/pt(0));

    double ct = std::cos(theta);
    double c3t = std::cos(3*theta);
    double st = std::sin(theta);
    double s3t = std::sin(3*theta);

    double mu = E/(1. - 2.*nu);
    double k =  (3.- nu)/(1. + nu);

    disp(0) = ((T*a)/(8.*mu)) * ( (r/a)*(k + 1.)*ct + ((2.*a)/r)*( (1. + k)*ct + c3t ) - ((2.*a*a*a)/(r*r*r))*c3t );
    disp(1) = ((T*a)/(8.*mu)) * ( (r/a)*(k - 3.)*st + ((2.*a)/r)*( (1. - k)*st + s3t ) - ((2.*a*a*a)/(r*r*r))*s3t );

//    disp(0) = ((8.*a*(1 + nu)*T)/E) * ((r/a)*(2./(1+nu))*ct + (a/r)*((4./(1+nu))*ct + c3t) - (a/r)*(a/r)*(a/r)*c3t);
//    disp(1) = ((8.*a*(1 + nu)*T)/E) * ((r/a)*((-2.*nu)/(1+nu))*st + (a/r)*(((-4.*nu)/(1+nu))*st + s3t) - (a/r)*(a/r)*(a/r)*s3t);

    return disp;
}

Eigen::VectorXd exact_plate_hole(Eigen::VectorXd pt)
{
    double a = 1.0;
    Eigen::VectorXd stress(3);
    double r     = std::sqrt(pt(0)*pt(0) + pt(1)*pt(1));
    double theta = std::atan(pt(1)/pt(0));

    double c2t = std::cos(2*theta);
    double c4t = std::cos(4*theta);
    double s2t = std::sin(2*theta);
    double s4t = std::sin(4*theta);

    double fac1 = (a/r)*(a/r);
    double fac2 = fac1*fac1;
    stress(0) = 1. - fac1*(1.5*c2t+c4t)+1.5*fac2*c4t; // xx
    stress(1) =    - fac1*(0.5*c2t-c4t)-1.5*fac2*c4t; // yy
    stress(2) =    - fac1*(0.5*s2t+s4t)+1.5*fac2*s4t; // xy

    return stress;
}


Eigen::Vector2d exact_plate_hole_left(Eigen::VectorXd pt)
{
    Eigen::VectorXd stress = exact_plate_hole(pt);
    Eigen::Vector2d stressRet(-stress(0), -stress(2));
    std::cout << "Load_Left: " << stressRet << std::endl;
    return stressRet;
}

Eigen::Vector2d exact_plate_hole_upper(Eigen::VectorXd pt)
{
    Eigen::VectorXd stress = exact_plate_hole(pt);
    Eigen::Vector2d stressRet(-stress(2), -stress(1));
    std::cout << "Load_Upper: " << stressRet << std::endl;
    return stressRet;
}


NuTo::Structure* buildPlateWithHole2DNeumann(const std::string &resultDir, int refine, int BC)
{
    NuTo::Structure* myStructure = new NuTo::Structure(2);

    /*********/
    // Mesh  //
    /*********/

    int noPtsX = 3;
    int noPtsY = 4;

    Eigen::Matrix<Eigen::VectorXd, Eigen::Dynamic, Eigen::Dynamic> controlPts(noPtsX, noPtsY);

    controlPts(0,0) = Eigen::Vector2d(-1, 0);
    controlPts(0,1) = Eigen::Vector2d(-1, 0.4142135623730951);
    controlPts(0,2) = Eigen::Vector2d(-0.4142135623730951, 1);
    controlPts(0,3) = Eigen::Vector2d(0,1);

    controlPts(1,0) = Eigen::Vector2d(-2.5,  0);
    controlPts(1,1) = Eigen::Vector2d(-2.5,  0.75);
    controlPts(1,2) = Eigen::Vector2d(-0.75, 2.5);
    controlPts(1,3) = Eigen::Vector2d( 0,    2.5);

    controlPts(2,0) = Eigen::Vector2d(-4,0);
    controlPts(2,1) = Eigen::Vector2d(-4,4);
    controlPts(2,2) = Eigen::Vector2d(-4,4);
    controlPts(2,3) = Eigen::Vector2d(0,4);

    Eigen::MatrixXd weights(noPtsX, noPtsY);
    weights.setOnes(noPtsX, noPtsY);
    weights(0,1) = (1. + 1./sqrt(2))/2.;
    weights(0,2) = (1. + 1./sqrt(2))/2.;


    Eigen::VectorXd knotsX(7);
    knotsX << 0, 0, 0, 0.5, 1, 1, 1;
    Eigen::VectorXd knotsY(6);
    knotsY << 0, 0, 0, 1, 1, 1;

    Eigen::Vector2i degree(2,2);

    NuTo::BSplineSurface surface(degree, knotsX, knotsY, controlPts, weights);

    int numX = 100;
    int numY = 1;
    Eigen::MatrixXd parameters(numX*numY, 2);
    int count = 0;
    double incrX = 1./(numX-1);
    for(int i = 0; i < numX; i++)
    {
            parameters(count, 0) = incrX*i;
            parameters(count, 1) = 1.;
            count++;
    }

    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> plateWithHole = surface.SurfacePoints(parameters);
    boost::filesystem::path resultFileName(resultDir);
    std::string ident("plateWithHole");
    resultFileName /= ident+".dat";
    plateWithHole.WriteToFile(resultFileName.string(), "  ");

    for(int i = 0; i < refine; i++)
    {
        surface.DuplicateKnots(0);
        surface.DuplicateKnots(1);
    }

    for(int i = 0; i < surface.GetNumControlPoints(1); i++)
    {
        for(int j = 0; j < surface.GetNumControlPoints(0); j++)
        {
            if (PRINTRESULT) std::cout << "(" << surface.GetControlPoint(i, j)(0) << "," <<  surface.GetControlPoint(i, j)(1) << ")" << "  ";
        }
    }

    plateWithHole = surface.SurfacePoints(parameters);
    ident = "plateWithHoleRefined";
    resultFileName = resultDir;
    resultFileName /= ident+".dat";
    plateWithHole.WriteToFile(resultFileName.string(), "  ");

    /**************/
    // Structure  //
    /**************/

    std::set<NuTo::Node::eDof> setOfDOFS;
    setOfDOFS.insert(NuTo::Node::COORDINATES);
    setOfDOFS.insert(NuTo::Node::DISPLACEMENTS);

    if (PRINTRESULT) std::cout << "Nodes:\n";
    /** create nodes **/
    count = 0;
    for(int i = 0; i < surface.GetNumControlPoints(1); i++)
    {
        for(int j = 0; j < surface.GetNumControlPoints(0); j++)
        {
            myStructure->NodeCreateDOFs(count, setOfDOFS, surface.GetControlPoint(i,j));
            if (PRINTRESULT) std::cout << "ID: " << count << ", Coordinate: "  << "(" << surface.GetControlPoint(i, j)(0) << "," <<  surface.GetControlPoint(i, j)(1) << "), weight: " << surface.GetWeights()(i,j) <<std::endl;
            count++;
        }
    }

    std::cout << "Knots X: \n" << surface.GetKnotVector(0).transpose() << std::endl;
    std::cout << "Knots Y: \n" << surface.GetKnotVector(1).transpose() << std::endl;

    std::cout << std::flush;

    if (PRINTRESULT)    std::cout << std::endl;

    int interpolationType = myStructure->InterpolationTypeCreate("IGA2D");

    std::vector<Eigen::VectorXd> vecKnots;
    vecKnots.push_back(surface.GetKnotVector(0));
    vecKnots.push_back(surface.GetKnotVector(1));

    myStructure->InterpolationTypeAdd(interpolationType,
                                      NuTo::Node::COORDINATES,
                                      NuTo::Interpolation::eTypeOrder::SPLINE,
                                      degree,
                                      vecKnots,
                                      surface.GetWeights());

    myStructure->InterpolationTypeAdd(interpolationType,
                                      NuTo::Node::DISPLACEMENTS,
                                      NuTo::Interpolation::eTypeOrder::SPLINE,
                                      degree,
                                      vecKnots,
                                      surface.GetWeights());

    /** create elements **/
    Eigen::VectorXi elementIncidence((degree(0)+1)*(degree(1)+1));
    int numElementsX = surface.GetNumIGAElements(0);
    int numElementsY = surface.GetNumIGAElements(1);
    for(int elementY = 0; elementY < numElementsY; elementY++)
    {
        for(int elementX = 0; elementX < numElementsX; elementX++)
        {
            elementIncidence = surface.GetElementControlPointIDs(elementX, elementY);
            if (PRINTRESULT) std::cout << "Incidence: \n" << elementIncidence << std::endl;
            myStructure->ElementCreate(interpolationType, elementIncidence, surface.GetElementKnots(elementX, elementY), surface.GetElementKnotIDs(elementX, elementY));
        }
    }

    /** create section **/
    double Thickness = 1.;
    int mySection = myStructure->SectionCreate("Plane_Stress");
    myStructure->SectionSetThickness(mySection,Thickness);

    /** create constitutive law **/
    double YoungsModulus = 1.e5;
    double PoissonRatio  = 0.3;
    int myMatLin = myStructure->ConstitutiveLawCreate("LINEAR_ELASTIC_ENGINEERING_STRESS");
    myStructure->ConstitutiveLawSetParameterDouble(myMatLin, NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, YoungsModulus);
    myStructure->ConstitutiveLawSetParameterDouble(myMatLin, NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, PoissonRatio);

    myStructure->ElementTotalSetConstitutiveLaw(myMatLin);
    myStructure->ElementTotalSetSection(mySection);

    /**********************/
    // Boundary condition //
    /**********************/

    int groupElementsLeft  = myStructure->GroupCreate("ELEMENTS");
    int groupElementsUpper = myStructure->GroupCreate("ELEMENTS");

    if (PRINTRESULT) std::cout << "Neumann elements left: ";
    int start = (surface.GetNumIGAElements(1)-1) * surface.GetNumIGAElements(0);
    for(int i = start; i < start + surface.GetNumIGAElements(0)/2; i++)
    {
        if (PRINTRESULT) std::cout << i << ", ";
        myStructure->GroupAddElement(groupElementsLeft, i);
    }

    if (PRINTRESULT) std::cout << "\n";


    if (PRINTRESULT) std::cout << "Neumann elements upper: ";
    start += surface.GetNumIGAElements(0)/2;
    for(int i = start; i < surface.GetNumIGAElements(); i++)
    {
        if (PRINTRESULT) std::cout << i << ", ";
        myStructure->GroupAddElement(groupElementsUpper, i);
    }

    if (PRINTRESULT) std::cout << "\n";

    auto LambdaGetBoundaryNodesLeft = [](NuTo::NodeBase* rNodePtr) -> bool
                                {
                                    double Tol = 1.e-6;
                                    if (rNodePtr->GetNum(NuTo::Node::COORDINATES)>0)
                                    {
                                        double x = rNodePtr->Get(NuTo::Node::COORDINATES)[0];
                                        if ((x >= -4 - Tol && x <= -4 + Tol))
                                        {
                                            return true;
                                        }
                                    }
                                    return false;
                                };  // GetBoundaryNodesLambda

    auto LambdaGetBoundaryNodesUpper = [](NuTo::NodeBase* rNodePtr) -> bool
                                {
                                    double Tol = 1.e-6;
                                    if (rNodePtr->GetNum(NuTo::Node::COORDINATES)>0)
                                    {
                                        double y = rNodePtr->Get(NuTo::Node::COORDINATES)[1];
                                        if ((y >= 4 - Tol && y <= 4 + Tol))
                                        {
                                            return true;
                                        }
                                    }
                                    return false;
                                };  // GetBoundaryNodesLambda


    auto LambdaGetBoundaryNodesLower = [](NuTo::NodeBase* rNodePtr) -> bool
                                {
                                    double Tol = 1.e-6;
                                    if (rNodePtr->GetNum(NuTo::Node::COORDINATES)>0)
                                    {
                                        double y = rNodePtr->Get(NuTo::Node::COORDINATES)[1];
                                        if ((y >= 0 - Tol && y <= 0 + Tol))
                                        {
                                            return true;
                                        }
                                    }
                                    return false;
                                };  // LambdaGetBoundaryNodesLowerSymmetry

    auto LambdaGetBoundaryNodesRight = [](NuTo::NodeBase* rNodePtr) -> bool
                                {
                                    double Tol = 1.e-6;
                                    if (rNodePtr->GetNum(NuTo::Node::COORDINATES)>0)
                                    {
                                        double x = rNodePtr->Get(NuTo::Node::COORDINATES)[0];
                                        if ((x >= 0 - Tol && x <= 0 + Tol))
                                        {
                                            return true;
                                        }
                                    }
                                    return false;
                                };  // LambdaGetBoundaryNodesRight

    auto LambdaGetBoundaryNodesCorner = [](NuTo::NodeBase* rNodePtr) -> bool
                                {
                                    double Tol = 1.e-6;
                                    if (rNodePtr->GetNum(NuTo::Node::COORDINATES)>0)
                                    {
                                        double x = rNodePtr->Get(NuTo::Node::COORDINATES)[0];
                                        double y = rNodePtr->Get(NuTo::Node::COORDINATES)[1];
                                        if ((y >= 4 - Tol && y <= 4 + Tol) && (x >= -4 - Tol && x <= -4 + Tol) )
                                        {
                                            return true;
                                        }
                                    }
                                    return false;
                                };  // GetBoundaryNodesLambda

    int groupNodesLeft  = myStructure->GroupCreate(NuTo::Groups::eGroupId::Nodes);
    int groupNodesUpper = myStructure->GroupCreate(NuTo::Groups::eGroupId::Nodes);

    int groupNodesLower = myStructure->GroupCreate(NuTo::Groups::eGroupId::Nodes);
    int groupNodesRight = myStructure->GroupCreate(NuTo::Groups::eGroupId::Nodes);

    int groupNodesCorner = myStructure->GroupCreate(NuTo::Groups::eGroupId::Nodes);

    myStructure->GroupAddNodeFunction(groupNodesLower, LambdaGetBoundaryNodesLower);
    myStructure->GroupAddNodeFunction(groupNodesRight, LambdaGetBoundaryNodesRight);

    myStructure->GroupAddNodeFunction(groupNodesLeft, LambdaGetBoundaryNodesLeft);
    myStructure->GroupAddNodeFunction(groupNodesUpper, LambdaGetBoundaryNodesUpper);

    myStructure->GroupAddNodeFunction(groupNodesCorner, LambdaGetBoundaryNodesCorner);

    NuTo::FullVector<int,Eigen::Dynamic> rMembers;
    myStructure->NodeGroupGetMembers(groupNodesLeft, rMembers);
    std::cout << "Node Members group Left: \n" << rMembers << std::endl;

    myStructure->NodeGroupGetMembers(groupNodesUpper, rMembers);
    std::cout << "Node Members group Upper: \n" << rMembers << std::endl;

    myStructure->ElementGroupGetMembers(groupElementsUpper, rMembers);
    std::cout << "Element Members group Upper: \n" << rMembers << std::endl;

    myStructure->ElementGroupGetMembers(groupElementsLeft, rMembers);
    std::cout << "Element Members group Left: \n" << rMembers << std::endl;


    myStructure->NodeGroupGetMembers(groupNodesLower, rMembers);
    std::cout << "Node Members group Lower: \n" << rMembers << std::endl;

    int countDBC = 0;
    NuTo::FullVector<double,Eigen::Dynamic> direction(2);
    direction << 0 ,1;
    if (PRINTRESULT) std::cout << "Dirichlet nodes in y dir : ";
    for(int i = 0; i < rMembers.rows(); i++)
    {
        if (PRINTRESULT) std::cout <<  rMembers(i) << ", ";
        countDBC = myStructure->ConstraintLinearSetDisplacementNode(rMembers(i), direction, 0.0);
    }
    std::cout << "\n";

    myStructure->NodeGroupGetMembers(groupNodesRight, rMembers);
    std::cout << "Node Members group Right: \n" << rMembers << std::endl;

    if (PRINTRESULT) std::cout << "Dirichlet nodes in x dir : ";
    direction << 1, 0;
    for(int i = 0; i < rMembers.rows(); i++)
    {
        if (PRINTRESULT) std::cout <<  rMembers(i) << ", ";
        countDBC = myStructure->ConstraintLinearSetDisplacementNode(rMembers(i), direction, 0.0);
    }

    std::cout << "\n";

    myStructure->SetNumLoadCases(1);

    if(BC == 0)
    {
        std::function<Eigen::Vector2d(Eigen::Vector2d)> stress_left  = exact_plate_hole_left;
        std::function<Eigen::Vector2d(Eigen::Vector2d)> stress_upper = exact_plate_hole_upper;

        myStructure->LoadSurfacePressureFunctionCreate2D(0, groupElementsLeft,  groupNodesLeft, stress_left);
        myStructure->LoadSurfacePressureFunctionCreate2D(0, groupElementsUpper, groupNodesUpper, stress_upper);
    }
    else
    {
        double Stress = -10.;
        myStructure->LoadSurfacePressureCreate2D(0, groupElementsLeft, groupNodesLeft, Stress);
    }

    countDBC++;

//    myStructure->NodeGroupGetMembers(groupNodesCorner, rMembers);
//    std::cout << "Node Members group Corner: \n" << rMembers << std::endl;
//    for(int dof = 0; dof < 1; dof++)
//    {
//        myStructure->ConstraintLinearEquationCreate (countDBC, rMembers(0), NuTo::Node::DISPLACEMENTS, dof, 1., 0.);
//        myStructure->ConstraintLinearEquationAddTerm(countDBC, rMembers(1), NuTo::Node::DISPLACEMENTS, dof, -1.);
//        countDBC++;
//    }

    myStructure->CalculateMaximumIndependentSets();
    myStructure->NodeBuildGlobalDofs();

    myStructure->NodeInfo(10);

    return myStructure;
}




void solve(NuTo::Structure *myStructure, double solution, const std::string &resultDir, const std::string &name, bool exc, double tol = 1.e-6)
{
    myStructure->SolveGlobalSystemStaticElastic();

    double nodeDisp = myStructure->NodeGetNodePtr(myStructure->GetNumNodes()-1)->Get(NuTo::Node::DISPLACEMENTS)[0];

    if (PRINTRESULT)
    {
        std::cout << "Displacement node numerical: \n";
        for(int i = 0; i < myStructure->GetNumNodes(); i++)
            std::cout << i << ": \n" << myStructure->NodeGetNodePtr(i)->Get(NuTo::Node::DISPLACEMENTS) << std::endl;
    }

    std::cout << "Dimension: " << myStructure->GetDimension() << std::endl;
    std::cout << "Displacement last node numerical:\n" << nodeDisp << std::endl;
    std::cout << "Displacement last node analytical:\n" << solution << std::endl;

    int visualizationGroup = myStructure->GroupGetElementsTotal();
    myStructure->AddVisualizationComponent(visualizationGroup, NuTo::VisualizeBase::DISPLACEMENTS);
    myStructure->AddVisualizationComponent(visualizationGroup, NuTo::VisualizeBase::ENGINEERING_STRAIN);
    myStructure->AddVisualizationComponent(visualizationGroup, NuTo::VisualizeBase::ENGINEERING_STRESS);

    myStructure->ExportVtkDataFileElements(resultDir+"/Elements" + name + ".vtu", true);
    myStructure->ExportVtkDataFileNodes(resultDir+"/Nodes" + name + ".vtu", true);

    if(exc)
    {
        if (fabs(nodeDisp - solution)/fabs(solution) > tol)
        {
            throw NuTo::Exception("[LobattoIntegration] : displacement is not correct");
        }
    }

}


//! @brief Imports a mesh file and builds the hessian and the internal gradient.
void Neumann(const std::string &resultDir, const std::string &rMeshFile, const std::string &fileName, int BC)
{
    NuTo::Structure myStructure(2);
    NuTo::FullMatrix<int, Eigen::Dynamic, Eigen::Dynamic> groupIndices = myStructure.ImportFromGmsh(rMeshFile + fileName, NuTo::ElementData::CONSTITUTIVELAWIP, NuTo::IpData::NOIPDATA);

    std::cout << groupIndices.size() << std::endl;

    int interpolationType = groupIndices.GetValue(0, 1);
    myStructure.InterpolationTypeAdd(interpolationType, NuTo::Node::DISPLACEMENTS, NuTo::Interpolation::LOBATTO3);

    myStructure.SetVerboseLevel(10);
    myStructure.ElementConvertToInterpolationType(groupIndices.GetValue(0, 0));

//    myStructure.InterpolationTypeSetIntegrationType(interpolationType, NuTo::IntegrationType::IntegrationType2D3NGauss3Ip, NuTo::IpData::NOIPDATA);

    myStructure.InterpolationTypeInfo(0);

    myStructure.NodeBuildGlobalDofs();
    int section = myStructure.SectionCreate("PLANE_STRESS");

    double Thickness = 1.;
    myStructure.SectionSetThickness(section, Thickness);

    myStructure.ElementTotalSetSection(section);

    int constitutiveLaw = myStructure.ConstitutiveLawCreate("Linear_Elastic_Engineering_Stress");
    myStructure.ConstitutiveLawSetParameterDouble(constitutiveLaw, NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, 1.e5);
    myStructure.ConstitutiveLawSetParameterDouble(constitutiveLaw, NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, 0.3);
    myStructure.ElementTotalSetConstitutiveLaw(constitutiveLaw);

    myStructure.CalculateMaximumIndependentSets();

    auto LambdaGetBoundaryNodesLeft = [](NuTo::NodeBase* rNodePtr) -> bool
                                {
                                    double Tol = 1.e-6;
                                    if (rNodePtr->GetNum(NuTo::Node::COORDINATES)>0)
                                    {
                                        double x = rNodePtr->Get(NuTo::Node::COORDINATES)[0];
                                        if ((x >= -4 - Tol && x <= -4 + Tol))
                                        {
                                            return true;
                                        }
                                    }
                                    return false;
                                };  // GetBoundaryNodesLambda

    auto LambdaGetBoundaryNodesUpper = [](NuTo::NodeBase* rNodePtr) -> bool
                                {
                                    double Tol = 1.e-6;
                                    if (rNodePtr->GetNum(NuTo::Node::COORDINATES)>0)
                                    {
                                        double y = rNodePtr->Get(NuTo::Node::COORDINATES)[1];
                                        if ((y >= 4 - Tol && y <= 4 + Tol))
                                        {
                                            return true;
                                        }
                                    }
                                    return false;
                                };  // GetBoundaryNodesLambda

    auto LambdaGetBoundaryNodesLowerSymmetry = [](NuTo::NodeBase* rNodePtr) -> bool
                                {
                                    double Tol = 1.e-6;
                                    if (rNodePtr->GetNum(NuTo::Node::COORDINATES)>0)
                                    {
                                        double y = rNodePtr->Get(NuTo::Node::COORDINATES)[1];
                                        if ((y >= 0 - Tol && y <= 0 + Tol))
                                        {
                                            return true;
                                        }
                                    }
                                    return false;
                                };  // LambdaGetBoundaryNodesLowerSymmetry

    auto LambdaGetBoundaryNodesRightSymmetry = [](NuTo::NodeBase* rNodePtr) -> bool
                                {
                                    double Tol = 1.e-6;
                                    if (rNodePtr->GetNum(NuTo::Node::COORDINATES)>0)
                                    {
                                        double x = rNodePtr->Get(NuTo::Node::COORDINATES)[0];
                                        if ((x >= 0 - Tol && x <= 0 + Tol))
                                        {
                                            return true;
                                        }
                                    }
                                    return false;
                                };  // LambdaGetBoundaryNodesLowerSymmetry



    int groupNodeBCLeft  = myStructure.GroupCreate(NuTo::Groups::eGroupId::Nodes);
    int groupNodeBCUpper = myStructure.GroupCreate(NuTo::Groups::eGroupId::Nodes);
    int groupNodeBCRightSymmetry = myStructure.GroupCreate(NuTo::Groups::eGroupId::Nodes);
    int groupNodeBCLowerSymmetry  = myStructure.GroupCreate(NuTo::Groups::eGroupId::Nodes);

    myStructure.GroupAddNodeFunction(groupNodeBCLeft, LambdaGetBoundaryNodesLeft);
    myStructure.GroupAddNodeFunction(groupNodeBCUpper, LambdaGetBoundaryNodesUpper);
    myStructure.GroupAddNodeFunction(groupNodeBCLowerSymmetry, LambdaGetBoundaryNodesLowerSymmetry);
    myStructure.GroupAddNodeFunction(groupNodeBCRightSymmetry, LambdaGetBoundaryNodesRightSymmetry);

    // Dirichlet symmetry //
    NuTo::FullVector<double,Eigen::Dynamic> direction(2);
    NuTo::FullVector<int,Eigen::Dynamic> rMembers;

    direction << 0 ,1;
    myStructure.NodeGroupGetMembers(groupNodeBCLowerSymmetry, rMembers);
    for(int i = 0; i < rMembers.rows(); i++)
    {
        myStructure.ConstraintLinearSetDisplacementNode(rMembers(i), direction, 0.0);
    }

    direction << 1, 0;

    myStructure.NodeGroupGetMembers(groupNodeBCRightSymmetry, rMembers);
    for(int i = 0; i < rMembers.rows(); i++)
    {
        myStructure.ConstraintLinearSetDisplacementNode(rMembers(i), direction, 0.0);
    }

    int groupelementBCLeft = myStructure.GroupCreate("ELEMENTS");
    int groupelementBCUpper = myStructure.GroupCreate("ELEMENTS");

    myStructure.GroupAddElementsFromNodes(groupelementBCLeft, groupNodeBCLeft, false);
    myStructure.GroupAddElementsFromNodes(groupelementBCUpper, groupNodeBCUpper, false);

    myStructure.SetNumLoadCases(1);

    if(BC == 0)
    {
        std::function<Eigen::Vector2d(Eigen::Vector2d)> stress_left  = exact_plate_hole_left;
        std::function<Eigen::Vector2d(Eigen::Vector2d)> stress_upper = exact_plate_hole_upper;

        myStructure.LoadSurfacePressureFunctionCreate2D(0, groupelementBCLeft,  groupNodeBCLeft, stress_left);
        myStructure.LoadSurfacePressureFunctionCreate2D(0, groupelementBCUpper, groupNodeBCUpper, stress_upper);
    }
    else
    {
        double Stress = -10.;
        myStructure.LoadSurfacePressureCreate2D(0, groupelementBCLeft, groupNodeBCLeft, Stress);
    }

    myStructure.CalculateMaximumIndependentSets();
    myStructure.NodeBuildGlobalDofs();

    myStructure.SolveGlobalSystemStaticElastic();

    int visualizationGroup = myStructure.GroupGetElementsTotal();
    myStructure.AddVisualizationComponent(visualizationGroup, NuTo::VisualizeBase::DISPLACEMENTS);
    myStructure.AddVisualizationComponent(visualizationGroup, NuTo::VisualizeBase::ENGINEERING_STRAIN);
    myStructure.AddVisualizationComponent(visualizationGroup, NuTo::VisualizeBase::ENGINEERING_STRESS);

    myStructure.ExportVtkDataFileElements(resultDir + "/Elements" + fileName + ".vtu", true);
    myStructure.ExportVtkDataFileNodes(resultDir + "/Nodes" + fileName + ".vtu", true);
}


//! @brief Imports a mesh file and builds the hessian and the internal gradient.
void Dirichlet(const std::string &resultDir, const std::string &rMeshFile, const std::string &fileName)
{
    NuTo::Structure myStructure(2);
    NuTo::FullMatrix<int, Eigen::Dynamic, Eigen::Dynamic> groupIndices = myStructure.ImportFromGmsh(rMeshFile + fileName, NuTo::ElementData::CONSTITUTIVELAWIP, NuTo::IpData::NOIPDATA);

    std::cout << groupIndices.size() << std::endl;

    int interpolationType = groupIndices.GetValue(0, 1);
    myStructure.InterpolationTypeAdd(interpolationType, NuTo::Node::DISPLACEMENTS, NuTo::Interpolation::LOBATTO3);

    myStructure.SetVerboseLevel(10);
    myStructure.ElementConvertToInterpolationType(groupIndices.GetValue(0, 0));

//    myStructure.InterpolationTypeSetIntegrationType(interpolationType, NuTo::IntegrationType::IntegrationType2D3NGauss3Ip, NuTo::IpData::NOIPDATA);

    myStructure.InterpolationTypeInfo(0);

    myStructure.NodeBuildGlobalDofs();
    int section = myStructure.SectionCreate("PLANE_STRESS");
    myStructure.SectionSetThickness(section, 1.);

    myStructure.ElementTotalSetSection(section);

    int constitutiveLaw = myStructure.ConstitutiveLawCreate("Linear_Elastic_Engineering_Stress");
    myStructure.ConstitutiveLawSetParameterDouble(constitutiveLaw, NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, 1.e5);
    myStructure.ConstitutiveLawSetParameterDouble(constitutiveLaw, NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, 0.3);
    myStructure.ElementTotalSetConstitutiveLaw(constitutiveLaw);

    myStructure.CalculateMaximumIndependentSets();

    auto LambdaGetBoundaryNodesRight = [](NuTo::NodeBase* rNodePtr) -> bool
                                {
                                    double Tol = 1.e-6;
                                    if (rNodePtr->GetNum(NuTo::Node::COORDINATES)>0)
                                    {
                                        double x = rNodePtr->Get(NuTo::Node::COORDINATES)[0];
                                        if ((x >= 4.0 - Tol && x <= 4.0 + Tol))
                                        {
                                            return true;
                                        }
                                    }
                                    return false;
                                };  // GetBoundaryNodesLambda

    auto LambdaGetBoundaryNodesLeft = [](NuTo::NodeBase* rNodePtr) -> bool
                                {
                                    double Tol = 1.e-6;
                                    if (rNodePtr->GetNum(NuTo::Node::COORDINATES)>0)
                                    {
                                        double x = rNodePtr->Get(NuTo::Node::COORDINATES)[0];
                                        if ((x >= -4 - Tol && x <= -4 + Tol))
                                        {
                                            return true;
                                        }
                                    }
                                    return false;
                                };  // GetBoundaryNodesLambda



    int groupNodeBCRight = myStructure.GroupCreate(NuTo::Groups::eGroupId::Nodes);
    int groupNodeBCLeft  = myStructure.GroupCreate(NuTo::Groups::eGroupId::Nodes);

    myStructure.GroupAddNodeFunction(groupNodeBCRight, LambdaGetBoundaryNodesRight);
    myStructure.GroupAddNodeFunction(groupNodeBCLeft, LambdaGetBoundaryNodesLeft);

    NuTo::FullVector<double,Eigen::Dynamic> direction(2);
    direction << 1., 0.;

    myStructure.ConstraintLinearSetDisplacementNodeGroup(groupNodeBCLeft,  direction, -0.01);
    myStructure.ConstraintLinearSetDisplacementNodeGroup(groupNodeBCRight, direction,  0.01);

    direction << 0., 1.;
    myStructure.ConstraintLinearSetDisplacementNode(0, direction, 0.);

    myStructure.CalculateMaximumIndependentSets();
    myStructure.NodeBuildGlobalDofs();

    myStructure.SolveGlobalSystemStaticElastic();

    int visualizationGroup = myStructure.GroupGetElementsTotal();
    myStructure.AddVisualizationComponent(visualizationGroup, NuTo::VisualizeBase::DISPLACEMENTS);
    myStructure.AddVisualizationComponent(visualizationGroup, NuTo::VisualizeBase::ENGINEERING_STRAIN);
    myStructure.AddVisualizationComponent(visualizationGroup, NuTo::VisualizeBase::ENGINEERING_STRESS);

    myStructure.ExportVtkDataFileElements(resultDir + "/Elements" + fileName + ".vtu", true);
    myStructure.ExportVtkDataFileNodes(resultDir + "/Nodes" + fileName + ".vtu", true);
}


int main()
{
    double solution = 0;
    NuTo::Structure* myStructure = NULL;

    std::string resultDir = "./ResultsIGA";

    if (boost::filesystem::exists(resultDir))
    {
        if (boost::filesystem::is_directory(resultDir))
        {
            boost::filesystem::remove_all(resultDir);
        }
    }

    // create result directory
    boost::filesystem::create_directory(resultDir);

    // 1D
    myStructure = buildStructure1D(solution, 1);
    solve(myStructure, solution, resultDir, "1d", true);

    // 2D constant stress

    myStructure = constantStress(solution, 0, resultDir);
    solve(myStructure, solution, resultDir, "Rectangle0", true);

    myStructure = constantStress(solution, 1, resultDir);
    solve(myStructure, solution, resultDir, "Rectangle1", true);

    contactPatchTest();

    // plate with hole
//    myStructure = buildPlateWithHole2DNeumann(resultDir, 2, 1);
//    solve(myStructure, solution, resultDir, "Hole2", false);

//    myStructure = buildPlateWithHole2DNeumann(resultDir, 5, 1);
//    solve(myStructure, solution, resultDir, "Hole5", false);

//    std::string meshFile2D = "/home/potto1/igafem/igafem/gmshFiles/";

//    std::string fileName = "11-voidedF2.msh";
////    Dirichlet(resultDir, meshFile2D, fileName);
//    Neumann(resultDir, meshFile2D, fileName, 1);

    return 0;
}
