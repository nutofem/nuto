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
NuTo::Structure* buildStructure1D(double& DisplacementCorrect)
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
    int numThreads = 4;
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

    if (PRINTRESULT) std::cout << "Control points: \n" << curve.GetControlPoints() << std::endl;
    if (PRINTRESULT) std::cout << "Knots: \n" << curve.GetKnotVector() << std::endl;

    int num = 100;
    NuTo::FullVector<double, Eigen::Dynamic> parameters(num);
    double inc = (1./num);
    for(int i = 1; i < num-1; i++) parameters[i] = i*inc;
    parameters[num-1] = 1.;

    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> curvepoints = curve.CurvePoints(parameters);

    // write the data
    std::string resultDir = std::string("./IGAResults");
    if (boost::filesystem::exists(resultDir))
    {
        if (boost::filesystem::is_directory(resultDir))
        {
            boost::filesystem::remove_all(resultDir);
        }
    }
    // create result directory
    boost::filesystem::create_directory(resultDir);

    boost::filesystem::path resultFileName(resultDir);
    std::string ident("curve");
    resultFileName /= ident+".dat";

    curvepoints.WriteToFile(resultFileName.string(), "  ");

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
    myStructure->InterpolationTypeAdd(interpolationType, NuTo::Node::COORDINATES, NuTo::Interpolation::eTypeOrder::SPLINE, vecDegree, vecKnots);
    myStructure->InterpolationTypeAdd(interpolationType, NuTo::Node::DISPLACEMENTS, NuTo::Interpolation::eTypeOrder::SPLINE, vecDegree, vecKnots);

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
NuTo::Structure* buildRect2D(double& DisplacementCorrect)
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

    /** Structure 2D **/
    NuTo::Structure* myStructure = new NuTo::Structure(2);

#ifdef _OPENMP
    int numThreads = 4;
    myStructure->SetNumProcessors(numThreads);
#endif

    /** create section **/
    int mySection = myStructure->SectionCreate("Plane_Stress");
    myStructure->SectionSetThickness(mySection,Thickness);

    /** create constitutive law **/
    int myMatLin = myStructure->ConstitutiveLawCreate("LINEAR_ELASTIC_ENGINEERING_STRESS");
    myStructure->ConstitutiveLawSetParameterDouble(myMatLin,NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS,YoungsModulus);
    myStructure->ConstitutiveLawSetParameterDouble(myMatLin,NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO,PoissonRatio);

    /** Knots and control points **/
    int numElementsX = 3;
    int numElementsY = 2;

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

    if (PRINTRESULT) std::cout << "KnotsX: \n" << knotsX << std::endl;

    if (PRINTRESULT) std::cout << "KnotsY: \n" << knotsY << std::endl;

    Eigen::Matrix<Eigen::VectorXd, Eigen::Dynamic, Eigen::Dynamic> controlPoints(numControlPointsX, numControlPointsY);

    for(int i = 0; i < numControlPointsY; i++)
    {
        for(int j = 0; j < numControlPointsX; j++)
        {
            double yCoord = (Height/(numControlPointsY-1))*i;
            double xCoord = (Length/(numControlPointsX-1))*j;
            Eigen::VectorXd temp(2);
            temp << xCoord, yCoord;
            controlPoints(j, i) = temp;
            if (PRINTRESULT) std::cout << std::endl << controlPoints(j,i) << std::endl << std::endl;
        }
    }

    NuTo::BSplineSurface surface(degree, knotsX, knotsY, controlPoints);

    std::set<NuTo::Node::eDof> setOfDOFS;
    setOfDOFS.insert(NuTo::Node::COORDINATES);
    setOfDOFS.insert(NuTo::Node::DISPLACEMENTS);

    if (PRINTRESULT) std::cout << "Nodes:\n";
    /** create nodes **/
    int count = 0;
    for(int i = 0; i < surface.GetNumControlPoints(1); i++)
    {
        for(int j = 0; j < surface.GetNumControlPoints(0); j++)
        {
            myStructure->NodeCreateDOFs(count, setOfDOFS, surface.GetControlPoint(j,i));
            if (PRINTRESULT) std::cout << "ID: " << count << ", Coordinate: \n" << surface.GetControlPoint(j,i) << std::endl;
            count++;
        }
    }

    if (PRINTRESULT)    std::cout << std::endl;

    int interpolationType = myStructure->InterpolationTypeCreate("IGA2D");

    std::vector<Eigen::VectorXd> vecKnots;
    vecKnots.push_back(surface.GetKnotVector(0));
    vecKnots.push_back(surface.GetKnotVector(1));

    myStructure->InterpolationTypeAdd(interpolationType,
                                      NuTo::Node::COORDINATES,
                                      NuTo::Interpolation::eTypeOrder::SPLINE,
                                      degree,
                                      vecKnots);

    myStructure->InterpolationTypeAdd(interpolationType,
                                      NuTo::Node::DISPLACEMENTS,
                                      NuTo::Interpolation::eTypeOrder::SPLINE,
                                      degree,
                                      vecKnots);

    /** create elements **/
    Eigen::VectorXi elementIncidence((degree(0)+1)*(degree(1)+1));
    numElementsX = surface.GetNumIGAElements(0);
    numElementsY = surface.GetNumIGAElements(1);
    for(int elementY = 0; elementY < numElementsY; elementY++)
    {
        for(int elementX = 0; elementX < numElementsX; elementX++)
        {
            elementIncidence = surface.GetElementControlPointIDs(elementX, elementY);
            if (PRINTRESULT) std::cout << "Incidence: \n" << elementIncidence << std::endl;
            myStructure->ElementCreate(interpolationType, elementIncidence, surface.GetElementKnots(elementX, elementY), surface.GetElementKnotIDs(elementX, elementY));
        }
    }

    myStructure->Info();

    /** assign constitutive law **/
    myStructure->ElementTotalSetConstitutiveLaw(myMatLin);
    myStructure->ElementTotalSetSection(mySection);

    /** Boundary condition **/

    // Dirichlet
    NuTo::FullVector<double,Eigen::Dynamic> direction(2);
    direction << 1 ,0;
    myStructure->ConstraintLinearSetDisplacementNode(0, direction, 0.0);
    direction << 0 ,1;
    myStructure->ConstraintLinearSetDisplacementNode(0, direction, 0.0);

    if (PRINTRESULT) std::cout <<  "Nodes fixed: " << 0 << ", ";

    direction << 1 ,0;
    for(int i = 1; i < numControlPointsY; i++)
    {
        myStructure->ConstraintLinearSetDisplacementNode(i*numControlPointsX, direction, 0.0);
        if (PRINTRESULT) std::cout << i*numControlPointsY << ", ";
    }

    if (PRINTRESULT) std::cout << std::endl;

    // right boundary
    int groupNumberNodesRight = myStructure->GroupCreate("NODES");

    if (PRINTRESULT) std::cout <<  "Pressure on nodes: ";

    for(int i = 1; i <= numControlPointsY; i++)
    {
        myStructure->GroupAddNode(groupNumberNodesRight, i*numControlPointsX - 1);
        if (PRINTRESULT) std::cout << i*numControlPointsY - 1 << ", ";
    }
    if (PRINTRESULT) std::cout << std::endl;

    if (PRINTRESULT) std::cout <<  "Pressure on elements: ";

    int groupNumberElementsRight = myStructure->GroupCreate("ELEMENTS");
    for(int i = 1; i <= numElementsY; i++)
    {
        myStructure->GroupAddElement(groupNumberElementsRight, i*numElementsX - 1);
        if (PRINTRESULT) std::cout << i*numElementsY - 1 << ", ";
    }
    if (PRINTRESULT) std::cout << std::endl;

    //NuTo::FullVector<double, 2> ForceVecRight({Force,0.});
    myStructure->SetNumLoadCases(1);
    //myStructure->LoadSurfaceConstDirectionCreate2D(0, groupNumberElementsRight, groupNumberNodesRight, ForceVecRight);
    myStructure->LoadSurfacePressureCreate2D(0, groupNumberElementsRight, groupNumberNodesRight, -Stress);

    myStructure->CalculateMaximumIndependentSets();
    myStructure->NodeBuildGlobalDofs();

    myStructure->NodeInfo(10);

    return myStructure;
}

NuTo::Structure* buildPlateWithHole2D(double& DisplacementCorrect)
{
     NuTo::Structure* myStructure = new NuTo::Structure(2);

     /*********/
     // Mesh  //
     /*********/

     int noPtsX = 4;
     int noPtsY = 3;

     Eigen::Matrix<Eigen::VectorXd, Eigen::Dynamic, Eigen::Dynamic> controlPts(noPtsX, noPtsY);

     controlPts(0,0) = Eigen::Vector2d(-1,0);
     controlPts(1,0) = Eigen::Vector2d(-1,0.4142135623730951);
     controlPts(2,0) = Eigen::Vector2d(-0.4142135623730951, 1);
     controlPts(3,0) = Eigen::Vector2d(0,1);

     controlPts(0,1) = Eigen::Vector2d(-2.5,0);
     controlPts(1,1) = Eigen::Vector2d(-2.5, 0.75);
     controlPts(2,1) = Eigen::Vector2d(-0.75, 2.5);
     controlPts(3,1) = Eigen::Vector2d(0, 2.5);

     controlPts(0,2) = Eigen::Vector2d(-4,0);
     controlPts(1,2) = Eigen::Vector2d(-4,4);
     controlPts(2,2) = Eigen::Vector2d(-4,4);
     controlPts(3,2) = Eigen::Vector2d(0,4);

     Eigen::VectorXd knotsX(7);
     knotsX << 0, 0, 0, 0.5, 1, 1, 1;
     Eigen::VectorXd knotsY(6);
     knotsY << 0, 0, 0, 1, 1, 1;

     Eigen::Vector2i degree(2,2);

     NuTo::BSplineSurface surface(degree, knotsX, knotsY, controlPts);


     /**************/
     // Structure  //
     /**************/

     std::set<NuTo::Node::eDof> setOfDOFS;
     setOfDOFS.insert(NuTo::Node::COORDINATES);
     setOfDOFS.insert(NuTo::Node::DISPLACEMENTS);

     if (PRINTRESULT) std::cout << "Nodes:\n";
     /** create nodes **/
     int count = 0;
     for(int i = 0; i < surface.GetNumControlPoints(1); i++)
     {
         for(int j = 0; j < surface.GetNumControlPoints(0); j++)
         {
             myStructure->NodeCreateDOFs(count, setOfDOFS, surface.GetControlPoint(j,i));
             if (PRINTRESULT) std::cout << "ID: " << count << ", Coordinate: \n" << surface.GetControlPoint(j,i) << std::endl;
             count++;
         }
     }

     if (PRINTRESULT)    std::cout << std::endl;

     int interpolationType = myStructure->InterpolationTypeCreate("IGA2D");

     std::vector<Eigen::VectorXd> vecKnots;
     vecKnots.push_back(surface.GetKnotVector(0));
     vecKnots.push_back(surface.GetKnotVector(1));

     myStructure->InterpolationTypeAdd(interpolationType,
                                       NuTo::Node::COORDINATES,
                                       NuTo::Interpolation::eTypeOrder::SPLINE,
                                       degree,
                                       vecKnots);

     myStructure->InterpolationTypeAdd(interpolationType,
                                       NuTo::Node::DISPLACEMENTS,
                                       NuTo::Interpolation::eTypeOrder::SPLINE,
                                       degree,
                                       vecKnots);

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
     double Thickness = 2.2;
     int mySection = myStructure->SectionCreate("Plane_Stress");
     myStructure->SectionSetThickness(mySection,Thickness);

     /** create constitutive law **/
     double YoungsModulus = 10.e5;
     double PoissonRatio  = 0.3;
     int myMatLin = myStructure->ConstitutiveLawCreate("LINEAR_ELASTIC_ENGINEERING_STRESS");
     myStructure->ConstitutiveLawSetParameterDouble(myMatLin,NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS,YoungsModulus);
     myStructure->ConstitutiveLawSetParameterDouble(myMatLin,NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO,PoissonRatio);

     myStructure->ElementTotalSetConstitutiveLaw(myMatLin);
     myStructure->ElementTotalSetSection(mySection);

     /**********************/
     // Boundary condition //
     /**********************/

     double Stress = -10;

     //--- Dirichlet ---//

     NuTo::FullVector<double,Eigen::Dynamic> direction(2);
     direction << 0 ,1;

     myStructure->ConstraintLinearSetDisplacementNode(0, direction, 0.0);
     myStructure->ConstraintLinearSetDisplacementNode(4, direction, 0.0);
     myStructure->ConstraintLinearSetDisplacementNode(8, direction, 0.0);

     direction << 1, 0;

     myStructure->ConstraintLinearSetDisplacementNode(3, direction, 0.0);
     myStructure->ConstraintLinearSetDisplacementNode(7, direction, 0.0);
     myStructure->ConstraintLinearSetDisplacementNode(11, direction, 0.0);

     //--- Neumann ---//

     int groupNumberNodesRight = myStructure->GroupCreate("NODES");

     myStructure->GroupAddNode(groupNumberNodesRight, 8);
     myStructure->GroupAddNode(groupNumberNodesRight, 9);
     myStructure->GroupAddNode(groupNumberNodesRight, 10);


     int groupNumberElementsRight = myStructure->GroupCreate("ELEMENTS");
     myStructure->GroupAddElement(groupNumberElementsRight, 0);

     myStructure->SetNumLoadCases(1);
     myStructure->LoadSurfacePressureCreate2D(0, groupNumberElementsRight, groupNumberNodesRight, -Stress);

     myStructure->CalculateMaximumIndependentSets();
     myStructure->NodeBuildGlobalDofs();

     myStructure->NodeInfo(10);

     return myStructure;
}

void solve(NuTo::Structure *myStructure, double solution, double tol = 1.e-6)
{
    myStructure->SolveGlobalSystemStaticElastic();

    int numNodes = myStructure->GetNumNodes();
    double nodeDisp = myStructure->NodeGetNodePtr(numNodes - 1)->GetDisplacement(0);

    if (PRINTRESULT)
    {
        std::cout << "Displacement node numerical: \n";
        for(int i = 0; i < numNodes; i++)  std::cout << i << ": " << myStructure->NodeGetNodePtr(i)->GetDisplacement(0) << std::endl;
    }

    std::cout << "Dimension: " << myStructure->GetDimension() << std::endl;
    std::cout << "Displacement last node numerical:\n" << nodeDisp << std::endl;
    std::cout << "Displacement last node analytical:\n" << solution << std::endl;

    int visualizationGroup = myStructure->GroupGetElementsTotal();
    myStructure->AddVisualizationComponent(visualizationGroup, NuTo::VisualizeBase::DISPLACEMENTS);
    myStructure->AddVisualizationComponent(visualizationGroup, NuTo::VisualizeBase::ENGINEERING_STRAIN);
    myStructure->AddVisualizationComponent(visualizationGroup, NuTo::VisualizeBase::ENGINEERING_STRESS);

    std::string resultDir = "./ResultsIGA";
    boost::filesystem::create_directory(resultDir);

    myStructure->ExportVtkDataFileElements(resultDir+"/Elements.vtu", true);
    myStructure->ExportVtkDataFileNodes(resultDir+"/Nodes.vtu", true);


    if (fabs(nodeDisp - solution)/fabs(solution) > tol)
    {
        throw NuTo::Exception("[LobattoIntegration] : displacement is not correct");
    }

}


int main()
{
    double solution = 0;
    NuTo::Structure* myStructure = NULL;

    // 1D
    myStructure = buildStructure1D(solution);
    solve(myStructure, solution);

    // 2D
    myStructure = buildRect2D(solution);
    solve(myStructure, solution);

    // 2D
//    myStructure = buildPlateWithHole2D(solution);
//    solve(myStructure, solution);

    return 0;
}
