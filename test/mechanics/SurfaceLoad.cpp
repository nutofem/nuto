#include <iomanip>

#include "mechanics/structures/StructureOutputBlockVector.h"
#include "mechanics/structures/unstructured/Structure.h"
#include "mechanics/MechanicsException.h"
#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/sections/SectionEnum.h"
#include "mechanics/interpolationtypes/InterpolationTypeEnum.h"

#include "math/SparseDirectSolverMUMPS.h"
#include "math/SparseMatrixCSRVector2General.h"

#define PRINTRESULT false


double pressure = 42.;

void CheckHydrostaticPressure(NuTo::Structure& rStructure)
{
    auto extLoadVector = rStructure.BuildGlobalExternalLoadVector(0).J[NuTo::Node::eDof::DISPLACEMENTS];

    int dimension = rStructure.GetDimension();
    int numNodes = rStructure.GetNumNodes();
    if (extLoadVector.rows() != numNodes*dimension)
        throw NuTo::MechanicsException("[SurfaceLoad::CheckHydrostaticPressure] F_ext.GetNumRows() != dimension * numNodes. Maybe you (mistakenly) applied constraints.");

    Eigen::VectorXd resForce = Eigen::VectorXd::Zero(dimension);
    for (int iNode = 0; iNode < numNodes; ++iNode)
        resForce += extLoadVector.block(dimension*iNode,0, dimension,1);

    if (extLoadVector.cwiseAbs().maxCoeff() < 1.e-6)
        throw NuTo::MechanicsException("[SurfaceLoad::CheckHydrostaticPressure] No external load at all!");

    double relativeError = resForce.cwiseAbs().maxCoeff() / extLoadVector.cwiseAbs().maxCoeff();
    std::cout << "Relative error: " << std::setw(10) << relativeError << "\t";
    if (relativeError > 1.e-8)
    {
        std::cout << "Result Load != 0:" << std::endl;
        std::cout << resForce << std::endl;
        std::cout << "Total external load vector:" << std::endl;
        std::cout << extLoadVector << std::endl;
        throw NuTo::MechanicsException("[SurfaceLoad::CheckHydrostaticPressure] The resulting forces should cancel out. But they don't.");
    }
}

void HydrostaticPressureTriangle2D(NuTo::Interpolation::eTypeOrder rInterpolationDisp)
{
    NuTo::Structure myStructure(2);
    myStructure.SetShowTime(false);

    Eigen::MatrixXd nodeCoords(2,3), nodeCoordsRot(2,3);

    nodeCoords << 0,2,2,
                  0,0,3;

    double angle = M_PI / 6;
    auto rotationMatrix = Eigen::Rotation2D<double>(angle);
    nodeCoordsRot = rotationMatrix.toRotationMatrix() * nodeCoords;

    std::vector<int> nodeIds = myStructure.NodesCreate(nodeCoordsRot);

    std::vector<int> surfaces(3);

    surfaces[0] = myStructure.GroupCreate("Nodes");
    myStructure.GroupAddNode(surfaces[0], nodeIds[0]);
    myStructure.GroupAddNode(surfaces[0], nodeIds[1]);

    surfaces[1] = myStructure.GroupCreate("Nodes");
    myStructure.GroupAddNode(surfaces[1], nodeIds[1]);
    myStructure.GroupAddNode(surfaces[1], nodeIds[2]);

    surfaces[2] = myStructure.GroupCreate("Nodes");
    myStructure.GroupAddNode(surfaces[2], nodeIds[2]);
    myStructure.GroupAddNode(surfaces[2], nodeIds[0]);

    myStructure.InterpolationTypeCreate(0, NuTo::Interpolation::eShapeType::TRIANGLE2D);

    myStructure.InterpolationTypeAdd(0, NuTo::Node::eDof::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    myStructure.InterpolationTypeAdd(0, NuTo::Node::eDof::DISPLACEMENTS, rInterpolationDisp);

    int elementId = myStructure.ElementCreate(0, nodeIds);
    int elementGroup = myStructure.GroupCreate("Elements");
    myStructure.GroupAddElement(elementGroup, elementId);
    myStructure.ElementConvertToInterpolationType(elementGroup);

    for (int iSurface = 0; iSurface < 3; ++iSurface)
        myStructure.LoadSurfacePressureCreate2D(0, elementGroup, surfaces[iSurface], pressure);

    int section = myStructure.SectionCreate(NuTo::eSectionType::PLANE_STRESS);
    myStructure.SectionSetThickness(section, 13.);
    myStructure.ElementTotalSetSection(section);

    CheckHydrostaticPressure(myStructure);
    std::cout << "[SurfaceLoad::HydrostaticPressureTriangle2D] " + NuTo::Interpolation::TypeOrderToString(rInterpolationDisp) + " done." << std::endl;
}


void HydrostaticPressureQuad2D(NuTo::Interpolation::eTypeOrder rInterpolationDisp)
{
    NuTo::Structure myStructure(2);
    myStructure.SetShowTime(false);

    Eigen::MatrixXd nodeCoords(2,4), nodeCoordsRot(2,4);

    nodeCoords << 0,2,2,0,
                  0,0,3,2;

    double angle = M_PI / 6;
    auto rotationMatrix = Eigen::Rotation2D<double>(angle);
    nodeCoordsRot = rotationMatrix.toRotationMatrix() * nodeCoords;

    std::vector<int> nodeIds = myStructure.NodesCreate(nodeCoordsRot);

    std::vector<int> surfaces(4);

    surfaces[0] = myStructure.GroupCreate("Nodes");
    myStructure.GroupAddNode(surfaces[0], nodeIds[0]);
    myStructure.GroupAddNode(surfaces[0], nodeIds[1]);

    surfaces[1] = myStructure.GroupCreate("Nodes");
    myStructure.GroupAddNode(surfaces[1], nodeIds[1]);
    myStructure.GroupAddNode(surfaces[1], nodeIds[2]);

    surfaces[2] = myStructure.GroupCreate("Nodes");
    myStructure.GroupAddNode(surfaces[2], nodeIds[2]);
    myStructure.GroupAddNode(surfaces[2], nodeIds[3]);

    surfaces[3] = myStructure.GroupCreate("Nodes");
    myStructure.GroupAddNode(surfaces[3], nodeIds[3]);
    myStructure.GroupAddNode(surfaces[3], nodeIds[0]);

    myStructure.InterpolationTypeCreate(0, NuTo::Interpolation::eShapeType::QUAD2D);
    myStructure.InterpolationTypeAdd(0, NuTo::Node::eDof::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    myStructure.InterpolationTypeAdd(0, NuTo::Node::eDof::DISPLACEMENTS, rInterpolationDisp);

    int elementId = myStructure.ElementCreate(0, nodeIds);
    int elementGroup = myStructure.GroupCreate("Elements");
    myStructure.GroupAddElement(elementGroup, elementId);
    myStructure.ElementConvertToInterpolationType(elementGroup);

    for (int iSurface = 0; iSurface < 4; ++iSurface)
        myStructure.LoadSurfacePressureCreate2D(0, elementGroup, surfaces[iSurface], pressure);

    int section = myStructure.SectionCreate(NuTo::eSectionType::PLANE_STRESS);
    myStructure.SectionSetThickness(section, 13.);
    myStructure.ElementTotalSetSection(section);

    CheckHydrostaticPressure(myStructure);
    std::cout << "[SurfaceLoad::HydrostaticPressureQuad2D] " + NuTo::Interpolation::TypeOrderToString(rInterpolationDisp) + " done." << std::endl;
}

void HydrostaticPressureTetrahedron3D(NuTo::Interpolation::eTypeOrder rInterpolationDisp)
{
    NuTo::Structure myStructure(3);
    myStructure.SetShowTime(false);

    Eigen::MatrixXd nodeCoords(3,4), nodeCoordsRot(3,4);

    nodeCoords << 0,2,0,0,
                  1,0,3,2,
                  0,0,0,4;

    Eigen::Matrix3d rotationMatrix(
            Eigen::AngleAxisd( M_PI/3., Eigen::Vector3d::UnitX())
          * Eigen::AngleAxisd(-M_PI/5., Eigen::Vector3d::UnitY())
          * Eigen::AngleAxisd( M_PI/7., Eigen::Vector3d::UnitZ()));

    nodeCoordsRot = rotationMatrix * nodeCoords;

    std::vector<int> nodeIds = myStructure.NodesCreate(nodeCoordsRot);

    std::vector<int> surfaces(4);

    surfaces[0] = myStructure.GroupCreate("Nodes");
    myStructure.GroupAddNode(surfaces[0], nodeIds[0]);
    myStructure.GroupAddNode(surfaces[0], nodeIds[1]);
    myStructure.GroupAddNode(surfaces[0], nodeIds[2]);

    surfaces[1] = myStructure.GroupCreate("Nodes");
    myStructure.GroupAddNode(surfaces[1], nodeIds[0]);
    myStructure.GroupAddNode(surfaces[1], nodeIds[1]);
    myStructure.GroupAddNode(surfaces[1], nodeIds[3]);


    surfaces[2] = myStructure.GroupCreate("Nodes");
    myStructure.GroupAddNode(surfaces[2], nodeIds[0]);
    myStructure.GroupAddNode(surfaces[2], nodeIds[2]);
    myStructure.GroupAddNode(surfaces[2], nodeIds[3]);


    surfaces[3] = myStructure.GroupCreate("Nodes");
    myStructure.GroupAddNode(surfaces[3], nodeIds[1]);
    myStructure.GroupAddNode(surfaces[3], nodeIds[2]);
    myStructure.GroupAddNode(surfaces[3], nodeIds[3]);


    myStructure.InterpolationTypeCreate(0, NuTo::Interpolation::eShapeType::TETRAHEDRON3D);
    myStructure.InterpolationTypeAdd(0, NuTo::Node::eDof::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    myStructure.InterpolationTypeAdd(0, NuTo::Node::eDof::DISPLACEMENTS, rInterpolationDisp);

    int elementId = myStructure.ElementCreate(0, nodeIds);
    int elementGroup = myStructure.GroupCreate("Elements");
    myStructure.GroupAddElement(elementGroup, elementId);
    myStructure.ElementConvertToInterpolationType(elementGroup);

    for (int iSurface = 0; iSurface < 4; ++iSurface)
        myStructure.LoadSurfacePressureCreate3D(0, elementGroup, surfaces[iSurface], pressure);

    CheckHydrostaticPressure(myStructure);
    std::cout << "[SurfaceLoad::HydrostaticPressureTetrahedron3D] " + NuTo::Interpolation::TypeOrderToString(rInterpolationDisp) + " done." << std::endl;
}

void HydrostaticPressureBrick3D(NuTo::Interpolation::eTypeOrder rInterpolationDisp)
{
    NuTo::Structure myStructure(3);
    myStructure.SetShowTime(false);

    Eigen::MatrixXd nodeCoords(3,8), nodeCoordsRot(3,8);

    nodeCoords <<
            1,3,3,0, 1,3,3,0,
            0,0,4,4, 0,0,4,4,
            1,0,1,0, 5,5,5,5;

    Eigen::Matrix3d rotationMatrix(
            Eigen::AngleAxisd( M_PI/3., Eigen::Vector3d::UnitX())
          * Eigen::AngleAxisd(-M_PI/5., Eigen::Vector3d::UnitY())
          * Eigen::AngleAxisd( M_PI/7., Eigen::Vector3d::UnitZ()));

    nodeCoordsRot = rotationMatrix * nodeCoords;

    std::vector<int> nodeIds = myStructure.NodesCreate(nodeCoordsRot);

    std::vector<int> surfaces(6);

    surfaces[0] = myStructure.GroupCreate("Nodes");
    myStructure.GroupAddNode(surfaces[0], nodeIds[0]);
    myStructure.GroupAddNode(surfaces[0], nodeIds[1]);
    myStructure.GroupAddNode(surfaces[0], nodeIds[2]);
    myStructure.GroupAddNode(surfaces[0], nodeIds[3]);

    surfaces[1] = myStructure.GroupCreate("Nodes");
    myStructure.GroupAddNode(surfaces[1], nodeIds[4]);
    myStructure.GroupAddNode(surfaces[1], nodeIds[5]);
    myStructure.GroupAddNode(surfaces[1], nodeIds[6]);
    myStructure.GroupAddNode(surfaces[1], nodeIds[7]);

    surfaces[2] = myStructure.GroupCreate("Nodes");
    myStructure.GroupAddNode(surfaces[2], nodeIds[0]);
    myStructure.GroupAddNode(surfaces[2], nodeIds[1]);
    myStructure.GroupAddNode(surfaces[2], nodeIds[5]);
    myStructure.GroupAddNode(surfaces[2], nodeIds[4]);

    surfaces[3] = myStructure.GroupCreate("Nodes");
    myStructure.GroupAddNode(surfaces[3], nodeIds[2]);
    myStructure.GroupAddNode(surfaces[3], nodeIds[3]);
    myStructure.GroupAddNode(surfaces[3], nodeIds[7]);
    myStructure.GroupAddNode(surfaces[3], nodeIds[6]);

    surfaces[4] = myStructure.GroupCreate("Nodes");
    myStructure.GroupAddNode(surfaces[4], nodeIds[1]);
    myStructure.GroupAddNode(surfaces[4], nodeIds[2]);
    myStructure.GroupAddNode(surfaces[4], nodeIds[6]);
    myStructure.GroupAddNode(surfaces[4], nodeIds[5]);

    surfaces[5] = myStructure.GroupCreate("Nodes");
    myStructure.GroupAddNode(surfaces[5], nodeIds[0]);
    myStructure.GroupAddNode(surfaces[5], nodeIds[3]);
    myStructure.GroupAddNode(surfaces[5], nodeIds[7]);
    myStructure.GroupAddNode(surfaces[5], nodeIds[4]);



    myStructure.InterpolationTypeCreate(0, NuTo::Interpolation::eShapeType::BRICK3D);
    myStructure.InterpolationTypeAdd(0, NuTo::Node::eDof::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    myStructure.InterpolationTypeAdd(0, NuTo::Node::eDof::DISPLACEMENTS, rInterpolationDisp);

    int elementId = myStructure.ElementCreate(0, nodeIds);
    int elementGroup = myStructure.GroupCreate("Elements");
    myStructure.GroupAddElement(elementGroup, elementId);
    myStructure.ElementConvertToInterpolationType(elementGroup);

    for (int iSurface = 0; iSurface < 6; ++iSurface)
        myStructure.LoadSurfacePressureCreate3D(0, elementGroup, surfaces[iSurface], pressure);


    CheckHydrostaticPressure(myStructure);
    std::cout << "[SurfaceLoad::HydrostaticPressureBrick3D] " + NuTo::Interpolation::TypeOrderToString(rInterpolationDisp) + " done." << std::endl;
}

void HydrostaticPressurePrism3D(NuTo::Interpolation::eTypeOrder rInterpolationDisp)
{
    NuTo::Structure myStructure(3);
    myStructure.SetShowTime(false);

    Eigen::MatrixXd nodeCoords(3,6), nodeCoordsRot(3,6);

    nodeCoords <<
        0,3,0,0,3,0,
        0,0,4,0,0,4,
        0,0,0,5,5,5;

    Eigen::Matrix3d rotationMatrix(
        Eigen::AngleAxisd( M_PI/3., Eigen::Vector3d::UnitX())
            * Eigen::AngleAxisd(-M_PI/5., Eigen::Vector3d::UnitY())
            * Eigen::AngleAxisd( M_PI/7., Eigen::Vector3d::UnitZ()));

    nodeCoordsRot = rotationMatrix * nodeCoords;

    std::vector<int> nodeIds = myStructure.NodesCreate(nodeCoordsRot);

    std::vector<int> surfaces(5);

    surfaces[0] = myStructure.GroupCreate("Nodes");
    myStructure.GroupAddNode(surfaces[0], nodeIds[0]);
    myStructure.GroupAddNode(surfaces[0], nodeIds[1]);
    myStructure.GroupAddNode(surfaces[0], nodeIds[2]);

    surfaces[1] = myStructure.GroupCreate("Nodes");
    myStructure.GroupAddNode(surfaces[1], nodeIds[3]);
    myStructure.GroupAddNode(surfaces[1], nodeIds[4]);
    myStructure.GroupAddNode(surfaces[1], nodeIds[5]);

    surfaces[2] = myStructure.GroupCreate("Nodes");
    myStructure.GroupAddNode(surfaces[2], nodeIds[1]);
    myStructure.GroupAddNode(surfaces[2], nodeIds[2]);
    myStructure.GroupAddNode(surfaces[2], nodeIds[5]);
    myStructure.GroupAddNode(surfaces[2], nodeIds[4]);

    surfaces[3] = myStructure.GroupCreate("Nodes");
    myStructure.GroupAddNode(surfaces[3], nodeIds[0]);
    myStructure.GroupAddNode(surfaces[3], nodeIds[2]);
    myStructure.GroupAddNode(surfaces[3], nodeIds[5]);
    myStructure.GroupAddNode(surfaces[3], nodeIds[3]);

    surfaces[4] = myStructure.GroupCreate("Nodes");
    myStructure.GroupAddNode(surfaces[4], nodeIds[1]);
    myStructure.GroupAddNode(surfaces[4], nodeIds[0]);
    myStructure.GroupAddNode(surfaces[4], nodeIds[3]);
    myStructure.GroupAddNode(surfaces[4], nodeIds[4]);

    myStructure.InterpolationTypeCreate(0, NuTo::Interpolation::eShapeType::PRISM3D);
    myStructure.InterpolationTypeAdd(0, NuTo::Node::eDof::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    myStructure.InterpolationTypeAdd(0, NuTo::Node::eDof::DISPLACEMENTS, rInterpolationDisp);

    int elementId = myStructure.ElementCreate(0, nodeIds);
    int elementGroup = myStructure.GroupCreate("Elements");
    myStructure.GroupAddElement(elementGroup, elementId);
    myStructure.ElementConvertToInterpolationType(elementGroup);

    for (int surface : surfaces)
        myStructure.LoadSurfacePressureCreate3D(0, elementGroup, surface, pressure);


    CheckHydrostaticPressure(myStructure);
    std::cout << "[SurfaceLoad::HydrostaticPressurePrism3D] " + NuTo::Interpolation::TypeOrderToString(rInterpolationDisp) + " done." << std::endl;
}

void CheckSurfaceLoad(NuTo::Structure& rStructure, const Eigen::VectorXd& rLoad, double rSurfaceArea)
{
    auto extLoadVector = rStructure.BuildGlobalExternalLoadVector(0).J[NuTo::Node::eDof::DISPLACEMENTS];

    int dimension = rStructure.GetDimension();
    int numNodes = rStructure.GetNumNodes();
    if (extLoadVector.rows() != numNodes*dimension)
        throw NuTo::MechanicsException("[SurfaceLoad::CheckSurfaceLoad] F_ext.GetNumRows() != dimension * numNodes. Maybe you (mistakenly) applied constraints.");

    Eigen::VectorXd resForce = Eigen::VectorXd::Zero(dimension);

    for (int iNode = 0; iNode < numNodes; ++iNode)
        resForce += extLoadVector.block(dimension*iNode,0, dimension,1);


    Eigen::VectorXd resForceCorrect = rLoad * rSurfaceArea;

    double relativeError = (resForce - resForceCorrect).cwiseAbs().maxCoeff() / extLoadVector.cwiseAbs().maxCoeff();
    std::cout << "Relative error: " << std::setw(10) << relativeError << "\t";
    if (relativeError > 1.e-6)
    {
        std::cout << "Result Load:" << std::endl;
        std::cout << resForce << std::endl;
        std::cout << "Result Load correct:" << std::endl;
        std::cout << resForceCorrect << std::endl;
        std::cout << "Total external load vector:" << std::endl;
        std::cout << extLoadVector << std::endl;
        throw NuTo::MechanicsException("[SurfaceLoad::CheckSurfaceLoad] Wrong external load calculation!");
    }
}

void SurfaceLoadTriangle2D(NuTo::Interpolation::eTypeOrder rInterpolationDisp)
{
    NuTo::Structure myStructure(2);
    myStructure.SetShowTime(false);

    Eigen::MatrixXd nodeCoords(2,3);

    double lx = 2;
    double ly = 5;
    double thickness = 13.;

    nodeCoords <<
            0, lx, 0,
           -ly, 0,ly;

    std::vector<int> nodeIds = myStructure.NodesCreate(nodeCoords);

    int s = myStructure.GroupCreate("Nodes");
    myStructure.GroupAddNode(s, nodeIds[0]);
    myStructure.GroupAddNode(s, nodeIds[1]);

    double surfaceArea = thickness*std::sqrt(lx*lx + ly*ly);

    myStructure.InterpolationTypeCreate(0, NuTo::Interpolation::eShapeType::TRIANGLE2D);
    myStructure.InterpolationTypeAdd(0, NuTo::Node::eDof::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    myStructure.InterpolationTypeAdd(0, NuTo::Node::eDof::DISPLACEMENTS, rInterpolationDisp);

    int elementId = myStructure.ElementCreate(0, nodeIds);
    int elementGroup = myStructure.GroupCreate("Elements");
    myStructure.GroupAddElement(elementGroup, elementId);
    myStructure.ElementConvertToInterpolationType(elementGroup);

    int section = myStructure.SectionCreate(NuTo::eSectionType::PLANE_STRESS);
    myStructure.SectionSetThickness(section, thickness);
    myStructure.ElementTotalSetSection(section);

    Eigen::Vector2d load({42., -M_PI});

    myStructure.LoadSurfaceConstDirectionCreate2D(0, elementGroup, s, load);

    CheckSurfaceLoad(myStructure, load, surfaceArea);
    std::cout << "[SurfaceLoad::SurfaceLoadTriangle2D] " + NuTo::Interpolation::TypeOrderToString(rInterpolationDisp) + " done." << std::endl;
}


void SurfaceLoadQuad2D(NuTo::Interpolation::eTypeOrder rInterpolationDisp)
{
    NuTo::Structure myStructure(2);
    myStructure.SetShowTime(false);

    Eigen::MatrixXd nodeCoords(2,4);

    double lx = 2;
    double ly = 5;
    double thickness = 13.;

    nodeCoords <<
            0, lx, 0, -lx,
           -ly, 0,ly, 0;

    std::vector<int> nodeIds = myStructure.NodesCreate(nodeCoords);

    int s = myStructure.GroupCreate("Nodes");
    myStructure.GroupAddNode(s, nodeIds[0]);
    myStructure.GroupAddNode(s, nodeIds[1]);

    double surfaceArea = thickness*std::sqrt(lx*lx + ly*ly);

    myStructure.InterpolationTypeCreate(0, NuTo::Interpolation::eShapeType::QUAD2D);
    myStructure.InterpolationTypeAdd(0, NuTo::Node::eDof::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    myStructure.InterpolationTypeAdd(0, NuTo::Node::eDof::DISPLACEMENTS, rInterpolationDisp);

    int elementId = myStructure.ElementCreate(0, nodeIds);
    int elementGroup = myStructure.GroupCreate("Elements");
    myStructure.GroupAddElement(elementGroup, elementId);
    myStructure.ElementConvertToInterpolationType(elementGroup);

    int section = myStructure.SectionCreate(NuTo::eSectionType::PLANE_STRESS);
    myStructure.SectionSetThickness(section, thickness);
    myStructure.ElementTotalSetSection(section);

    Eigen::Vector2d load({42., -M_PI});

    myStructure.LoadSurfaceConstDirectionCreate2D(0, elementGroup, s, load);

    CheckSurfaceLoad(myStructure, load, surfaceArea);
    std::cout << "[SurfaceLoad::SurfaceLoadQuad2D] " + NuTo::Interpolation::TypeOrderToString(rInterpolationDisp) + " done." << std::endl;
}

void SurfaceLoadQuad2DIGA(NuTo::Interpolation::eTypeOrder rInterpolationDisp)
{
    NuTo::Structure myStructure(2);
    myStructure.SetShowTime(false);

    int interpolationType = myStructure.InterpolationTypeCreate("IGA2D");

    Eigen::MatrixXd nodeCoords(2,4);

    double lx = 2;
    double ly = 5;
    double thickness = 13.;

    nodeCoords <<
            0, lx, 0, -lx,
           -ly, 0,ly, 0;


    Eigen::Vector2i degree(1,1);

    Eigen::MatrixXd weights(2,2);
    weights << 1, 1, 1, 1;

    Eigen::VectorXd knots(4);
    knots << 0,0,1,1;

    std::vector<Eigen::VectorXd> vecKnots;
    vecKnots.push_back(knots);
    vecKnots.push_back(knots);

    myStructure.InterpolationTypeAdd(interpolationType,
                                     NuTo::Node::eDof::COORDINATES,
                                     NuTo::Interpolation::eTypeOrder::SPLINE,
                                     degree,
                                     vecKnots,
                                     weights);

    myStructure.InterpolationTypeAdd(interpolationType,
                                     NuTo::Node::eDof::DISPLACEMENTS,
                                     NuTo::Interpolation::eTypeOrder::SPLINE,
                                     degree,
                                     vecKnots,
                                     weights);

    std::set<NuTo::Node::eDof> setOfDOFS;
    setOfDOFS.insert(NuTo::Node::eDof::COORDINATES);
    setOfDOFS.insert(NuTo::Node::eDof::DISPLACEMENTS);

    for(int i = 0; i < nodeCoords.cols(); i++)
        myStructure.NodeCreateDOFs(0, setOfDOFS, nodeCoords.col(i));

    Eigen::Matrix2d elementKnots(2,2);
    elementKnots << 0, 1, 0, 1;

    Eigen::Vector2i elementKnotIDs(2,2);

    Eigen::VectorXi elementIncidence(4);
    elementIncidence << 0, 1, 2, 3;

    myStructure.ElementCreate(interpolationType,
                               elementIncidence,
                               elementKnots,
                               elementKnotIDs);

    int section = myStructure.SectionCreate(NuTo::eSectionType::PLANE_STRESS);
    myStructure.SectionSetThickness(section, thickness);
    myStructure.ElementTotalSetSection(section);
}

void SurfaceLoadTetrahedron3D(NuTo::Interpolation::eTypeOrder rInterpolationDisp)
{
    NuTo::Structure myStructure(3);
    myStructure.SetShowTime(false);

    Eigen::MatrixXd nodeCoords(3,4), nodeCoordsRot(3,4);

    double lx = 2;
    double ly = 5;
    double lz = 12;

    nodeCoords <<
            0,lx, 0, 0,
            0, 0,ly, 0,
            0, 0, 0,lz;

    // rotate randomly
    Eigen::Matrix3d rotationMatrix(
            Eigen::AngleAxisd( M_PI/2.9, Eigen::Vector3d::UnitX())
          * Eigen::AngleAxisd(-M_PI/15., Eigen::Vector3d::UnitY())
          * Eigen::AngleAxisd( M_PI/7.5, Eigen::Vector3d::UnitZ()));

    nodeCoordsRot = rotationMatrix * nodeCoords;

    std::vector<int> nodeIds = myStructure.NodesCreate(nodeCoordsRot);

    int s = myStructure.GroupCreate("Nodes");
    myStructure.GroupAddNode(s, nodeIds[0]);
    myStructure.GroupAddNode(s, nodeIds[1]);
    myStructure.GroupAddNode(s, nodeIds[3]);

    double surfaceArea = .5*lx*lz;

    myStructure.InterpolationTypeCreate(0, NuTo::Interpolation::eShapeType::TETRAHEDRON3D);
    myStructure.InterpolationTypeAdd(0, NuTo::Node::eDof::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    myStructure.InterpolationTypeAdd(0, NuTo::Node::eDof::DISPLACEMENTS, rInterpolationDisp);

    int elementId = myStructure.ElementCreate(0, nodeIds);
    int elementGroup = myStructure.GroupCreate("Elements");
    myStructure.GroupAddElement(elementGroup, elementId);
    myStructure.ElementConvertToInterpolationType(elementGroup);

    Eigen::Vector3d load({42., -M_PI, 1.e3});

    myStructure.LoadSurfaceConstDirectionCreate3D(0, elementGroup, s, load);

    CheckSurfaceLoad(myStructure, load, surfaceArea);
    std::cout << "[SurfaceLoad::SurfaceLoadTetrahedron3D] " + NuTo::Interpolation::TypeOrderToString(rInterpolationDisp) + " done." << std::endl;
}

void SurfaceLoadBrick3D(NuTo::Interpolation::eTypeOrder rInterpolationDisp)
{
    NuTo::Structure myStructure(3);
    myStructure.SetShowTime(false);

    Eigen::MatrixXd nodeCoords(3,8), nodeCoordsRot(3,8);

    double lx = 2;
    double ly = 5;
    double lz = 12;

    nodeCoords <<
            0,lx,lx, 0,  0,lx,lx, 0,
            0, 0,ly,ly,  0, 0,ly,ly,
            0, 0, 0, 0, lz,lz,lz,lz;

    // rotate randomly
    Eigen::Matrix3d rotationMatrix(
            Eigen::AngleAxisd( M_PI/2.9, Eigen::Vector3d::UnitX())
          * Eigen::AngleAxisd(-M_PI/15., Eigen::Vector3d::UnitY())
          * Eigen::AngleAxisd( M_PI/7.5, Eigen::Vector3d::UnitZ()));

    nodeCoordsRot = rotationMatrix * nodeCoords;

    std::vector<int> nodeIds = myStructure.NodesCreate(nodeCoordsRot);

    int s = myStructure.GroupCreate("Nodes");
    myStructure.GroupAddNode(s, nodeIds[0]);
    myStructure.GroupAddNode(s, nodeIds[1]);
    myStructure.GroupAddNode(s, nodeIds[5]);
    myStructure.GroupAddNode(s, nodeIds[4]);

    double surfaceArea = lx*lz;

    myStructure.InterpolationTypeCreate(0, NuTo::Interpolation::eShapeType::BRICK3D);
    myStructure.InterpolationTypeAdd(0, NuTo::Node::eDof::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    myStructure.InterpolationTypeAdd(0, NuTo::Node::eDof::DISPLACEMENTS, rInterpolationDisp);

    int elementId = myStructure.ElementCreate(0, nodeIds);
    int elementGroup = myStructure.GroupCreate("Elements");
    myStructure.GroupAddElement(elementGroup, elementId);
    myStructure.ElementConvertToInterpolationType(elementGroup);

    Eigen::Vector3d load({42., -M_PI, 4.2e6});

    myStructure.LoadSurfaceConstDirectionCreate3D(0, elementGroup, s, load);

    CheckSurfaceLoad(myStructure, load, surfaceArea);
    std::cout << "[SurfaceLoad::SurfaceLoadBrick3D] " + NuTo::Interpolation::TypeOrderToString(rInterpolationDisp) + " done." << std::endl;
}

int main()
{

    // 1) all surfaces need testing!
    // 2) pressure (normal) and surface load (certain direction) to be tested
    // 3) not axis aligned
    // 4) analytical solution as reference
    try
    {
        HydrostaticPressureTriangle2D(      NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
        HydrostaticPressureTriangle2D(      NuTo::Interpolation::eTypeOrder::EQUIDISTANT2);
        HydrostaticPressureTriangle2D(      NuTo::Interpolation::eTypeOrder::EQUIDISTANT3);
        HydrostaticPressureTriangle2D(      NuTo::Interpolation::eTypeOrder::EQUIDISTANT4);

        HydrostaticPressureQuad2D(          NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
        HydrostaticPressureQuad2D(          NuTo::Interpolation::eTypeOrder::EQUIDISTANT2);
        HydrostaticPressureQuad2D(          NuTo::Interpolation::eTypeOrder::LOBATTO2);
        HydrostaticPressureQuad2D(          NuTo::Interpolation::eTypeOrder::LOBATTO3);
        HydrostaticPressureQuad2D(          NuTo::Interpolation::eTypeOrder::LOBATTO4);

        HydrostaticPressureTetrahedron3D(   NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
        HydrostaticPressureTetrahedron3D(   NuTo::Interpolation::eTypeOrder::EQUIDISTANT2);

        //HydrostaticPressurePrism3D(NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
        //HydrostaticPressurePrism3D(         NuTo::Interpolation::eTypeOrder::EQUIDISTANT2);

        HydrostaticPressureBrick3D(         NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
        HydrostaticPressureBrick3D(         NuTo::Interpolation::eTypeOrder::EQUIDISTANT2);
        HydrostaticPressureBrick3D(         NuTo::Interpolation::eTypeOrder::LOBATTO2);
        HydrostaticPressureBrick3D(         NuTo::Interpolation::eTypeOrder::LOBATTO3);
        HydrostaticPressureBrick3D(         NuTo::Interpolation::eTypeOrder::LOBATTO4);

        SurfaceLoadTriangle2D(              NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
        SurfaceLoadTriangle2D(              NuTo::Interpolation::eTypeOrder::EQUIDISTANT2);
        SurfaceLoadTriangle2D(              NuTo::Interpolation::eTypeOrder::EQUIDISTANT3);
        SurfaceLoadTriangle2D(              NuTo::Interpolation::eTypeOrder::EQUIDISTANT4);

        SurfaceLoadQuad2D(                  NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
        SurfaceLoadQuad2D(                  NuTo::Interpolation::eTypeOrder::EQUIDISTANT2);
        SurfaceLoadQuad2D(                  NuTo::Interpolation::eTypeOrder::LOBATTO2);
        SurfaceLoadQuad2D(                  NuTo::Interpolation::eTypeOrder::LOBATTO3);
        SurfaceLoadQuad2D(                  NuTo::Interpolation::eTypeOrder::LOBATTO4);

        SurfaceLoadTetrahedron3D(           NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
        SurfaceLoadTetrahedron3D(           NuTo::Interpolation::eTypeOrder::EQUIDISTANT2);

        SurfaceLoadBrick3D(                 NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
        SurfaceLoadBrick3D(                 NuTo::Interpolation::eTypeOrder::EQUIDISTANT2);
        SurfaceLoadBrick3D(                 NuTo::Interpolation::eTypeOrder::LOBATTO2);
        SurfaceLoadBrick3D(                 NuTo::Interpolation::eTypeOrder::LOBATTO3);
        SurfaceLoadBrick3D(                 NuTo::Interpolation::eTypeOrder::LOBATTO4);
    }
    catch (NuTo::Exception& e)
    {
        std::cout << e.ErrorMessage() << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
