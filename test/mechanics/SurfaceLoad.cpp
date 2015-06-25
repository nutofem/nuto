#include "nuto/math/FullMatrix.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/mechanics/MechanicsException.h"

#include "nuto/math/SparseDirectSolverMUMPS.h"
#include "nuto/math/SparseMatrixCSRVector2General.h"

#define PRINTRESULT false


double pressure = 42.;

void CheckHydrostaticPressure(NuTo::Structure& rStructure)
{
    NuTo::FullVector<double, Eigen::Dynamic> extLoadVector;
    rStructure.BuildGlobalExternalLoadVector(0, extLoadVector);

    int dimension = rStructure.GetDimension();
    int numNodes = rStructure.GetNumNodes();
    if (extLoadVector.GetNumRows() != numNodes*dimension)
        throw NuTo::MechanicsException("[SurfaceLoad::CheckHydrostaticPressure] F_ext.GetNumRows() != dimension * numNodes. Maybe you (mistakenly) applied constraints.");

    NuTo::FullVector<double, Eigen::Dynamic> resForce(dimension);

    resForce.setZero();
    for (int iNode = 0; iNode < numNodes; ++iNode)
        resForce += extLoadVector.block(dimension*iNode,0, dimension,1);

    if (extLoadVector.Abs().Max() < 1.e-8)
        throw NuTo::MechanicsException("[SurfaceLoad::CheckHydrostaticPressure] No external load at all!");

    if (resForce.Abs().Max() > 1.e-8)
    {
        std::cout << "Resultant Load != 0:" << std::endl;
        resForce.Info();
        std::cout << "Total external load vector:" << std::endl;
        extLoadVector.Info();
        throw NuTo::MechanicsException("[SurfaceLoad::CheckHydrostaticPressure] The resulting forces should cancel out. But they don't.");
    }
}

void HydrostaticPressureTriangle2D(NuTo::Interpolation::eTypeOrder rInterpolationCoords,
                                   NuTo::Interpolation::eTypeOrder rInterpolationDisp)
{
    NuTo::Structure myStructure(2);
    myStructure.SetShowTime(false);

    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> nodeCoords(2,3), nodeCoordsRot(2,3);

    nodeCoords << 0,2,2,
                  0,0,3;

    double angle = M_PI / 6;
    auto rotationMatrix = Eigen::Rotation2D<double>(angle);
    nodeCoordsRot = rotationMatrix.toRotationMatrix() * nodeCoords;

    NuTo::FullMatrix<int, Eigen::Dynamic> nodeIds = myStructure.NodesCreate(nodeCoordsRot);

    NuTo::FullVector<int, 3> s; // surfaces

    s(0) = myStructure.GroupCreate("Nodes");
    myStructure.GroupAddNode(s.GetValue(0), nodeIds.GetValue(0));
    myStructure.GroupAddNode(s.GetValue(0), nodeIds.GetValue(1));

    s(1) = myStructure.GroupCreate("Nodes");
    myStructure.GroupAddNode(s.GetValue(1), nodeIds.GetValue(1));
    myStructure.GroupAddNode(s.GetValue(1), nodeIds.GetValue(2));

    s(2) = myStructure.GroupCreate("Nodes");
    myStructure.GroupAddNode(s.GetValue(2), nodeIds.GetValue(2));
    myStructure.GroupAddNode(s.GetValue(2), nodeIds.GetValue(0));

    myStructure.InterpolationTypeCreate(0, NuTo::Interpolation::TRIANGLE2D);

    myStructure.InterpolationTypeAdd(0, NuTo::Node::COORDINATES, rInterpolationCoords);
    myStructure.InterpolationTypeAdd(0, NuTo::Node::DISPLACEMENTS, rInterpolationDisp);

    int elementGroup = myStructure.ElementsCreate(0,nodeIds);
    myStructure.ElementConvertToInterpolationType(elementGroup);

    for (int iSurface = 0; iSurface < 3; ++iSurface)
        myStructure.LoadSurfacePressureCreate2D(0, elementGroup, s.GetValue(iSurface), pressure);

    int section = myStructure.SectionCreate(NuTo::Section::PLANE_STRESS);
    myStructure.SectionSetThickness(section, 13.);
    myStructure.ElementTotalSetSection(section);

    CheckHydrostaticPressure(myStructure);
    std::cout << "[SurfaceLoad::HydrostaticPressureTriangle2D] done." << std::endl;
}


void HydrostaticPressureQuad2D(NuTo::Interpolation::eTypeOrder rInterpolationCoords,
                               NuTo::Interpolation::eTypeOrder rInterpolationDisp)
{
    NuTo::Structure myStructure(2);
    myStructure.SetShowTime(false);

    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> nodeCoords(2,4), nodeCoordsRot(2,4);

    nodeCoords << 0,2,2,0,
                  0,0,3,2;

    double angle = M_PI / 6;
    auto rotationMatrix = Eigen::Rotation2D<double>(angle);
    nodeCoordsRot = rotationMatrix.toRotationMatrix() * nodeCoords;

    NuTo::FullMatrix<int, Eigen::Dynamic> nodeIds = myStructure.NodesCreate(nodeCoordsRot);

    NuTo::FullVector<int, 4> s; // surfaces

    s(0) = myStructure.GroupCreate("Nodes");
    myStructure.GroupAddNode(s.GetValue(0), nodeIds.GetValue(0));
    myStructure.GroupAddNode(s.GetValue(0), nodeIds.GetValue(1));

    s(1) = myStructure.GroupCreate("Nodes");
    myStructure.GroupAddNode(s.GetValue(1), nodeIds.GetValue(1));
    myStructure.GroupAddNode(s.GetValue(1), nodeIds.GetValue(2));

    s(2) = myStructure.GroupCreate("Nodes");
    myStructure.GroupAddNode(s.GetValue(2), nodeIds.GetValue(2));
    myStructure.GroupAddNode(s.GetValue(2), nodeIds.GetValue(3));

    s(3) = myStructure.GroupCreate("Nodes");
    myStructure.GroupAddNode(s.GetValue(3), nodeIds.GetValue(3));
    myStructure.GroupAddNode(s.GetValue(3), nodeIds.GetValue(0));

    myStructure.InterpolationTypeCreate(0, NuTo::Interpolation::QUAD2D);
    myStructure.InterpolationTypeAdd(0, NuTo::Node::COORDINATES, rInterpolationCoords);
    myStructure.InterpolationTypeAdd(0, NuTo::Node::DISPLACEMENTS, rInterpolationDisp);

    int elementGroup = myStructure.ElementsCreate(0,nodeIds);
    myStructure.ElementConvertToInterpolationType(elementGroup);

    for (int iSurface = 0; iSurface < 4; ++iSurface)
        myStructure.LoadSurfacePressureCreate2D(0, elementGroup, s.GetValue(iSurface), pressure);

    int section = myStructure.SectionCreate(NuTo::Section::PLANE_STRESS);
    myStructure.SectionSetThickness(section, 13.);
    myStructure.ElementTotalSetSection(section);

    CheckHydrostaticPressure(myStructure);
    std::cout << "[SurfaceLoad::HydrostaticPressureQuad2D] done." << std::endl;
}

void HydrostaticPressureTetrahedron3D(NuTo::Interpolation::eTypeOrder rInterpolationCoords,
                                      NuTo::Interpolation::eTypeOrder rInterpolationDisp)
{
    NuTo::Structure myStructure(3);
    myStructure.SetShowTime(false);

    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> nodeCoords(3,4), nodeCoordsRot(3,4);

    nodeCoords << 0,2,0,0,
                  1,0,3,2,
                  0,0,0,4;

    Eigen::Matrix3d rotationMatrix(
            Eigen::AngleAxisd( M_PI/3., Eigen::Vector3d::UnitX())
          * Eigen::AngleAxisd(-M_PI/5., Eigen::Vector3d::UnitY())
          * Eigen::AngleAxisd( M_PI/7., Eigen::Vector3d::UnitZ()));

    nodeCoordsRot = rotationMatrix * nodeCoords;

    NuTo::FullMatrix<int, Eigen::Dynamic> nodeIds = myStructure.NodesCreate(nodeCoordsRot);

    NuTo::FullVector<int, 4> s; // surfaces

    s(0) = myStructure.GroupCreate("Nodes");
    myStructure.GroupAddNode(s.GetValue(0), nodeIds.GetValue(0));
    myStructure.GroupAddNode(s.GetValue(0), nodeIds.GetValue(1));
    myStructure.GroupAddNode(s.GetValue(0), nodeIds.GetValue(2));

    s(1) = myStructure.GroupCreate("Nodes");
    myStructure.GroupAddNode(s.GetValue(1), nodeIds.GetValue(0));
    myStructure.GroupAddNode(s.GetValue(1), nodeIds.GetValue(1));
    myStructure.GroupAddNode(s.GetValue(1), nodeIds.GetValue(3));


    s(2) = myStructure.GroupCreate("Nodes");
    myStructure.GroupAddNode(s.GetValue(2), nodeIds.GetValue(0));
    myStructure.GroupAddNode(s.GetValue(2), nodeIds.GetValue(2));
    myStructure.GroupAddNode(s.GetValue(2), nodeIds.GetValue(3));


    s(3) = myStructure.GroupCreate("Nodes");
    myStructure.GroupAddNode(s.GetValue(3), nodeIds.GetValue(1));
    myStructure.GroupAddNode(s.GetValue(3), nodeIds.GetValue(2));
    myStructure.GroupAddNode(s.GetValue(3), nodeIds.GetValue(3));


    myStructure.InterpolationTypeCreate(0, NuTo::Interpolation::TETRAHEDRON3D);
    myStructure.InterpolationTypeAdd(0, NuTo::Node::COORDINATES, rInterpolationCoords);
    myStructure.InterpolationTypeAdd(0, NuTo::Node::DISPLACEMENTS, rInterpolationDisp);

    int elementGroup = myStructure.ElementsCreate(0,nodeIds);
    myStructure.ElementConvertToInterpolationType(elementGroup);

    for (int iSurface = 0; iSurface < 4; ++iSurface)
        myStructure.LoadSurfacePressureCreate3D(0, elementGroup, s.GetValue(iSurface), pressure);

    CheckHydrostaticPressure(myStructure);
    std::cout << "[SurfaceLoad::HydrostaticPressureTetrahedron3D] done." << std::endl;
}

void HydrostaticPressureBrick3D(NuTo::Interpolation::eTypeOrder rInterpolationCoords,
                                NuTo::Interpolation::eTypeOrder rInterpolationDisp)
{
    NuTo::Structure myStructure(3);
    myStructure.SetShowTime(false);

    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> nodeCoords(3,8), nodeCoordsRot(3,8);

    nodeCoords <<
            1,3,3,0, 1,3,3,0,
            0,0,4,4, 0,0,4,4,
            1,0,1,0, 5,5,5,5;

    Eigen::Matrix3d rotationMatrix(
            Eigen::AngleAxisd( M_PI/3., Eigen::Vector3d::UnitX())
          * Eigen::AngleAxisd(-M_PI/5., Eigen::Vector3d::UnitY())
          * Eigen::AngleAxisd( M_PI/7., Eigen::Vector3d::UnitZ()));

    nodeCoordsRot = rotationMatrix * nodeCoords;

    NuTo::FullMatrix<int, Eigen::Dynamic> nodeIds = myStructure.NodesCreate(nodeCoordsRot);

    NuTo::FullVector<int, 6> s; // surfaces

    s(0) = myStructure.GroupCreate("Nodes");
    myStructure.GroupAddNode(s.GetValue(0), nodeIds.GetValue(0));
    myStructure.GroupAddNode(s.GetValue(0), nodeIds.GetValue(1));
    myStructure.GroupAddNode(s.GetValue(0), nodeIds.GetValue(2));
    myStructure.GroupAddNode(s.GetValue(0), nodeIds.GetValue(3));

    s(1) = myStructure.GroupCreate("Nodes");
    myStructure.GroupAddNode(s.GetValue(1), nodeIds.GetValue(4));
    myStructure.GroupAddNode(s.GetValue(1), nodeIds.GetValue(5));
    myStructure.GroupAddNode(s.GetValue(1), nodeIds.GetValue(6));
    myStructure.GroupAddNode(s.GetValue(1), nodeIds.GetValue(7));

    s(2) = myStructure.GroupCreate("Nodes");
    myStructure.GroupAddNode(s.GetValue(2), nodeIds.GetValue(0));
    myStructure.GroupAddNode(s.GetValue(2), nodeIds.GetValue(1));
    myStructure.GroupAddNode(s.GetValue(2), nodeIds.GetValue(5));
    myStructure.GroupAddNode(s.GetValue(2), nodeIds.GetValue(4));

    s(3) = myStructure.GroupCreate("Nodes");
    myStructure.GroupAddNode(s.GetValue(3), nodeIds.GetValue(2));
    myStructure.GroupAddNode(s.GetValue(3), nodeIds.GetValue(3));
    myStructure.GroupAddNode(s.GetValue(3), nodeIds.GetValue(7));
    myStructure.GroupAddNode(s.GetValue(3), nodeIds.GetValue(6));

    s(4) = myStructure.GroupCreate("Nodes");
    myStructure.GroupAddNode(s.GetValue(4), nodeIds.GetValue(1));
    myStructure.GroupAddNode(s.GetValue(4), nodeIds.GetValue(2));
    myStructure.GroupAddNode(s.GetValue(4), nodeIds.GetValue(6));
    myStructure.GroupAddNode(s.GetValue(4), nodeIds.GetValue(5));

    s(5) = myStructure.GroupCreate("Nodes");
    myStructure.GroupAddNode(s.GetValue(5), nodeIds.GetValue(0));
    myStructure.GroupAddNode(s.GetValue(5), nodeIds.GetValue(3));
    myStructure.GroupAddNode(s.GetValue(5), nodeIds.GetValue(7));
    myStructure.GroupAddNode(s.GetValue(5), nodeIds.GetValue(4));



    myStructure.InterpolationTypeCreate(0, NuTo::Interpolation::BRICK3D);
    myStructure.InterpolationTypeAdd(0, NuTo::Node::COORDINATES, rInterpolationCoords);
    myStructure.InterpolationTypeAdd(0, NuTo::Node::DISPLACEMENTS, rInterpolationDisp);

    int elementGroup = myStructure.ElementsCreate(0,nodeIds);
    myStructure.ElementConvertToInterpolationType(elementGroup);

    for (int iSurface = 0; iSurface < 6; ++iSurface)
        myStructure.LoadSurfacePressureCreate3D(0, elementGroup, s.GetValue(iSurface), pressure);


    CheckHydrostaticPressure(myStructure);
    std::cout << "[SurfaceLoad::HydrostaticPressureBrick3D] done." << std::endl;
}

void CheckSurfaceLoad(NuTo::Structure& rStructure, const NuTo::FullVector<double, Eigen::Dynamic>& rLoad, double rSurfaceArea)
{
    NuTo::FullVector<double, Eigen::Dynamic> extLoadVector;
    rStructure.BuildGlobalExternalLoadVector(0, extLoadVector);

    int dimension = rStructure.GetDimension();
    int numNodes = rStructure.GetNumNodes();
    if (extLoadVector.GetNumRows() != numNodes*dimension)
        throw NuTo::MechanicsException("[SurfaceLoad::CheckSurfaceLoad] F_ext.GetNumRows() != dimension * numNodes. Maybe you (mistakenly) applied constraints.");

    NuTo::FullVector<double, Eigen::Dynamic> resForce(dimension), resForceCorrect(dimension);

    resForce.setZero();
    for (int iNode = 0; iNode < numNodes; ++iNode)
        resForce += extLoadVector.block(dimension*iNode,0, dimension,1);


    resForceCorrect = rLoad * rSurfaceArea;
    if ((resForce - resForceCorrect).cwiseAbs().maxCoeff() > 1.e-8)
    {
        std::cout << "Resultant Load:" << std::endl;
        resForce.Info();
        std::cout << "Resultant Load correct:" << std::endl;
        resForceCorrect.Info();
        std::cout << "Total external load vector:" << std::endl;
        extLoadVector.Info();
        throw NuTo::MechanicsException("[SurfaceLoad::CheckSurfaceLoad] Wrong external load calculation!");
    }
}

void SurfaceLoadTriangle2D(NuTo::Interpolation::eTypeOrder rInterpolationCoords,
                           NuTo::Interpolation::eTypeOrder rInterpolationDisp)
{
    NuTo::Structure myStructure(2);
    myStructure.SetShowTime(false);

    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> nodeCoords(2,3);

    double lx = 2;
    double ly = 5;
    double thickness = 13.;

    nodeCoords <<
            0, lx, 0,
           -ly, 0,ly;

    NuTo::FullMatrix<int, Eigen::Dynamic> nodeIds = myStructure.NodesCreate(nodeCoords);

    int s = myStructure.GroupCreate("Nodes");
    myStructure.GroupAddNode(s, nodeIds.GetValue(0));
    myStructure.GroupAddNode(s, nodeIds.GetValue(1));

    double surfaceArea = thickness*std::sqrt(lx*lx + ly*ly);

    myStructure.InterpolationTypeCreate(0, NuTo::Interpolation::TRIANGLE2D);
    myStructure.InterpolationTypeAdd(0, NuTo::Node::COORDINATES, rInterpolationCoords);
    myStructure.InterpolationTypeAdd(0, NuTo::Node::DISPLACEMENTS, rInterpolationDisp);

    int elementGroup = myStructure.ElementsCreate(0,nodeIds);
    myStructure.ElementConvertToInterpolationType(elementGroup);

    int section = myStructure.SectionCreate(NuTo::Section::PLANE_STRESS);
    myStructure.SectionSetThickness(section, thickness);
    myStructure.ElementTotalSetSection(section);

    NuTo::FullVector<double, Eigen::Dynamic> load(2);
    load << 42.,-M_PI;

    myStructure.LoadSurfaceConstDirectionCreate2D(0, elementGroup, s, load);

    CheckSurfaceLoad(myStructure, load, surfaceArea);
    std::cout << "[SurfaceLoad::SurfaceLoadTriangle2D] done." << std::endl;
}


void SurfaceLoadQuad2D(NuTo::Interpolation::eTypeOrder rInterpolationCoords,
                       NuTo::Interpolation::eTypeOrder rInterpolationDisp)
{
    NuTo::Structure myStructure(2);
    myStructure.SetShowTime(false);

    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> nodeCoords(2,4);

    double lx = 2;
    double ly = 5;
    double thickness = 13.;

    nodeCoords <<
            0, lx, 0, -lx,
           -ly, 0,ly, 0;

    NuTo::FullMatrix<int, Eigen::Dynamic> nodeIds = myStructure.NodesCreate(nodeCoords);

    int s = myStructure.GroupCreate("Nodes");
    myStructure.GroupAddNode(s, nodeIds.GetValue(0));
    myStructure.GroupAddNode(s, nodeIds.GetValue(1));

    double surfaceArea = thickness*std::sqrt(lx*lx + ly*ly);

    myStructure.InterpolationTypeCreate(0, NuTo::Interpolation::QUAD2D);
    myStructure.InterpolationTypeAdd(0, NuTo::Node::COORDINATES, rInterpolationCoords);
    myStructure.InterpolationTypeAdd(0, NuTo::Node::DISPLACEMENTS, rInterpolationDisp);

    int elementGroup = myStructure.ElementsCreate(0,nodeIds);
    myStructure.ElementConvertToInterpolationType(elementGroup);

    int section = myStructure.SectionCreate(NuTo::Section::PLANE_STRESS);
    myStructure.SectionSetThickness(section, thickness);
    myStructure.ElementTotalSetSection(section);

    NuTo::FullVector<double, Eigen::Dynamic> load(2);
    load << 42.,-M_PI;

    myStructure.LoadSurfaceConstDirectionCreate2D(0, elementGroup, s, load);

    CheckSurfaceLoad(myStructure, load, surfaceArea);
    std::cout << "[SurfaceLoad::SurfaceLoadQuad2D] done." << std::endl;
}

void SurfaceLoadTetrahedron3D(NuTo::Interpolation::eTypeOrder rInterpolationCoords,
                              NuTo::Interpolation::eTypeOrder rInterpolationDisp)
{
    NuTo::Structure myStructure(3);
    myStructure.SetShowTime(false);

    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> nodeCoords(3,4), nodeCoordsRot(3,4);

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

    NuTo::FullMatrix<int, Eigen::Dynamic> nodeIds = myStructure.NodesCreate(nodeCoordsRot);

    int s = myStructure.GroupCreate("Nodes");
    myStructure.GroupAddNode(s, nodeIds.GetValue(0));
    myStructure.GroupAddNode(s, nodeIds.GetValue(1));
    myStructure.GroupAddNode(s, nodeIds.GetValue(3));

    double surfaceArea = .5*lx*lz;

    myStructure.InterpolationTypeCreate(0, NuTo::Interpolation::TETRAHEDRON3D);
    myStructure.InterpolationTypeAdd(0, NuTo::Node::COORDINATES, rInterpolationCoords);
    myStructure.InterpolationTypeAdd(0, NuTo::Node::DISPLACEMENTS, rInterpolationDisp);

    int elementGroup = myStructure.ElementsCreate(0,nodeIds);
    myStructure.ElementConvertToInterpolationType(elementGroup);

    NuTo::FullVector<double, Eigen::Dynamic> load(3);
    load << 42.,-M_PI,4.2e6;

    myStructure.LoadSurfaceConstDirectionCreate3D(0, elementGroup, s, load);

    CheckSurfaceLoad(myStructure, load, surfaceArea);
    std::cout << "[SurfaceLoad::SurfaceLoadTetrahedron3D] done." << std::endl;
}

void SurfaceLoadBrick3D(NuTo::Interpolation::eTypeOrder rInterpolationCoords,
                        NuTo::Interpolation::eTypeOrder rInterpolationDisp)
{
    NuTo::Structure myStructure(3);
    myStructure.SetShowTime(false);

    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> nodeCoords(3,8), nodeCoordsRot(3,8);

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

    NuTo::FullMatrix<int, Eigen::Dynamic> nodeIds = myStructure.NodesCreate(nodeCoordsRot);

    int s = myStructure.GroupCreate("Nodes");
    myStructure.GroupAddNode(s, nodeIds.GetValue(0));
    myStructure.GroupAddNode(s, nodeIds.GetValue(1));
    myStructure.GroupAddNode(s, nodeIds.GetValue(5));
    myStructure.GroupAddNode(s, nodeIds.GetValue(4));

    double surfaceArea = lx*lz;

    myStructure.InterpolationTypeCreate(0, NuTo::Interpolation::BRICK3D);
    myStructure.InterpolationTypeAdd(0, NuTo::Node::COORDINATES, rInterpolationCoords);
    myStructure.InterpolationTypeAdd(0, NuTo::Node::DISPLACEMENTS, rInterpolationDisp);

    int elementGroup = myStructure.ElementsCreate(0,nodeIds);
    myStructure.ElementConvertToInterpolationType(elementGroup);

    NuTo::FullVector<double, Eigen::Dynamic> load(3);
    load << 42.,-M_PI,4.2e6;

    myStructure.LoadSurfaceConstDirectionCreate3D(0, elementGroup, s, load);

    CheckSurfaceLoad(myStructure, load, surfaceArea);
    std::cout << "[SurfaceLoad::SurfaceLoadBrick3D] done." << std::endl;
}


int main()
{

    // 1) all surfaces need testing!
    // 2) pressure (normal) and surface load (certain direction) to be tested
    // 3) not axis aligned
    // 4) analytical solution as reference
//    try
//    {

        HydrostaticPressureTriangle2D(NuTo::Interpolation::EQUIDISTANT1, NuTo::Interpolation::EQUIDISTANT3);
        HydrostaticPressureQuad2D(NuTo::Interpolation::EQUIDISTANT1, NuTo::Interpolation::EQUIDISTANT1);
        HydrostaticPressureBrick3D(NuTo::Interpolation::EQUIDISTANT1, NuTo::Interpolation::EQUIDISTANT1);
        HydrostaticPressureTetrahedron3D(NuTo::Interpolation::EQUIDISTANT1, NuTo::Interpolation::EQUIDISTANT2);

        SurfaceLoadTriangle2D(NuTo::Interpolation::EQUIDISTANT1, NuTo::Interpolation::EQUIDISTANT1);

        SurfaceLoadQuad2D(NuTo::Interpolation::EQUIDISTANT1,NuTo::Interpolation::EQUIDISTANT1);
        SurfaceLoadQuad2D(NuTo::Interpolation::EQUIDISTANT1,NuTo::Interpolation::EQUIDISTANT2);

        SurfaceLoadTetrahedron3D(NuTo::Interpolation::EQUIDISTANT1,NuTo::Interpolation::EQUIDISTANT2);
        SurfaceLoadBrick3D(NuTo::Interpolation::EQUIDISTANT1,NuTo::Interpolation::EQUIDISTANT1);

//    }
//    catch (NuTo::Exception& e)
//    {
//        std::cout << e.ErrorMessage() << std::endl;
//        return -1;
//    }

    return 0;
}
