#include "MeshGenerator.h"


#include "math/FullVector.h"
#include "mechanics/structures/unstructured/Structure.h"
#include "boost/progress.hpp"


//! @brief Creates a 1D mesh based on an equidistant line segment mesh. The final shape is created by a mapping function.
//! @param rStructure ... reference to structure
//! @param rSection ... section ID
//! @param rInterpolationType ... interpolation type ID
//! @param rNumElements ... vector with number of elements for each dimension
//! @param rLength ... vector with length in each dimension
//! @param rMappingFunction ... mapping function
void NuTo::MeshGenerator::GenerateMesh1D(NuTo::Structure &rStructure,
                                         int rSection,
                                         int rConstitutiveLaw,
                                         int rInterpolationType,
                                         std::array<int, 1> rNumElements,
                                         std::function<Eigen::VectorXd (double)> rMappingFunction)
{
    std::vector<int> NodeIDs = MeshGenerator::GetNodeCoordinatesLineSegmentMapped(rStructure,
                                                                                  rNumElements,
                                                                                  rMappingFunction);

    std::cout << std::endl << "Creating Elements of Mesh" << std::endl << std::endl;

    std::vector<int> Nodes(2);
    for(int x_count = 0; x_count < rNumElements[0]; x_count++)
    {
        Nodes[0] = NodeIDs[x_count    ];
        Nodes[1] = NodeIDs[x_count + 1];

        int elementNumber = rStructure.ElementCreate(rInterpolationType, Nodes);

        // set element section
        rStructure.ElementSetSection(elementNumber,rSection);

        rStructure.ElementSetConstitutiveLaw(elementNumber,rConstitutiveLaw);
    }

}



//! @brief Calculates the coordinates for a equidistant line segment [0,1] and maps them to new coordinates with help of a given function
//! @param rNumElements ... vector with number of elements for each dimension
//! @param rMappingFunction ... mapping function
std::vector<int> NuTo::MeshGenerator::GetNodeCoordinatesLineSegmentMapped(NuTo::Structure &rStructure,
                                                                          std::array<int,1> rNumElements,
                                                                          std::function<Eigen::VectorXd (double)> rMappingFunction)
{
    assert(rNumElements[0]>0);

    std::vector<int> NodeIDs;


    double delta_x = 1.0/ static_cast<double>(rNumElements[0]);

    int NumNodesX = rNumElements[0]+1;

    NodeIDs.resize(NumNodesX);

    std::cout << std::endl << "Creating Nodes of Mesh" << std::endl << std::endl;

//    boost::progress_display show_progress( NodeIDs.size() );


    int NodeNum = 0;

    Eigen::VectorXd Coordinates(1);


        for(int x_count = 0; x_count < NumNodesX; x_count++)
        {
            NodeNum = x_count;

            Coordinates = rMappingFunction(static_cast<double>(x_count) * delta_x);

            NodeIDs[NodeNum] = rStructure.NodeCreate(Coordinates);

//            ++show_progress;
        }

        return NodeIDs;
}




//! @brief Creates a line segment mesh
//! @param rStructure ... reference to structure
//! @param rSection ... section ID
//! @param rInterpolationType ... interpolation type ID
//! @param rNumElements ... vector with number of elements for each dimension
//! @param rLength ... vector with length in each dimension
void NuTo::MeshGenerator::MeshLineSegment(NuTo::Structure &rStructure,
                                          int rSection,
                                          int rConstitutiveLaw,
                                          int rInterpolationType,
                                          std::array<int, 1> rNumElements,
                                          std::array<double, 1> rLength)
{
    assert(rLength[0]>0.0);

    auto MappingFunction = [&rLength](double rX) -> Eigen::VectorXd
                            {
                                Eigen::VectorXd CoordVec(1);
                                CoordVec[0] = rX * rLength[0];
                                return CoordVec;
                            };

    GenerateMesh1D(rStructure,
                   rSection,
                   rConstitutiveLaw,
                   rInterpolationType,
                   rNumElements,
                   MappingFunction);
}






//! @brief Creates a 2D mesh based on an equidistant square mesh. The final shape is created by a mapping function.
//! @param rStructure ... reference to structure
//! @param rSection ... section ID
//! @param rInterpolationType ... interpolation type ID
//! @param rNumElements ... vector with number of elements for each dimension
//! @param rLength ... vector with length in each dimension

//! @param rMappingFunction ... mapping function
void NuTo::MeshGenerator::GenerateMesh2D(NuTo::Structure &rStructure,
                                         int rSection,
                                         int rConstitutiveLaw,
                                         int rInterpolationType,
                                         std::array<int, 2> rNumElements,
                                         std::function<Eigen::VectorXd (double, double)> rMappingFunction)
{

    std::vector<int> NodeIDs = MeshGenerator::GetNodeCoordinatesSquarePlaneMapped(rStructure,
                                                                                  rNumElements,
                                                                                  rMappingFunction);
    std::cout << std::endl << "Creating Elements of Mesh" << std::endl << std::endl;

    int NumNodesX = rNumElements[0]+1;
    std::vector<int> Nodes(4);
    for(int y_count = 0; y_count < rNumElements[1]; y_count++)
    {
        for(int x_count = 0; x_count < rNumElements[0]; x_count++)
        {
            Nodes[0] = NodeIDs[x_count       +  y_count       * NumNodesX];
            Nodes[1] = NodeIDs[x_count + 1   +  y_count       * NumNodesX];
            Nodes[2] = NodeIDs[x_count + 1   + (y_count + 1)  * NumNodesX];
            Nodes[3] = NodeIDs[x_count       + (y_count + 1)  * NumNodesX];

            int elementNumber = rStructure.ElementCreate(rInterpolationType, Nodes);

            // set element section
            rStructure.ElementSetSection(elementNumber,rSection);

            rStructure.ElementSetConstitutiveLaw(elementNumber,rConstitutiveLaw);
        }
    }
}



//! @brief Calculates the coordinates for a equidistant, square plane mesh [-1,1] and maps them to new coordinates with help of a given function
//! @param rNumElements ... vector with number of elements for each dimension
//! @param rMappingFunction ... mapping function
std::vector<int> NuTo::MeshGenerator::GetNodeCoordinatesSquarePlaneMapped(NuTo::Structure &rStructure,
                                                                          std::array<int, 2> rNumElements,
                                                                          std::function<Eigen::VectorXd (double, double)> rMappingFunction)
{
    assert(rNumElements[0]>0);
    assert(rNumElements[1]>0);

    std::vector<int> NodeIDs;


    double delta_x = 1.0/ static_cast<double>(rNumElements[0]);
    double delta_y = 1.0/ static_cast<double>(rNumElements[1]);

    int NumNodesX = rNumElements[0]+1;
    int NumNodesY = rNumElements[1]+1;

    NodeIDs.resize(NumNodesX * NumNodesY);

    std::cout << std::endl << "Creating Nodes of Mesh" << std::endl << std::endl;

//    boost::progress_display show_progress( NodeIDs.size() );


    int NodeNum = 0;

    Eigen::VectorXd Coordinates(2);


    for(int y_count = 0; y_count < NumNodesY; y_count++)
    {
        for(int x_count = 0; x_count < NumNodesX; x_count++)
        {
            NodeNum = x_count + y_count * NumNodesX;

            Coordinates = rMappingFunction(static_cast<double>(x_count) * delta_x,
                                           static_cast<double>(y_count) * delta_y);

            NodeIDs[NodeNum] = rStructure.NodeCreate(Coordinates);

//            ++show_progress;
        }
    }

    return NodeIDs;
}


//! @brief Creates a rectangular mesh
//! @param rStructure ... reference to structure
//! @param rSection ... section ID
//! @param rInterpolationType ... interpolation type ID
//! @param rNumElements ... vector with number of elements for each dimension
//! @param rLength ... vector with length in each dimension
void NuTo::MeshGenerator::MeshRectangularPlane(NuTo::Structure &rStructure,
                                               int rSection,
                                               int rConstitutiveLaw,
                                               int rInterpolationType,
                                               std::array<int, 2> rNumElements,
                                               std::array<double, 2> rLength)
{

    assert(rLength[0]>0.0);
    assert(rLength[1]>0.0);

    auto MappingFunction = [&rLength](double rX, double rY) -> Eigen::VectorXd
                            {
                                Eigen::VectorXd CoordVec(2);
                                CoordVec <<  rX * rLength[0], rY * rLength[1];
                                return CoordVec;
                            };

    GenerateMesh2D(rStructure,
                   rSection,
                   rConstitutiveLaw,
                   rInterpolationType,
                   rNumElements,
                   MappingFunction);
}




//! @brief Creates a 3D mesh based on an equidistant cubic mesh. The final shape is created by a mapping function.
//! @param rStructure ... reference to structure
//! @param rSection ... section ID
//! @param rInterpolationType ... interpolation type ID
//! @param rNumElements ... vector with number of elements for each dimension
//! @param rLength ... vector with length in each dimension
//! @param rMappingFunction ... mapping function
void NuTo::MeshGenerator::GenerateMesh3D(NuTo::Structure &rStructure,
                                         int rSection,
                                         int rConstitutiveLaw,
                                         int rInterpolationType,
                                         std::array<int, 3> rNumElements,
                                         std::function<Eigen::VectorXd (double, double, double)> rMappingFunction)
{
    std::vector<int> NodeIDs = MeshGenerator::GetNodeCoordinatesCuboidMapped(rStructure,
                                                                             rNumElements,
                                                                             rMappingFunction);
    std::cout << std::endl << "Creating Elements of Mesh" << std::endl << std::endl;

    int NumNodesX = rNumElements[0]+1;
    int NumNodesY = rNumElements[1]+1;
    std::vector<int> Nodes(8);
    for(int z_count = 0; z_count < rNumElements[2]; z_count++)
    {
        for(int y_count = 0; y_count < rNumElements[1]; y_count++)
        {
            for(int x_count = 0; x_count < rNumElements[0]; x_count++)
            {
                Nodes[0] = NodeIDs[x_count       +  y_count       * NumNodesX   +  z_count       * NumNodesX * NumNodesY];
                Nodes[1] = NodeIDs[x_count + 1   +  y_count       * NumNodesX   +  z_count       * NumNodesX * NumNodesY];
                Nodes[2] = NodeIDs[x_count + 1   + (y_count + 1)  * NumNodesX   +  z_count       * NumNodesX * NumNodesY];
                Nodes[3] = NodeIDs[x_count       + (y_count + 1)  * NumNodesX   +  z_count       * NumNodesX * NumNodesY];
                Nodes[4] = NodeIDs[x_count       + y_count        * NumNodesX   + (z_count + 1)  * NumNodesX * NumNodesY];
                Nodes[5] = NodeIDs[x_count + 1   + y_count        * NumNodesX   + (z_count + 1)  * NumNodesX * NumNodesY];
                Nodes[6] = NodeIDs[x_count + 1   + (y_count + 1)  * NumNodesX   + (z_count + 1)  * NumNodesX * NumNodesY];
                Nodes[7] = NodeIDs[x_count       + (y_count + 1)  * NumNodesX   + (z_count + 1)  * NumNodesX * NumNodesY];

                int elementNumber = rStructure.ElementCreate(rInterpolationType, Nodes);

                // set element section
                rStructure.ElementSetSection(elementNumber,rSection);

                rStructure.ElementSetConstitutiveLaw(elementNumber,rConstitutiveLaw);
            }
        }
    }
}








//! @brief Calculates the coordinates for a equidistant cuboid mesh and maps them to new coordinates with help of a given function (e.g. cylinder mesh based on cuboid coordinates)
//! @param rNumElements ... vector with number of elements for each dimension
//! @param rLength ... vector with length in each dimension
//! @param rMappingFunction ... mapping function
std::vector<int> NuTo::MeshGenerator::GetNodeCoordinatesCuboidMapped(NuTo::Structure &rStructure,
                                                                     std::array<int, 3> rNumElements,
                                                                     std::function<Eigen::VectorXd (double, double, double)> rMappingFunction)
{
    assert(rNumElements[0]>0);
    assert(rNumElements[1]>0);
    assert(rNumElements[2]>0);

    std::vector<int> NodeIDs;


    double delta_x = 1.0/ static_cast<double>(rNumElements[0]);
    double delta_y = 1.0/ static_cast<double>(rNumElements[1]);
    double delta_z = 1.0/ static_cast<double>(rNumElements[2]);

    int NumNodesX = rNumElements[0]+1;
    int NumNodesY = rNumElements[1]+1;
    int NumNodesZ = rNumElements[2]+1;

    NodeIDs.resize(NumNodesX * NumNodesY * NumNodesZ);

    std::cout << std::endl << "Creating Nodes of Mesh" << std::endl << std::endl;

//    boost::progress_display show_progress( NodeIDs.size() );


    int NodeNum = 0;

    Eigen::VectorXd Coordinates(3);

    for(int z_count = 0; z_count < NumNodesZ; z_count++)
    {
        for(int y_count = 0; y_count < NumNodesY; y_count++)
        {
            for(int x_count = 0; x_count < NumNodesX; x_count++)
            {
                NodeNum = x_count + y_count* NumNodesX + z_count * NumNodesX * NumNodesY;

                Coordinates = rMappingFunction(static_cast<double>(x_count) * delta_x,
                                               static_cast<double>(y_count) * delta_y,
                                               static_cast<double>(z_count) * delta_z);

                NodeIDs[NodeNum] = rStructure.NodeCreate(Coordinates);

//                ++show_progress;
            }
        }
    }

    return NodeIDs;
}



//! @brief Creates a cuboid mesh
//! @param rStructure ... reference to structure
//! @param rSection ... section ID
//! @param rInterpolationType ... interpolation type ID
//! @param rNumElements ... vector with number of elements for each dimension
//! @param rLength ... vector with length in each dimension
void NuTo::MeshGenerator::MeshCuboid(NuTo::Structure &rStructure,
                                     int rSection,
                                     int rConstitutiveLaw,
                                     int rInterpolationType,
                                     std::array<int, 3> rNumElements,
                                     std::array<double, 3> rLength)
{

    assert(rLength[0]>0.0);
    assert(rLength[1]>0.0);
    assert(rLength[2]>0.0);

    auto MappingFunction = [&rLength](double rX, double rY, double rZ) -> Eigen::VectorXd
                            {
                                Eigen::VectorXd CoordVec(3);
                                CoordVec << rX * rLength[0], rY * rLength[1], rZ * rLength[2];
                                return CoordVec;
                            };

    GenerateMesh3D(rStructure,
                   rSection,
                   rConstitutiveLaw,
                   rInterpolationType,
                   rNumElements,
                   MappingFunction);
}


//! @brief Creates a cylinder mesh
//! @param rStructure ... reference to structure
//! @param rSection ... section ID
//! @param rInterpolationType ... interpolation type ID
//! @param rNumElements ... vector with number of elements for each dimension
//! @param rRadius ... radius of the cylinder
//! @param rHeight ... height of the cylinder
void NuTo::MeshGenerator::MeshCylinder(NuTo::Structure &rStructure,
                                       int rSection,
                                       int rConstitutiveLaw,
                                       int rInterpolationType,
                                       std::array<int, 3> rNumElements,
                                       double rRadius,
                                       double rHeight, std::array<double, 2> rDensitiyFactor)
{
    assert(rRadius>0.0);
    assert(rHeight>0.0);

    auto MappingFunction = [&rRadius,rHeight](double rX, double rY, double rZ) -> Eigen::VectorXd
                            {
                                rX = rX *2 -1;
                                rY = rY *2 -1;
                                rZ = rZ *2 -1;
                                rX *= 1.+ (1.-std::abs(rX))/2.;
                                rY *= 1.+ (1.-std::abs(rY))/2.;
                                rZ *= 1.+ (1.-std::abs(rZ))/2.;
                                Eigen::VectorXd CoordVec(3);
                                CoordVec <<  rX * sqrt(1 - (rY * rY) / 2.0 ) * rRadius / 2.0,
                                             rY * sqrt(1 - (rX * rX) / 2.0 ) * rRadius / 2.0,
                                             rZ * rHeight / 2.0;
//                                double radius = sqrt(CoordVec[0]*CoordVec[0] + CoordVec[1]*CoordVec[1]);
//                                double coordCorrectionRadius = 1.0 + (rRadius/2.0 - radius) / rRadius/2.0;
//                                CoordVec[0]*=coordCorrectionRadius;
//                                CoordVec[1]*=coordCorrectionRadius;
                                return CoordVec;
                            };

    GenerateMesh3D(rStructure,
                   rSection,
                   rConstitutiveLaw,
                   rInterpolationType,
                   rNumElements,
                   MappingFunction);

}
