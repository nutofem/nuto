#pragma once

#include "nuto/math/FullVector.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"

#include <array>
#include <vector>


//VHIRTHAMTODO Add dimensions to specific function names like "MeshCylinder" +"3D"
namespace NuTo
{

class MeshGenerator
{
    //! @brief constructor -> deleted, because class only uses static functions
    MeshGenerator() = delete;



    /*---------------------------------------------*\
    |*             1D mesh generation              *|
    \*---------------------------------------------*/

private:

    //! @brief Creates a 1D mesh based on an equidistant line segment mesh. The final shape is created by a mapping function.
    //! @param rStructure ... reference to structure
    //! @param rSection ... section ID
    //! @param rInterpolationType ... interpolation type ID
    //! @param rNumElements ... vector with number of elements for each dimension
    //! @param rLength ... vector with length in each dimension
    //! @param rMappingFunction ... mapping function
    static void GenerateMesh1D(NuTo::Structure &rStructure,
                                        int rSection,
                                        int rConstitutiveLaw,
                                        int rInterpolationType,
                                        std::array<int, 1> rNumElements,
                                        std::function<NuTo::FullVector<double,Eigen::Dynamic>(double)> rMappingFunction);


    //! @brief Calculates the coordinates for a equidistant line segment [0,1] and maps them to new coordinates with help of a given function
    //! @param rNumElements ... vector with number of elements for each dimension
    //! @param rMappingFunction ... mapping function
    static std::vector<int> GetNodeCoordinatesLineSegmentMapped(NuTo::Structure &rStructure,
                                                                std::array<int,1> rNumElements,
                                                                std::function<NuTo::FullVector<double,Eigen::Dynamic>(double)> rMappingFunction);

public:

    //! @brief Creates a line segment mesh
    //! @param rStructure ... reference to structure
    //! @param rSection ... section ID
    //! @param rInterpolationType ... interpolation type ID
    //! @param rNumElements ... vector with number of elements for each dimension
    //! @param rLength ... vector with length in each dimension
    static void MeshLineSegment(NuTo::Structure &rStructure,
                                int rSection,
                                int rConstitutiveLaw,
                                int rInterpolationType,
                                std::array<int, 1> rNumElements,
                                std::array<double,1> rLength);


    /*---------------------------------------------*\
    |*             2D mesh generation              *|
    \*---------------------------------------------*/

private:

    //! @brief Creates a 2D mesh based on an equidistant square mesh. The final shape is created by a mapping function.
    //! @param rStructure ... reference to structure
    //! @param rSection ... section ID
    //! @param rInterpolationType ... interpolation type ID
    //! @param rNumElements ... vector with number of elements for each dimension
    //! @param rLength ... vector with length in each dimension
    //! @param rMappingFunction ... mapping function
    static void GenerateMesh2D(NuTo::Structure &rStructure,
                                        int rSection,
                                        int rConstitutiveLaw,
                                        int rInterpolationType,
                                        std::array<int, 2> rNumElements,
                                        std::function<NuTo::FullVector<double,Eigen::Dynamic>(double, double)> rMappingFunction);


    //! @brief Calculates the coordinates for a equidistant, square plane mesh [0,1] and maps them to new coordinates with help of a given function
    //! @param rNumElements ... vector with number of elements for each dimension
    //! @param rMappingFunction ... mapping function
    static std::vector<int> GetNodeCoordinatesSquarePlaneMapped(NuTo::Structure &rStructure,
                                                                std::array<int,2> rNumElements,
                                                                std::function<NuTo::FullVector<double,Eigen::Dynamic>(double, double)> rMappingFunction);
public:
    //! @brief Creates a rectangular mesh
    //! @param rStructure ... reference to structure
    //! @param rSection ... section ID
    //! @param rInterpolationType ... interpolation type ID
    //! @param rNumElements ... vector with number of elements for each dimension
    //! @param rLength ... vector with length in each dimension
    static void MeshRectangularPlane(NuTo::Structure &rStructure,
                                        int rSection,
                                        int rConstitutiveLaw,
                                        int rInterpolationType,
                                        std::array<int, 2> rNumElements,
                                        std::array<double,2> rLength);


    /*---------------------------------------------*\
    |*             3D mesh generation              *|
    \*---------------------------------------------*/

private:

    //! @brief Creates a 3D mesh based on an equidistant cubic mesh. The final shape is created by a mapping function.
    //! @param rStructure ... reference to structure
    //! @param rSection ... section ID
    //! @param rInterpolationType ... interpolation type ID
    //! @param rNumElements ... vector with number of elements for each dimension
    //! @param rLength ... vector with length in each dimension
    //! @param rMappingFunction ... mapping function
    static void GenerateMesh3D(NuTo::Structure &rStructure,
                                        int rSection,
                                        int rConstitutiveLaw,
                                        int rInterpolationType,
                                        std::array<int, 3> rNumElements,
                                        std::function<NuTo::FullVector<double, Eigen::Dynamic> (double, double, double)> rMappingFunction);


    //! @brief Calculates the coordinates for a equidistant cuboid mesh
    //! @param rNumElements ... vector with number of elements for each dimension
    //! @param rLength ... vector with length in each dimension
    static std::vector<int> GetNodeCoordinatesCuboid(NuTo::Structure &rStructure,
                                                     std::array<int,3> rNumElements,
                                                     std::array<double,3> rLength);

    //! @brief Calculates the coordinates for a equidistant cuboid mesh and maps them to new coordinates with help of a given function (e.g. cylinder mesh based on cuboid coordinates)
    //! @param rNumElements ... vector with number of elements for each dimension
    //! @param rLength ... vector with length in each dimension
    //! @param rMappingFunction ... mapping function
    static std::vector<int> GetNodeCoordinatesCuboidMapped(NuTo::Structure &rStructure,
                                                           std::array<int,3> rNumElements,
                                                           std::function<NuTo::FullVector<double,Eigen::Dynamic>(double, double, double)> rMappingFunction);

public:

    //! @brief Creates a cuboid mesh
    //! @param rStructure ... reference to structure
    //! @param rSection ... section ID
    //! @param rInterpolationType ... interpolation type ID
    //! @param rNumElements ... vector with number of elements for each dimension
    //! @param rLength ... vector with length in each dimension
    static void MeshCuboid(NuTo::Structure &rStructure,
                           int rSection,
                           int rConstitutiveLaw,
                           int rInterpolationType,
                           std::array<int, 3> rNumElements,
                           std::array<double,3> rLength);

    //! @brief Creates a cylinder mesh
    //! @param rStructure ... reference to structure
    //! @param rSection ... section ID
    //! @param rInterpolationType ... interpolation type ID
    //! @param rNumElements ... vector with number of elements for each dimension
    //! @param rRadius ... radius of the cylinder
    //! @param rHeight ... height of the cylinder
    static void MeshCylinder(NuTo::Structure &rStructure,
                             int rSection,
                             int rConstitutiveLaw,
                             int rInterpolationType,
                             std::array<int,3> rNumElements,
                             double rRadius,
                             double rHeight);

};

}
