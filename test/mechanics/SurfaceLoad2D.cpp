#include "nuto/math/FullMatrix.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/mechanics/MechanicsException.h"

#include "nuto/math/SparseDirectSolverMUMPS.h"
#include "nuto/math/SparseMatrixCSRVector2General.h"

#define PRINTRESULT false

int main()
{
    try
    {
        /**
         *  Test of the surface integration in 2D.
         *  The corresponding 2D node methods are currently (only)
         *  used in the calculation of the external force vector.
         *  Thus, this is tested here.
         *
         *  Setup (with poorly drawn load)
         *   ^
         *   |
         *
         *
         *   x/\
         * ly|\/\
         *   | \/\
         *   |  \/\
         *   x   x/\
         *   |    \/\
         *   |     \/
         *   x--x---x    ----->
         *          lx
         */

        //create structure
        NuTo::Structure myStructure(2);

        double lx = 42;
        double ly = 6174;
        double thickness = 13.37;
        double pressure = 3.141529;

        //create nodes
        NuTo::FullVector<double, Eigen::Dynamic> node0Coords({ 0, 0});
        NuTo::FullVector<double, Eigen::Dynamic> node1Coords({lx, 0});
        NuTo::FullVector<double, Eigen::Dynamic> node2Coords({ 0,ly});

        NuTo::FullVector<double, Eigen::Dynamic> node3Coords({lx/2.,    0.});
        NuTo::FullVector<double, Eigen::Dynamic> node4Coords({lx/2., ly/2.});
        NuTo::FullVector<double, Eigen::Dynamic> node5Coords({   0., ly/2.});


        int node0Index = myStructure.NodeCreate("displacements", node0Coords);
        int node1Index = myStructure.NodeCreate("displacements", node1Coords);
        int node2Index = myStructure.NodeCreate("displacements", node2Coords);
        int node3Index = myStructure.NodeCreate("displacements", node3Coords);
        int node4Index = myStructure.NodeCreate("displacements", node4Coords);
        int node5Index = myStructure.NodeCreate("displacements", node5Coords);


        //create element
        NuTo::FullMatrix<int, Eigen::Dynamic, Eigen::Dynamic> nodeNumbers(6,1);
        nodeNumbers << node0Index, node1Index, node2Index, node3Index, node4Index, node5Index;

        myStructure.ElementsCreate("PLANE2D6N", nodeNumbers);

        //Calculate maximum independent sets for parallelization (openmp)
        myStructure.CalculateMaximumIndependentSets();

        int nodeGroup = myStructure.GroupCreate("NODES");
        int elementGroup = myStructure.GroupCreate("ELEMENTS");

        myStructure.GroupAddNode(nodeGroup, node1Index);
        myStructure.GroupAddNode(nodeGroup, node2Index);
        myStructure.GroupAddNode(nodeGroup, node4Index);

        myStructure.GroupAddElementsFromNodes(elementGroup, nodeGroup, false);



        //create section
        int mySection = myStructure.SectionCreate("Plane_Strain");
        myStructure.SectionSetThickness(mySection,thickness);

        myStructure.ElementTotalSetSection(mySection);

        myStructure.LoadSurfacePressureCreate2D(0,elementGroup, nodeGroup, pressure);


        myStructure.SetVerboseLevel(10);
        myStructure.Info();

        // build global external load vector
        NuTo::FullVector<double,Eigen::Dynamic> extForceVector;
        myStructure.BuildGlobalExternalLoadVector(0,extForceVector);

        std::cout << "external force vector" << std::endl;
        extForceVector.Info();

        int numNodes = myStructure.GetNumNodes();

        // sum up all external node forces
        double Fx = 0;
        double Fy = 0;
        for (int iNode = 0; iNode < numNodes; ++iNode)
        {
            Fx += extForceVector(2*iNode);
            Fy += extForceVector(2*iNode+1);
        }

        /**
         * correct result:
         *   F = - Int_A (n * p) dA
         *     = - t * Int_S (n * p) dS
         *     = -t * n * p * sqrt(lx^2+ly^2)
         *
         */

        double abs_n = sqrt(lx*lx + ly*ly);
        double nx = ly / abs_n;
        double ny = lx / abs_n;

        double Fx_analytical = -thickness * nx * pressure * abs_n;
        double Fy_analytical = -thickness * ny * pressure * abs_n;

        std::cout << "Fx_analytical = " << Fx_analytical << std::endl;
        std::cout << "Fx_numerical  = " << Fx << std::endl;
        std::cout << "Fy_analytical = " << Fy_analytical << std::endl;
        std::cout << "Fy_numerical  = " << Fy << std::endl;

        if (std::abs(Fx_analytical - Fx) > 1.e-8 or std::abs(Fy_analytical - Fy) > 1.e-8)
        {

            throw NuTo::Exception("[SurfaceLoad2D] wrong forces.");
        }

    }
    catch (NuTo::Exception& e)
    {
        std::cout << e.ErrorMessage() << std::endl;
        return -1;
    }

    return 0;
}
