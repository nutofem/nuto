#include <iostream>
#include "base/Exception.h"
#include "math/FullMatrix.h"
#include "math/FullVector.h"
#include "math/SparseMatrixCSRGeneral.h"
#include "math/SparseMatrixCSRVector2General.h"
#include "math/SparseDirectSolverMUMPS.h"
#include "mechanics/structures/unstructured/Structure.h"
#include "mechanics/structures/StructureOutputBlockVector.h"
#include "mechanics/dofSubMatrixStorage/BlockScalar.h"
#include "mechanics/MechanicsEnums.h"
#include "visualize/VisualizeEnum.h"

int main()
{
    // definitions
    double YoungsModulus = 20000.;
    double PoissonsRatio = 0.2;
    double Width = 1000.;
    double Height = 1000.;
    double Length = 1000.;
    int NumElementsX = 6;
    int NumElementsY = 6;
    int NumElementsZ = 6;
    double Force = 1.;
    bool EnableDisplacementControl = true;
    double BoundaryDisplacement = 0.1;

    // create structure
    NuTo::Structure myStructure(3);

    // create material law
    int Material1 = myStructure.ConstitutiveLawCreate("Linear_Elastic_Engineering_Stress");
    myStructure.ConstitutiveLawSetParameterDouble(Material1,
            NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, YoungsModulus);
    myStructure.ConstitutiveLawSetParameterDouble(Material1,
            NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, PoissonsRatio);

    int mySection1 = myStructure.SectionCreate("VOLUME");

    // create nodes
    NuTo::FullVector<double,Eigen::Dynamic> nodeCoordinates(3);
    int node = 0;
    for(int zCount = 0; zCount < NumElementsZ + 1; zCount++)
    {
        nodeCoordinates(2) = (double)zCount * Height/(double)NumElementsZ;
        for(int yCount = 0; yCount < NumElementsY + 1; yCount++)
        {
            nodeCoordinates(1) = (double)yCount * Width/(double)NumElementsY;
            for(int xCount = 0; xCount < NumElementsX + 1; xCount++)
            {
                nodeCoordinates(0) = (double)xCount * Length/(double)NumElementsX;
                myStructure.NodeCreate(node, nodeCoordinates);
                node++;
            }
        }
    }

    int InterpolationType = myStructure.InterpolationTypeCreate("BRICK3D");
    myStructure.InterpolationTypeAdd(InterpolationType, NuTo::Node::eDof::COORDINATES,
            NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    myStructure.InterpolationTypeAdd(InterpolationType, NuTo::Node::eDof::DISPLACEMENTS,
            NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);


    // create elements
    std::vector<int> elementIncidence(8);
    int element = 0;
    for(int zCount = 0; zCount < NumElementsZ; zCount++)
    {
        for(int yCount = 0; yCount < NumElementsY; yCount++)
        {
            for(int xCount = 0; xCount < NumElementsX; xCount++)
            {
                int node1 = zCount * (NumElementsX + 1) * (NumElementsY + 1) + yCount * (NumElementsX + 1) + xCount;
                elementIncidence[0] = node1;
                elementIncidence[1] = node1 + 1;
                elementIncidence[2] = node1 + NumElementsX + 2;
                elementIncidence[3] = node1 + NumElementsX + 1;
                elementIncidence[4] = node1 + (NumElementsX + 1) * (NumElementsY + 1);
                elementIncidence[5] = node1 + (NumElementsX + 1) * (NumElementsY + 1) + 1;
                elementIncidence[6] = node1 + (NumElementsX + 1) * (NumElementsY + 1) + NumElementsX + 2;
                elementIncidence[7] = node1 + (NumElementsX + 1) * (NumElementsY + 1) + NumElementsX + 1;

                myStructure.ElementCreate(InterpolationType, elementIncidence);
                myStructure.ElementSetConstitutiveLaw(element,Material1);
                myStructure.ElementSetSection(element,mySection1);
                element ++;
            }
        }
    }

    myStructure.ElementTotalConvertToInterpolationType();

    // boundary conditions
    NuTo::FullVector<double,Eigen::Dynamic> direction(3);
    direction(0)= 1;
    direction(1)= 0;
    direction(2)= 0;
    for(int zCount = 0; zCount < NumElementsZ + 1; zCount++)
    {
        for(int yCount = 0; yCount < NumElementsY + 1; yCount++)
        {
            int node = zCount * (NumElementsX + 1) * (NumElementsY + 1) + yCount * (NumElementsX + 1);
            myStructure.ConstraintLinearSetDisplacementNode(node, direction, 0.0);
        }
    }
    direction(0)= 0;
    direction(1)= 0;
    direction(2)= 1;
    myStructure.ConstraintLinearSetDisplacementNode(0, direction, 0.0);
    myStructure.ConstraintLinearSetDisplacementNode(NumElementsY * (NumElementsX + 1), direction, 0.0);
    direction(0)= 0;
    direction(1)= 1;
    direction(2)= 0;
    myStructure.ConstraintLinearSetDisplacementNode(0, direction, 0.0);

    // apply nodes
    if(EnableDisplacementControl)
    {
        std::cout << "Displacement control" << std::endl;
        // boundary displacments
        direction(0)= 1;
        direction(1)= 0;
        direction(2)= 0;
        for(int zCount = 0; zCount < NumElementsZ + 1; zCount++)
        {
            for(int yCount = 0; yCount < NumElementsY + 1; yCount++)
            {
                int node = zCount * (NumElementsX + 1) * (NumElementsY + 1) + yCount * (NumElementsX + 1) + NumElementsX;
                myStructure.ConstraintLinearSetDisplacementNode(node, direction, BoundaryDisplacement);
            }
        }
    }
    else
    {
        std::cout << "Load control" << std::endl;
        // apply load to nodes
        direction(0)= 1;
        direction(1)= 0;
        direction(2)= 0;
        for(int zCount = 0; zCount < NumElementsZ + 1; zCount++)
        {
            double nodeForce;
            if(zCount == 0 || zCount == NumElementsZ)
            {
                nodeForce = Force / (4 *NumElementsY * NumElementsZ);
            }
            else
            {
                nodeForce = Force / (2 *NumElementsY * NumElementsZ);
            }
            int node = zCount * (NumElementsX + 1) * (NumElementsY + 1) + NumElementsX;
            myStructure.LoadCreateNodeForce(1, node, direction, nodeForce);
            for(int yCount = 1; yCount < NumElementsY; yCount++)
            {
                node = zCount * (NumElementsX + 1) * (NumElementsY + 1) + yCount * (NumElementsX + 1) + NumElementsX;
                myStructure.LoadCreateNodeForce(1, node, direction, 2 * nodeForce);
            }
            node = (zCount + 1) * (NumElementsX + 1) * (NumElementsY + 1) - 1;
            myStructure.LoadCreateNodeForce(1 ,node, direction, nodeForce);
        }
    }

    // start analysis
    myStructure.CalculateMaximumIndependentSets();
    myStructure.SolveGlobalSystemStaticElastic(1);
    auto residual = myStructure.BuildGlobalInternalGradient() - myStructure.BuildGlobalExternalLoadVector(1);

    std::cout << "residual: " << residual.J.CalculateNormL2() << std::endl;

    // visualize results
    int visualizationGroup = myStructure.GroupCreate(NuTo::eGroupId::Elements);
    myStructure.GroupAddElementsTotal(visualizationGroup);

    myStructure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::DISPLACEMENTS);
    myStructure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::ENGINEERING_STRAIN);
    myStructure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::ENGINEERING_STRESS);

    myStructure.ExportVtkDataFileElements("Brick8N.vtk");
    return 0;
}

