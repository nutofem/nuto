#include <iostream>
#include "nuto/math/FullMatrix.h"
#include "nuto/math/SparseMatrixCSRGeneral.h"
#include "nuto/math/SparseDirectSolverMUMPS.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/mechanics/structures/StructureOutputBlockVector.h"
#include "nuto/mechanics/constitutive/laws/AdditiveInputExplicit.h"
#include "nuto/mechanics/dofSubMatrixStorage/BlockScalar.h"
#include "nuto/mechanics/MechanicsEnums.h"
#include "nuto/visualize/VisualizeEnum.h"

int main()
{
    // definitions
    double YoungsModulus = 20000.;
    double Area = 100. * 100.;
    double Length = 1000.;
    int NumElements = 10;
    double Force = 1.;
    bool EnableDisplacementControl = false;
    double BoundaryDisplacement = 0.1;

    // create one-dimensional structure
    NuTo::Structure structure(1);

    // create section
    auto Section1 = structure.SectionCreate("Truss");
    structure.SectionSetArea(Section1, Area);

    // create material law
    auto Material1 = structure.ConstitutiveLawCreate("Linear_Elastic_Engineering_Stress");
    structure.ConstitutiveLawSetParameterDouble(
            Material1, NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, YoungsModulus);

    // create nodes
    NuTo::FullVector<double, Eigen::Dynamic> nodeCoordinates(1);
    for (int node = 0; node < NumElements + 1; node++)
    {
        std::cout << "create node: " << node << " coordinates: " << node * Length / NumElements << std::endl;
        nodeCoordinates(0) = node * Length / NumElements;
        structure.NodeCreate(node, nodeCoordinates);
    }

    int InterpolationType = structure.InterpolationTypeCreate("Truss1D");
    structure.InterpolationTypeAdd(InterpolationType, NuTo::Node::eDof::COORDINATES,
            NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    structure.InterpolationTypeAdd(InterpolationType, NuTo::Node::eDof::DISPLACEMENTS,
            NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);

    // create elements
    NuTo::FullVector<int, Eigen::Dynamic> elementIncidence(2);
    for (int element = 0; element < NumElements; element++)
    {
        std::cout << "create element: " << element << " nodes: " << element << "," << element + 1 << std::endl;
        elementIncidence(0) = element;
        elementIncidence(1) = element + 1;
        structure.ElementCreate(InterpolationType, elementIncidence);
        structure.ElementSetSection(element, Section1);
        structure.ElementSetConstitutiveLaw(element, Material1);
    }

    structure.ElementTotalConvertToInterpolationType();

    // set boundary conditions and loads
    NuTo::FullVector<double, Eigen::Dynamic> direction(1);
    direction(0) = 1;
    structure.ConstraintLinearSetDisplacementNode(0, direction, 0.0);
    structure.SetNumLoadCases(1);
    if (EnableDisplacementControl)
    {
        std::cout << "Displacement control" << std::endl;
        structure.ConstraintLinearSetDisplacementNode(NumElements, direction, BoundaryDisplacement);
    }
    else
    {
        std::cout << "Load control" << std::endl;
        structure.LoadCreateNodeForce(0, NumElements, direction, Force);
    }

    // start analysis
    structure.SolveGlobalSystemStaticElastic(1);
    auto residual = structure.BuildGlobalInternalGradient() - structure.BuildGlobalExternalLoadVector(1);

    std::cout << "residual: " << residual.J.CalculateNormL2() << std::endl;

    // visualize results
    int visualizationGroup = structure.GroupCreate(NuTo::eGroupId::Elements);
    structure.GroupAddElementsTotal(visualizationGroup);

    structure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::DISPLACEMENTS);
    structure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::ENGINEERING_STRAIN);
    structure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::ENGINEERING_STRESS);
    structure.ExportVtkDataFileElements("Truss1D2N.vtk");

    return 0;
}
