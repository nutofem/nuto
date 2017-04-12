#include <iostream>

#include "math/SparseMatrixCSRGeneral.h"
#include "math/SparseDirectSolverMUMPS.h"

#include "mechanics/MechanicsEnums.h"
#include "mechanics/constitutive/laws/AdditiveInputExplicit.h"
#include "mechanics/constraints/ConstraintCompanion.h"
#include "mechanics/dofSubMatrixStorage/BlockScalar.h"
#include "mechanics/sections/SectionTruss.h"
#include "mechanics/structures/unstructured/Structure.h"
#include "mechanics/structures/StructureOutputBlockVector.h"
#include "visualize/VisualizeEnum.h"

using namespace NuTo;

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
    Structure structure(1);

    // create section
    auto Section1 = SectionTruss::Create(Area);

    // create material law
    auto Material1 = structure.ConstitutiveLawCreate("Linear_Elastic_Engineering_Stress");
    structure.ConstitutiveLawSetParameterDouble(Material1, Constitutive::eConstitutiveParameter::YOUNGS_MODULUS,
                                                YoungsModulus);

    // create nodes
    Eigen::VectorXd nodeCoordinates(1);
    for (int node = 0; node < NumElements + 1; node++)
    {
        std::cout << "create node: " << node << " coordinates: " << node * Length / NumElements << std::endl;
        nodeCoordinates(0) = node * Length / NumElements;
        structure.NodeCreate(node, nodeCoordinates);
    }

    int InterpolationType = structure.InterpolationTypeCreate("Truss1D");
    structure.InterpolationTypeAdd(InterpolationType, Node::eDof::COORDINATES, Interpolation::eTypeOrder::EQUIDISTANT1);
    structure.InterpolationTypeAdd(InterpolationType, Node::eDof::DISPLACEMENTS,
                                   Interpolation::eTypeOrder::EQUIDISTANT1);

    // create elements
    std::vector<int> elementIncidence(2);
    for (int element = 0; element < NumElements; element++)
    {
        std::cout << "create element: " << element << " nodes: " << element << "," << element + 1 << std::endl;
        elementIncidence[0] = element;
        elementIncidence[1] = element + 1;
        structure.ElementCreate(InterpolationType, elementIncidence);
        structure.ElementSetSection(element, Section1);
        structure.ElementSetConstitutiveLaw(element, Material1);
    }

    structure.ElementTotalConvertToInterpolationType();

    // set boundary conditions and loads
    auto& origin = structure.NodeGetAtCoordinate(Eigen::Matrix<double, 1, 1>::Zero());
    structure.Constraints().Add(Node::eDof::DISPLACEMENTS, Constraint::Value(origin));
    if (EnableDisplacementControl)
    {
        std::cout << "Displacement control" << std::endl;
        auto& lastNode = *structure.NodeGetNodePtr(NumElements);
        structure.Constraints().Add(Node::eDof::DISPLACEMENTS, Constraint::Value(lastNode, BoundaryDisplacement));
    }
    else
    {
        std::cout << "Load control" << std::endl;
        Eigen::VectorXd direction(1);
        direction(0) = 1;
        structure.LoadCreateNodeForce(NumElements, direction, Force);
    }

    // start analysis
    structure.SolveGlobalSystemStaticElastic();
    auto residual = structure.BuildGlobalInternalGradient() - structure.BuildGlobalExternalLoadVector();

    std::cout << "residual: " << residual.J.CalculateNormL2() << std::endl;

    // visualize results
    int visualizationGroup = structure.GroupCreate(eGroupId::Elements);
    structure.GroupAddElementsTotal(visualizationGroup);

    structure.AddVisualizationComponent(visualizationGroup, eVisualizeWhat::DISPLACEMENTS);
    structure.AddVisualizationComponent(visualizationGroup, eVisualizeWhat::ENGINEERING_STRAIN);
    structure.AddVisualizationComponent(visualizationGroup, eVisualizeWhat::ENGINEERING_STRESS);
    structure.ExportVtkDataFileElements("Truss1D2N.vtk");

    return 0;
}
