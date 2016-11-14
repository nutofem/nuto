#include <iostream>
#include "nuto/math/FullMatrix.h"
#include "nuto/math/SparseMatrixCSRGeneral.h"
#include "nuto/math/SparseDirectSolverMUMPS.h"
#include "nuto/mechanics/dofSubMatrixStorage/BlockScalar.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/mechanics/structures/StructureOutputBlockVector.h"
#include "nuto/mechanics/MechanicsEnums.h"
#include "nuto/visualize/VisualizeEnum.h"

int main()
{
    // definitions
    double YoungsModulus = 10000.;
    double YoungsModulusDiff = 20000.; // difference from node 0 to lats node, +/-
    double Area = 10. * 10.;
    int NumElements = 4;
    double Force = 1.;
    bool EnableDisplacementControl = true;
    double BoundaryDisplacement = 1;
    std::ofstream output;

    // create one-dimensional structure
    NuTo::Structure structure(1);

    // create section
    int Section1 = structure.SectionCreate("Truss");
    structure.SectionSetArea(Section1, Area);

    // create material law
    int Material1 = structure.ConstitutiveLawCreate("LINEAR_ELASTIC_ENGINEERING_STRESS");
    int Material2 = structure.ConstitutiveLawCreate("LINEAR_ELASTIC_ENGINEERING_STRESS");
    structure.ConstitutiveLawSetParameterDouble(
            Material1, NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, YoungsModulus);
    structure.ConstitutiveLawSetParameterDouble(
            Material2, NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, YoungsModulus + YoungsModulusDiff);

    std::string filename = "Truss1D2NReference-nodes";
    output.open(filename.c_str());
    if (output)
    {
        output << "#" << filename << "\n";
        output << "nodes\n";
    }
    else
    {
        std::cout << "Truss1DMultiphase] error output file.\n";
    }

    // create nodes
    NuTo::FullVector<double, Eigen::Dynamic> nodeCoordinates(1);
    nodeCoordinates(0) = 0.;
    output << nodeCoordinates(0) << "\n";
    structure.NodeCreate(0, nodeCoordinates);
    nodeCoordinates(0) = 2000.;
    output << nodeCoordinates(0) << "\n";
    structure.NodeCreate(1, nodeCoordinates);
    nodeCoordinates(0) = 3000.0;
    output << nodeCoordinates(0) << "\n";
    structure.NodeCreate(2, nodeCoordinates);
    nodeCoordinates(0) = 4000.0;
    output << nodeCoordinates(0) << "\n";
    structure.NodeCreate(3, nodeCoordinates);
    nodeCoordinates(0) = 6000.0;
    output << nodeCoordinates(0) << "\n";
    structure.NodeCreate(4, nodeCoordinates);
    output.close();

    int interpolationType = structure.InterpolationTypeCreate("TRUSS1D");
    structure.InterpolationTypeAdd(
            interpolationType, NuTo::Node::eDof::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    structure.InterpolationTypeAdd(
            interpolationType, NuTo::Node::eDof::DISPLACEMENTS, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);

    filename = "Truss1D2NReference-E";
    output.open(filename.c_str());
    if (output)
    {
        output << "#" << filename << "\n";
        output << "modul\n";
    }
    else
        std::cout << __LINE__ << " Truss1D2NReference] error output file.\n";
    // create elements
    NuTo::FullVector<int, Eigen::Dynamic> elementIncidence(2);
    for (int element = 0; element < NumElements; element++)
    {
        std::cout << "create element: " << element << " nodes: " << element << "," << element + 1 << std::endl;
        elementIncidence(0) = element;
        elementIncidence(1) = element + 1;
        structure.ElementCreate(interpolationType, elementIncidence);
        structure.ElementSetSection(element, Section1);
        //			structure.ElementSetConstitutiveLaw(element,Material1);
        if (element < 2)
        {
            structure.ElementSetConstitutiveLaw(element, Material1);
            output << structure.ConstitutiveLawGetParameterDouble(
                              Material1, NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS) << "\n";
        }
        else
        {
            structure.ElementSetConstitutiveLaw(element, Material2);
            output << structure.ConstitutiveLawGetParameterDouble(
                              Material2, NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS) << "\n";
        }
    }
    output.close();

    structure.ElementTotalConvertToInterpolationType(1e-6, 3);

    filename = "Truss1D2NReference-ipCoor";
    output.open(filename.c_str());
    if (output)
    {
        output << "#" << filename << "\n";
        output << "ipCoor\n";
        output << "1000\n2500\n3500\n5000\n";
        output.close();
    }
    else
    {
        std::cout << __LINE__ << " Truss1D2NReference] error output file.\n";
    }

    // set boundary conditions and loads
    NuTo::FullVector<double, Eigen::Dynamic> direction(1);
    direction(0) = 1;
    structure.ConstraintLinearSetDisplacementNode(0, direction, 0.0);
    if (EnableDisplacementControl)
    {
        std::cout << "Displacement control" << std::endl;
        structure.ConstraintLinearSetDisplacementNode(NumElements, direction, BoundaryDisplacement);
    }
    else
    {
        std::cout << "Load control" << std::endl;
        structure.LoadCreateNodeForce(1, NumElements, direction, Force);
    }
    // start analysis
    structure.CalculateMaximumIndependentSets();
    structure.SolveGlobalSystemStaticElastic(1);

    auto displacementVector = structure.NodeExtractDofValues(0);

    NuTo::FullVector<double, Eigen::Dynamic> rDisplacements;
    filename = "Truss1D2NReference-disp";
    output.open(filename.c_str());
    if (output)
    {
        output << "#" << filename << "\n";
        output << "disp\n 0.\n";
        if (EnableDisplacementControl)
        {
            for (int element = 1; element < NumElements; element++)
            {
                structure.NodeGetDisplacements(element, rDisplacements);
                output << rDisplacements << "\n";
            }
            output << "1.\n";
        }
        else
        {
            output << displacementVector;
        }
        output.close();
    }

    // calculate residual
    auto residual = structure.BuildGlobalInternalGradient() - structure.BuildGlobalExternalLoadVector(1);

    std::cout << "residual: " << residual.J.CalculateNormL2() << std::endl;

    std::cout << "element 0: strain\n" << structure.ElementGetEngineeringStrain(0) << std::endl;
    std::cout << "           stress\n" << structure.ElementGetEngineeringStress(0) << std::endl;

    std::cout << "element 0: strain\n" << structure.ElementGetEngineeringStrain(1) << std::endl;
    std::cout << "           stress\n" << structure.ElementGetEngineeringStress(1) << std::endl;

    std::cout << "element 0: strain\n" << structure.ElementGetEngineeringStrain(2) << std::endl;
    std::cout << "           stress\n" << structure.ElementGetEngineeringStress(2) << std::endl;

    std::cout << "element 0: strain\n" << structure.ElementGetEngineeringStrain(3) << std::endl;
    std::cout << "           stress\n" << structure.ElementGetEngineeringStress(3) << std::endl;

    filename = "Truss1D2NReference-strain00";
    output.open(filename.c_str());
    if (output)
    {
        output << "#" << filename << "\n";
        output << "strain00\n";
        for (int element = 0; element < NumElements; element++)
        {
            output << structure.ElementGetEngineeringStrain(element).GetRow(0) << "\n";
        }
        output.close();
    }
    else
    {
        std::cout << __LINE__ << " Truss1D2NReference] error output file.\n";
    }

    filename = "Truss1D2NReference-stress00";
    output.open(filename.c_str());
    if (output)
    {
        output << "#" << filename << "\n";
        output << "stress00\n";
        for (int element = 0; element < NumElements; element++)
        {
            output << structure.ElementGetEngineeringStress(element).GetRow(0) << "\n";
        }
        output.close();
    }
    else
    {
        std::cout << __LINE__ << " Truss1D2NReference] error output file.\n";
    }

    // visualize results
    int visualizationGroup = structure.GroupCreate(NuTo::eGroupId::Elements);
    structure.GroupAddElementsTotal(visualizationGroup);

    structure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::DISPLACEMENTS);
    structure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::ENGINEERING_STRAIN);
    structure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::ENGINEERING_STRESS);
    structure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::SECTION);
    structure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::CONSTITUTIVE);

    structure.ExportVtkDataFileElements("Truss1D2NReference.vtk");

    return 0;
}
