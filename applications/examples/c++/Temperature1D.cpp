/*! @file Temperature1D.cpp
 * Heat Conduction in 1D, with one Dirichlet and one Neumann boundary condition.
 *
 * Solve the heat equation in 1D, \f$ ρ c_T \dot{u} - k u,_{xx} = 0 \f$, under
 * the boundary conditions \f$u|_{x = 0} = 20\f$ and \f$q|_{x=l} = 10 \f$.
 * @image html Temperature1D.png
 */

#include <iostream>
#include "math/SparseMatrixCSRGeneral.h"
#include "math/SparseDirectSolverMUMPS.h"
#include "mechanics/structures/unstructured/Structure.h"
#include "mechanics/structures/StructureOutputBlockVector.h"
#include "mechanics/dofSubMatrixStorage/BlockScalar.h"
#include "mechanics/MechanicsEnums.h"
#include "visualize/VisualizeEnum.h"
#include "mechanics/mesh/MeshGenerator.h"

int main()
{
    // geometry/mesh
    double area = 1.0;
    double length = 1.0;
    int num_elements = 10;

    // boundaries
    double boundary_temperature = 20.0;
    double boundary_flux = 10.0;

    // material
    double conductivity = 1.0;

    //! create one-dimensional structure
    NuTo::Structure myStructure(1);

    int interpolationType = NuTo::MeshGenerator::Grid(myStructure, {length}, {num_elements}).second;

    myStructure.InterpolationTypeAdd(interpolationType, NuTo::Node::eDof::TEMPERATURE,
                                     NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    myStructure.InterpolationTypeSetIntegrationType(interpolationType, "IntegrationType1D2NGauss2Ip");
    myStructure.ElementTotalConvertToInterpolationType();

    // create section
    auto truss = myStructure.SectionCreate("Truss");
    myStructure.SectionSetArea(truss, area);
    myStructure.ElementTotalSetSection(truss);

    // create material law
    auto material = myStructure.ConstitutiveLawCreate("Heat_Conduction");
    myStructure.ConstitutiveLawSetParameterDouble(material, 
        NuTo::Constitutive::eConstitutiveParameter::THERMAL_CONDUCTIVITY, conductivity);
    myStructure.ElementTotalSetConstitutiveLaw(material);

    // set boundary conditions and loads
    Eigen::VectorXd direction(1);
    direction(0) = 1.0;
    myStructure.ConstraintLinearSetTemperatureNode(0, boundary_temperature);
    myStructure.SetNumLoadCases(1);
    myStructure.LoadCreateNodeHeatFlux(0, num_elements, direction, boundary_flux);

    // start analysis
    myStructure.SolveGlobalSystemStaticElastic(0);
    auto residual = myStructure.BuildGlobalInternalGradient() - myStructure.BuildGlobalExternalLoadVector(0);
    std::cout << "residual: " << residual.J.CalculateNormL2() << std::endl;

    // visualize results
    int visualizationGroup = myStructure.GroupCreate(NuTo::eGroupId::Elements);
    myStructure.GroupAddElementsTotal(visualizationGroup);

    myStructure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::TEMPERATURE);
    myStructure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::HEAT_FLUX);
    myStructure.ExportVtkDataFileElements("Temperature1D.vtk");

    return 0;
}
