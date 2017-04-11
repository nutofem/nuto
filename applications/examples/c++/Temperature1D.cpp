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
#include "mechanics/MechanicsEnums.h"
#include "mechanics/constraints/ConstraintCompanion.h"
#include "mechanics/dofSubMatrixStorage/BlockScalar.h"
#include "mechanics/mesh/MeshGenerator.h"
#include "mechanics/structures/unstructured/Structure.h"
#include "mechanics/sections/SectionTruss.h"
#include "mechanics/structures/StructureOutputBlockVector.h"
#include "visualize/VisualizeEnum.h"

using namespace NuTo;

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
    Structure structure(1);

    int interpolationType = MeshGenerator::Grid(structure, {length}, {num_elements}).second;

    structure.InterpolationTypeAdd(interpolationType, Node::eDof::TEMPERATURE,
                                     Interpolation::eTypeOrder::EQUIDISTANT1);
    structure.InterpolationTypeSetIntegrationType(interpolationType, "IntegrationType1D2NGauss2Ip");
    structure.ElementTotalConvertToInterpolationType();

    // create section
    auto truss = SectionTruss::Create(area);
    structure.ElementTotalSetSection(truss);

    // create material law
    auto material = structure.ConstitutiveLawCreate("Heat_Conduction");
    structure.ConstitutiveLawSetParameterDouble(material, 
        Constitutive::eConstitutiveParameter::THERMAL_CONDUCTIVITY, conductivity);
    structure.ElementTotalSetConstitutiveLaw(material);

    // set boundary conditions and loads
    auto& origin = structure.NodeGetAtCoordinate(Eigen::Matrix<double, 1, 1>::Zero());
    structure.Constraints().Add(Node::eDof::TEMPERATURE, Constraint::Value(origin, boundary_temperature));

    structure.SetNumLoadCases(1);
    Eigen::VectorXd direction(1);
    direction(0) = 1.0;
    structure.LoadCreateNodeHeatFlux(0, num_elements, direction, boundary_flux);

    // start analysis
    structure.SolveGlobalSystemStaticElastic(0);
    auto residual = structure.BuildGlobalInternalGradient() - structure.BuildGlobalExternalLoadVector(0);
    std::cout << "residual: " << residual.J.CalculateNormL2() << std::endl;

    // visualize results
    int visualizationGroup = structure.GroupCreate(eGroupId::Elements);
    structure.GroupAddElementsTotal(visualizationGroup);

    structure.AddVisualizationComponent(visualizationGroup, eVisualizeWhat::TEMPERATURE);
    structure.AddVisualizationComponent(visualizationGroup, eVisualizeWhat::HEAT_FLUX);
    structure.ExportVtkDataFileElements("Temperature1D.vtk");

    return 0;
}
