/*! @file Temperature1D.cpp
 * Heat Conduction in 1D, with one Dirichlet and one Neumann boundary condition.
 *
 * Solve the heat equation in 1D, \f$ ρ c_T \dot{u} - k u,_{xx} = 0 \f$, under
 * the boundary conditions \f$u|_{x = 0} = 20\f$ and \f$q|_{x=l} = 10 \f$.
 * @image html Temperature1D.png
 */

#include <iostream>
#include "nuto/math/SparseMatrixCSRGeneral.h"
#include "nuto/math/SparseDirectSolverMUMPS.h"
#include "nuto/mechanics/MechanicsEnums.h"
#include "nuto/mechanics/constraints/ConstraintCompanion.h"
#include "nuto/mechanics/dofSubMatrixStorage/BlockScalar.h"
#include "nuto/mechanics/mesh/MeshGenerator.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/mechanics/sections/SectionTruss.h"
#include "nuto/mechanics/structures/StructureOutputBlockVector.h"
#include "nuto/visualize/VisualizeEnum.h"

using namespace NuTo;

int main()
{
    // geometry/mesh
    double area = 1.0;
    double length = 1.0;
    int num_elements = 10;

    // boundaries
    double boundary_temperature = 20.0;

    // material
    double conductivity = 1.0;

    //! create one-dimensional structure
    Structure structure(1);

    int interpolationType = MeshGenerator::Grid(structure, {length}, {num_elements}).second;

    structure.InterpolationTypeAdd(interpolationType, Node::eDof::TEMPERATURE, Interpolation::eTypeOrder::EQUIDISTANT1);
    structure.InterpolationTypeSetIntegrationType(interpolationType, "IntegrationType1D2NGauss2Ip");
    structure.ElementTotalConvertToInterpolationType();

    // create section
    auto truss = SectionTruss::Create(area);
    structure.ElementTotalSetSection(truss);

    // create material law
    auto material = structure.ConstitutiveLawCreate("Heat_Conduction");
    structure.ConstitutiveLawSetParameterDouble(material, Constitutive::eConstitutiveParameter::THERMAL_CONDUCTIVITY,
                                                conductivity);
    structure.ElementTotalSetConstitutiveLaw(material);

    // set boundary conditions and loads
    auto& origin = structure.NodeGetAtCoordinate(0);
    structure.Constraints().Add(Node::eDof::TEMPERATURE, Constraint::Value(origin, boundary_temperature));

    auto& nodeRight = structure.NodeGetAtCoordinate(length);
    structure.Constraints().Add(Node::eDof::TEMPERATURE, Constraint::Value(nodeRight, boundary_temperature + 20));

    // start analysis
    structure.SolveGlobalSystemStaticElastic();
    auto residual = structure.BuildGlobalInternalGradient() - structure.BuildGlobalExternalLoadVector();
    std::cout << "residual: " << residual.J.CalculateNormL2() << std::endl;

    // visualize results
    int visualizationGroup = structure.GroupCreate(eGroupId::Elements);
    structure.GroupAddElementsTotal(visualizationGroup);

    structure.AddVisualizationComponent(visualizationGroup, eVisualizeWhat::TEMPERATURE);
    structure.AddVisualizationComponent(visualizationGroup, eVisualizeWhat::HEAT_FLUX);
    structure.ExportVtkDataFileElements("Temperature1D.vtk");

    return 0;
}
