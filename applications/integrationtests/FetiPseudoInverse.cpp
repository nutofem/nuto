//============================================================================
// Name        : FetiPseudoInverse.cpp
// Author      : Philip Huschke
// Version     : 21 Nov 2016
// Copyright   :
// Description : Modified LU factorization of a singular matrix
//               to obtain rigid body modes
//
//
//============================================================================


#include "nuto/math/FullVector.h"
#include "nuto/mechanics/constitutive/ConstitutiveEnum.h"
#include "nuto/mechanics/elements/ElementBase.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeEnum.h"
#include "nuto/mechanics/interpolationtypes/InterpolationTypeEnum.h"
#include "nuto/mechanics/groups/GroupEnum.h"
#include "nuto/mechanics/nodes/NodeEnum.h"
#include "nuto/mechanics/sections/SectionEnum.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/mechanics/timeIntegration/NewmarkDirect.h"
#include "nuto/visualize/VisualizeEnum.h"
#include "nuto/mechanics/structures/StructureOutputBlockMatrix.h"
#include "eigen3/Eigen/Dense"

#include <boost/filesystem.hpp>
#include <iostream>
#include <fstream>
#include <string>

const Eigen::Vector2d directionX = Eigen::Vector2d::UnitX();
const Eigen::Vector2d directionY = Eigen::Vector2d::UnitY();

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main()
{
    NuTo::Structure structure(2);

    Eigen::Vector2d coords;

    // create node 0
    coords(0) = 0.0;
    coords(1) = 0.0;
    structure.NodeCreate(0, coords);

    // create node 1
    coords(0) = 1.0;
    coords(1) = 0.0;
    structure.NodeCreate(1, coords);

    // create node 2
    coords(0) = 1.0;
    coords(1) = 1.0;
    structure.NodeCreate(2, coords);

    // create node 3
    coords(0) = 0.0;
    coords(1) = 1.0;
    structure.NodeCreate(3, coords);

    int interpolationType = structure.InterpolationTypeCreate("Quad2D");
    structure.InterpolationTypeAdd(interpolationType, NuTo::Node::eDof::COORDINATES,   NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    structure.InterpolationTypeAdd(interpolationType, NuTo::Node::eDof::DISPLACEMENTS, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);

    Eigen::Vector4i nodes;
    nodes(0) = 0;
    nodes(1) = 1;
    nodes(2) = 2;
    nodes(3) = 3;
    structure.ElementCreate(interpolationType, nodes);
    structure.ElementTotalConvertToInterpolationType();

    // section
    int sectionId = structure.SectionCreate(NuTo::eSectionType::PLANE_STRESS);
    structure.SectionSetThickness(sectionId, 1.0);
    structure.ElementTotalSetSection(sectionId);

    // material
    int materialId = structure.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
    structure.ConstitutiveLawSetParameterDouble(materialId, NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, 8.0);
    structure.ConstitutiveLawSetParameterDouble(materialId, NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, 0);
    structure.ElementTotalSetConstitutiveLaw(materialId);

    //    // constraints
    //    int groupNodesLeftBoundary = structure.GroupCreate(NuTo::eGroupId::Nodes);

    //    structure.GroupAddNodeCoordinateRange(groupNodesLeftBoundary,0,-1.e-6,+1.e-6);

    //    structure.ConstraintLinearSetDisplacementNodeGroup(groupNodesLeftBoundary, directionX, 0);
    //    structure.ConstraintLinearSetDisplacementNodeGroup(groupNodesLeftBoundary, directionY, 0);


    NuTo::StructureOutputBlockMatrix stiffnessMatrix = structure.BuildGlobalHessian0();


    Eigen::MatrixXd  mat = stiffnessMatrix.JJ.ExportToFullMatrix();


    structure.NodeInfo(10);




    std::cout << "stiffnessMatrix: \n" << mat << std::endl;


}
