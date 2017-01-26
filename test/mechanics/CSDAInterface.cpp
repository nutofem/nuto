/*
 * CSDAInterface.cpp
 *
 *  Created on: 28 October 2016
 *      Author: Thomas Titscher
 */

#include <cmath>
#include <boost/filesystem/operations.hpp>
#include <mechanics/interpolationtypes/InterpolationType.h>
#include "visualize/VisualizeEnum.h"
#include "mechanics/constitutive/ConstitutiveEnum.h"
#include "mechanics/groups/GroupEnum.h"
#include "mechanics/structures/unstructured/Structure.h"
#include "mechanics/interpolationtypes/InterpolationTypeEnum.h"
#include "mechanics/interpolationtypes/InterpolationBase.h"
#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/nodes/NodeBase.h"
#include "mechanics/structures/StructureOutputBlockVector.h"
#include "mechanics/tools/GlobalFractureEnergyIntegrator.h"
#include "mechanics/timeIntegration/NewmarkDirect.h"
#include "mechanics/elements/ElementBase.h"
#include "mechanics/elements/ContinuumElement.h"

/*               3   2
 *   /|          /  /           \
 *   /|         /e0/          ---\
 *   /|        /  /           ---/
 *   /|      _/  /              /
 *           0   1
 *
 *  + interface angle & thickness
 *  - angle 90 = vertical
 *  - (0,0) is at the middle of the structure
 *
 */

int FindLocalElementIndex(int rGlobalDof, Eigen::VectorXi rGlobalDofs)
{
    for (int i = 0; i < rGlobalDofs.rows(); ++i)
    {
        if (rGlobalDofs[i] == rGlobalDof)
            return i;
    }
    throw;
}

void CheckFractureEnergy2D(int rAngleDegree, double rInterfaceThickness)
{
    NuTo::Structure s(2);
    s.SetShowTime(false);
    s.GetLogger().SetQuiet(true);

    const double ly2 = 2.;  // half of Length_y
    const double lz = 6.;

    const double angleRad = M_PI / 180. * rAngleDegree;

    const double projectedThickness = rInterfaceThickness / std::sin(angleRad);
    const double xInterfaceOffset =  ly2 / std::tan(angleRad);
    const double interfaceLength = 2 * ly2 / std::sin(angleRad);

    s.GetLogger() << "projectedThickness " << projectedThickness << '\n';
    s.GetLogger() << "xInterfaceOffset " << xInterfaceOffset << '\n';
    s.GetLogger() << "interfaceLength " << interfaceLength << '\n';

    // lower nodes
    s.NodeCreate(0, Eigen::Vector2d({-xInterfaceOffset - projectedThickness / 2., -ly2}));
    s.NodeCreate(1, Eigen::Vector2d({-xInterfaceOffset + projectedThickness / 2., -ly2}));

    // upper nodes
    s.NodeCreate(2, Eigen::Vector2d({+xInterfaceOffset + projectedThickness / 2.,  ly2}));
    s.NodeCreate(3, Eigen::Vector2d({+xInterfaceOffset - projectedThickness / 2.,  ly2}));


    int it = s.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::QUAD2D);
    s.InterpolationTypeAdd(it, NuTo::Node::eDof::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    s.InterpolationTypeAdd(it, NuTo::Node::eDof::DISPLACEMENTS, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);


    std::vector<int> ids({0, 1, 2, 3});
    s.ElementCreate(0, it, ids);

    s.ElementTotalConvertToInterpolationType();

    int mySection = s.SectionCreate("Plane_Stress");
    s.SectionSetThickness(mySection, lz);
    s.ElementTotalSetSection(mySection);
    s.ElementTotalConvertToInterpolationType();

    using namespace NuTo::Constitutive;
    s.ConstitutiveLawCreate(0, eConstitutiveType::LOCAL_DAMAGE_MODEL);

    constexpr double fractureEnergy         = 0.1;

    s.ConstitutiveLawSetParameterDouble(0, eConstitutiveParameter::YOUNGS_MODULUS,       20000.);
    s.ConstitutiveLawSetParameterDouble(0, eConstitutiveParameter::POISSONS_RATIO,       0.0);
    s.ConstitutiveLawSetParameterDouble(0, eConstitutiveParameter::TENSILE_STRENGTH,     4.);
    s.ConstitutiveLawSetParameterDouble(0, eConstitutiveParameter::COMPRESSIVE_STRENGTH, 4.);
    s.ConstitutiveLawSetParameterDouble(0, eConstitutiveParameter::FRACTURE_ENERGY,      fractureEnergy / rInterfaceThickness);
    s.ConstitutiveLawSetDamageLaw(0, eDamageLawType::ISOTROPIC_EXPONENTIAL_SOFTENING);

    s.ElementSetConstitutiveLaw(0, 0);

    s.NodeBuildGlobalDofs();
    int dofBC1 = s.NodeGetNodePtr(1)->GetDof(NuTo::Node::eDof::DISPLACEMENTS, 0);
    int dofBC2 = s.NodeGetNodePtr(2)->GetDof(NuTo::Node::eDof::DISPLACEMENTS, 0);

    int localDofIndex1 = FindLocalElementIndex(dofBC1, s.ElementBuildGlobalDofsRow(0)[NuTo::Node::eDof::DISPLACEMENTS]);
    int localDofIndex2 = FindLocalElementIndex(dofBC2, s.ElementBuildGlobalDofsRow(0)[NuTo::Node::eDof::DISPLACEMENTS]);

    int numLoadSteps = 200;
    double bcEnd = 0.4;

    Eigen::VectorXd displ(numLoadSteps+1);
    Eigen::VectorXd force(numLoadSteps+1);

    auto globalDofs = s.NodeExtractDofValues(0);
    auto& globalDisplacementDofs = globalDofs.J[NuTo::Node::eDof::DISPLACEMENTS];

    for (int i = 0; i < numLoadSteps+1; ++i)
    {
        double bc = bcEnd * i / (numLoadSteps);
        displ[i] = bc;
        globalDisplacementDofs[dofBC1] = bc;
        globalDisplacementDofs[dofBC2] = bc;

        s.StructureBase::NodeMergeDofValues(globalDofs);
        auto internalForces = s.ElementBuildInternalGradient(0)[NuTo::Node::eDof::DISPLACEMENTS];
        force[i] = (internalForces[localDofIndex1] + internalForces[localDofIndex2]);
    }

//    std::cout << force << std::endl;

    NuTo::Tools::GlobalFractureEnergyIntegrator integrator(force, displ);
    double crackArea = lz * interfaceLength;
    double globalFractureEnergy = integrator.IntegrateSofteningCurve(crackArea, 0.01);
    double error = std::abs(fractureEnergy - globalFractureEnergy);
    double tolerance = fractureEnergy / 10.;

    std::cout << "angle: " << rAngleDegree << "\t thickness: " << rInterfaceThickness << "\t GF: " << globalFractureEnergy << "\t Error: " << error << std::endl;
    if (error > tolerance)
    {
        throw;
    }
}


void CSDA2D()
{

    /*        \/
     * 3    2  7     6
     *
     *
     * 0    1  4     5
     *
     */


    constexpr double thickness2 = 0.1;
    constexpr double lx2 = 10;
    constexpr double ly = 5;
    constexpr double lz = 2;



    NuTo::Structure s(2);

    s.NodeCreate(0, Eigen::Vector2d({-lx2, 0}));
    s.NodeCreate(1, Eigen::Vector2d({-thickness2, 0}));
    s.NodeCreate(2, Eigen::Vector2d({-thickness2, ly}));
    s.NodeCreate(3, Eigen::Vector2d({-lx2, ly}));

    s.NodeCreate(4, Eigen::Vector2d({thickness2, 0}));
    s.NodeCreate(5, Eigen::Vector2d({lx2, 0}));
    s.NodeCreate(6, Eigen::Vector2d({lx2, ly}));
    s.NodeCreate(7, Eigen::Vector2d({thickness2, ly}));


    int it = s.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::QUAD2D);
    s.InterpolationTypeAdd(it, NuTo::Node::eDof::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    s.InterpolationTypeAdd(it, NuTo::Node::eDof::DISPLACEMENTS, NuTo::Interpolation::eTypeOrder::EQUIDISTANT2);

    int it2 = s.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::TRIANGLE2D);
    s.InterpolationTypeAdd(it2, NuTo::Node::eDof::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    s.InterpolationTypeAdd(it2, NuTo::Node::eDof::DISPLACEMENTS, NuTo::Interpolation::eTypeOrder::EQUIDISTANT2);

    s.ElementCreate(0, it, {0, 1, 2, 3});
    s.ElementCreate(1, it, {4, 5, 6, 7});

     s.ElementCreate(2, it, {1,4,7,2});
//    s.ElementCreate(2, it2, {1, 4, 7});
//    s.ElementCreate(3, it2, {1, 7, 2});

    using namespace NuTo::Constitutive;
    int LIN = 0;
    int CSDA = 1;
    s.ConstitutiveLawCreate(LIN, eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
    s.ConstitutiveLawSetParameterDouble(LIN, eConstitutiveParameter::YOUNGS_MODULUS,       20000.);
    s.ConstitutiveLawSetParameterDouble(LIN, eConstitutiveParameter::POISSONS_RATIO,       0.0);

    s.ConstitutiveLawCreate(CSDA, eConstitutiveType::LOCAL_DAMAGE_MODEL);
    constexpr double fractureEnergy         = 0.1;
    s.ConstitutiveLawSetParameterDouble(CSDA, eConstitutiveParameter::YOUNGS_MODULUS,       20000.);
    s.ConstitutiveLawSetParameterDouble(CSDA, eConstitutiveParameter::POISSONS_RATIO,       0.0);
    s.ConstitutiveLawSetParameterDouble(CSDA, eConstitutiveParameter::TENSILE_STRENGTH,     4.);
    s.ConstitutiveLawSetParameterDouble(CSDA, eConstitutiveParameter::COMPRESSIVE_STRENGTH, 40.);
    s.ConstitutiveLawSetParameterDouble(CSDA, eConstitutiveParameter::FRACTURE_ENERGY,      fractureEnergy / (thickness2*2.));
    s.ConstitutiveLawSetDamageLaw(CSDA, eDamageLawType::ISOTROPIC_EXPONENTIAL_SOFTENING);

    s.ElementSetConstitutiveLaw(0, LIN);
    s.ElementSetConstitutiveLaw(1, LIN);
    s.ElementSetConstitutiveLaw(2, CSDA);
//    s.ElementSetConstitutiveLaw(3, CSDA);

    int mySection = s.SectionCreate("Plane_Strain");
    s.SectionSetThickness(mySection, lz);
    s.ElementTotalSetSection(mySection);
    s.ElementTotalConvertToInterpolationType();


    int nodeFixXY = s.NodeGetIdAtCoordinate(Eigen::Vector2d({-lx2, 0}), 1.e-5);
    int nodeFixY = s.NodeGetIdAtCoordinate(Eigen::Vector2d({lx2, 0}), 1.e-5);
    int nodeBC = s.NodeGetIdAtCoordinate(Eigen::Vector2d({thickness2, ly}), 1.e-5);

    s.ConstraintLinearSetDisplacementNode(nodeFixXY, Eigen::Vector2d({1,0}), 0);
    s.ConstraintLinearSetDisplacementNode(nodeFixXY, Eigen::Vector2d({0,1}), 0);
    s.ConstraintLinearSetDisplacementNode(nodeFixY, Eigen::Vector2d({0,1}), 0);

    int BC = s.ConstraintLinearSetDisplacementNode(nodeBC, Eigen::Vector2d({0,1}), 0);

    s.NodeBuildGlobalDofs();
    std::cout << s.GetNumTotalActiveDofs() << std::endl;
    std::cout << s.GetNumTotalDependentDofs() << std::endl;

    double deltaD = .5;

    Eigen::Matrix2d dispRHS;
    dispRHS << 0, 0, 1, -deltaD;

//    s.AddVisualizationComponent(s.GroupGetElementsTotal(), NuTo::eVisualizeWhat::DISPLACEMENTS);

    NuTo::NewmarkDirect newmark(&s);

    s.SetShowTime(false);
    newmark.SetShowTime(false);

    newmark.SetTimeStep(0.1);
    newmark.SetMinTimeStep(0.001);
    newmark.SetMaxTimeStep(0.1);
    newmark.SetToleranceForce(1e-6);
    newmark.SetAutomaticTimeStepping(true);
    newmark.SetPerformLineSearch(true);
    newmark.SetMaxNumIterations(20);

    bool deleteDirectory = true;
    newmark.SetResultDirectory("./CSDA2D", deleteDirectory);

    newmark.AddTimeDependentConstraint(BC, dispRHS);

//    newmark.AddResultNodeDisplacements("Displ", nodeBC);
//    int groupNodeBC = s.GroupCreate(NuTo::eGroupId::Nodes);
//    s.GroupAddNode(groupNodeBC, nodeBC);
//    newmark.AddResultGroupNodeForce("Force", groupNodeBC);


    newmark.Solve(1);
}


int GetNodeGroupFromElements(NuTo::Structure& rS, int rG)
{
    int n = rS.GroupCreate(NuTo::eGroupId::Nodes);
    rS.GroupAddNodesFromElements(n, rG);
    return n;
}

std::vector<NuTo::NodeBase*> GetPrismNodes(NuTo::Structure &rS, int rGroupMaster, int rGroupSlave)
{
    int nMaster = GetNodeGroupFromElements(rS, rGroupMaster);
    int nSlave = GetNodeGroupFromElements(rS, rGroupSlave);
    int gNodesPrism = rS.GroupIntersection(nMaster, nSlave);
    rS.GroupDelete(nMaster);
    rS.GroupDelete(nSlave);

    std::vector<NuTo::NodeBase*> v;
    v.reserve(rS.GroupGetNumMembers(gNodesPrism));
    for (int nodeId : rS.GroupGetMemberIds(gNodesPrism))
        v.push_back(rS.NodeGetNodePtr(nodeId));

    return v;
}

struct ElementSurface
{
    NuTo::ElementBase* elem;
    int surface;
};

std::vector<NuTo::NodeBase*> GetSurfaceNodes(ElementSurface rElementSurface)
{
    Eigen::VectorXi surfaceNodeIndices = rElementSurface.elem->GetInterpolationType().GetSurfaceNodeIndices(rElementSurface.surface);
    int numSurfaceNodes = surfaceNodeIndices.rows();
    std::vector<NuTo::NodeBase *> surfaceNodes(numSurfaceNodes);

    for (int iSurfaceNode = 0; iSurfaceNode < numSurfaceNodes; ++iSurfaceNode)
    {
        surfaceNodes[iSurfaceNode] = rElementSurface.elem->GetNode(surfaceNodeIndices(iSurfaceNode, 0));
    }
    return surfaceNodes;
}

int FindSurfaceId(NuTo::ElementBase* rElement, std::vector<NuTo::NodeBase*> rNodes)
{
    const auto& it = rElement->GetInterpolationType();
    for (int iSurface = 0; iSurface < it.GetNumSurfaces(); ++iSurface)
    {
        bool elementSurfaceNodesMatchBoundaryNodes = true;
        Eigen::VectorXi surfaceNodeIndices = it.GetSurfaceNodeIndices(iSurface);

        std::vector<NuTo::NodeBase *> surfaceNodes = GetSurfaceNodes({rElement, iSurface});

        //check, if all surface nodes are in the node group
        for (unsigned int countNode = 0; countNode < surfaceNodes.size(); countNode++)
        {
            if (std::find(rNodes.begin(), rNodes.end(), surfaceNodes[countNode]) == rNodes.end())
            {
                //this surface has at least one node that is not in the list, continue
                elementSurfaceNodesMatchBoundaryNodes = false;
                break;
            }
        }
        if (elementSurfaceNodesMatchBoundaryNodes)
            return iSurface;
    }
    return -42; // error code.
}

std::vector<ElementSurface> GetElementSurfaceVector(NuTo::Structure& rS, int rG, std::vector<NuTo::NodeBase*>& rNodes)
{
    std::vector<ElementSurface> v;
    v.reserve(rS.GroupGetNumMembers(rG));
    for (int elementId : rS.GroupGetMemberIds(rG))
    {
        NuTo::ElementBase* e = rS.ElementGetElementPtr(elementId);
        int surfaceId = FindSurfaceId(e, rNodes);
        if (surfaceId >= 0)
            v.push_back({e, surfaceId});
        // else: element is not part of the surface
    }
    return v;
}

bool AreEqual(const std::vector<NuTo::NodeBase*>& r1, const std::vector<NuTo::NodeBase*>& r2)
{
    for (const auto& n1 : r1)
    {
        bool n1IsInr2 = false;
        for (const auto& n2 : r2)
        {
            if (n1 == n2)
                n1IsInr2 = true;
        }
        if (not n1IsInr2)
            return false;
    }
    return true;
}

std::vector<std::pair<ElementSurface, ElementSurface>> FindMatchingElements(NuTo::Structure& rS, int rGroupMaster, int rGroupSlave)
{
    std::vector<NuTo::NodeBase*> gNodesPrism = GetPrismNodes(rS, rGroupMaster, rGroupSlave);

    std::vector<ElementSurface> eMaster = GetElementSurfaceVector(rS, rGroupMaster, gNodesPrism);
    std::vector<ElementSurface> eSlave = GetElementSurfaceVector(rS, rGroupSlave, gNodesPrism);

    std::vector<std::pair<ElementSurface, ElementSurface>> pairs;
    pairs.reserve(eMaster.size());

    for (auto& esMaster : eMaster)
    {
        bool surfaceFound = false;
        std::vector<NuTo::NodeBase*> surfaceNodesMaster = GetSurfaceNodes(esMaster);
        for (auto& esSlave : eSlave)
        {
            std::vector<NuTo::NodeBase*> surfaceNodesSlave = GetSurfaceNodes(esSlave);
            if (AreEqual(surfaceNodesSlave, surfaceNodesMaster))
            {
                surfaceFound = true;
                pairs.push_back(std::make_pair(esMaster, esSlave));
            }
        }
        if (not surfaceFound)
            throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "No matching Slaveate surface found.");
    }
    return pairs;
}

bool HasOnlyCoordinateInterpolation(NuTo::Structure& rS, int gElement)
{
    for (int elementId : rS.GroupGetMemberIds(gElement))
    {
        auto* e = rS.ElementGetElementPtr(elementId);
        auto dofs = e->GetInterpolationType().GetDofs();
        if (dofs.size() != 1 or *dofs.begin() != NuTo::Node::eDof::COORDINATES)
            return false;
    }
    return true;
}

std::vector<std::pair<ElementSurface, ElementSurface>> FindElementPairsContainingTheNode(const NuTo::NodeBase* rNode, const std::vector<std::pair<ElementSurface, ElementSurface>>& rPairs)
{
    std::vector<std::pair<ElementSurface, ElementSurface>> matchingPairs;
    for (const auto& pair : rPairs)
    {
        NuTo::ElementBase* element = pair.first.elem;
        for (int i = 0; i < element->GetNumNodes(); ++i)
        {
            if (element->GetNode(i) == rNode)
            {
                matchingPairs.push_back(pair);
                continue;
            }
        }
    }
    return matchingPairs;
}

int GetNodeCoordinatesIndex(const NuTo::NodeBase* rNode, const NuTo::ElementBase* rElement)
{
    for (int i = 0; i < rElement->GetNumNodes(NuTo::Node::eDof::COORDINATES); ++i)
    {
        if (rElement->GetNode(i, NuTo::Node::eDof::COORDINATES) == rNode)
            return i;
    }
    throw;
}

Eigen::VectorXd GetLocalSurfaceCoordinates(int rNodeIndex, const ElementSurface& rElementSurface)
{
    const auto& it = rElementSurface.elem->GetInterpolationType().Get(NuTo::Node::eDof::COORDINATES);

    // find surface parameters
    auto localNodeCoordinates = it.GetNaturalNodeCoordinates(rNodeIndex);
    Eigen::VectorXd R = it.CalculateNaturalSurfaceCoordinates(Eigen::VectorXd::Zero(2), rElementSurface.surface) - localNodeCoordinates;
    Eigen::MatrixXd dRdS = it.CalculateDerivativeNaturalSurfaceCoordinates(Eigen::VectorXd::Zero(2), rElementSurface.surface);

    Eigen::JacobiSVD<Eigen::MatrixXd> svd(dRdS, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::VectorXd surfaceParameters = -svd.solve(R);
    assert((localNodeCoordinates - it.CalculateNaturalSurfaceCoordinates(surfaceParameters, rElementSurface.surface)).norm() < 1.e-10);
//    std::cout << "Surface parameters " << surfaceParameters.transpose() << std::endl;
    return surfaceParameters;
}

Eigen::VectorXd GetLocalSurfaceCoordinates(const NuTo::NodeBase* rNode, const ElementSurface& rElementSurface)
{
    return GetLocalSurfaceCoordinates(GetNodeCoordinatesIndex(rNode, rElementSurface.elem), rElementSurface);
}




Eigen::VectorXd CalculateNormalAtNode(const NuTo::NodeBase* rNode, const ElementSurface& rElementSurface)
{
    Eigen::MatrixXd nodeCoordinates = rElementSurface.elem->ExtractNodeValues(0, NuTo::Node::eDof::COORDINATES);
    const auto& it = rElementSurface.elem->GetInterpolationType().Get(NuTo::Node::eDof::COORDINATES);
    Eigen::VectorXd ipCoordsSurface = GetLocalSurfaceCoordinates(rNode, rElementSurface);
    Eigen::VectorXd ipCoordsNatural = it.CalculateNaturalSurfaceCoordinates(ipCoordsSurface, rElementSurface.surface);

    Eigen::MatrixXd derivativeShapeFunctionsNatural = it.CalculateDerivativeShapeFunctionsNatural(ipCoordsNatural);
    const Eigen::Matrix3d jacobian = rElementSurface.elem->AsContinuumElement3D().CalculateJacobian(derivativeShapeFunctionsNatural, nodeCoordinates);

    Eigen::MatrixXd derivativeNaturalSurfaceCoordinates = it.CalculateDerivativeNaturalSurfaceCoordinates(ipCoordsSurface, rElementSurface.surface); // = [dXi / dAlpha]
    Eigen::Vector3d dXdAlpha = jacobian * derivativeNaturalSurfaceCoordinates.col(0);
    Eigen::Vector3d dXdBeta  = jacobian * derivativeNaturalSurfaceCoordinates.col(1);

//    std::cout << "nodeCoordinates \n" << nodeCoordinates << std::endl;
//    std::cout << "jacobian \n" << jacobian << std::endl;
//    std::cout << "derivativeNaturalSurfaceCoordinates \n" << derivativeNaturalSurfaceCoordinates << std::endl;
//    std::cout << "dXdAlpha \n" << dXdAlpha.transpose() << std::endl;
//    std::cout << "dXdBeta \n" << dXdBeta.transpose() << std::endl;


    Eigen::Vector3d surfaceNormalVector = dXdAlpha.cross(dXdBeta); // = || [dX / dXi] * [dXi / dAlpha] ||
    surfaceNormalVector.normalize();
    return surfaceNormalVector;
}

NuTo::NodeBase* CloneNode(NuTo::Structure& rS, const NuTo::NodeBase* rNode)
{
    int nodeId = rS.NodeCreate(rNode->Get(NuTo::Node::eDof::COORDINATES));
    return rS.NodeGetNodePtr(nodeId);
}

NuTo::Interpolation::eTypeOrder GetCoordinateInterpolation(NuTo::Structure& rS, int rGroupMaster, int rGroupSlave)
{
    auto* firstElement = rS.ElementGetElementPtr(rS.GroupGetMemberIds(rGroupMaster)[0]);
    NuTo::Interpolation::eTypeOrder coordinateInterpolation = firstElement->GetInterpolationType().Get(NuTo::Node::eDof::COORDINATES).GetTypeOrder();
    for (int elementId : rS.GroupGetMemberIds(rGroupMaster))
    {
        auto* e = rS.ElementGetElementPtr(elementId);
        auto type = e->GetInterpolationType().Get(NuTo::Node::eDof::COORDINATES).GetTypeOrder();
        if (type != coordinateInterpolation)
            throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "All elements in the groups must have the same coordinate interpolation.");
    }
    return coordinateInterpolation;
}

NuTo::NodeBase* InterpolateNode(
    NuTo::Structure& rS,
    NuTo::NodeBase* rMaster,
    NuTo::NodeBase* rSlave,
    std::map<NuTo::NodeBase*, NuTo::NodeBase*>& rMasterNodeToInterpolatedNode)
{
    auto it = rMasterNodeToInterpolatedNode.find(rMaster);
    if (it != rMasterNodeToInterpolatedNode.end())
        return it->second;


    Eigen::Vector3d coords = (rMaster->Get(NuTo::Node::eDof::COORDINATES) + rSlave->Get(NuTo::Node::eDof::COORDINATES)) / 2.;
    int nodeId = rS.NodeCreate(coords);

    NuTo::NodeBase* newNode = rS.NodeGetNodePtr(nodeId);
    rMasterNodeToInterpolatedNode[rMaster] = newNode;
    return newNode;
}

int CreatePrisms(NuTo::Structure& rS, int rGroupMaster, int rGroupSlave, double rThickness)
{
    if (not HasOnlyCoordinateInterpolation(rS, rGroupMaster) or not HasOnlyCoordinateInterpolation(rS, rGroupSlave))
        throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "Elements must only have COORDINATES interpolation.");

    NuTo::Interpolation::eTypeOrder coordinateInterpolation = GetCoordinateInterpolation(rS, rGroupMaster, rGroupSlave);

    int it = rS.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::PRISM3D);
    rS.InterpolationTypeAdd(it, NuTo::Node::eDof::COORDINATES, coordinateInterpolation);
    rS.InterpolationTypeAdd(it, NuTo::Node::eDof::DISPLACEMENTS, NuTo::Interpolation::eTypeOrder::EQUIDISTANT2);

    auto pairs = FindMatchingElements(rS, rGroupMaster, rGroupSlave);
    std::cout << pairs.size() << std::endl;

    // UPDATE NODES:
    std::vector<NuTo::NodeBase*> gNodesPrism = GetPrismNodes(rS, rGroupMaster, rGroupSlave);

    std::map<NuTo::NodeBase*, NuTo::NodeBase*> clonedNodesMapping;
    std::map<NuTo::NodeBase*, Eigen::Vector3d> normals;


    for (NuTo::NodeBase* node : gNodesPrism)
    {
        auto elementPairsContainingTheNode = FindElementPairsContainingTheNode(node, pairs);
        normals[node] = CalculateNormalAtNode(node, elementPairsContainingTheNode[0].first);
    }

    for (NuTo::NodeBase* node : gNodesPrism)
    {
        NuTo::NodeBase* clonedNode = CloneNode(rS, node);
        Eigen::Vector3d nodeCoordinates = node->Get(NuTo::Node::eDof::COORDINATES);
        const Eigen::Vector3d& normal = normals[node];

//        std::cout << "normal\n" << normal << std::endl;

        node->Set(NuTo::Node::eDof::COORDINATES, nodeCoordinates - rThickness * 0.5 * normal);
        clonedNode->Set(NuTo::Node::eDof::COORDINATES, nodeCoordinates + rThickness * 0.5 * normal);

        clonedNodesMapping[node] = clonedNode;

        for (int slaveElementId : rS.GroupGetMemberIds(rGroupSlave))
            rS.ElementGetElementPtr(slaveElementId)->ExchangeNodePtr(node, clonedNode);
    }


    int gPrism = rS.GroupCreate(NuTo::eGroupId::Elements);

    std::map<NuTo::NodeBase*, NuTo::NodeBase*> masterNodeToInterpolatedNode;

    // CREATE PRISMS
    for (auto& pair : pairs)
    {
        if (coordinateInterpolation == NuTo::Interpolation::eTypeOrder::EQUIDISTANT1)
        {
            auto surfaceNodes = GetSurfaceNodes(pair.first);

            // linear elements
            assert(surfaceNodes.size() == 3);
            for (int i = 0; i < 3; ++i)
            {
                NuTo::NodeBase* masterNode = surfaceNodes[i];
                NuTo::NodeBase* slaveNode = clonedNodesMapping[masterNode];
                surfaceNodes.push_back(slaveNode);
            }

            int elementId = rS.ElementCreate(it, surfaceNodes);
            rS.GroupAddElement(gPrism, elementId);
        }
        else if (coordinateInterpolation == NuTo::Interpolation::eTypeOrder::EQUIDISTANT2)
        {
            const auto& itBase = pair.first.elem->GetInterpolationType().Get(NuTo::Node::eDof::COORDINATES);

            std::vector<NuTo::NodeBase*> sortedNodes(6);

            std::cout << itBase.GetNumSurfaceNodes(pair.first.surface) << std::endl;
            for (int i = 0; i < itBase.GetNumSurfaceNodes(pair.first.surface); ++i)
            {
                int nodeIndex = itBase.GetSurfaceNodeIndex(pair.first.surface, i);
                Eigen::Vector2d surfaceCoords = GetLocalSurfaceCoordinates(nodeIndex, pair.first);
//                std::cout << surfaceCoords << std::endl;

                if (surfaceCoords.isApprox(Eigen::Vector2d{0,0}))
                    sortedNodes[0] = pair.first.elem->GetNode(nodeIndex);
                if (surfaceCoords.isApprox(Eigen::Vector2d{1,0}))
                    sortedNodes[1] = pair.first.elem->GetNode(nodeIndex);
                if (surfaceCoords.isApprox(Eigen::Vector2d{0,1}))
                    sortedNodes[2] = pair.first.elem->GetNode(nodeIndex);
                if (surfaceCoords.isApprox(Eigen::Vector2d{.5,0}))
                    sortedNodes[3] = pair.first.elem->GetNode(nodeIndex);
                if (surfaceCoords.isApprox(Eigen::Vector2d{.5,.5}))
                    sortedNodes[4] = pair.first.elem->GetNode(nodeIndex);
                if (surfaceCoords.isApprox(Eigen::Vector2d{0,.5}))
                    sortedNodes[5] = pair.first.elem->GetNode(nodeIndex);
            }

            /*
                                     Prism18:

                                           w
                                           ^
                                           |
                                           3
                                         ,/|`\
                                       12  |  13
                                     ,/    |    `\
                                    4------14-----5
                                    |      8      |
                                    |    ,/|`\    |
                                    |  15  |  16  |
                                    |,/    |    `\|
                                   ,10-----17-----11
                                 ,/ |      0      | `\
                                u   |    ,/ `\    |   v
                                    |  ,6     `7  |
                                    |,/         `\|
                                    1------9------2

             */

            std::vector<NuTo::NodeBase*> nodeVector(18);
            // lower corners
            nodeVector[0] = sortedNodes[0];
            nodeVector[1] = sortedNodes[1];
            nodeVector[2] = sortedNodes[2];

            // upper corners
            nodeVector[3] = clonedNodesMapping[sortedNodes[0]];
            nodeVector[4] = clonedNodesMapping[sortedNodes[1]];
            nodeVector[5] = clonedNodesMapping[sortedNodes[2]];

            // lower mids
            nodeVector[6] = sortedNodes[3];
            nodeVector[9] = sortedNodes[4];
            nodeVector[7] = sortedNodes[5];

            // upper mids
            nodeVector[12] = clonedNodesMapping[sortedNodes[3]];
            nodeVector[14] = clonedNodesMapping[sortedNodes[4]];
            nodeVector[13] = clonedNodesMapping[sortedNodes[5]];

            // mid corners
            nodeVector[8] =  InterpolateNode(rS, nodeVector[0], nodeVector[3], masterNodeToInterpolatedNode);
            nodeVector[10] = InterpolateNode(rS, nodeVector[1], nodeVector[4], masterNodeToInterpolatedNode);
            nodeVector[11] = InterpolateNode(rS, nodeVector[2], nodeVector[5], masterNodeToInterpolatedNode);

            // mid mids
            nodeVector[15] = InterpolateNode(rS, nodeVector[6], nodeVector[12], masterNodeToInterpolatedNode);
            nodeVector[16] = InterpolateNode(rS, nodeVector[7], nodeVector[13], masterNodeToInterpolatedNode);
            nodeVector[17] = InterpolateNode(rS, nodeVector[9], nodeVector[14], masterNodeToInterpolatedNode);

//            int i = 0;
//            for (auto* n : nodeVector)
//            {
////                std::cout << i << ": \t" << n->Get(NuTo::Node::eDof::COORDINATES).transpose() << std::endl;
////                ++i;
//            }

            int elementId = rS.ElementCreate(it, nodeVector);
            rS.GroupAddElement(gPrism, elementId);
        }



        // find the corner nodes master (3)
        // find the nodes in between the corner nodes, maybe via natural surface coordinates (3)
        // sort them

        // find the matching slave nodes (6)
        // interpolate between master and slave (another 6) BOOM.

        // create
    }
    return gPrism;
}



void PrismCreate(NuTo::Interpolation::eTypeOrder rCoordinateInterpolation)
{
    constexpr double thickness = .1;
    constexpr double lx = 10;
    constexpr double ly = 2;
    constexpr double lz = 5;

    NuTo::Structure s(3);

    int it = s.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::TETRAHEDRON3D);
    s.InterpolationTypeAdd(it, NuTo::Node::eDof::COORDINATES, rCoordinateInterpolation);

    int gMatrix = s.GroupCreate(NuTo::eGroupId::Elements);
    int gAggreg = s.GroupCreate(NuTo::eGroupId::Elements);

    if (rCoordinateInterpolation == NuTo::Interpolation::eTypeOrder::EQUIDISTANT1)
    {
        s.NodeCreate(0, Eigen::Vector3d({-lx, 0, 0}));
        s.NodeCreate(1, Eigen::Vector3d({0,-ly, 0}));
        s.NodeCreate(2, Eigen::Vector3d({0, 0, 0}));
        s.NodeCreate(3, Eigen::Vector3d({0, ly, 0}));
        s.NodeCreate(4, Eigen::Vector3d({0, -ly/2., lz/2.}));
        s.NodeCreate(5, Eigen::Vector3d({0, ly/2., lz/2.}));
        s.NodeCreate(6, Eigen::Vector3d({0, 0, lz}));
        s.NodeCreate(7, Eigen::Vector3d({lx, 0, 0}));

        s.ElementCreate(1, it, {0, 1, 2, 4});
        s.ElementCreate(2, it, {0, 2, 3, 5});
        s.ElementCreate(3, it, {0, 4, 2, 5});
        s.ElementCreate(4, it, {0, 4, 5, 6});
        s.ElementCreate(5, it, {7, 1, 2, 4});
        s.ElementCreate(6, it, {7, 2, 3, 5});
        s.ElementCreate(7, it, {7, 4, 2, 5});
        s.ElementCreate(8, it, {7, 4, 5, 6});


        s.GroupAddElement(gMatrix, 1);
        s.GroupAddElement(gMatrix, 2);
        s.GroupAddElement(gMatrix, 3);
        s.GroupAddElement(gMatrix, 4);

        s.GroupAddElement(gAggreg, 5);
        s.GroupAddElement(gAggreg, 6);
        s.GroupAddElement(gAggreg, 7);
        s.GroupAddElement(gAggreg, 8);
        CreatePrisms(s, gMatrix, gAggreg, thickness);
    }
    else
    {
        Eigen::Vector3d(0.5, 0.0, 0.0);
        Eigen::Vector3d(0.5, 0.5, 0.0);
        Eigen::Vector3d(0.0, 0.5, 0.0);
        Eigen::Vector3d(0.0, 0.0, 0.5);
        Eigen::Vector3d(0.0, 0.5, 0.5);
        Eigen::Vector3d(0.5, 0.0, 0.5);



        s.NodeCreate(0, Eigen::Vector3d({0,   -ly/2, 0}));
        s.NodeCreate(1, Eigen::Vector3d({lx,      0, 0}));
        s.NodeCreate(2, Eigen::Vector3d({0,    ly/2, 0}));
        s.NodeCreate(3, Eigen::Vector3d({0,       0, lz}));

        s.NodeCreate(4, Eigen::Vector3d({lx/2.,-ly/4, 0}));
        s.NodeCreate(5, Eigen::Vector3d({lx/2., ly/4, 0}));
        s.NodeCreate(6, Eigen::Vector3d({0,        0, 0}));

        s.NodeCreate(7, Eigen::Vector3d({0,    -ly/4, lz/2}));
        s.NodeCreate(8, Eigen::Vector3d({0,     ly/4, lz/2}));
        s.NodeCreate(9, Eigen::Vector3d({lx/2.,   0, lz/2}));


        s.NodeCreate(10, Eigen::Vector3d({-lx,      0, 0}));
        s.NodeCreate(11, Eigen::Vector3d({-lx/2.,-ly/4, 0}));
        s.NodeCreate(12, Eigen::Vector3d({-lx/2., ly/4, 0}));
        s.NodeCreate(13, Eigen::Vector3d({-lx/2.,   0, lz/2}));



        s.ElementCreate(1, it, {0,  1, 2, 3,  4,  5, 6, 7, 8,  9});
        s.ElementCreate(2, it, {0, 10, 2, 3, 11, 12, 6, 7, 8, 13});
//            s.ElementCreate(2, it, {1, 2, 3, 7});

        s.GroupAddElement(gMatrix, 1);
        s.GroupAddElement(gAggreg, 2);

        CreatePrisms(s, gMatrix, gAggreg, thickness);
    }
    s.InterpolationTypeAdd(it, NuTo::Node::eDof::DISPLACEMENTS, NuTo::Interpolation::eTypeOrder::EQUIDISTANT2);
    s.ElementTotalConvertToInterpolationType();
}



void CSDA3D()
{

    NuTo::Structure s(3);
    auto ids = s.ImportFromGmsh("CSDAMesh.msh");

    int gMatrix = ids[0].first;
    int gAggreg = ids[1].first;


    using namespace NuTo::Constitutive;
    int LIN = 0;
    int CSDA = 1;
    s.ConstitutiveLawCreate(LIN, eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
    s.ConstitutiveLawSetParameterDouble(LIN, eConstitutiveParameter::YOUNGS_MODULUS,       20000.);
    s.ConstitutiveLawSetParameterDouble(LIN, eConstitutiveParameter::POISSONS_RATIO,       0.0);

    const double thickness = 0.1;

    s.ConstitutiveLawCreate(CSDA, eConstitutiveType::LOCAL_DAMAGE_MODEL);
    constexpr double fractureEnergy         = 0.1;
    s.ConstitutiveLawSetParameterDouble(CSDA, eConstitutiveParameter::YOUNGS_MODULUS,       200.);
    s.ConstitutiveLawSetParameterDouble(CSDA, eConstitutiveParameter::POISSONS_RATIO,       0.0);
    s.ConstitutiveLawSetParameterDouble(CSDA, eConstitutiveParameter::TENSILE_STRENGTH,     4.);
    s.ConstitutiveLawSetParameterDouble(CSDA, eConstitutiveParameter::COMPRESSIVE_STRENGTH, 40.);
    s.ConstitutiveLawSetParameterDouble(CSDA, eConstitutiveParameter::FRACTURE_ENERGY,      fractureEnergy / thickness);
    s.ConstitutiveLawSetDamageLaw(CSDA, eDamageLawType::ISOTROPIC_EXPONENTIAL_SOFTENING);

    std::cout << "GetNumNodes() \n" << s.GetNumNodes() << std::endl;

    int gPrisms = CreatePrisms(s, gMatrix, gAggreg, thickness);
    s.ElementTotalSetConstitutiveLaw(LIN);
    s.ElementGroupSetConstitutiveLaw(gPrisms, CSDA);

    std::cout << "GetNumNodes() \n" << s.GetNumNodes() << std::endl;


    s.InterpolationTypeAdd(ids[0].second, NuTo::Node::eDof::DISPLACEMENTS, NuTo::Interpolation::eTypeOrder::EQUIDISTANT2);
    s.InterpolationTypeAdd(ids[1].second, NuTo::Node::eDof::DISPLACEMENTS, NuTo::Interpolation::eTypeOrder::EQUIDISTANT2);

    int mySection = s.SectionCreate("Volume");
    s.ElementTotalSetSection(mySection);
    s.ElementTotalConvertToInterpolationType();

    std::cout << "GetNumNodes() \n" << s.GetNumNodes() << std::endl;

    s.ElementInfo(10);
    s.NodeInfo(10);

    int nodeFixXYZ = s.NodeGetIdAtCoordinate(Eigen::Vector3d({-5, 0, 0}), 1.e-5);
    s.ConstraintLinearSetDisplacementNode(nodeFixXYZ, Eigen::Vector3d::UnitX(), 0);
    s.ConstraintLinearSetDisplacementNode(nodeFixXYZ, Eigen::Vector3d::UnitY(), 0);
    s.ConstraintLinearSetDisplacementNode(nodeFixXYZ, Eigen::Vector3d::UnitZ(), 0);

    int nodeFixYZ = s.NodeGetIdAtCoordinate(Eigen::Vector3d({5, 0, 0}), 1.e-5);
    s.ConstraintLinearSetDisplacementNode(nodeFixYZ, Eigen::Vector3d::UnitY(), 0);
    s.ConstraintLinearSetDisplacementNode(nodeFixYZ, Eigen::Vector3d::UnitZ(), 0);

    int groupNodeFixZ = s.GroupCreate(NuTo::eGroupId::Nodes);
    s.GroupAddNodeRadiusRange(groupNodeFixZ, Eigen::Vector3d({0, 0, 2}), 0, 2*thickness);

    s.ConstraintLinearSetDisplacementNodeGroup(groupNodeFixZ, Eigen::Vector3d::UnitY(), 0);
//    s.ConstraintLinearSetDisplacementNodeGroup(groupNodeFixZ, Eigen::Vector3d::UnitX(), 0);

//    int nodeL = s.NodeGetIdAtCoordinate(Eigen::Vector3d({- thickness/2., 0, lz}), 1.e-6);
//    int nodeR = s.NodeGetIdAtCoordinate(Eigen::Vector3d({+ thickness/2., 0, lz}), 1.e-6);
//    int constraintId = s.ConstraintLinearEquationCreate(nodeL, "X_DISPLACEMENT", 1, 0);
//    s.ConstraintLinearEquationAddTerm(constraintId, nodeR, "X_DISPLACEMENT", -1);


    int BC = s.ConstraintLinearSetDisplacementNodeGroup(groupNodeFixZ, Eigen::Vector3d::UnitZ(), 0);

    s.NodeBuildGlobalDofs();
    std::cout << s.GetNumTotalActiveDofs() << std::endl;
    std::cout << s.GetNumTotalDependentDofs() << std::endl;

    double deltaD = .5;

    Eigen::Matrix2d dispRHS;
    dispRHS << 0, 0, 1, -deltaD;

    s.AddVisualizationComponent(s.GroupGetElementsTotal(), NuTo::eVisualizeWhat::DISPLACEMENTS);

    NuTo::NewmarkDirect newmark(&s);

    s.SetShowTime(false);
    newmark.SetShowTime(false);

    newmark.SetTimeStep(0.1);
    newmark.SetMinTimeStep(1.e-12);
    newmark.SetMaxTimeStep(0.1);
    newmark.SetToleranceForce(1.e-06);
    newmark.SetAutomaticTimeStepping(true);
    newmark.SetPerformLineSearch(true);
    newmark.SetMaxNumIterations(100);

    bool deleteDirectory = true;
    newmark.SetResultDirectory("./CSDA3D", deleteDirectory);

    newmark.AddTimeDependentConstraint(BC, dispRHS);

//    newmark.AddResultNodeDisplacements("Displ", nodeBC);
//    int groupNodeBC = s.GroupCreate(NuTo::eGroupId::Nodes);
//    s.GroupAddNode(groupNodeBC, nodeBC);
//    newmark.AddResultGroupNodeForce("Force", groupNodeBC);


    newmark.Solve(1);
}

int main()
{
    PrismCreate(NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    PrismCreate(NuTo::Interpolation::eTypeOrder::EQUIDISTANT2);


//    CSDA2D();
    CSDA3D();

//
//    CheckFractureEnergy2D(90, .1);
//    CheckFractureEnergy2D(90, .01);
//    CheckFractureEnergy2D(90, .001);
//
//    CheckFractureEnergy2D(75, .001);

    return EXIT_SUCCESS;
}
