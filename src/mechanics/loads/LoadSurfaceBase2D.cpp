#include <set>

#include "mechanics/groups/Group.h"
#include "mechanics/nodes/NodeBase.h"
#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/sections/Section.h"
#include "mechanics/structures/StructureBase.h"
#include "mechanics/integrationtypes/IntegrationTypeBase.h"
#include "mechanics/loads/LoadSurfaceBase2D.h"
#include "mechanics/loads/LoadSurfacePressureFunction2D.h"
#include "mechanics/elements/ContinuumElement.h"
#include "mechanics/integrationtypes/IntegrationTypeEnum.h"
#include "mechanics/interpolationtypes/InterpolationBase.h"
#include "mechanics/interpolationtypes/InterpolationType.h"
#include "mechanics/interpolationtypes/InterpolationTypeEnum.h"

using namespace NuTo;

LoadSurfaceBase2D::LoadSurfaceBase2D(StructureBase* rStructure, int rElementGroupId, int rNodeGroupId)
    : LoadBase()
{
    // get element group
    const Group<ElementBase>* elementGroup = rStructure->GroupGetGroupPtr(rElementGroupId)->AsGroupElement();

    // get node group
    const Group<NodeBase>* nodeGroup = rStructure->GroupGetGroupPtr(rNodeGroupId)->AsGroupNode();

    // since the search is done via the id's, the surface nodes are ptr, so make another set with the node ptrs
    std::set<const NodeBase*> nodePtrSet;
    for (auto node : *nodeGroup)
    {
        nodePtrSet.insert(node.second);
    }

    // loop over all elements
    Eigen::VectorXi surfaceNodeIndices;
    std::vector<const NodeBase*> surfaceNodes;

    for (auto itElement : *elementGroup)
    {
        // check if plane element
        auto& elementPtr = *dynamic_cast<ContinuumElement<2>*>(itElement.second);
        const InterpolationType& interpolationType = elementPtr.GetInterpolationType();

        // loop over all surfaces
        for (int iSurface = 0; iSurface < interpolationType.GetNumSurfaces(); iSurface++)
        {
            bool addSurface = true;
            surfaceNodeIndices = interpolationType.GetSurfaceNodeIndices(iSurface);
            int numSurfaceNodes = surfaceNodeIndices.rows();
            surfaceNodes.resize(numSurfaceNodes);

            for (int iSurfaceNode = 0; iSurfaceNode < numSurfaceNodes; ++iSurfaceNode)
            {
                surfaceNodes[iSurfaceNode] = elementPtr.GetNode(surfaceNodeIndices(iSurfaceNode, 0));
            }

            // check, if all surface nodes are in the node group
            for (auto& surfaceNode : surfaceNodes)
            {
                if (nodePtrSet.find(surfaceNode) == nodePtrSet.end())
                {
                    // this surface has at least one node that is not in the list, continue
                    addSurface = false;
                }
            }

            if (addSurface)
            {
                mElements2D.push_back(std::make_pair(&elementPtr, iSurface));
            }
        }
    }

    // set standard integration types for triangles and quads, this can be modified according to the needs
    mIntegrationType2NPtr = rStructure->GetPtrIntegrationType(eIntegrationType::IntegrationType1D2NGauss1Ip);
    mIntegrationType3NPtr = rStructure->GetPtrIntegrationType(eIntegrationType::IntegrationType1D2NGauss2Ip);
    mIntegrationType4NPtr = rStructure->GetPtrIntegrationType(eIntegrationType::IntegrationType1D2NGauss3Ip);
    mIntegrationType5NPtr = rStructure->GetPtrIntegrationType(eIntegrationType::IntegrationType1D2NGauss4Ip);

    mIntegrationType3NPtrLobatto = rStructure->GetPtrIntegrationType(eIntegrationType::IntegrationType1D2NLobatto3Ip);
    mIntegrationType4NPtrLobatto = rStructure->GetPtrIntegrationType(eIntegrationType::IntegrationType1D2NLobatto4Ip);
    mIntegrationType5NPtrLobatto = rStructure->GetPtrIntegrationType(eIntegrationType::IntegrationType1D2NLobatto5Ip);
}


void LoadSurfaceBase2D::AddLoadToGlobalSubVectors(StructureOutputBlockVector& externalLoad) const
{
    for (auto it : mElements2D)
    {
        const auto* elementPtr = it.first;
        int surface = it.second;

        const InterpolationBase& interpolationTypeDisps =
                elementPtr->GetInterpolationType().Get(Node::eDof::DISPLACEMENTS);
        const InterpolationBase& interpolationTypeCoords =
                elementPtr->GetInterpolationType().Get(Node::eDof::COORDINATES);

        IntegrationTypeBase* integrationType(nullptr);
        switch (interpolationTypeDisps.GetTypeOrder())
        {
        case Interpolation::eTypeOrder::EQUIDISTANT1:
            integrationType = mIntegrationType2NPtr;
            break;
        case Interpolation::eTypeOrder::EQUIDISTANT2:
            integrationType = mIntegrationType3NPtr;
            break;
        case Interpolation::eTypeOrder::EQUIDISTANT3:
            integrationType = mIntegrationType4NPtr;
            break;
        case Interpolation::eTypeOrder::EQUIDISTANT4:
            integrationType = mIntegrationType5NPtr;
            break;
        case Interpolation::eTypeOrder::LOBATTO2:
            integrationType = mIntegrationType3NPtrLobatto;
            break;
        case Interpolation::eTypeOrder::LOBATTO3:
            integrationType = mIntegrationType4NPtrLobatto;
            break;
        case Interpolation::eTypeOrder::LOBATTO4:
            integrationType = mIntegrationType5NPtrLobatto;
            break;
        case Interpolation::eTypeOrder::SPLINE:
        {
            int degree = interpolationTypeDisps.GetSurfaceDegree(surface);
            switch (degree)
            {
            case 0: // (0+1)/2 = 0.5 ips
                integrationType = mIntegrationType2NPtr;
                break;
            case 1: // 2/2 = 1 ips or (1+3)/2 = 2 ips lobatto
                integrationType = mIntegrationType2NPtr;
                break;
            case 2: // (2+1)/2 = 1.5 ips or (2+3)/2 = 2.5 ips lobatto
                integrationType = mIntegrationType3NPtr;
                break;
            case 3: // (3+1)/2 = 2 ips or (3+3)/2 = 3 ips lobatto
                integrationType = mIntegrationType3NPtr;
                break;
            case 4: // (4+1)/2 = 2.5 ips or (4+3)/2 = 3.5 ips lobatto
                integrationType = mIntegrationType4NPtr;
                break;
            default:
                throw MechanicsException(__PRETTY_FUNCTION__, "Interpolation for exact integration of " +
                                                                      std::to_string(degree) + " IGA not implemented");
            }
            break;
        }
        default:
            throw MechanicsException(__PRETTY_FUNCTION__,
                                     "Integration types only for 2, 3, 4 and 5 nodes (on the surface) implemented.");
        }

        Eigen::VectorXd nodeCoordinates = elementPtr->ExtractNodeValues(0, Node::eDof::COORDINATES);

        Eigen::Matrix<double, 1, 1> ipCoordsSurface;
        Eigen::Vector2d ipCoordsNatural;
        Eigen::Vector2d ipCoordsGlobal;

        Eigen::MatrixXd derivativeNaturalSurfaceCoordinates;
        Eigen::MatrixXd derivativeShapeFunctionsNatural;
        Eigen::VectorXd shapeFunctions;

        // loop over surface integration points
        for (int theIp = 0; theIp < integrationType->GetNumIntegrationPoints(); theIp++)
        {
            // #######################################
            // ##  Calculate IP coordinates in
            // ## surface CS, natural CS and global CS
            // #######################################
            ipCoordsSurface = integrationType->GetLocalIntegrationPointCoordinates(theIp);


            if (interpolationTypeDisps.GetTypeOrder() == Interpolation::eTypeOrder::SPLINE)
            {
                ipCoordsNatural = interpolationTypeCoords.CalculateNaturalSurfaceCoordinatesIGA(
                        ipCoordsSurface, surface, elementPtr->GetKnots());
                derivativeShapeFunctionsNatural = interpolationTypeCoords.DerivativeShapeFunctionsNaturalIGA(
                        ipCoordsNatural, elementPtr->GetKnotIDs());
                shapeFunctions = interpolationTypeDisps.ShapeFunctionsIGA(ipCoordsNatural, elementPtr->GetKnotIDs());
            }
            else
            {
                ipCoordsNatural = interpolationTypeCoords.CalculateNaturalSurfaceCoordinates(ipCoordsSurface, surface);
                derivativeShapeFunctionsNatural =
                        interpolationTypeCoords.DerivativeShapeFunctionsNatural(ipCoordsNatural);
                shapeFunctions = interpolationTypeDisps.ShapeFunctions(ipCoordsNatural);
            }

            // Eigen::MatrixXd N = interpolationTypeCoords.MatrixN(ipCoordsNatural);
            ipCoordsGlobal = elementPtr->InterpolateDofGlobal(ipCoordsNatural, NuTo::Node::eDof::COORDINATES);
            // #######################################
            // ##  Calculate the surface jacobian
            // ## = || [dX / dXi] * [dXi / dAlpha] ||
            // #######################################


            const Eigen::Matrix2d jacobianStd =
                    elementPtr->CalculateJacobian(derivativeShapeFunctionsNatural, nodeCoordinates); // = [dX / dXi]

            // in case of non IGA - just an identity matrix
            const Eigen::Matrix2d jacobianIGA =
                    elementPtr->CalculateJacobianParametricSpaceIGA(); // = [dXi / d\tilde{Xi}]


            derivativeNaturalSurfaceCoordinates = interpolationTypeCoords.CalculateDerivativeNaturalSurfaceCoordinates(
                    ipCoordsSurface, surface); // = [dXi / dAlpha]

            double detJacobian = (jacobianStd * jacobianIGA * derivativeNaturalSurfaceCoordinates)
                                         .norm(); // = || [dX / dXi] * [dXi / dAlpha] ||

            // #######################################
            // ##  Calculate surface normal vector
            // ## = ( dY / dAlpha     -dX/dAlpha).T
            // #######################################
            // dXdAlpha :  [2 x NumNodes*2] * [NumNodes*2 x 1] * [2 x 1] = [2 x 1]
            Eigen::Vector2d surfaceTangentVector = jacobianStd * derivativeNaturalSurfaceCoordinates;
            surfaceTangentVector.normalize();
            Eigen::Vector2d surfaceNormalVector;
            surfaceNormalVector(0) = surfaceTangentVector(1, 0);
            surfaceNormalVector(1) = -surfaceTangentVector(0, 0);


            // calculate weighting factor
            double thickness = elementPtr->GetSection()->GetThickness();
            double factor = thickness * (integrationType->GetIntegrationPointWeight(theIp)) * detJacobian;

            // calculate surface load
            Eigen::Vector2d loadVector;
            CalculateSurfaceLoad(ipCoordsGlobal, surfaceNormalVector, loadVector);
            loadVector *= factor;

            // add load vector to global vector
            for (int iNode = 0; iNode < shapeFunctions.rows(); iNode++)
            {
                const NodeBase* node = elementPtr->GetNode(interpolationTypeDisps.GetNodeIndex(iNode));
                assert(node->GetNum(Node::eDof::DISPLACEMENTS) == 2);
                for (int iDispDof = 0; iDispDof < 2; iDispDof++)
                {
                    int theDof = node->GetDof(Node::eDof::DISPLACEMENTS, iDispDof);
                    double theLoad = shapeFunctions[iNode] * loadVector(iDispDof);
                    if (theDof < externalLoad.J[Node::eDof::DISPLACEMENTS].rows())
                    {
                        externalLoad.J[Node::eDof::DISPLACEMENTS](theDof) += theLoad;
                    }
                    else
                    {
                        externalLoad.K[Node::eDof::DISPLACEMENTS](
                                theDof - externalLoad.J[Node::eDof::DISPLACEMENTS].rows()) += theLoad;
                    }
                }
            }
        }
    }
}
