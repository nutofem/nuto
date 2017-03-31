#include <set>

#include "mechanics/groups/Group.h"
#include "mechanics/nodes/NodeBase.h"
#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/structures/StructureBase.h"
#include "mechanics/integrationtypes/IntegrationTypeEnum.h"
#include "mechanics/integrationtypes/IntegrationTypeBase.h"
#include "mechanics/interpolationtypes/InterpolationBase.h"
#include "mechanics/interpolationtypes/InterpolationType.h"
#include "mechanics/interpolationtypes/InterpolationTypeEnum.h"
#include "mechanics/loads/LoadSurfaceBase3D.h"
#include "mechanics/elements/ContinuumElement.h"

using namespace NuTo;

LoadSurfaceBase3D::LoadSurfaceBase3D(StructureBase* rStructure, int rElementGroupId, int rNodeGroupId)
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
        // check if solid element
        auto& elementPtr = *dynamic_cast<ContinuumElement<3>*>(itElement.second);
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
                    // this surface has at least on node that is not in the list, continue
                    addSurface = false;
                }
            }

            if (addSurface)
            {
                mVolumeElements.push_back(std::make_pair(&elementPtr, iSurface));
            }
        }
    }

    // set standard integration types for triangles and quads, this can be modified according to the needs
    mIntegrationTypeTriangleGauss1 = rStructure->GetPtrIntegrationType(eIntegrationType::IntegrationType2D3NGauss1Ip);
    mIntegrationTypeTriangleGauss2 = rStructure->GetPtrIntegrationType(eIntegrationType::IntegrationType2D3NGauss3Ip);

    mIntegrationTypeQuadGauss1 = rStructure->GetPtrIntegrationType(eIntegrationType::IntegrationType2D4NGauss4Ip);
    mIntegrationTypeQuadGauss2 = rStructure->GetPtrIntegrationType(eIntegrationType::IntegrationType2D4NGauss4Ip);

    mIntegrationTypeQuadLobatto2 = rStructure->GetPtrIntegrationType(eIntegrationType::IntegrationType2D4NLobatto9Ip);
    mIntegrationTypeQuadLobatto3 = rStructure->GetPtrIntegrationType(eIntegrationType::IntegrationType2D4NLobatto16Ip);
    mIntegrationTypeQuadLobatto4 = rStructure->GetPtrIntegrationType(eIntegrationType::IntegrationType2D4NLobatto25Ip);
}


void LoadSurfaceBase3D::AddLoadToGlobalSubVectors(Eigen::VectorXd& rActiceDofsLoadVector,
                                                  Eigen::VectorXd& rDependentDofsLoadVector) const
{
    for (auto it : mVolumeElements)
    {

        const auto* elementPtr = it.first;
        int surface = it.second;

        const InterpolationBase& interpolationTypeDisps =
                elementPtr->GetInterpolationType().Get(Node::eDof::DISPLACEMENTS);
        const InterpolationBase& interpolationTypeCoords =
                elementPtr->GetInterpolationType().Get(Node::eDof::COORDINATES);

        IntegrationTypeBase* integrationType(nullptr);
        switch (elementPtr->GetInterpolationType().GetShapeType())
        {
        case Interpolation::eShapeType::TETRAHEDRON3D:
        {
            switch (interpolationTypeDisps.GetTypeOrder())
            {
            case Interpolation::eTypeOrder::EQUIDISTANT1:
                integrationType = mIntegrationTypeTriangleGauss1;
                break;
            case Interpolation::eTypeOrder::EQUIDISTANT2:
                integrationType = mIntegrationTypeTriangleGauss2;
                break;
            default:
                break;
            }
        }
        break;
        case Interpolation::eShapeType::BRICK3D:
        {
            switch (interpolationTypeDisps.GetTypeOrder())
            {
            case Interpolation::eTypeOrder::EQUIDISTANT1:
                integrationType = mIntegrationTypeQuadGauss1;
                break;

            case Interpolation::eTypeOrder::EQUIDISTANT2:
                integrationType = mIntegrationTypeQuadGauss2;
                break;

            case Interpolation::eTypeOrder::LOBATTO2:
                integrationType = mIntegrationTypeQuadLobatto2;
                break;

            case Interpolation::eTypeOrder::LOBATTO3:
                integrationType = mIntegrationTypeQuadLobatto3;
                break;

            case Interpolation::eTypeOrder::LOBATTO4:
                integrationType = mIntegrationTypeQuadLobatto4;
                break;

            default:
                break;
            }
        }
        break;
        default:
            break;
        }
        if (integrationType == nullptr)
            throw MechanicsException(__PRETTY_FUNCTION__,
                                     "Integration types only for 2, 3, 4 and 5 nodes (on the surface) implemented.");

        Eigen::MatrixXd nodeCoordinates = elementPtr->ExtractNodeValues(0, Node::eDof::COORDINATES);

        Eigen::Vector2d ipCoordsSurface;
        Eigen::Vector3d ipCoordsNatural;
        Eigen::Vector3d ipCoordsGlobal;

        Eigen::MatrixXd derivativeNaturalSurfaceCoordinates;
        Eigen::Vector3d dXdAlpha, dXdBeta;
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
            ipCoordsNatural = interpolationTypeCoords.CalculateNaturalSurfaceCoordinates(ipCoordsSurface, surface);
            ipCoordsGlobal = interpolationTypeCoords.CalculateMatrixN(ipCoordsNatural) * nodeCoordinates;

            // #######################################
            // ##  Calculate the surface jacobian
            // ##     normal vector
            // ## = [dX / dAlpha] x [dX / dBeta]
            // #######################################
            derivativeShapeFunctionsNatural =
                    interpolationTypeCoords.CalculateDerivativeShapeFunctionsNatural(ipCoordsNatural);
            const Eigen::Matrix3d jacobian =
                    elementPtr->CalculateJacobian(derivativeShapeFunctionsNatural, nodeCoordinates); // = [dX / dXi]

            derivativeNaturalSurfaceCoordinates = interpolationTypeCoords.CalculateDerivativeNaturalSurfaceCoordinates(
                    ipCoordsSurface, surface); // = [dXi / dAlpha]
            dXdAlpha = jacobian * derivativeNaturalSurfaceCoordinates.col(0);
            dXdBeta = jacobian * derivativeNaturalSurfaceCoordinates.col(1);

            Eigen::Vector3d surfaceNormalVector = dXdAlpha.cross(dXdBeta); // = || [dX / dXi] * [dXi / dAlpha] ||

            double detJacobian = surfaceNormalVector.norm();

            surfaceNormalVector.normalize();


            // calculate weighting factor
            double factor(detJacobian * (integrationType->GetIntegrationPointWeight(theIp)));

            // calculate surface load
            Eigen::Vector3d loadVector;
            CalculateSurfaceLoad(ipCoordsGlobal, surfaceNormalVector, loadVector);
            loadVector *= factor;

            // calculate 3D shape functions
            shapeFunctions = interpolationTypeDisps.CalculateShapeFunctions(ipCoordsNatural);

            // add load vector to global vector
            for (int iNode = 0; iNode < shapeFunctions.rows(); iNode++)
            {
                const NodeBase* node = elementPtr->GetNode(interpolationTypeDisps.GetNodeIndex(iNode));
                assert(node->GetNum(Node::eDof::DISPLACEMENTS) == 3);
                for (int iDispDof = 0; iDispDof < 3; iDispDof++)
                {
                    int theDof = node->GetDof(Node::eDof::DISPLACEMENTS, iDispDof);
                    double theLoad = shapeFunctions[iNode] * loadVector(iDispDof);
                    if (theDof < rActiceDofsLoadVector.rows())
                    {
                        rActiceDofsLoadVector(theDof) += theLoad;
                    }
                    else
                    {
                        rDependentDofsLoadVector(theDof - rActiceDofsLoadVector.rows()) += theLoad;
                    }
                }
            }
        }
    }
}


#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_IMPLEMENT(LoadSurfaceBase3D)
#endif
