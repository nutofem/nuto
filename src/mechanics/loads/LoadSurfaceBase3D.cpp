// $Id: LoadLoadSurfaceBase3D.cpp 178 2009-12-11 20:53:12Z eckardt4 $

#include <set>
#include "math/FullVector.h"
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

//! @brief constructor
NuTo::LoadSurfaceBase3D::LoadSurfaceBase3D(int rLoadCase, StructureBase* rStructure, int rElementGroupId, int rNodeGroupId) :
        LoadBase(rLoadCase)
{
    //get element group
    const Group<ElementBase> *elementGroup = rStructure->GroupGetGroupPtr(rElementGroupId)->AsGroupElement();

    //get node group
    const Group<NodeBase> *nodeGroup = rStructure->GroupGetGroupPtr(rNodeGroupId)->AsGroupNode();

    //since the search is done via the id's, the surface nodes are ptr, so make another set with the node ptrs
    std::set<const NodeBase*> nodePtrSet;
    for (Group<NodeBase>::const_iterator itNode = nodeGroup->begin(); itNode != nodeGroup->end(); itNode++)
    {
        nodePtrSet.insert(itNode->second);
    }

    //loop over all elements
    Eigen::VectorXi surfaceNodeIndices;
    std::vector<const NodeBase*> surfaceNodes;
    for (auto itElement : *elementGroup)
    {
        try
        {
            //check if solid element
            ContinuumElement<3>& elementPtr = itElement.second->AsContinuumElement3D();
            const InterpolationType& interpolationType = elementPtr.GetInterpolationType();

            //loop over all surfaces
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

//                for (const NodeBase* node : nodePtrSet)
//                {
//                    std::cout << elementPtr->GetStructure()->NodeGetId(node)  << " ";
//                }
//                std::cout << std::endl;
//                std::cout << surfaceNodeIndices.transpose() << std::endl << std::endl;

                //check, if all surface nodes are in the node group
                for (unsigned int countNode = 0; countNode < surfaceNodes.size(); countNode++)
                {
                    if (nodePtrSet.find(surfaceNodes[countNode]) == nodePtrSet.end())
                    {
                        //this surface has at least on node that is not in the list, continue
                        addSurface = false;
                    }
                }

                if (addSurface)
                {
//                    std::cout << "Surface added!" << std::endl;
                    mVolumeElements.push_back(std::make_pair(&elementPtr, iSurface));
                }
            }
        } catch (NuTo::MechanicsException &e)
        {
            std::stringstream ss;
            assert(rStructure->ElementGetId(itElement.second) == itElement.first);
            ss << itElement.first;
            e.AddMessage("[NuTo::LoadSurfaceBase3D::LoadSurfaceBase3D] Error calculating surfaces for surface loads in element " + ss.str() + "(Maybe not a solid element?).");
            throw;
        } catch (...)
        {
            std::stringstream ss;
            assert(rStructure->ElementGetId(itElement.second) == itElement.first);
            ss << itElement.first;
            throw NuTo::MechanicsException("[NuTo::LoadSurfaceBase3D::LoadSurfaceBase3D] Error calculating surfaces for surface loads in element " + ss.str() + "(Maybe not a solid element?).");
        }

    }

    //set standard integration types for triangles and quads, this can be modified according to the needs
    mIntegrationTypeTriangleGauss1 = rStructure->GetPtrIntegrationType(eIntegrationType::IntegrationType2D3NGauss1Ip);
    mIntegrationTypeTriangleGauss2 = rStructure->GetPtrIntegrationType(eIntegrationType::IntegrationType2D3NGauss3Ip);


    mIntegrationTypeQuadGauss1 = rStructure->GetPtrIntegrationType(eIntegrationType::IntegrationType2D4NGauss4Ip);
    mIntegrationTypeQuadGauss2 = rStructure->GetPtrIntegrationType(eIntegrationType::IntegrationType2D4NGauss4Ip);
    mIntegrationTypeQuadLobatto2 = rStructure->GetPtrIntegrationType(eIntegrationType::IntegrationType2D4NLobatto9Ip);
    mIntegrationTypeQuadLobatto3 = rStructure->GetPtrIntegrationType(eIntegrationType::IntegrationType2D4NLobatto16Ip);
    mIntegrationTypeQuadLobatto4 = rStructure->GetPtrIntegrationType(eIntegrationType::IntegrationType2D4NLobatto25Ip);

}

//! @brief adds the load to global sub-vectors
//! @param rActiceDofsLoadVector ... global load vector which correspond to the active dofs
//! @param rDependentDofsLoadVector ... global load vector which correspond to the dependent dofs
void NuTo::LoadSurfaceBase3D::AddLoadToGlobalSubVectors(int rLoadCase, NuTo::FullVector<double, Eigen::Dynamic>& rActiceDofsLoadVector, NuTo::FullVector<double, Eigen::Dynamic>& rDependentDofsLoadVector) const
{
    if (rLoadCase != mLoadCase)
        return;
    for (auto it : mVolumeElements)
    {

        const auto* elementPtr = it.first;
        int surface = it.second;

//        std::cout << "Surface: " << surface << std::endl;

        const InterpolationBase& interpolationTypeDisps = elementPtr->GetInterpolationType().Get(Node::eDof::DISPLACEMENTS);
        const InterpolationBase& interpolationTypeCoords = elementPtr->GetInterpolationType().Get(Node::eDof::COORDINATES);

        IntegrationTypeBase* integrationType(0);
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
            throw MechanicsException("[NuTo::LoadSurfaceBase2D::LoadSurfaceBase2D] integration types only for 2, 3, 4 and 5 nodes (on the surface) implemented.");

        Eigen::MatrixXd nodeCoordinates = elementPtr->ExtractNodeValues(0, Node::eDof::COORDINATES);

        Eigen::Matrix<double, 2, 1> ipCoordsSurface;
        Eigen::Matrix<double, 3, 1> ipCoordsNatural;
        NuTo::FullVector<double, 3> ipCoordsGlobal;

        Eigen::MatrixXd derivativeNaturalSurfaceCoordinates;
        Eigen::Vector3d dXdAlpha, dXdBeta;
        Eigen::MatrixXd derivativeShapeFunctionsNatural;
        Eigen::VectorXd shapeFunctions;

        //loop over surface integration points
        for (int theIp = 0; theIp < integrationType->GetNumIntegrationPoints(); theIp++)
        {
            // #######################################
            // ##  Calculate IP coordinates in
            // ## surface CS, natural CS and global CS
            // #######################################
            double tmp[2];
            integrationType->GetLocalIntegrationPointCoordinates2D(theIp, tmp);
            ipCoordsSurface(0) = tmp[0];
            ipCoordsSurface(1) = tmp[1];
            ipCoordsNatural = interpolationTypeCoords.CalculateNaturalSurfaceCoordinates(ipCoordsSurface, surface);
            ipCoordsGlobal = interpolationTypeCoords.CalculateMatrixN(ipCoordsNatural) * nodeCoordinates;

            // #######################################
            // ##  Calculate the surface jacobian
            // ##     normal vector
            // ## = [dX / dAlpha] x [dX / dBeta]
            // #######################################
            derivativeShapeFunctionsNatural = interpolationTypeCoords.CalculateDerivativeShapeFunctionsNatural(ipCoordsNatural);
            const Eigen::Matrix3d jacobian = elementPtr->CalculateJacobian(derivativeShapeFunctionsNatural, nodeCoordinates);                     // = [dX / dXi]

            derivativeNaturalSurfaceCoordinates = interpolationTypeCoords.CalculateDerivativeNaturalSurfaceCoordinates(ipCoordsSurface, surface); // = [dXi / dAlpha]
            dXdAlpha = jacobian * derivativeNaturalSurfaceCoordinates.col(0);
            dXdBeta  = jacobian * derivativeNaturalSurfaceCoordinates.col(1);

            NuTo::FullVector<double, 3> surfaceNormalVector = dXdAlpha.cross(dXdBeta); // = || [dX / dXi] * [dXi / dAlpha] ||

            double detJacobian = surfaceNormalVector.Norm();

            surfaceNormalVector.normalize();


            //calculate weighting factor
            double factor(detJacobian * (integrationType->GetIntegrationPointWeight(theIp)));

//            std::cout << "DetJacobian: " << detJacobian << std::endl;

            //calculate surface load
            FullVector<double, 3> loadVector;
            CalculateSurfaceLoad(ipCoordsGlobal, surfaceNormalVector, loadVector);
            loadVector *= factor;

            // calculate 3D shape functions
            shapeFunctions = interpolationTypeDisps.CalculateShapeFunctions(ipCoordsNatural);
//            std::cout << shapeFunctions.transpose() << std::endl;

            //add load vector to global vector
            for (int iNode = 0; iNode < shapeFunctions.rows(); iNode++)
            {
                const NodeBase* node = elementPtr->GetNode(interpolationTypeDisps.GetNodeIndex(iNode));
                assert(node->GetNum(Node::eDof::DISPLACEMENTS) == 3);
                for (int iDispDof = 0; iDispDof < 3; iDispDof++)
                {
                    int theDof = node->GetDof(Node::eDof::DISPLACEMENTS, iDispDof);
                    double theLoad = shapeFunctions[iNode] * loadVector(iDispDof);
                    if (theDof < rActiceDofsLoadVector.GetNumRows())
                    {
                        //std::cout << "add to dof " << theDof << " " << shapeFunctions[iNode]*loadVector(countDispDof) << std::endl;
                        rActiceDofsLoadVector(theDof) += theLoad;
                    }
                    else
                    {
                        rDependentDofsLoadVector(theDof - rActiceDofsLoadVector.GetNumRows()) += theLoad;
                    }
                }
            }
        }
    }
}


#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::LoadSurfaceBase3D)
#endif
