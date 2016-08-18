// $Id: LoadLoadSurfaceBase2D.cpp 178 2009-12-11 20:53:12Z eckardt4 $
#include <set>
#include "nuto/mechanics/groups/Group.h"
#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/structures/StructureBase.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeEnum.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeBase.h"
#include "nuto/mechanics/loads/LoadSurfaceBase2D.h"
#include "nuto/mechanics/elements/ContinuumElement.h"

//! @brief constructor
NuTo::LoadSurfaceBase2D::LoadSurfaceBase2D(int rLoadCase, StructureBase* rStructure, int rElementGroupId, int rNodeGroupId) :
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

//    std::cout << "number of loaded nodes " << nodeGroup->GetNumMembers() << std::endl;

//loop over all elements
    Eigen::VectorXi surfaceNodeIndices;
    std::vector<const NodeBase*> surfaceNodes;

    for (auto itElement : *elementGroup)
    {
        try
        {
            //check if plane element
            ContinuumElement<2>& elementPtr = itElement.second->AsContinuumElement2D();
            const InterpolationType* InterpolationType = elementPtr.GetInterpolationType();

            //loop over all surfaces
            for (int iSurface = 0; iSurface < InterpolationType->GetNumSurfaces(); iSurface++)
            {
                bool addSurface = true;
                surfaceNodeIndices = InterpolationType->GetSurfaceNodeIndices(iSurface);
                int numSurfaceNodes = surfaceNodeIndices.rows();
                surfaceNodes.resize(numSurfaceNodes);

                for (int iSurfaceNode = 0; iSurfaceNode < numSurfaceNodes; ++iSurfaceNode)
                {
                    surfaceNodes[iSurfaceNode] = elementPtr.GetNode(surfaceNodeIndices.at(iSurfaceNode, 0));
                }

                //check, if all surface nodes are in the node group
                for (unsigned int countNode = 0; countNode < surfaceNodes.size(); countNode++)
                {
                    if (nodePtrSet.find(surfaceNodes[countNode]) == nodePtrSet.end())
                    {
                        //this surface has at least one node that is not in the list, continue
                        addSurface = false;
                    }
                }

                if (addSurface)
                {
                    mElements2D.push_back(std::make_pair(&elementPtr, iSurface));
//            		double nodeCoordinates[2];
//            		surfaceNodes[0]->GetCoordinates2D(nodeCoordinates);
//            		surfaceNodes[1]->GetCoordinates2D(nodeCoordinates);
                }
            }
        } catch (NuTo::MechanicsException &e)
        {
            std::stringstream ss;
            assert(rStructure->ElementGetId(itElement.second) == itElement.first);
            ss << itElement.first;
            e.AddMessage("[NuTo::LoadSurfaceBase2D::LoadSurfaceBase2D] Error calculating surfaces for surface loads in element " + ss.str() + "(Maybe not a solid element?).");
            throw e;
        } catch (...)
        {
            std::stringstream ss;
            assert(rStructure->ElementGetId(itElement.second) == itElement.first);
            ss << itElement.first;
            throw NuTo::MechanicsException("[NuTo::LoadSurfaceBase2D::LoadSurfaceBase2D] Error calculating surfaces for surface loads in element " + ss.str() + "(Maybe not a solid element?).");
        }
    }

    //set standard integration types for triangles and quads, this can be modified according to the needs
    mIntegrationType2NPtr = rStructure->GetPtrIntegrationType(IntegrationType::IntegrationType1D2NGauss1Ip);
    mIntegrationType3NPtr = rStructure->GetPtrIntegrationType(IntegrationType::IntegrationType1D2NGauss2Ip);
    mIntegrationType4NPtr = rStructure->GetPtrIntegrationType(IntegrationType::IntegrationType1D2NGauss3Ip);
    mIntegrationType5NPtr = rStructure->GetPtrIntegrationType(IntegrationType::IntegrationType1D2NGauss4Ip);

    mIntegrationType3NPtrLobatto = rStructure->GetPtrIntegrationType(IntegrationType::IntegrationType1D2NLobatto3Ip);
    mIntegrationType4NPtrLobatto = rStructure->GetPtrIntegrationType(IntegrationType::IntegrationType1D2NLobatto4Ip);
    mIntegrationType5NPtrLobatto = rStructure->GetPtrIntegrationType(IntegrationType::IntegrationType1D2NLobatto5Ip);

}

//! @brief adds the load to global sub-vectors
//! @param rActiveDofsLoadVector ... global load vector which correspond to the active dofs
//! @param rDependentDofsLoadVector ... global load vector which correspond to the dependent dofs
void NuTo::LoadSurfaceBase2D::AddLoadToGlobalSubVectors(int rLoadCase, NuTo::FullVector<double, Eigen::Dynamic>& rActiveDofsLoadVector, NuTo::FullVector<double, Eigen::Dynamic>& rDependentDofsLoadVector) const
{
    if (rLoadCase != mLoadCase)
        return;
    for (auto it : mElements2D)
    {
        const auto* elementPtr = it.first;
        int surface = it.second;

        const InterpolationBase& interpolationTypeDisps = elementPtr->GetInterpolationType()->Get(Node::DISPLACEMENTS);
        const InterpolationBase& interpolationTypeCoords = elementPtr->GetInterpolationType()->Get(Node::COORDINATES);

        IntegrationTypeBase* integrationType(0);
        switch (interpolationTypeDisps.GetTypeOrder())
        {
        case Interpolation::EQUIDISTANT1:
            integrationType = mIntegrationType2NPtr;
            break;
        case Interpolation::EQUIDISTANT2:
            integrationType = mIntegrationType3NPtr;
            break;
        case Interpolation::EQUIDISTANT3:
            integrationType = mIntegrationType4NPtr;
            break;
        case Interpolation::EQUIDISTANT4:
            integrationType = mIntegrationType5NPtr;
            break;
        case Interpolation::LOBATTO2:
            integrationType = mIntegrationType3NPtrLobatto;
            break;
        case Interpolation::LOBATTO3:
            integrationType = mIntegrationType4NPtrLobatto;
            break;
        case Interpolation::LOBATTO4:
            integrationType = mIntegrationType5NPtrLobatto;
            break;
        case Interpolation::SPLINE:
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
                throw MechanicsException(__PRETTY_FUNCTION__, "Interpolation for exact integration of " + std::to_string(degree) + " IGA not implemented");
            }
            break;
        }
        default:
            throw MechanicsException("[NuTo::LoadSurfaceBase2D::LoadSurfaceBase2D] integration types only for 2, 3, 4 and 5 nodes (on the surface) implemented.");
        }

        Eigen::VectorXd nodeCoordinates = elementPtr->ExtractNodeValues(0, Node::COORDINATES);

        Eigen::Matrix<double, 1, 1> ipCoordsSurface;
        Eigen::Matrix<double, 2, 1> ipCoordsNatural;
        NuTo::FullVector<double, 2> ipCoordsGlobal;

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
            double tmp;
            integrationType->GetLocalIntegrationPointCoordinates1D(theIp, tmp);
            ipCoordsSurface(0) = tmp;

            if (interpolationTypeDisps.GetTypeOrder() == Interpolation::SPLINE)
                ipCoordsNatural = interpolationTypeCoords.CalculateNaturalSurfaceCoordinates(ipCoordsSurface, surface, elementPtr->GetKnots());
            else
                ipCoordsNatural = interpolationTypeCoords.CalculateNaturalSurfaceCoordinates(ipCoordsSurface, surface);

            ipCoordsGlobal = interpolationTypeCoords.CalculateMatrixN(ipCoordsNatural) * nodeCoordinates;

            // #######################################
            // ##  Calculate the surface jacobian
            // ## = || [dX / dXi] * [dXi / dAlpha] ||
            // #######################################

            derivativeShapeFunctionsNatural = interpolationTypeCoords.CalculateDerivativeShapeFunctionsNatural(ipCoordsNatural);

            const Eigen::Matrix2d jacobianStd = elementPtr->CalculateJacobian(derivativeShapeFunctionsNatural, nodeCoordinates);// = [dX / dXi]

            // in case of non IGA - just an identity matrix
            const Eigen::Matrix2d jacobianIGA = elementPtr->CalculateJacobianParametricSpaceIGA();// = [dXi / d\tilde{Xi}]

            const Eigen::Matrix2d jacobian = jacobianStd*jacobianIGA;

            derivativeNaturalSurfaceCoordinates = interpolationTypeCoords.CalculateDerivativeNaturalSurfaceCoordinates(ipCoordsSurface, surface); // = [dXi / dAlpha]
            double detJacobian = (jacobian * derivativeNaturalSurfaceCoordinates).norm();                           // = || [dX / dXi] * [dXi / dAlpha] ||

            // #######################################
            // ##  Calculate surface normal vector
            // ## = ( dY / dAlpha     -dX/dAlpha).T
            // #######################################
            // dXdAlpha :  [2 x NumNodes*2] * [NumNodes*2 x 1] * [2 x 1] = [2 x 1]
            Eigen::Vector2d surfaceTangentVector = jacobian * derivativeNaturalSurfaceCoordinates;
            surfaceTangentVector.normalize();
            NuTo::FullVector<double, 2> surfaceNormalVector;
            surfaceNormalVector(0) = surfaceTangentVector.at(1, 0);
            surfaceNormalVector(1) = -surfaceTangentVector.at(0, 0);

//            std::cout << "Global IP coordinate:        " << ipCoordsGlobal.transpose()      << std::endl;
//            std::cout << "Surface tangent vector @ IP: " << surfaceTangentVector.transpose()<< std::endl;
//            std::cout << "Surface normal vector  @ IP: " << surfaceNormalVector.transpose() << std::endl;

            // calculate 2D shape functions
            shapeFunctions = interpolationTypeDisps.CalculateShapeFunctions(ipCoordsNatural);

            //calculate weighting factor
            double thickness = elementPtr->GetSection()->GetThickness();
            double factor = thickness * (integrationType->GetIntegrationPointWeight(theIp)) * detJacobian;

            //calculate surface load
            FullVector<double, 2> loadVector;
            CalculateSurfaceLoad(ipCoordsGlobal, surfaceNormalVector, loadVector);
            loadVector *= factor;

//			std::cout << "load vector with weights \n" << loadVector << std::endl;
//			std::cout << "detJ " << detJacobian << " weight " << integrationType->GetIntegrationPointWeight(theIp) << std::endl;

            //add load vector to global vector
            for (int iNode = 0; iNode < shapeFunctions.rows(); iNode++)
            {
                const NodeBase* node = elementPtr->GetNode(interpolationTypeDisps.GetNodeIndex(iNode));
                assert(node->GetNum(Node::DISPLACEMENTS) == 2);
                for (int iDispDof = 0; iDispDof < 2; iDispDof++)
                {
                    int theDof = node->GetDof(Node::DISPLACEMENTS, iDispDof);
                    double theLoad = shapeFunctions[iNode] * loadVector(iDispDof);
                    if (theDof < rActiveDofsLoadVector.GetNumRows())
                    {
//						std::cout << "add to dof " << theDof << " " << theLoad << std::endl;
                        rActiveDofsLoadVector(theDof) += theLoad;
                    }
                    else
                    {
                        rDependentDofsLoadVector(theDof - rActiveDofsLoadVector.GetNumRows()) += theLoad;
                    }
                }
            }
        }
    }
}

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::LoadSurfaceBase2D)
#endif
