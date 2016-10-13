#include "nuto/base/ErrorEnum.h"
#include "nuto/math/FullMatrix.h"
#include "nuto/mechanics/constitutive/ConstitutiveBase.h"
#include "nuto/mechanics/dofSubMatrixStorage/BlockFullMatrix.h"
#include "nuto/mechanics/elements/ContinuumContactElement.h"
#include "nuto/mechanics/elements/ContinuumElement.h"
#include "nuto/mechanics/elements/ElementDataBase.h"
#include "nuto/mechanics/elements/ElementEnum.h"
#include "nuto/mechanics/elements/ElementOutputBase.h"
#include "nuto/mechanics/elements/ElementOutputIpData.h"
#include "nuto/mechanics/elements/EvaluateDataContinuumBoundary.h"
#include "nuto/mechanics/elements/IpDataEnum.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeBase.h"
#include "nuto/mechanics/interpolationtypes/InterpolationBase.h"
#include "nuto/mechanics/interpolationtypes/InterpolationType.h"
#include "nuto/mechanics/nodes/NodeEnum.h"
#include "nuto/mechanics/sections/SectionTruss.h"
#include "nuto/mechanics/sections/SectionPlane.h"
#include "nuto/mechanics/structures/StructureBase.h"
#include "nuto/mechanics/constitutive/ConstitutiveEnum.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveIOMap.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveScalar.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveVector.h"
#include "nuto/mechanics/constitutive/inputoutput/EngineeringStrain.h"
#include "nuto/mechanics/constitutive/inputoutput/EngineeringStress.h"
#include "nuto/mechanics/groups/Group.h"
#include "nuto/mechanics/nodes/NodeBase.h"

template <int TDim>
NuTo::ContinuumContactElement<TDim>::ContinuumContactElement(const ContinuumElement<TDim> *rSlaveElement,
                                                             int rSurfaceId,
                                                             int rElementGroupId,
                                                             int rNodeGroupId,
                                                             const IntegrationTypeBase *rIntegrationType)
    : ContinuumBoundaryElement<TDim>(rSlaveElement, rSurfaceId), mIntegrationType(rIntegrationType)
{
    //get element group
    const Group<ElementBase> *elementGroup = this->mStructure->GroupGetGroupPtr(rElementGroupId)->AsGroupElement();

    //get node group
    const Group<NodeBase> *nodeGroup = this->mStructure->GroupGetGroupPtr(rNodeGroupId)->AsGroupNode();

    //since the search is done via the id's, the surface nodes are ptr, so make another set with the node ptrs
    std::set<const NodeBase*> nodePtrSet;
    for (Group<NodeBase>::const_iterator itNode = nodeGroup->begin(); itNode != nodeGroup->end(); itNode++)
    {
        nodePtrSet.insert(itNode->second);
    }

    // std::cout << "number of loaded nodes " << nodeGroup->GetNumMembers() << std::endl;

    //loop over all elements
    Eigen::VectorXi surfaceNodeIndices;
    std::vector<const NodeBase*> surfaceNodes;

    for (auto &itElement : *elementGroup)
    {
        try
        {
            ElementBase* base = itElement.second;
            //check if plane element
            ContinuumElement<TDim>* elementPtr = dynamic_cast<ContinuumElement<TDim>*>(base);
            const InterpolationType* interpolationType = elementPtr->GetInterpolationType();

            //loop over all surfaces
            for (int iSurface = 0; iSurface < interpolationType->GetNumSurfaces(); iSurface++)
            {
                bool addSurface = true;
                surfaceNodeIndices = interpolationType->GetSurfaceNodeIndices(iSurface);
                int numSurfaceNodes = surfaceNodeIndices.rows();
                surfaceNodes.resize(numSurfaceNodes);

                for (int iSurfaceNode = 0; iSurfaceNode < numSurfaceNodes; ++iSurfaceNode)
                {
                    surfaceNodes[iSurfaceNode] = elementPtr->GetNode(surfaceNodeIndices(iSurfaceNode, 0));
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
                    mElementsMaster.push_back(std::make_pair(elementPtr, iSurface));
                }
            }
        } catch (NuTo::MechanicsException &e)
        {
            std::stringstream ss;
            assert(this->mStructure->ElementGetId(itElement.second) == itElement.first);
            ss << itElement.first;
            e.AddMessage("[NuTo::ContinuumContactElement] Error calculating surfaces for surface loads in element " + ss.str() + "(Maybe not a solid element?).");
            throw e;
        } catch (...)
        {
            std::stringstream ss;
            assert(this->mStructure->ElementGetId(itElement.second) == itElement.first);
            ss << itElement.first;
            throw NuTo::MechanicsException("[NuTo::ContinuumContactElement] Error calculating surfaces for surface loads in element " + ss.str() + "(Maybe not a solid element?).");
        }
    }
}

template<int TDim>
void NuTo::ContinuumContactElement<TDim>::CalculateElementOutputs(std::map<Element::eOutput, std::shared_ptr<ElementOutputBase>>& rElementOutput,
                                                                  EvaluateDataContinuumBoundary<TDim> &rData,
                                                                  int rTheIP,
                                                                  const ConstitutiveInputMap& constitutiveInput,
                                                                  const ConstitutiveOutputMap& constitutiveOutput) const
{
    rData.mDetJxWeightIPxSection = this->CalculateDetJxWeightIPxSection(rData.mDetJacobian, rTheIP); // formerly known as "factor"

    for (auto it : rElementOutput)
    {
        switch (it.first)
        {
        case Element::eOutput::GAP_MATRIX_MORTAR:
            CalculateElementOutputGapMatrixMortar(it.second->GetBlockFullMatrixDouble(), rData, constitutiveOutput, rTheIP);
            break;
        case Element::eOutput::INTERNAL_GRADIENT:
        case Element::eOutput::HESSIAN_0_TIME_DERIVATIVE:
        case Element::eOutput::HESSIAN_1_TIME_DERIVATIVE:
        case Element::eOutput::HESSIAN_2_TIME_DERIVATIVE:
        case Element::eOutput::LUMPED_HESSIAN_2_TIME_DERIVATIVE:
        case Element::eOutput::UPDATE_STATIC_DATA:
        case Element::eOutput::UPDATE_TMP_STATIC_DATA:
            break;
        case Element::eOutput::IP_DATA:
            this->CalculateElementOutputIpData(it.second->GetIpData(), constitutiveOutput, rTheIP);
            break;
        case Element::eOutput::GLOBAL_ROW_DOF:
        case Element::eOutput::GLOBAL_COLUMN_DOF:
            break;
        default:
            throw MechanicsException(__PRETTY_FUNCTION__, "element output not implemented.");
        }
    }
}

template<int TDim>
void NuTo::ContinuumContactElement<TDim>::ProjectIntegrationPointOnMaster()
{

}

template<int TDim>
void NuTo::ContinuumContactElement<TDim>::CalculateElementOutputGapMatrixMortar(BlockFullMatrix<double>& rGapMatrix,
                                                                                EvaluateDataContinuumBoundary<TDim> &rData,
                                                                                const ConstitutiveOutputMap& constitutiveOutput,
                                                                                int rTheIP) const
{
    // ===> Projection of the rTheIP on the master element => \xi^s_{IP}, \xi^m_*, n^m_*

    // ===> Get the position \xi^s_{IP}
    Eigen::VectorXd coordinatedIPSlave;

    double temp[2];
    switch (TDim)
    {
    case 1:
    {
        mIntegrationType->GetLocalIntegrationPointCoordinates1D(rTheIP, temp[0]);
        coordinatedIPSlave.resize(1);
        coordinatedIPSlave(0) = temp[0];
        break;
    }
    case 2:
    {
        mIntegrationType->GetLocalIntegrationPointCoordinates2D(rTheIP, temp);
        coordinatedIPSlave.resize(2);
        coordinatedIPSlave << temp[0], temp[1];
        break;
    }
    default:
        throw MechanicsException(__PRETTY_FUNCTION__, "Contact between 2D and 3D bodies implemented.");
        break;
    }

    // ===> Get the starting point for iteration
    double minDistance = std::numeric_limits<double>::infinity();
    Eigen::VectorXd parameterMin;
    for(auto &it : mElementsMaster)
    {
        const auto* elementPtr = it.first;
        int surfaceId = it.second;

        // ===> Get the position on the master curve/surface
        // ===> Compare and set to minimum if though

        const InterpolationBase& interpolationTypeCoords = elementPtr->GetInterpolationType()->Get(Node::eDof::COORDINATES);
        Eigen::VectorXd referenceCoordinates(1);
        referenceCoordinates(0);
        Eigen::VectorXd parameter = interpolationTypeCoords.CalculateNaturalSurfaceCoordinates(referenceCoordinates, surfaceId, elementPtr->GetKnots());
        Eigen::VectorXd coordinatesMaster = elementPtr->InterpolateDofGlobalSurfaceDerivative(0, parameter, 0, 0);
        double distance = (coordinatesMaster - coordinatedIPSlave).norm();
        if(minDistance < distance)
        {
            minDistance = distance;
            parameterMin = parameter;
        }
    }

    // ===> Newton for projection (get \xi^m_*, n^m_*)


    // ===> Assemble the element gap matrix => \int_{\Gamma_e} F(ShapeFunctionsSlave(\xi^s), ShapeFunctionsMaster(\xi^*), n^*) d\Gamma
}

template class NuTo::ContinuumContactElement<1>;
template class NuTo::ContinuumContactElement<2>;

