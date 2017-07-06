#include "mechanics/constitutive/ConstitutiveBase.h"
#include "mechanics/dofSubMatrixStorage/BlockFullMatrix.h"
#include "mechanics/elements/ContinuumContactElement.h"
#include "mechanics/elements/ContinuumElement.h"
#include "mechanics/elements/ElementEnum.h"
#include "mechanics/elements/ElementOutputBase.h"
#include "mechanics/elements/ElementOutputIpData.h"
#include "mechanics/elements/EvaluateDataContinuumBoundary.h"
#include "mechanics/integrationtypes/IntegrationTypeBase.h"
#include "mechanics/interpolationtypes/InterpolationBase.h"
#include "mechanics/interpolationtypes/InterpolationType.h"
#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/sections/SectionTruss.h"
#include "mechanics/sections/SectionPlane.h"
#include "mechanics/constitutive/ConstitutiveEnum.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveIOMap.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveScalar.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveVector.h"
#include "mechanics/constitutive/inputoutput/EngineeringStrain.h"
#include "mechanics/constitutive/inputoutput/EngineeringStress.h"
#include "mechanics/groups/Group.h"
#include "mechanics/nodes/NodeBase.h"

template <int TDim>
NuTo::ContinuumContactElement<TDim>::ContinuumContactElement(const ContinuumElement<TDim>& rSlaveElement,
                                                             int rSurfaceId, const Group<ElementBase>* elementGroup,
                                                             const Group<NodeBase>* nodeGroup,
                                                             const IntegrationTypeBase& integrationType)
    : ContinuumBoundaryElement<TDim>(rSlaveElement, integrationType, rSurfaceId)
    , mIntegrationType(&integrationType)
{
    // since the search is done via the id's, the surface nodes are ptr, so make another set with the node ptrs
    std::set<const NodeBase*> nodePtrSet;
    for (auto node : *nodeGroup)
    {
        nodePtrSet.insert(node.second);
    }

    // std::cout << "number of loaded nodes " << nodeGroup->GetNumMembers() << std::endl;

    // loop over all elements
    Eigen::VectorXi surfaceNodeIndices;
    std::vector<const NodeBase*> surfaceNodes;

    for (auto& itElement : *elementGroup)
    {
        try
        {
            ElementBase* base = itElement.second;
            // check if plane element
            ContinuumElement<TDim>* elementPtr = dynamic_cast<ContinuumElement<TDim>*>(base);
            const InterpolationType& interpolationType = elementPtr->GetInterpolationType();

            // loop over all surfaces
            for (int iSurface = 0; iSurface < interpolationType.GetNumSurfaces(); iSurface++)
            {
                bool addSurface = true;
                surfaceNodeIndices = interpolationType.GetSurfaceNodeIndices(iSurface);
                int numSurfaceNodes = surfaceNodeIndices.rows();
                surfaceNodes.resize(numSurfaceNodes);

                for (int iSurfaceNode = 0; iSurfaceNode < numSurfaceNodes; ++iSurfaceNode)
                {
                    surfaceNodes[iSurfaceNode] = elementPtr->GetNode(surfaceNodeIndices(iSurfaceNode, 0));
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
                    mElementsMaster.push_back(std::make_pair(elementPtr, iSurface));
                }
            }
        }
        catch (NuTo::MechanicsException& e)
        {
            std::stringstream ss;
            ss << itElement.first;
            e.AddMessage("[NuTo::ContinuumContactElement] Error calculating surfaces for surface loads in element " +
                         ss.str() + "(Maybe not a solid element?).");
            throw;
        }
        catch (...)
        {
            std::stringstream ss;
            ss << itElement.first;
            throw NuTo::MechanicsException(
                    "[NuTo::ContinuumContactElement] Error calculating surfaces for surface loads in element " +
                    ss.str() + "(Maybe not a solid element?).");
        }
    }
}

template <int TDim>
void NuTo::ContinuumContactElement<TDim>::CalculateElementOutputs(
        std::map<Element::eOutput, std::shared_ptr<ElementOutputBase>>& rElementOutput,
        EvaluateDataContinuumBoundary<TDim>& rData, int rTheIP, const ConstitutiveInputMap&,
        const ConstitutiveOutputMap& constitutiveOutput) const
{
    rData.mDetJxWeightIPxSection =
            this->CalculateDetJxWeightIPxSection(rData.mDetJacobian, rTheIP); // formerly known as "factor"

    for (auto it : rElementOutput)
    {
        switch (it.first)
        {
        case Element::eOutput::GAP_MATRIX_MORTAR:
            CalculateElementOutputGapMatrixMortar(it.second->GetBlockFullMatrixDouble(), rData, constitutiveOutput,
                                                  rTheIP);
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

template <int TDim>
void NuTo::ContinuumContactElement<TDim>::ProjectIntegrationPointOnMaster()
{
}

template <int TDim>
void NuTo::ContinuumContactElement<TDim>::CalculateElementOutputGapMatrixMortar(
        BlockFullMatrix<double>&, EvaluateDataContinuumBoundary<TDim>&,
        const ConstitutiveOutputMap&, int rTheIP) const
{
    // ===> Projection of the rTheIP on the master element => \xi^s_{IP}, \xi^m_*, n^m_*

    // ===> Get the position \xi^s_{IP}
    Eigen::VectorXd coordinatedIPSlave = mIntegrationType->GetLocalIntegrationPointCoordinates(rTheIP);

    // ===> Get the starting point for iteration
    double minDistance = std::numeric_limits<double>::infinity();
    Eigen::VectorXd parameterMin;
    for (auto& it : mElementsMaster)
    {
        const auto* elementPtr = it.first;
        int surfaceId = it.second;

        // ===> Get the position on the master curve/surface
        // ===> Compare and set to minimum if though

        const InterpolationBase& interpolationTypeCoords =
                elementPtr->GetInterpolationType().Get(Node::eDof::COORDINATES);
        Eigen::VectorXd referenceCoordinates(1);
        referenceCoordinates(0);
        Eigen::VectorXd parameter = interpolationTypeCoords.CalculateNaturalSurfaceCoordinatesIGA(
                referenceCoordinates, surfaceId, elementPtr->GetKnots());
        Eigen::VectorXd coordinatesMaster = elementPtr->InterpolateDofGlobalSurfaceDerivative(0, parameter, 0, 0);
        double distance = (coordinatesMaster - coordinatedIPSlave).norm();
        if (minDistance < distance)
        {
            minDistance = distance;
            parameterMin = parameter;
        }
    }

    // ===> Newton for projection (get \xi^m_*, n^m_*)


    // ===> Assemble the element gap matrix => \int_{\Gamma_e} F(ShapeFunctionsSlave(\xi^s),
    // ShapeFunctionsMaster(\xi^*), n^*) d\Gamma
}

template class NuTo::ContinuumContactElement<1>;
template class NuTo::ContinuumContactElement<2>;
