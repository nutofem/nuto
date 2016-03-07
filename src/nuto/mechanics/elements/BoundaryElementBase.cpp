/*
 * BoundaryElementBase.cpp
 *
 *  Created on: 5 Jun 2015
 *      Author: ttitsche
 */

#include "nuto/mechanics/elements/BoundaryElementBase.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeBase.h"
#include "nuto/mechanics/structures/StructureBase.h"

NuTo::BoundaryElementBase::BoundaryElementBase(const ElementBase* rBaseElement, int rSurfaceId)
: ElementBase::ElementBase(rBaseElement->GetStructure(), rBaseElement->GetElementDataType(), rBaseElement->GetIpDataType(0), rBaseElement->GetInterpolationType()),
  mBaseElement(rBaseElement), mSurfaceId(rSurfaceId), mBoundaryConditionType(BoundaryType::NOT_SET)
{

}
void NuTo::BoundaryElementBase::SetSection(const SectionBase* rSection)
{
}

const NuTo::SectionBase* NuTo::BoundaryElementBase::GetSection() const
{
    return mBaseElement->GetSection();
}

const Eigen::VectorXd NuTo::BoundaryElementBase::GetIntegrationPointVolume() const
{
    throw MechanicsException("[NuTo::BoundaryElementBase::GetIntegrationPointVolume] Not implemented.");
}

const Eigen::Vector3d NuTo::BoundaryElementBase::GetGlobalIntegrationPointCoordinates(int rIpNum) const
{
    Eigen::VectorXd naturalSurfaceIpCoordinates;
    switch (GetStructure()->GetDimension())
    {
        case 1:
        {
            double ipCoordinate;
            GetIntegrationType()->GetLocalIntegrationPointCoordinates1D(rIpNum, ipCoordinate);
            naturalSurfaceIpCoordinates.resize(1);
            naturalSurfaceIpCoordinates(0) = ipCoordinate;
            break;
        }
        case 2:
        {
            double ipCoordinates[2];
            GetIntegrationType()->GetLocalIntegrationPointCoordinates2D(rIpNum, ipCoordinates);
            naturalSurfaceIpCoordinates.resize(2);
            naturalSurfaceIpCoordinates(0) = ipCoordinates[0];
            naturalSurfaceIpCoordinates(1) = ipCoordinates[1];
            break;
        }
        case 3:
        {
            double ipCoordinates[3];
            GetIntegrationType()->GetLocalIntegrationPointCoordinates3D(rIpNum, ipCoordinates);
            naturalSurfaceIpCoordinates.resize(3);
            naturalSurfaceIpCoordinates(0) = ipCoordinates[0];
            naturalSurfaceIpCoordinates(1) = ipCoordinates[1];
            naturalSurfaceIpCoordinates(2) = ipCoordinates[2];
            break;
        }
        default:
            break;
    }


    Eigen::VectorXd naturalIpCoordinates = mInterpolationType->Get(Node::COORDINATES).CalculateNaturalSurfaceCoordinates(naturalSurfaceIpCoordinates, mSurfaceId);

    Eigen::VectorXd shapeFunctions = mInterpolationType->Get(Node::COORDINATES).CalculateShapeFunctions(naturalIpCoordinates);
    Eigen::MatrixXd nodeCoordinates = ExtractNodeValues(0, Node::COORDINATES);
    assert(shapeFunctions.rows() == nodeCoordinates.cols());

    Eigen::Vector3d globalIntegrationPointCoordinates = Eigen::Vector3d::Zero();
    for (int iRow = 0; iRow < nodeCoordinates.rows(); iRow++)
        globalIntegrationPointCoordinates[iRow] = nodeCoordinates.row(iRow) * shapeFunctions;

    return globalIntegrationPointCoordinates;
}

const Eigen::MatrixXd NuTo::BoundaryElementBase::ExtractNodeValues(int rTimeDerivative, Node::eAttributes rAttribute) const
{
    return mBaseElement->ExtractNodeValues(rTimeDerivative, rAttribute);
}

NuTo::BoundaryType::eType NuTo::BoundaryElementBase::GetBoundaryConditionType() const
{
    return mBoundaryConditionType;
}

void NuTo::BoundaryElementBase::SetBoundaryConditionType(BoundaryType::eType rBoundaryConditionType)
{
    mBoundaryConditionType = rBoundaryConditionType;
}

#ifdef ENABLE_VISUALIZE
void NuTo::BoundaryElementBase::Visualize(VisualizeUnstructuredGrid& rVisualize, const std::list<std::shared_ptr<NuTo::VisualizeComponent>>& rVisualizationList)
{
    if (GetStructure()->GetVerboseLevel() > 10)
        std::cout << "[NuTo::BoundaryElementBase::Visualize] Pleeeaaase, implement the visualization for me!!!" << std::endl;
}
#endif
