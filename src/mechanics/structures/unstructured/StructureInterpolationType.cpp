/*
 * StructureInterpolationType.cpp
 *
 *  Created on: 31 Mar 2015
 *      Author: ttitsche
 */

#include "mechanics/groups/GroupEnum.h"
#include "mechanics/elements/ElementBase.h"
#include "mechanics/elements/IpDataEnum.h"

#include "mechanics/structures/unstructured/Structure.h"
#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/integrationtypes/IntegrationTypeBase.h"
#include "mechanics/interpolationtypes/InterpolationType.h"
#include "mechanics/interpolationtypes/InterpolationTypeEnum.h"

int NuTo::Structure::InterpolationTypeCreate(const std::string& rShape)
{
    return InterpolationTypeCreate(Interpolation::ShapeTypeToEnum(rShape));
}

int NuTo::Structure::InterpolationTypeCreate(NuTo::Interpolation::eShapeType rShape)
{
    int interpolationTypeNumber = GetUnusedId(mInterpolationTypeMap);
    this->InterpolationTypeCreate(interpolationTypeNumber, rShape);
    return interpolationTypeNumber;
}


void NuTo::Structure::InterpolationTypeSetIntegrationType(int rInterpolationTypeId, const std::string& rIntegrationType)
{
    InterpolationTypeSetIntegrationType(rInterpolationTypeId, GetPtrIntegrationType(rIntegrationType));
}

void NuTo::Structure::InterpolationTypeSetIntegrationType(int rInterpolationTypeId, eIntegrationType rIntegrationType)
{
    InterpolationTypeSetIntegrationType(rInterpolationTypeId, GetPtrIntegrationType(rIntegrationType));
}

void NuTo::Structure::InterpolationTypeSetIntegrationType(int rInterpolationTypeId, IntegrationTypeBase* rIntegrationType)
{
    InterpolationType* interpolationType = InterpolationTypeGet(rInterpolationTypeId);
    interpolationType->ClearCache();

    // update all elements
    // disable show time

    bool showTime = GetShowTime();
    SetShowTime(false);

    int elementGroupId = GroupCreate("Elements");
    GroupAddElementFromType(elementGroupId, rInterpolationTypeId);

    for (int elementId : GroupGetMemberIds(elementGroupId))
    {
        ElementBase* element = ElementGetElementPtr(elementId);
        element->SetIntegrationType(*rIntegrationType);
    }

    GroupDelete(elementGroupId);
    SetShowTime(showTime);
}


void NuTo::Structure::InterpolationTypeInfo(int rInterpolationTypeId) const
{
    mLogger << InterpolationTypeGet(rInterpolationTypeId)->Info();
}

void NuTo::Structure::InterpolationTypeCreate(int rInterpolationTypeId, NuTo::Interpolation::eShapeType rShape)
{
    // check if constitutive law identifier exists
    if (mInterpolationTypeMap.find(rInterpolationTypeId) != mInterpolationTypeMap.end())
        throw NuTo::MechanicsException("[NuTo::Structure::InterpolationTypeCreate] Interpolation type already exists.");

    mInterpolationTypeMap.insert(rInterpolationTypeId, new InterpolationType(rShape, GetDimension()));

}

void NuTo::Structure::InterpolationTypeAdd(int rInterpolationTypeId, const std::string& rDofType, const std::string& rTypeOrder)
{
    InterpolationTypeAdd(rInterpolationTypeId,Node::DofToEnum(rDofType), Interpolation::TypeOrderToEnum(rTypeOrder));
}

void NuTo::Structure::InterpolationTypeAdd(int rInterpolationTypeId,
                                           NuTo::Node::eDof rDofType,
                                           NuTo::Interpolation::eTypeOrder rTypeOrder,
                                           const Eigen::VectorXi &rDegree,
                                           const std::vector<Eigen::VectorXd> &rKnots,
                                           const Eigen::MatrixXd &rWeights)
{
    InterpolationType* interpolationType = InterpolationTypeGet(rInterpolationTypeId);
    interpolationType->AddDofInterpolation(rDofType, rTypeOrder, rDegree, rKnots, rWeights);

    eIntegrationType integrationTypeEnum = interpolationType->GetStandardIntegrationType();
    const IntegrationTypeBase& integrationType = *this->GetPtrIntegrationType(integrationTypeEnum);

    interpolationType->ClearCache();
    if (mVerboseLevel > 2)
        mLogger << "[NuTo::Structure::InterpolationTypeAdd] Updated IntegrationType to " << IntegrationTypeToString(integrationType.GetEnumType()) << ".\n";

    // update all elements
    // disable show time

    bool showTime = GetShowTime();
    SetShowTime(false);

    int elementGroupId = GroupCreate("Elements");
    GroupAddElementFromType(elementGroupId, rInterpolationTypeId);

    for (int elementId : GroupGetMemberIds(elementGroupId))
    {
        ElementBase* element = ElementGetElementPtr(elementId);
        element->SetIntegrationType(integrationType);
    }

    GroupDelete(elementGroupId);
    UpdateDofStatus();
    SetShowTime(showTime);
}

void NuTo::Structure::InterpolationTypeAdd(int rInterpolationTypeId, NuTo::Node::eDof rDofType, NuTo::Interpolation::eTypeOrder rTypeOrder)
{
    InterpolationType* interpolationType = InterpolationTypeGet(rInterpolationTypeId);
    interpolationType->AddDofInterpolation(rDofType, rTypeOrder);

    eIntegrationType integrationTypeEnum = interpolationType->GetStandardIntegrationType();
    const IntegrationTypeBase& integrationType = *this->GetPtrIntegrationType(integrationTypeEnum);

    interpolationType->ClearCache();
    if (mVerboseLevel > 2)
        mLogger << "[NuTo::Structure::InterpolationTypeAdd] Updated IntegrationType to " << IntegrationTypeToString(integrationType.GetEnumType()) << ".\n";

    // update all elements
    // disable show time

    bool showTime = GetShowTime();
    SetShowTime(false);

    int elementGroupId = GroupCreate(eGroupId::Elements);
    GroupAddElementFromType(elementGroupId, rInterpolationTypeId);

    for (int elementId : GroupGetMemberIds(elementGroupId))
    {
        ElementBase* element = ElementGetElementPtr(elementId);
        element->SetIntegrationType(integrationType);
    }

    GroupDelete(elementGroupId);
    UpdateDofStatus();
    SetShowTime(showTime);
}

NuTo::InterpolationType* NuTo::Structure::InterpolationTypeGet(int rInterpolationTypeId)
{
    auto itIterator = mInterpolationTypeMap.find(rInterpolationTypeId);
    // check if identifier exists
    if (itIterator == mInterpolationTypeMap.end())
        throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "Interpolation type does not exist.");
    return itIterator->second;
}

const NuTo::InterpolationType* NuTo::Structure::InterpolationTypeGet(int rInterpolationTypeId) const
{
    auto itIterator = mInterpolationTypeMap.find(rInterpolationTypeId);
    // check if identifier exists
    if (itIterator == mInterpolationTypeMap.end())
        throw NuTo::MechanicsException(__PRETTY_FUNCTION__, "Interpolation type does not exist.");
    return itIterator->second;
}
