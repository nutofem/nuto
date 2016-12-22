/*
 * StructureInterpolationType.cpp
 *
 *  Created on: 31 Mar 2015
 *      Author: ttitsche
 */

#include "mechanics/groups/GroupEnum.h"
#include "mechanics/elements/ElementBase.h"
#include "mechanics/elements/IpDataEnum.h"
#include "math/FullMatrix.h"
#include "mechanics/structures/unstructured/Structure.h"
#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/integrationtypes/IntegrationTypeBase.h"
#include "mechanics/interpolationtypes/InterpolationType.h"
#include "mechanics/interpolationtypes/InterpolationTypeEnum.h"

//! @brief creates a new interpolation type, calls the enum method
//! @param rShape ... element shape "TRUSS", "TRIANGLE", "QUAD", "TET", "BRICK", etc.
//! @return ... interpolation type id
int NuTo::Structure::InterpolationTypeCreate(const std::string& rShape)
{
    Interpolation::eShapeType shapeType = Interpolation::ShapeTypeToEnum(rShape);

    //find unused integer id
    int interpolationTypeNumber = mInterpolationTypeMap.size();
    boost::ptr_map<int,InterpolationType>::iterator it = mInterpolationTypeMap.find(interpolationTypeNumber);
    while (it!=mInterpolationTypeMap.end())
    {
        interpolationTypeNumber++;
        it = mInterpolationTypeMap.find(interpolationTypeNumber);
    }

    this->InterpolationTypeCreate(interpolationTypeNumber, shapeType);
    return interpolationTypeNumber;
}

int NuTo::Structure::InterpolationTypeCreate(NuTo::Interpolation::eShapeType rShape)
{

    //find unused integer id
    int interpolationTypeNumber = mInterpolationTypeMap.size();
    boost::ptr_map<int,InterpolationType>::iterator it = mInterpolationTypeMap.find(interpolationTypeNumber);
    while (it!=mInterpolationTypeMap.end())
    {
        interpolationTypeNumber++;
        it = mInterpolationTypeMap.find(interpolationTypeNumber);
    }

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
    auto iterator = mInterpolationTypeMap.find(rInterpolationTypeId);
    if (iterator == mInterpolationTypeMap.end())
        throw MechanicsException("[NuTo::Structure::InterpolationTypeSetIntegrationType] InterpolationType with ID"
                + std::to_string(iterator->first) + " does not exist.");

    InterpolationType* interpolationType = iterator->second;
    interpolationType->UpdateIntegrationType(*rIntegrationType);

    // update all elements
    // disable show time

    bool showTime = GetShowTime();
    SetShowTime(false);

    int elementGroupId = GroupCreate("Elements");
    GroupAddElementFromType(elementGroupId, rInterpolationTypeId);

    NuTo::FullVector<int, Eigen::Dynamic> elementIds = GroupGetMemberIds(elementGroupId);
    for (int iElement = 0; iElement < elementIds.GetNumRows(); ++iElement)
    {
        ElementBase* element = ElementGetElementPtr(elementIds.GetValue(iElement));
        element->SetIntegrationType(*rIntegrationType);
    }

    GroupDelete(elementGroupId);
    SetShowTime(showTime);
}



//! @brief prints the info to the interpolation type
//! @param rInterpolationTypeId ... interpolation type id
void NuTo::Structure::InterpolationTypeInfo(int rInterpolationTypeId) const
{
    boost::ptr_map<int,InterpolationType>::const_iterator itIterator = mInterpolationTypeMap.find(rInterpolationTypeId);
    // check if identifier exists
    if (itIterator == mInterpolationTypeMap.end())
        throw NuTo::MechanicsException("[NuTo::Structure::InterpolationTypeSetIsConstitutiveInput] Interpolation type does not exist.");


    mLogger << itIterator->second->Info();
}

//! @brief creates a new interpolation type
//! @param rInterpolationTypeId ... interpolation type id
//! @param rShape ... element shape "1DTRUSS", "2DTRIANGLE", "2DQUAD", "3DTET", "3DBRICK", etc.
//! @return ... index in the mInterpolationType map
void NuTo::Structure::InterpolationTypeCreate(int rInterpolationTypeId, NuTo::Interpolation::eShapeType rShape)
{
    // check if constitutive law identifier exists
    if (mInterpolationTypeMap.find(rInterpolationTypeId) != mInterpolationTypeMap.end())
        throw NuTo::MechanicsException("[NuTo::Structure::InterpolationTypeCreate] Interpolation type already exists.");

    mInterpolationTypeMap.insert(rInterpolationTypeId, new InterpolationType(rShape, GetDimension()));

}

//! @brief adds a dof to a interpolation type, calls the enum method
//! @param rInterpolationTypeId ... interpolation type id
//! @param rDofType ... dof type
//! @param rTypeOrder ... type and order of interpolation
void NuTo::Structure::InterpolationTypeAdd(int rInterpolationTypeId, const std::string& rDofType, const std::string& rTypeOrder)
{
    InterpolationTypeAdd(rInterpolationTypeId,Node::DofToEnum(rDofType), Interpolation::TypeOrderToEnum(rTypeOrder));
}

//! @brief adds a dof to IGA interpolation type
//! @param rInterpolationTypeId ... interpolation type id
//! @param rDofType ... dof type
//! @param rTypeOrder ... type and order of interpolation
void NuTo::Structure::InterpolationTypeAdd(int rInterpolationTypeId,
                                           NuTo::Node::eDof rDofType,
                                           NuTo::Interpolation::eTypeOrder rTypeOrder,
                                           const Eigen::VectorXi &rDegree,
                                           const std::vector<Eigen::VectorXd> &rKnots,
                                           const Eigen::MatrixXd &rWeights)
{
    boost::ptr_map<int,InterpolationType>::iterator itIterator = mInterpolationTypeMap.find(rInterpolationTypeId);
    // check if identifier exists
    if (itIterator == mInterpolationTypeMap.end())
        throw NuTo::MechanicsException("[NuTo::Structure::InterpolationTypeAdd] Interpolation type does not exist.");

    InterpolationType* interpolationType = itIterator->second;
    interpolationType->AddDofInterpolation(rDofType, rTypeOrder, rDegree, rKnots, rWeights);

    eIntegrationType integrationTypeEnum = interpolationType->GetStandardIntegrationType();
    const IntegrationTypeBase& integrationType = *this->GetPtrIntegrationType(integrationTypeEnum);

    interpolationType->UpdateIntegrationType(integrationType);
    if (mVerboseLevel > 2)
        mLogger << "[NuTo::Structure::InterpolationTypeAdd] Updated IntegrationType to " << integrationType.GetStrIdentifier() << ".\n";

    // update all elements
    // disable show time

    bool showTime = GetShowTime();
    SetShowTime(false);

    int elementGroupId = GroupCreate("Elements");
    GroupAddElementFromType(elementGroupId, rInterpolationTypeId);

    NuTo::FullVector<int, Eigen::Dynamic> elementIds = GroupGetMemberIds(elementGroupId);
    for (int iElement = 0; iElement < elementIds.GetNumRows(); ++iElement)
    {
        ElementBase* element = ElementGetElementPtr(elementIds.GetValue(iElement));
        element->SetIntegrationType(integrationType);
    }

    GroupDelete(elementGroupId);
    UpdateDofStatus();
    SetShowTime(showTime);
}

//! @brief adds a dof to a interpolation type
//! @param rInterpolationTypeId ... interpolation type id
//! @param rDofType ... dof type
//! @param rTypeOrder ... type and order of interpolation
void NuTo::Structure::InterpolationTypeAdd(int rInterpolationTypeId, NuTo::Node::eDof rDofType, NuTo::Interpolation::eTypeOrder rTypeOrder)
{
    boost::ptr_map<int,InterpolationType>::iterator itIterator = mInterpolationTypeMap.find(rInterpolationTypeId);
    // check if identifier exists
    if (itIterator == mInterpolationTypeMap.end())
        throw NuTo::MechanicsException("[NuTo::Structure::InterpolationTypeAdd] Interpolation type does not exist.");

    InterpolationType* interpolationType = itIterator->second;
    interpolationType->AddDofInterpolation(rDofType, rTypeOrder);

    eIntegrationType integrationTypeEnum = interpolationType->GetStandardIntegrationType();
    const IntegrationTypeBase& integrationType = *this->GetPtrIntegrationType(integrationTypeEnum);

    interpolationType->UpdateIntegrationType(integrationType);
    if (mVerboseLevel > 2)
        mLogger << "[NuTo::Structure::InterpolationTypeAdd] Updated IntegrationType to " << integrationType.GetStrIdentifier() << ".\n";

    // update all elements
    // disable show time

    bool showTime = GetShowTime();
    SetShowTime(false);

    int elementGroupId = GroupCreate(eGroupId::Elements);
    GroupAddElementFromType(elementGroupId, rInterpolationTypeId);

    NuTo::FullVector<int, Eigen::Dynamic> elementIds = GroupGetMemberIds(elementGroupId);
    for (int iElement = 0; iElement < elementIds.GetNumRows(); ++iElement)
    {
        ElementBase* element = ElementGetElementPtr(elementIds.GetValue(iElement));
        element->SetIntegrationType(integrationType);
    }

    GroupDelete(elementGroupId);
    UpdateDofStatus();
    SetShowTime(showTime);
}
