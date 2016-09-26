/*
 * StructureInterpolationType.cpp
 *
 *  Created on: 31 Mar 2015
 *      Author: ttitsche
 */

#include "nuto/mechanics/elements/ElementBase.h"
#include "nuto/mechanics/elements/IpDataEnum.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/mechanics/nodes/NodeEnum.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeBase.h"
#include "nuto/mechanics/interpolationtypes/InterpolationType.h"
#include "nuto/mechanics/interpolationtypes/InterpolationTypeEnum.h"

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

//! @brief creates a new interpolation type, calls the enum method
//! @param rShape ... element shape "TRUSS", "TRIANGLE", "QUAD", "TET", "BRICK", etc.
//! @return ... interpolation type id
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


//! @brief sets the integration type for a specific interpolation type
//! @param rInterpolationTypeId ... interpolation type id
//! @param rIntegrationType ... integration type string
void NuTo::Structure::InterpolationTypeSetIntegrationType(int rInterpolationTypeId, const std::string& rIntegrationType, const std::string& rIpDataType)
{
    InterpolationTypeSetIntegrationType(rInterpolationTypeId, GetPtrIntegrationType(rIntegrationType), IpData::IpDataTypeToEnum(rIpDataType));
}

//! @brief sets the integration type for a specific interpolation type
//! @param rInterpolationTypeId ... interpolation type id
//! @param rIntegrationType ... integration type enum
void NuTo::Structure::InterpolationTypeSetIntegrationType(int rInterpolationTypeId, eIntegrationType rIntegrationType, IpData::eIpDataType rIpDataType)
{
    InterpolationTypeSetIntegrationType(rInterpolationTypeId, GetPtrIntegrationType(rIntegrationType), rIpDataType);
}

//! @brief sets the integration type for a specific interpolation type
//! @param rInterpolationTypeId ... interpolation type id
//! @param rIntegrationType ... integration type pointer
//! @param rIpDataType ... ip data type enum
void NuTo::Structure::InterpolationTypeSetIntegrationType(int rInterpolationTypeId, IntegrationTypeBase* rIntegrationType, IpData::eIpDataType rIpDataType)
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
        element->SetIntegrationType(rIntegrationType, rIpDataType);
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
    const IntegrationTypeBase* integrationType = this->GetPtrIntegrationType(integrationTypeEnum);

    interpolationType->UpdateIntegrationType(*integrationType);
    if (mVerboseLevel > 2)
        mLogger << "[NuTo::Structure::InterpolationTypeAdd] Updated IntegrationType to " << integrationType->GetStrIdentifier() << ".\n";

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
        element->SetIntegrationType(integrationType, element->GetIpDataType(0));
    }

    GroupDelete(elementGroupId);
    UpdateDofStatus();
    SetShowTime(showTime);

}
