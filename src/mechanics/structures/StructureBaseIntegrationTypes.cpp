#include <sstream>
#include <iostream>
#include "mechanics/structures/StructureBase.h"
#include "mechanics/integrationtypes/IntegrationType0DBoundary.h"
#include "mechanics/integrationtypes/IntegrationType1D2NGauss1Ip.h"
#include "mechanics/integrationtypes/IntegrationType1D2NGauss2Ip.h"
#include "mechanics/integrationtypes/IntegrationType1D2NBoundaryGauss3Ip.h"
#include "mechanics/integrationtypes/IntegrationType1D2NGauss3Ip.h"
#include "mechanics/integrationtypes/IntegrationType1D2NGauss4Ip.h"
#include "mechanics/integrationtypes/IntegrationType1D2NGauss5Ip.h"
#include "mechanics/integrationtypes/IntegrationType1D2NLobatto3Ip.h"
#include "mechanics/integrationtypes/IntegrationType1D2NLobatto4Ip.h"
#include "mechanics/integrationtypes/IntegrationType1D2NLobatto5Ip.h"
#include "mechanics/integrationtypes/IntegrationType2D3NGauss13Ip.h"
#include "mechanics/integrationtypes/IntegrationType2D3NGauss16Ip.h"
#include "mechanics/integrationtypes/IntegrationType2D3NGauss1Ip.h"
#include "mechanics/integrationtypes/IntegrationType2D3NGauss3Ip.h"
#include "mechanics/integrationtypes/IntegrationType2D3NGauss4Ip.h"
#include "mechanics/integrationtypes/IntegrationType2D3NGauss6Ip.h"
#include "mechanics/integrationtypes/IntegrationType2D3NGauss12Ip.h"
#include "mechanics/integrationtypes/IntegrationType2D3NGauss12IpDetail.h"
#include "mechanics/integrationtypes/IntegrationType2D4NGauss1Ip.h"
#include "mechanics/integrationtypes/IntegrationType2D4NGauss4Ip.h"
#include "mechanics/integrationtypes/IntegrationType2D4NGauss9Ip.h"
#include "mechanics/integrationtypes/IntegrationType2D4NLobatto9Ip.h"
#include "mechanics/integrationtypes/IntegrationType2D4NLobatto16Ip.h"
#include "mechanics/integrationtypes/IntegrationType2D4NLobatto25Ip.h"
#include "mechanics/integrationtypes/IntegrationType3D4NGauss1Ip.h"
#include "mechanics/integrationtypes/IntegrationType3D4NGauss4Ip.h"
#include "mechanics/integrationtypes/IntegrationType3D6NGauss1Ip.h"
#include "mechanics/integrationtypes/IntegrationType3D6NGauss2x3Ip.h"
#include "mechanics/integrationtypes/IntegrationType3D8NGauss1Ip.h"
#include "mechanics/integrationtypes/IntegrationType3D8NGauss2x2x2Ip.h"
#include "mechanics/integrationtypes/IntegrationType3D8NLobatto.h"

//! @brief ... Info routine that prints general information about the allocated integration types
//! an integration type is only allocated if required (from created elements)
//! @param ... rVerboseLevel describes how detailed the information is
void NuTo::StructureBase::IntegrationTypeInfo(int rVerboseLevel) const
{
    std::cout << "number of integration types : " << mIntegrationTypeMap.size() << std::endl;
    for (boost::ptr_map<std::string, IntegrationTypeBase>::const_iterator it = mIntegrationTypeMap.begin(); it != mIntegrationTypeMap.end(); ++it)
    {
        it->second->Info(rVerboseLevel);
    }
}

//! @brief ... Returns a pointer to an integration type
//! if the integration type does not exist, the integration type is created
//! @param identIntegrationType Identifier for an integration type
NuTo::IntegrationTypeBase* NuTo::StructureBase::GetPtrIntegrationType(NuTo::eIntegrationType rEnumIntegrationType)
{
    std::string enumString = IntegrationTypeToString(rEnumIntegrationType);
    boost::ptr_map<std::string, IntegrationTypeBase>::iterator it = mIntegrationTypeMap.find(enumString);
    if (it != mIntegrationTypeMap.end())
        return it->second;
    else
    {
        //integration type does not exist, allocate the type
        NuTo::IntegrationTypeBase *ptrIntegrationType;
        switch (rEnumIntegrationType)
        {
        case NuTo::eIntegrationType::IntegrationType0DBoundary:
            ptrIntegrationType = new NuTo::IntegrationType0DBoundary();
            break;
        case NuTo::eIntegrationType::IntegrationType1D2NGauss1Ip:
            ptrIntegrationType = new NuTo::IntegrationType1D2NGauss1Ip();
            break;
        case NuTo::eIntegrationType::IntegrationType1D2NGauss2Ip:
            ptrIntegrationType = new NuTo::IntegrationType1D2NGauss2Ip();
            break;
        case NuTo::eIntegrationType::IntegrationType1D2NBoundaryGauss3Ip:
            ptrIntegrationType = new NuTo::IntegrationType1D2NBoundaryGauss3Ip();
            break;
        case NuTo::eIntegrationType::IntegrationType1D2NGauss3Ip:
            ptrIntegrationType = new NuTo::IntegrationType1D2NGauss3Ip();
            break;
        case NuTo::eIntegrationType::IntegrationType1D2NGauss4Ip:
            ptrIntegrationType = new NuTo::IntegrationType1D2NGauss4Ip();
            break;
        case NuTo::eIntegrationType::IntegrationType1D2NGauss5Ip:
            ptrIntegrationType = new NuTo::IntegrationType1D2NGauss5Ip();
            break;
        case NuTo::eIntegrationType::IntegrationType1D2NLobatto3Ip:
            ptrIntegrationType = new NuTo::IntegrationType1D2NLobatto3Ip();
            break;
        case NuTo::eIntegrationType::IntegrationType1D2NLobatto4Ip:
            ptrIntegrationType = new NuTo::IntegrationType1D2NLobatto4Ip();
            break;
        case NuTo::eIntegrationType::IntegrationType1D2NLobatto5Ip:
            ptrIntegrationType = new NuTo::IntegrationType1D2NLobatto5Ip();
            break;
        case NuTo::eIntegrationType::IntegrationType2D3NGauss13Ip:
            ptrIntegrationType = new NuTo::IntegrationType2D3NGauss13Ip();
            break;
        case NuTo::eIntegrationType::IntegrationType2D3NGauss16Ip:
            ptrIntegrationType = new NuTo::IntegrationType2D3NGauss16Ip();
            break;
        case NuTo::eIntegrationType::IntegrationType2D3NGauss1Ip:
            ptrIntegrationType = new NuTo::IntegrationType2D3NGauss1Ip();
            break;
        case NuTo::eIntegrationType::IntegrationType2D3NGauss3Ip:
            ptrIntegrationType = new NuTo::IntegrationType2D3NGauss3Ip();
            break;
        case NuTo::eIntegrationType::IntegrationType2D3NGauss4Ip:
            ptrIntegrationType = new NuTo::IntegrationType2D3NGauss4Ip();
            break;
        case NuTo::eIntegrationType::IntegrationType2D3NGauss6Ip:
            ptrIntegrationType = new NuTo::IntegrationType2D3NGauss6Ip();
            break;
        case NuTo::eIntegrationType::IntegrationType2D3NGauss12Ip:
            ptrIntegrationType = new NuTo::IntegrationType2D3NGauss12Ip();
            break;
        case NuTo::eIntegrationType::IntegrationType2D3NGauss12IpDetail:
            ptrIntegrationType = new NuTo::IntegrationType2D3NGauss12IpDetail();
            break;
        case NuTo::eIntegrationType::IntegrationType2D4NGauss1Ip:
            ptrIntegrationType = new NuTo::IntegrationType2D4NGauss1Ip();
            break;
        case NuTo::eIntegrationType::IntegrationType2D4NGauss4Ip:
            ptrIntegrationType = new NuTo::IntegrationType2D4NGauss4Ip();
            break;
        case NuTo::eIntegrationType::IntegrationType2D4NGauss9Ip:
            ptrIntegrationType = new NuTo::IntegrationType2D4NGauss9Ip();
            break;
        case NuTo::eIntegrationType::IntegrationType2D4NLobatto9Ip:
            ptrIntegrationType = new NuTo::IntegrationType2D4NLobatto9Ip();
            break;
        case NuTo::eIntegrationType::IntegrationType2D4NLobatto16Ip:
            ptrIntegrationType = new NuTo::IntegrationType2D4NLobatto16Ip();
            break;
        case NuTo::eIntegrationType::IntegrationType2D4NLobatto25Ip:
            ptrIntegrationType = new NuTo::IntegrationType2D4NLobatto25Ip();
            break;
        case NuTo::eIntegrationType::IntegrationType3D4NGauss1Ip:
            ptrIntegrationType = new NuTo::IntegrationType3D4NGauss1Ip();
            break;
        case NuTo::eIntegrationType::IntegrationType3D4NGauss4Ip:
            ptrIntegrationType = new NuTo::IntegrationType3D4NGauss4Ip();
            break;
        case NuTo::eIntegrationType::IntegrationType3D6NGauss2x3Ip:
            ptrIntegrationType = new NuTo::IntegrationType3D6NGauss2x3Ip();
            break;
        case NuTo::eIntegrationType::IntegrationType3D6NGauss1Ip:
            ptrIntegrationType = new NuTo::IntegrationType3D6NGauss1Ip();
            break;
        case NuTo::eIntegrationType::IntegrationType3D8NGauss1Ip:
            ptrIntegrationType = new NuTo::IntegrationType3D8NGauss1Ip();
            break;
        case NuTo::eIntegrationType::IntegrationType3D8NGauss2x2x2Ip:
            ptrIntegrationType = new NuTo::IntegrationType3D8NGauss2x2x2Ip();
            break;
        case  NuTo::eIntegrationType::IntegrationType3D8NLobatto3x3x3Ip:
            ptrIntegrationType = new NuTo::IntegrationType3D8NLobatto<3>();
            break;
        case  NuTo::eIntegrationType::IntegrationType3D8NLobatto4x4x4Ip:
            ptrIntegrationType = new NuTo::IntegrationType3D8NLobatto<4>();
            break;
        case  NuTo::eIntegrationType::IntegrationType3D8NLobatto5x5x5Ip:
            ptrIntegrationType = new NuTo::IntegrationType3D8NLobatto<5>();
            break;
        default:
            throw MechanicsException("[NuTo::StructureBase::GetPtrIntegrationType] Enum of integration type does not exist.");
        }
        it = mIntegrationTypeMap.insert(enumString, ptrIntegrationType).first;
        return it->second;
    }
}

//! @brief ... Returns a pointer to an integration type
//! if the integration type does not exist, the integration type is created
//! @param identIntegrationType Identifier for an integration type
NuTo::IntegrationTypeBase* NuTo::StructureBase::GetPtrIntegrationType(std::string rIdentIntegrationType)
{
    return GetPtrIntegrationType(IntegrationTypeToEnum(rIdentIntegrationType));
}
