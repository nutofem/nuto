// $Id$

#include <sstream>
#include <iostream>
#include "nuto/mechanics/structures/StructureBase.h"
#include "nuto/mechanics/integrationtypes/IntegrationType0DBoundary.h"
#include "nuto/mechanics/integrationtypes/IntegrationType1D2NConstVariableIp.h"
#include "nuto/mechanics/integrationtypes/IntegrationType1D2NGauss1Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType1D2NGauss2Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType1D2NBoundaryGauss3Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType1D2NGauss3Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType1D2NGauss4Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType1D2NGauss5Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType1D2NGaussNIp.h"
#include "nuto/mechanics/integrationtypes/IntegrationType1D2NLobatto3Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType1D2NLobatto4Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType1D2NLobatto5Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType1D2NLobatto6Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType2D3NGauss13Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType2D3NGauss16Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType2D3NGauss1Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType2D3NGauss3Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType2D3NGauss4Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType2D3NGauss6Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType2D3NGauss12Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType2D3NGauss12IpDetail.h"
#include "nuto/mechanics/integrationtypes/IntegrationType2D4NConstVariableIp.h"
#include "nuto/mechanics/integrationtypes/IntegrationType2D4NGauss1Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType2D4NGauss4Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType2D4NGauss9Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType2D4NLobatto9Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType2D4NLobatto16Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType2D4NLobatto25Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType2D4NModTriangle.h"
#include "nuto/mechanics/integrationtypes/IntegrationType2D4NModVariableIp.h"
#include "nuto/mechanics/integrationtypes/IntegrationType3D4NGauss1Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType3D4NGauss4Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType3D8NGauss1Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType3D8NGauss2x2x2Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType3D8NLobatto.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeTensorProduct.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeEnum.h"

//! @brief ... Info routine that prints general information about the allocated integration types
//! an integration type is only allocated if required (from created elements)
//! @param ... rVerboseLevel describes how detailed the information is
void NuTo::StructureBase::IntegrationTypeInfo(int rVerboseLevel) const
{
    std::cout << "number of integration types : " << mIntegrationTypeMap.size() << std::endl;
    for (boost::ptr_map<std::string, IntegrationTypeBase>::const_iterator it = mIntegrationTypeMap.begin();
         it != mIntegrationTypeMap.end(); it++)
    {
        it->second->Info(rVerboseLevel);
    }
}

//! @brief ... Returns a pointer to an integration type
//! if the integration type does not exist, the integration type is created
//! @param identIntegrationType Identifier for an integration type
NuTo::IntegrationTypeBase* NuTo::StructureBase::GetPtrIntegrationType(NuTo::eIntegrationType rEnumIntegrationType)
{
    boost::ptr_map<std::string, IntegrationTypeBase>::iterator it =
            mIntegrationTypeMap.find(mMappingIntEnum2String[static_cast<int>(rEnumIntegrationType)]);
    if (it != mIntegrationTypeMap.end())
        return it->second;
    else
    {
        // integration type does not exist, allocate the type
        NuTo::IntegrationTypeBase* ptrIntegrationType;
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
        case NuTo::eIntegrationType::IntegrationType1D2NGauss6Ip:
            ptrIntegrationType = new NuTo::IntegrationType1D2NGaussNIp<6>();
            break;
        case NuTo::eIntegrationType::IntegrationType1D2NGauss7Ip:
            ptrIntegrationType = new NuTo::IntegrationType1D2NGaussNIp<7>();
            break;
        case NuTo::eIntegrationType::IntegrationType1D2NGauss8Ip:
            ptrIntegrationType = new NuTo::IntegrationType1D2NGaussNIp<8>();
            break;
        case NuTo::eIntegrationType::IntegrationType1D2NGauss9Ip:
            ptrIntegrationType = new NuTo::IntegrationType1D2NGaussNIp<9>();
            break;
        case NuTo::eIntegrationType::IntegrationType1D2NGauss10Ip:
            ptrIntegrationType = new NuTo::IntegrationType1D2NGaussNIp<10>();
            break;
        case NuTo::eIntegrationType::IntegrationType1D2NGauss11Ip:
            ptrIntegrationType = new NuTo::IntegrationType1D2NGaussNIp<11>();
            break;
        case NuTo::eIntegrationType::IntegrationType1D2NGauss12Ip:
            ptrIntegrationType = new NuTo::IntegrationType1D2NGaussNIp<12>();
            break;
        case NuTo::eIntegrationType::IntegrationType1D2NGauss13Ip:
            ptrIntegrationType = new NuTo::IntegrationType1D2NGaussNIp<13>();
            break;
        case NuTo::eIntegrationType::IntegrationType1D2NGauss14Ip:
            ptrIntegrationType = new NuTo::IntegrationType1D2NGaussNIp<14>();
            break;
        case NuTo::eIntegrationType::IntegrationType1D2NGauss15Ip:
            ptrIntegrationType = new NuTo::IntegrationType1D2NGaussNIp<15>();
            break;
        case NuTo::eIntegrationType::IntegrationType1D2NGauss16Ip:
            ptrIntegrationType = new NuTo::IntegrationType1D2NGaussNIp<16>();
            break;
        case NuTo::eIntegrationType::IntegrationType1D2NGauss17Ip:
            ptrIntegrationType = new NuTo::IntegrationType1D2NGaussNIp<17>();
            break;
        case NuTo::eIntegrationType::IntegrationType1D2NGauss18Ip:
            ptrIntegrationType = new NuTo::IntegrationType1D2NGaussNIp<18>();
            break;
        case NuTo::eIntegrationType::IntegrationType1D2NGauss19Ip:
            ptrIntegrationType = new NuTo::IntegrationType1D2NGaussNIp<19>();
            break;
        case NuTo::eIntegrationType::IntegrationType1D2NGauss20Ip:
            ptrIntegrationType = new NuTo::IntegrationType1D2NGaussNIp<20>();
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
        case NuTo::eIntegrationType::IntegrationType1D2NLobatto6Ip:
            ptrIntegrationType = new NuTo::IntegrationType1D2NLobatto6Ip();
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
        case NuTo::eIntegrationType::IntegrationType3D8NGauss1Ip:
            ptrIntegrationType = new NuTo::IntegrationType3D8NGauss1Ip();
            break;
        case NuTo::eIntegrationType::IntegrationType3D8NGauss2x2x2Ip:
            ptrIntegrationType = new NuTo::IntegrationType3D8NGauss2x2x2Ip();
            break;
        case NuTo::eIntegrationType::IntegrationType3D8NLobatto3x3x3Ip:
            ptrIntegrationType = new NuTo::IntegrationType3D8NLobatto<3>();
            break;
        case NuTo::eIntegrationType::IntegrationType3D8NLobatto4x4x4Ip:
            ptrIntegrationType = new NuTo::IntegrationType3D8NLobatto<4>();
            break;
        case NuTo::eIntegrationType::IntegrationType3D8NLobatto5x5x5Ip:
            ptrIntegrationType = new NuTo::IntegrationType3D8NLobatto<5>();
            break;
        case NuTo::eIntegrationType::IntegrationTypeTensorProduct2D4N4:
            ptrIntegrationType = new NuTo::IntegrationTypeTensorProduct<2>(4, NuTo::eIntegrationMethod::GAUSS);
            break;
        case NuTo::eIntegrationType::IntegrationTypeTensorProduct2D4N5:
            ptrIntegrationType = new NuTo::IntegrationTypeTensorProduct<2>(5, NuTo::eIntegrationMethod::GAUSS);
            break;
        case NuTo::eIntegrationType::IntegrationTypeTensorProduct2D4N6:
            ptrIntegrationType = new NuTo::IntegrationTypeTensorProduct<2>(6, NuTo::eIntegrationMethod::GAUSS);
            break;
        case NuTo::eIntegrationType::IntegrationTypeTensorProduct2D4N7:
            ptrIntegrationType = new NuTo::IntegrationTypeTensorProduct<2>(7, NuTo::eIntegrationMethod::GAUSS);
            break;
        case NuTo::eIntegrationType::IntegrationTypeTensorProduct2D4N8:
            ptrIntegrationType = new NuTo::IntegrationTypeTensorProduct<2>(8, NuTo::eIntegrationMethod::GAUSS);
            break;
        case NuTo::eIntegrationType::IntegrationTypeTensorProduct2D4N9:
            ptrIntegrationType = new NuTo::IntegrationTypeTensorProduct<2>(9, NuTo::eIntegrationMethod::GAUSS);
            break;
        case NuTo::eIntegrationType::IntegrationTypeTensorProduct2D4N10:
            ptrIntegrationType = new NuTo::IntegrationTypeTensorProduct<2>(10, NuTo::eIntegrationMethod::GAUSS);
            break;

        default:
            throw MechanicsException(
                    "[NuTo::StructureBase::GetPtrIntegrationType] Enum of integration type does not exist.");
        }
        // insert new integration type with mMappingIntEnum2String[rEnumIntegrationType] being the string identifier
        if (mMappingIntEnum2String[static_cast<int>(rEnumIntegrationType)].length() == 0)
        {
            throw MechanicsException("[NuTo::StructureBase::GetPtrIntegrationType] Enum to String mapping of "
                                     "integration type does not exist, check constructor of StructureBase.");
        }
        it = mIntegrationTypeMap
                     .insert(mMappingIntEnum2String[static_cast<int>(rEnumIntegrationType)], ptrIntegrationType)
                     .first;
        return it->second;
    }
}

//! @brief ... Returns a pointer to an integration type
//! if the integration type does not exist, the integration type is created
//! @param identIntegrationType Identifier for an integration type
NuTo::IntegrationTypeBase* NuTo::StructureBase::GetPtrIntegrationType(const std::string& rIdentIntegrationType)
{
    std::string IntegrationTypeString;
    std::transform(rIdentIntegrationType.begin(), rIdentIntegrationType.end(),
                   std::back_inserter(IntegrationTypeString), (int (*)(int))toupper);

    boost::ptr_map<std::string, IntegrationTypeBase>::iterator it = mIntegrationTypeMap.find(IntegrationTypeString);
    if (it != mIntegrationTypeMap.end())
        return it->second;
    else
    {
        // check in the list of standard integration types
        for (int theIntegrationType = 0; theIntegrationType < (int)mMappingIntEnum2String.size(); theIntegrationType++)
        {
            if (mMappingIntEnum2String[theIntegrationType] == IntegrationTypeString)
            {
                return GetPtrIntegrationType((NuTo::eIntegrationType)theIntegrationType);
            }
        }

        // no standard integration type with a fixed number of integration points has been detected
        // try the integration types with variable number of integration points
        if (IntegrationTypeString.substr(0, 9) == "1D2NCONST")
        {
            // Allocate an integration type in 1D with variable number of integration points
            if (IntegrationTypeString.substr(IntegrationTypeString.length() - 2, 2) != "IP")
                throw MechanicsException("[NuTo::StructureBase::GetPtrIntegrationType] The name of an integration type\
 with variable number of Ips is e.g. 1D2NCONST100IP.");
            std::istringstream is(IntegrationTypeString.substr(9, IntegrationTypeString.length() - 11));
            int numIp;
            try
            {
                is >> numIp;
            }
            catch (...)
            {
                throw MechanicsException("[NuTo::StructureBase::GetPtrIntegrationType] Error converting number of "
                                         "integration points to integer.");
            }
            it = mIntegrationTypeMap
                         .insert(const_cast<std::string&>(IntegrationTypeString),
                                 new NuTo::IntegrationType1D2NConstVariableIp(numIp))
                         .first;
            return it->second;
        }
        else if (IntegrationTypeString.substr(0, 9) == "2D4NCONST")
        {
            // Allocate an integration type in 2D with variable number of integration points
            if (IntegrationTypeString.substr(IntegrationTypeString.length() - 2, 2) != "IP")
                throw MechanicsException("[NuTo::StructureBase::GetPtrIntegrationType] The name of an integration type\
 with variable number of Ips is e.g. 2D4NCONST100IP.");
            std::istringstream is(IntegrationTypeString.substr(9, IntegrationTypeString.length() - 11));
            int numIp;
            try
            {
                is >> numIp;
            }
            catch (...)
            {
                throw MechanicsException("[NuTo::StructureBase::GetPtrIntegrationType] Error converting number of "
                                         "integration points to integer.");
            }
            it = mIntegrationTypeMap
                         .insert(const_cast<std::string&>(IntegrationTypeString),
                                 new NuTo::IntegrationType2D4NConstVariableIp(numIp))
                         .first;
            return it->second;
        }
        else if (IntegrationTypeString.substr(0, 7) == "2D4NMOD")
        {
            // Allocate an integration type in 2D with variable number of integration points
            // find first incidence of 'IP'
            std::string::size_type posIP = IntegrationTypeString.find_first_of("IP");
            if (posIP == std::string::npos)
            {
                throw MechanicsException(
                        "[NuTo::StructureBase::GetPtrIntegrationType] The name of an modifiable integration type\
with variable number of Ips is e.g. 2D4NMOD100IP12.");
            }
            int numIp;
            try
            {
                std::istringstream is(IntegrationTypeString.substr(7, posIP - 7));
                is >> numIp;
            }
            catch (...)
            {
                throw MechanicsException("[NuTo::StructureBase::GetPtrIntegrationType] Error converting number of "
                                         "integration points to integer.");
            }
            int Id;
            try
            {
                std::istringstream is(
                        IntegrationTypeString.substr((posIP + 2), IntegrationTypeString.length() - (posIP + 2)));
                is >> Id;
            }
            catch (...)
            {
                throw MechanicsException(
                        "[NuTo::StructureBase::GetPtrIntegrationType] Error converting ID to integer.");
            }
            it = mIntegrationTypeMap
                         .insert(const_cast<std::string&>(IntegrationTypeString),
                                 new NuTo::IntegrationType2D4NModVariableIp(
                                         const_cast<std::string&>(IntegrationTypeString), numIp))
                         .first;
            return it->second;
        }
        else if (IntegrationTypeString.substr(0, 7) == "2D4NTRI")
        {
            // Allocate an integration type in 2D with variable number of integration points
            int Id;
            try
            {
                std::istringstream is(IntegrationTypeString.substr(7, IntegrationTypeString.length() - 7));
                is >> Id;
            }
            catch (...)
            {
                throw MechanicsException(
                        "[NuTo::StructureBase::GetPtrIntegrationType] Error converting ID to integer.");
            }
            it = mIntegrationTypeMap
                         .insert(const_cast<std::string&>(IntegrationTypeString),
                                 new NuTo::IntegrationType2D4NModTriangle(
                                         const_cast<std::string&>(IntegrationTypeString)))
                         .first;
            return it->second;
        }
        else
        {
            throw MechanicsException("[NuTo::StructureBase::GetPtrIntegrationType] Integration type " +
                                     IntegrationTypeString +
                                     " does not exist. Create the integration type before it is used.");
        }
    }
}
