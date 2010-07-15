#include <sstream>
#include <iostream>
#include "nuto/mechanics/structures/StructureBase.h"
#include "nuto/mechanics/integrationtypes/IntegrationType1D2NConstVariableIp.h"
#include "nuto/mechanics/integrationtypes/IntegrationType1D2NGauss1Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType1D2NGauss2Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType2D3NGauss1Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType2D3NGauss3Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType2D4NGauss1Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType2D4NGauss4Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType3D4NGauss1Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType3D8NGauss1Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType3D8NGauss2x2x2Ip.h"

//! @brief ... Info routine that prints general information about the allocated integration types
//! an integration type is only allocated if required (from created elements)
//! @param ... rVerboseLevel describes how detailed the information is
void NuTo::StructureBase::IntegrationTypeInfo(int rVerboseLevel)const
{
    std::cout<<"number of integration types : " << mIntegrationTypeMap.size() <<std::endl;
    for (boost::ptr_map<std::string,IntegrationTypeBase>::const_iterator it = mIntegrationTypeMap.begin(); it!=mIntegrationTypeMap.end(); it++ )
    {
        it->second->Info(rVerboseLevel);
    }
}

//! @brief ... Returns a pointer to an integration type
//! if the integration type does not exist, the integration type is created
//! @param identIntegrationType Identifier for an integration type
const NuTo::IntegrationTypeBase* NuTo::StructureBase::GetPtrIntegrationType
       (NuTo::IntegrationType::eIntegrationType rEnumIntegrationType)
{
    boost::ptr_map<std::string,IntegrationTypeBase>::iterator it = mIntegrationTypeMap.find(mMappingIntEnum2String[rEnumIntegrationType]);
    if (it!=mIntegrationTypeMap.end())
        return it->second;
    else
    {
        //integration type does not exist, allocate the type
        NuTo::IntegrationTypeBase *ptrIntegrationType;
        switch(rEnumIntegrationType)
        {
        case  NuTo::IntegrationType::IntegrationType1D2NGauss1Ip:
            ptrIntegrationType = new NuTo::IntegrationType1D2NGauss1Ip();
        break;
        case  NuTo::IntegrationType::IntegrationType1D2NGauss2Ip:
            ptrIntegrationType = new NuTo::IntegrationType1D2NGauss2Ip();
        break;
        case  NuTo::IntegrationType::IntegrationType2D3NGauss1Ip:
            ptrIntegrationType = new NuTo::IntegrationType2D3NGauss1Ip();
        break;
        case  NuTo::IntegrationType::IntegrationType2D3NGauss3Ip:
            ptrIntegrationType = new NuTo::IntegrationType2D3NGauss3Ip();
        break;
        case  NuTo::IntegrationType::IntegrationType2D4NGauss1Ip:
            ptrIntegrationType = new NuTo::IntegrationType2D4NGauss1Ip();
        break;
        case  NuTo::IntegrationType::IntegrationType2D4NGauss4Ip:
            ptrIntegrationType = new NuTo::IntegrationType2D4NGauss4Ip();
        break;
        case  NuTo::IntegrationType::IntegrationType3D4NGauss1Ip:
            ptrIntegrationType = new NuTo::IntegrationType3D4NGauss1Ip();
        break;
        case  NuTo::IntegrationType::IntegrationType3D8NGauss1Ip:
            ptrIntegrationType = new NuTo::IntegrationType3D8NGauss1Ip();
        break;
        case  NuTo::IntegrationType::IntegrationType3D8NGauss2x2x2Ip:
            ptrIntegrationType = new NuTo::IntegrationType3D8NGauss2x2x2Ip();
        break;
        default:
            throw MechanicsException("[NuTo::StructureBase::GetPtrIntegrationType] Enum of integration type does not exist.");
        }
        // insert new integration type with mMappingIntEnum2String[rEnumIntegrationType] being the string identifier
        it = mIntegrationTypeMap.insert(mMappingIntEnum2String[rEnumIntegrationType], ptrIntegrationType).first;
        return it->second;
    }
}

//! @brief ... Returns a pointer to an integration type
//! if the integration type does not exist, the integration type is created
//! @param identIntegrationType Identifier for an integration type
const NuTo::IntegrationTypeBase* NuTo::StructureBase::GetPtrIntegrationType(const std::string& rIdentIntegrationType)
{
    std::string IntegrationTypeString;
    std::transform(rIdentIntegrationType.begin(), rIdentIntegrationType.end(), std::back_inserter(IntegrationTypeString), (int(*)(int)) toupper);

    boost::ptr_map<std::string,IntegrationTypeBase>::const_iterator it = mIntegrationTypeMap.find(IntegrationTypeString);
    if (it!=mIntegrationTypeMap.end())
        return it->second;
    else
    {
        //check in the list of standard integration types
        for (int theIntegrationType =0; theIntegrationType<(int)mMappingIntEnum2String.size(); theIntegrationType++)
        {
            if (mMappingIntEnum2String[theIntegrationType]==IntegrationTypeString)
            {
                return GetPtrIntegrationType((NuTo::IntegrationType::eIntegrationType)theIntegrationType);
            }
        }

        //no standard integration type with a fixed number of integration points has been detected
        //try the integration types with variable number of integration points
        if (IntegrationTypeString.substr(0,9)=="1D2NCONST")
        {
            //Allocate an integration type in 1D with variable number of integration points
            if (IntegrationTypeString.substr(IntegrationTypeString.length()-2,2)!="IP")
                throw MechanicsException("[NuTo::StructureBase::GetPtrIntegrationType] The name of an integration type\
 with variable number of Ips is e.g. 1D2NCONST100IP.");
            std::istringstream is(IntegrationTypeString.substr(9,IntegrationTypeString.length()-11));
            int numIp;
            try
            {
                is >> numIp;
            }
            catch(...)
            {
                throw MechanicsException("[NuTo::StructureBase::GetPtrIntegrationType] Error converting number of integration points to integer.");
            }
            it = mIntegrationTypeMap.insert(const_cast<std::string&>(IntegrationTypeString), new NuTo::IntegrationType1D2NConstVariableIp(numIp)).first;
            return it->second;
        }
        else
        {
            throw MechanicsException("[NuTo::StructureBase::GetPtrIntegrationType] Integration type " + IntegrationTypeString +
                    " does not exist. Create the integration type before it is used.");

        }
    }
}
