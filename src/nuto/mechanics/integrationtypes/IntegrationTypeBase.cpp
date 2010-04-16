#include "nuto/mechanics/integrationtypes/IntegrationTypeBase.h"
//! @brief constructor
NuTo::IntegrationTypeBase::IntegrationTypeBase()
{}

//! @brief info about the integration type
//! @param rVerboseLevel determines how detailed the information is
void NuTo::IntegrationTypeBase::Info(int rVerboseLevel)const
{
    std::cout << GetStrIdentifier() << std::endl;
    if (rVerboseLevel>2)
    {
        double localCoord[3];
        for (int count=0; count<GetNumIntegrationPoints(); count++)
        {
            std::cout << "    IP " << count << " weight " << GetIntegrationPointWeight(count) << std::endl;
            std::cout << "        coordinates " ;
            switch (GetCoordinateDimension())
            {
            case 1:
                GetLocalIntegrationPointCoordinates1D(count,localCoord[0]);
                std::cout << localCoord[0] << std::endl;
                break;
            case 2:
                GetLocalIntegrationPointCoordinates2D(count,localCoord);
                std::cout << localCoord[0] << localCoord[1] << std::endl;
                break;
            case 3:
                GetLocalIntegrationPointCoordinates3D(count,localCoord);
                std::cout << localCoord[0] << localCoord[1] << localCoord[2] << std::endl;
                break;
            default:
                throw MechanicsException("[NuTo::IntegrationTypeBase::Info] Invalid dimension of integration point coordinates.");
            }
        }
    }
}

//! @brief returns the local coordinates of an integration point
//! @param rIpNum integration point (counting from zero)
//! @param rCoordinates (result)
void NuTo::IntegrationTypeBase::GetLocalIntegrationPointCoordinates1D(int rIpNum, double& rCoordinates)const
{
    throw MechanicsException("[NuTo::IntegrationTypeBase::GetLocalIntegrationPointCoordinates] This integration type does not support 1D coordinates.");
}

//! @brief returns the local coordinates of an integration point
//! @param rIpNum integration point (counting from zero)
//! @param rCoordinates (result)
void NuTo::IntegrationTypeBase::GetLocalIntegrationPointCoordinates2D(int rIpNum, double rCoordinates[2])const
{
    throw MechanicsException("[NuTo::IntegrationTypeBase::GetLocalIntegrationPointCoordinates] This integration type does not support 2D coordinates.");
}

//! @brief returns the local coordinates of an integration point
//! @param rIpNum integration point (counting from zero)
//! @param rCoordinates (result)
void NuTo::IntegrationTypeBase::GetLocalIntegrationPointCoordinates3D(int rIpNum, double rCoordinates[3])const
{
    throw MechanicsException("[NuTo::IntegrationTypeBase::GetLocalIntegrationPointCoordinates] This integration type does not support 3D coordinates.");
}

