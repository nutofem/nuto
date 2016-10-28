#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/integrationtypes/IntegrationTypeBase.h"
#include "nuto/mechanics/integrationtypes/IntegrationPointBase.h"

NuTo::IntegrationTypeBase::IntegrationTypeBase() {}

void NuTo::IntegrationTypeBase::Info(int rVerboseLevel) const
{
    std::cout << GetStrIdentifier() << std::endl;
    if (rVerboseLevel>2)
    {
        double localCoord[3];
        for (unsigned int count = 0; count<GetNumIntegrationPoints(); count++)
        {
            std::cout << "    IP " << count << " weight " << GetIntegrationPointWeight(count) << std::endl;
            std::cout << "        coordinates " ;
            switch (GetCoordinateDimension())
            {
            case 1:
                GetLocalIntegrationPointCoordinates1D(count,localCoord[0]);
                std::cout << "[ " << localCoord[0] << " ]" << std::endl;
                break;
            case 2:
                GetLocalIntegrationPointCoordinates2D(count,localCoord);
                std::cout << "[ " << localCoord[0] << " ; " << localCoord[1] << " ]" << std::endl;
                break;
            case 3:
                GetLocalIntegrationPointCoordinates3D(count,localCoord);
                std::cout << "[ " << localCoord[0] << " ; " << localCoord[1] << " ; " << localCoord[2] << " ]" << std::endl;
                break;
            default:
                throw MechanicsException(__PRETTY_FUNCTION__, "Invalid dimension of integration point coordinates.");
            }
        }
    }
}


void NuTo::IntegrationTypeBase::GetLocalIntegrationPointCoordinates1D(int, double&) const
{
    throw MechanicsException(__PRETTY_FUNCTION__, 
            "Integration type " + GetStrIdentifier() + " does not support 1D coordinates.");
}


void NuTo::IntegrationTypeBase::GetLocalIntegrationPointCoordinates2D(int, double*) const
{
    throw MechanicsException(__PRETTY_FUNCTION__,
            "Integration type " + GetStrIdentifier() + " does not support 1D coordinates.");
}


void NuTo::IntegrationTypeBase::GetLocalIntegrationPointCoordinates3D(int, double*) const
{
    throw MechanicsException(__PRETTY_FUNCTION__,
            "Integration type " + GetStrIdentifier() + " does not support 1D coordinates.");
}


void NuTo::IntegrationTypeBase::AddIntegrationPoints(std::vector<std::vector<double>>&, const unsigned short)
{
    throw MechanicsException(__PRETTY_FUNCTION__, "Cannot add an IP to integration type " + GetStrIdentifier() + ".");
}


void NuTo::IntegrationTypeBase::AddIntegrationPoint(const IntegrationPointBase&)
{
    throw MechanicsException(__PRETTY_FUNCTION__, "Cannot add an IP to integration type " + GetStrIdentifier() + ".");
}


void NuTo::IntegrationTypeBase::DeleteIntegrationPoint(const int)
{
    throw MechanicsException(__PRETTY_FUNCTION__, "Cannot delete an IP to integration type " + GetStrIdentifier() + ".");
}

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::IntegrationTypeBase)
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::IntegrationTypeBase)
#endif // ENABLE_SERIALIZATION
