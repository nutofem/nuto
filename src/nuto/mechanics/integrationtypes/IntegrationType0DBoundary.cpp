#include "nuto/mechanics/elements/ElementEnum.h"
#include "nuto/mechanics/integrationtypes/IntegrationType0DBoundary.h"
#include <assert.h>


NuTo::IntegrationType0DBoundary::IntegrationType0DBoundary() {}


void NuTo::IntegrationType0DBoundary::GetLocalIntegrationPointCoordinates1D(int, double&) const
{
    throw MechanicsException(__PRETTY_FUNCTION__, "Ip number out of range.");
}


unsigned int NuTo::IntegrationType0DBoundary::GetNumIntegrationPoints() const
{
    return 1;
}


double NuTo::IntegrationType0DBoundary::GetIntegrationPointWeight(int rIpNum) const
{
    switch (rIpNum)
    {
    case 0 :
        return 1;
    default:
        throw MechanicsException("Ip number out of range.");
    }
}


std::string NuTo::IntegrationType0DBoundary::GetStrIdentifier() const
{
    return GetStrIdentifierStatic();
}


std::string NuTo::IntegrationType0DBoundary::GetStrIdentifierStatic()
{
    return std::string("0DBOUNDARY");
}

#ifdef ENABLE_VISUALIZE
void NuTo::IntegrationType0DBoundary::GetVisualizationCells(
    unsigned int& NumVisualizationPoints,
    std::vector<double>& VisualizationPointLocalCoordinates,
    unsigned int& NumVisualizationCells,
    std::vector<NuTo::eCellTypes>& VisualizationCellType,
    std::vector<unsigned int>& VisualizationCellsIncidence,
    std::vector<unsigned int>& VisualizationCellsIP) const
{
    // no visualisation since its a 0D element
    NumVisualizationPoints = 0;
    NumVisualizationCells = 0;
}

#endif // ENABLE_VISUALIZE

bool NuTo::IntegrationType0DBoundary::CheckElementCompatibility(
        NuTo::Element::eElementType rElementType) const
{
    switch (rElementType)
    {
    case Element::eElementType::CONTINUUMBOUNDARYELEMENT:
    case Element::eElementType::CONTINUUMBOUNDARYELEMENTCONSTRAINEDCONTROLNODE:
        return true;

    default:
        return false;

    }
}

#ifdef ENABLE_SERIALIZATION
    template<class Archive>
    void NuTo::IntegrationType0DBoundary::serialize(Archive & ar, const unsigned int version)
    {
#ifdef DEBUG_SERIALIZATION
        std::cout << "start serialize IntegrationType0DBoundary" << std::endl;
#endif
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NuTo::IntegrationTypeBase);
#ifdef DEBUG_SERIALIZATION
        std::cout << "finish serialize IntegrationType0DBoundary" << std::endl;
#endif
    }
#endif // ENABLE_SERIALIZATION

