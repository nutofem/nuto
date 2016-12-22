
#include "mechanics/elements/ElementEnum.h"
#include "mechanics/integrationtypes/IntegrationType0DBoundary.h"
#include <assert.h>


// constructor
NuTo::IntegrationType0DBoundary::IntegrationType0DBoundary()
{
}

//! @brief returns the local coordinates of an integration point
//! @param rIpNum integration point (counting from zero)
//! @param rCoordinates (result)
void NuTo::IntegrationType0DBoundary::GetLocalIntegrationPointCoordinates1D(int rIpNum, double& rCoordinates)const
{
//    switch (rIpNum)
//    {
//    case 0 :
//        rCoordinates = -1; // this is the one located on the boundary of the real boundary element
//        break;
//    default:
        throw MechanicsException("[NuTo::IntegrationType0DBoundary::GetLocalIntegrationPointCoordinates] Ip number out of range.");
//    }
}


//! @brief returns the total number of integration points for this integration type
//! @return number of integration points
int NuTo::IntegrationType0DBoundary::GetNumIntegrationPoints()const
{
    return 1;
}

//! @brief returns the weight of an integration point
//! @param rIpNum integration point (counting from zero)
//! @return weight of integration points
double NuTo::IntegrationType0DBoundary::GetIntegrationPointWeight(int rIpNum)const
{
    switch (rIpNum)
    {
    case 0 :
        return 1;
    default:
        throw MechanicsException("[NuTo::IntegrationType0DBoundary::GetIntegrationPointWeight] Ip number out of range.");
    }
}

//! @brief returns a string with the identifier of the integration type
//! @return identifier
std::string NuTo::IntegrationType0DBoundary::GetStrIdentifier()const
{
    return GetStrIdentifierStatic();
}

//! @brief returns a string with the identifier of the integration type
//! @return identifier
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
//    VisualizationPointLocalCoordinates.push_back(-1);
//    VisualizationPointLocalCoordinates.push_back(-0.3873);
//    VisualizationPointLocalCoordinates.push_back(0.3873);
//    VisualizationPointLocalCoordinates.push_back(1);
    NumVisualizationCells = 0;
//    VisualizationCellType.push_back(NuTo::eCellTypes::LINE);
//    VisualizationCellType.push_back(NuTo::eCellTypes::LINE);
//    VisualizationCellType.push_back(NuTo::eCellTypes::LINE);
//    VisualizationCellsIncidence.push_back(0);
//    VisualizationCellsIncidence.push_back(1);
//    VisualizationCellsIncidence.push_back(1);
//    VisualizationCellsIncidence.push_back(2);
//    VisualizationCellsIncidence.push_back(2);
//    VisualizationCellsIncidence.push_back(3);
//    VisualizationCellsIP.push_back(1);
//    VisualizationCellsIP.push_back(2);
//    VisualizationCellsIP.push_back(3);
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


