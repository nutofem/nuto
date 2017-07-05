#include "mechanics/elements/ElementEnum.h"
#include "mechanics/integrationtypes/IntegrationType0DBoundary.h"

using namespace NuTo;

Eigen::VectorXd NuTo::IntegrationType0DBoundary::GetLocalIntegrationPointCoordinates(int rIpNum) const
{
    throw MechanicsException(
            "[NuTo::IntegrationType0DBoundary::GetLocalIntegrationPointCoordinates] Ip number out of range.");
}


//! @brief returns the total number of integration points for this integration type
//! @return number of integration points
int NuTo::IntegrationType0DBoundary::GetNumIntegrationPoints() const
{
    return 1;
}

//! @brief returns the weight of an integration point
//! @param rIpNum integration point (counting from zero)
//! @return weight of integration points
double NuTo::IntegrationType0DBoundary::GetIntegrationPointWeight(int rIpNum) const
{
    switch (rIpNum)
    {
    case 0:
        return 1;
    default:
        throw MechanicsException(
                "[NuTo::IntegrationType0DBoundary::GetIntegrationPointWeight] Ip number out of range.");
    }
}


int IntegrationType0DBoundary::GetDimension() const
{
    return 0;
}

#ifdef ENABLE_VISUALIZE
void NuTo::IntegrationType0DBoundary::GetVisualizationCells(unsigned int& NumVisualizationPoints,
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
