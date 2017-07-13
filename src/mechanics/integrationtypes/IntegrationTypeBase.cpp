#include <iostream>

#include "mechanics/integrationtypes/IntegrationTypeBase.h"

void NuTo::IntegrationTypeBase::Info(int rVerboseLevel) const
{
    if (rVerboseLevel > 2)
    {
        for (int count = 0; count < GetNumIntegrationPoints(); count++)
        {
            std::cout << "    IP " << count << " weight " << GetIntegrationPointWeight(count) << std::endl;
            std::cout << "        coordinates ";
            std::cout << GetLocalIntegrationPointCoordinates(count);
        }
    }
}

Eigen::MatrixXd NuTo::IntegrationTypeBase::GetNaturalIntegrationPointCoordinates() const
{
    throw Exception(__PRETTY_FUNCTION__, "Not implemented in base class.");
}

NuTo::IntegrationTypeBase::IpCellInfo NuTo::IntegrationTypeBase::GetVisualizationCells() const
{
    unsigned int NumVisualizationPoints;
    std::vector<double> VisualizationPointLocalCoordinates;
    unsigned int NumVisualizationCells;
    std::vector<NuTo::eCellTypes> VisualizationCellType;
    std::vector<unsigned int> VisualizationCellsIncidence;
    std::vector<unsigned int> VisualizationCellsIP;

    GetVisualizationCells(NumVisualizationPoints, VisualizationPointLocalCoordinates, NumVisualizationCells,
                          VisualizationCellType, VisualizationCellsIncidence, VisualizationCellsIP);
    IpCellInfo ipCellInfo;
    if (NumVisualizationCells == 0)
        return ipCellInfo;


    // transform cell vertex coordinates
    const int dim = VisualizationPointLocalCoordinates.size() / NumVisualizationPoints;
    Eigen::MatrixXd visualizationPointNaturalCoordinates =
            Eigen::MatrixXd::Map(VisualizationPointLocalCoordinates.data(), dim, NumVisualizationPoints);

    ipCellInfo.vertices.resize(NumVisualizationPoints);
    for (unsigned i = 0; i < NumVisualizationPoints; ++i)
        ipCellInfo.vertices[i].localCoords = visualizationPointNaturalCoordinates.col(i);

    // transform cells
    ipCellInfo.cells.resize(NumVisualizationCells);
    auto it = VisualizationCellsIncidence.begin();
    for (unsigned i = 0; i < NumVisualizationCells; ++i)
    {
        CellInfo& cellInfo = ipCellInfo.cells[i];
        cellInfo.cellType = VisualizationCellType[i];
        cellInfo.ipId = VisualizationCellsIP[i];

        const int numPoints = Visualize::GetNumPoints(cellInfo.cellType);
        cellInfo.pointIds = std::vector<int>(it, it + numPoints);
        std::advance(it, numPoints);
    }

    return ipCellInfo;
}
