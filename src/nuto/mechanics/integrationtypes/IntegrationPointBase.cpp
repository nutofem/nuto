// $Id$

#include "nuto/mechanics/integrationtypes/IntegrationPointBase.h"

// constructor
NuTo::IntegrationPointBase::IntegrationPointBase()
{
	mWeight=0;
	mCoords=std::vector<double>(0);
}

// constructor
NuTo::IntegrationPointBase::IntegrationPointBase(const std::vector<double>& rCoords, const double& rWeight, const std::vector<double>& rBoundingBox):
	mWeight(rWeight) , mCoords(rCoords)
{
#ifdef ENABLE_VISUALIZE
    mNumVisualizationPoints=4;
    mVisualizationCellType=NuTo::CellBase::QUAD;
    mVisualizationPointLocalCoordinates = rBoundingBox;
    mVisualizationCellsIncidence=std::vector<unsigned int>({0,1,2,3});
#endif // ENABLE_VISUALIZE
}

// destructor
NuTo::IntegrationPointBase::~IntegrationPointBase()
{
}

// Setter
void NuTo::IntegrationPointBase::SetWeight(const double &rWeight)
{
	mWeight=rWeight;
}

// Getter
void NuTo::IntegrationPointBase::GetWeight(double& rWeight)
{
	rWeight=mWeight;
}

double NuTo::IntegrationPointBase::GetWeight() const
{
	return mWeight;
}

std::vector<double> NuTo::IntegrationPointBase::GetLocalCoords() const
{
	return mCoords;
}

void NuTo::IntegrationPointBase::GetLocalCoords(std::vector<double>& rCoords)
{
	rCoords=mCoords;
}

#ifdef ENABLE_VISUALIZE
void NuTo::IntegrationPointBase::GetVisualizationCell( 	unsigned int& rNumVisualizationPoints,
														NuTo::CellBase::eCellTypes& rVisualizationCellType,
														std::vector<double>& rVisualizationPointLocalCoordinates,
														std::vector<unsigned int>& rVisualizationCellsIncidence ) const
{
	rNumVisualizationPoints=mNumVisualizationPoints;
	rVisualizationCellType=mVisualizationCellType;
	rVisualizationPointLocalCoordinates=mVisualizationPointLocalCoordinates;
	rVisualizationCellsIncidence=mVisualizationCellsIncidence;
}
#endif // ENABLE_VISUALIZE
