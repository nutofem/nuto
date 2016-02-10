// $Id$

#include "nuto/mechanics/integrationtypes/IntegrationPointBase.h"

#include "nuto/mechanics/MechanicsException.h"

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
	switch (rBoundingBox.size()){
	case 6:
		mNumVisualizationPoints=3;
		mVisualizationCellType=NuTo::CellBase::TRIANGLE;
		mVisualizationPointLocalCoordinates = rBoundingBox;
		mVisualizationCellsIncidence=std::vector<unsigned int>({0,1,2});
		break;
	case 8:
		mNumVisualizationPoints=4;
		mVisualizationCellType=NuTo::CellBase::QUAD;
		mVisualizationPointLocalCoordinates = rBoundingBox;
		mVisualizationCellsIncidence=std::vector<unsigned int>({0,1,2,3});
		break;
	default:
		throw NuTo::MechanicsException("[" + std::string(__PRETTY_FUNCTION__) + "] Invalid size of bounding box's coordinates for this integration point!");
	}
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

#ifdef ENABLE_SERIALIZATION
template<class Archive> void NuTo::IntegrationPointBase::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize NodeBase" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_NVP(mWeight);
    ar & boost::serialization::make_array(mCoords.data(), mCoords.size());
#ifdef ENABLE_VISUALIZE
    ar & BOOST_SERIALIZATION_NVP(mNumVisualizationPoints);
    ar & BOOST_SERIALIZATION_NVP(mVisualizationCellType);
    ar & boost::serialization::make_array(mVisualizationPointLocalCoordinates.data(), mVisualizationPointLocalCoordinates.size());
    ar & boost::serialization::make_array(mVisualizationCellsIncidence.data(), mVisualizationCellsIncidence.size());
#endif // ENABLE_VISUALIZE
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize NodeBase" << std::endl;
#endif
}

#endif // ENABLE_SERIALIZATION

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
