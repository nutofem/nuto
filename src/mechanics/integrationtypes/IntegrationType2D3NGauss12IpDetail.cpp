// $Id$

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif //ENABLE_SERIALIZATION

#ifdef ENABLE_VISUALIZE
#include "visualize/VisualizeEnum.h"
#endif // ENABLE_VISUALIZE

#include "mechanics/integrationtypes/IntegrationType2D3NGauss12IpDetail.h"
#include "math/DelaunayVoronoi.h"
#include <assert.h>

//! @brief constructor
NuTo::IntegrationType2D3NGauss12IpDetail::IntegrationType2D3NGauss12IpDetail()
{
    const double a = 0.063089104491502;
    const double b = 0.249286745170910;
    const double c = 0.310352451033785;
    const double d = 0.053145049844816;

    mIntegrationPointCoordinates.clear();
    mIntegrationPointCoordinates.push_back(Eigen::Vector2d(    a    ,     a    ));
    mIntegrationPointCoordinates.push_back(Eigen::Vector2d(1 - 2 * a,     a    ));
    mIntegrationPointCoordinates.push_back(Eigen::Vector2d(    a    , 1 - 2 * a));
    mIntegrationPointCoordinates.push_back(Eigen::Vector2d(    b    ,     b    ));
    mIntegrationPointCoordinates.push_back(Eigen::Vector2d(1 - 2 * b,     b    ));
    mIntegrationPointCoordinates.push_back(Eigen::Vector2d(    b    , 1 - 2 * b));
    mIntegrationPointCoordinates.push_back(Eigen::Vector2d(    c    ,     d    ));
    mIntegrationPointCoordinates.push_back(Eigen::Vector2d(    d    ,     c    ));
    mIntegrationPointCoordinates.push_back(Eigen::Vector2d(1 - c - d,     c    ));
    mIntegrationPointCoordinates.push_back(Eigen::Vector2d(1 - c - d,     d    ));
    mIntegrationPointCoordinates.push_back(Eigen::Vector2d(    c    , 1 - c - d));
    mIntegrationPointCoordinates.push_back(Eigen::Vector2d(    d    , 1 - c - d));

    const double e = 0.025422453185103;
    const double f = 0.058393137863189;
    const double g = 0.041425537809187;

    mIntegrationPointWeights.clear();
    mIntegrationPointWeights.push_back(e);
    mIntegrationPointWeights.push_back(e);
    mIntegrationPointWeights.push_back(e);
    mIntegrationPointWeights.push_back(f);
    mIntegrationPointWeights.push_back(f);
    mIntegrationPointWeights.push_back(f);
    mIntegrationPointWeights.push_back(g);
    mIntegrationPointWeights.push_back(g);
    mIntegrationPointWeights.push_back(g);
    mIntegrationPointWeights.push_back(g);
    mIntegrationPointWeights.push_back(g);
    mIntegrationPointWeights.push_back(g);

#ifdef ENABLE_VISUALIZE
    DelaunayVoronoi voronoi(mIntegrationPointCoordinates, true);

    std::vector<Eigen::Vector2d> boundaryPoints(3);
    boundaryPoints[0] = Eigen::Vector2d(0.,0.);
    boundaryPoints[1] = Eigen::Vector2d(1.,0.);
    boundaryPoints[2] = Eigen::Vector2d(0.,1.);

    voronoi.SetBoundary(boundaryPoints);



    std::vector<Eigen::Vector2d> points;
    std::vector<std::array<unsigned int,3>> cells;

    voronoi.CalculateVisualizationCellsTriangle(points, cells, mVisualizationCellIPIndices);

    unsigned int numPoints = points.size();
    unsigned int numCells = cells.size();

    mVisualizationPointCoordinates.resize(numPoints*2);
    mVisualizationCellIndices.resize(numCells*3);
    mVisualizationCellTypes.resize(numCells);


    for (unsigned int i = 0; i < numPoints; ++i)
    {
        mVisualizationPointCoordinates[2*i  ] = points[i].x();
        mVisualizationPointCoordinates[2*i+1] = points[i].y();
    }

    for (unsigned int i = 0; i < numCells; ++i)
    {
        mVisualizationCellTypes[i] = NuTo::eCellTypes::TRIANGLE;
        mVisualizationCellIndices[3*i  ] = cells[i][0];
        mVisualizationCellIndices[3*i+1] = cells[i][1];
        mVisualizationCellIndices[3*i+2] = cells[i][2];
    }




#endif // ENABLE_VISUALIZE

}

//! @brief returns the local coordinates of an integration point
//! @param rIpNum integration point (counting from zero)
//! @param rCoordinates (result)
void NuTo::IntegrationType2D3NGauss12IpDetail::GetLocalIntegrationPointCoordinates2D(int rIpNum, double rCoordinates[2]) const
{
    if (rIpNum < 0 || rIpNum >= 12)
        throw MechanicsException("[NuTo::IntegrationType2D3NGauss12Ip::GetLocalIntegrationPointCoordinates] Ip number out of range.");

    rCoordinates[0] = mIntegrationPointCoordinates[rIpNum][0];
    rCoordinates[1] = mIntegrationPointCoordinates[rIpNum][1];
}

//! @brief returns the total number of integration points for this integration type
//! @return number of integration points
int NuTo::IntegrationType2D3NGauss12IpDetail::GetNumIntegrationPoints() const
{
    return 12;
}

//! @brief returns the weight of an integration point
//! @param rIpNum integration point (counting from zero)
//! @return weight of integration points
double NuTo::IntegrationType2D3NGauss12IpDetail::GetIntegrationPointWeight(int rIpNum) const
{
    if (rIpNum < 0 || rIpNum >= 12)
        throw MechanicsException("[NuTo::IntegrationType2D3NGauss12Ip::GetIntegrationPointWeight] Ip number out of range.");

    return mIntegrationPointWeights[rIpNum];
}

//! @brief returns a string with the identifier of the integration type
//! @return identifier
std::string NuTo::IntegrationType2D3NGauss12IpDetail::GetStrIdentifier() const
{
    return GetStrIdentifierStatic();
}

//! @brief returns a string with the identifier of the integration type
//! @return identifier
std::string NuTo::IntegrationType2D3NGauss12IpDetail::GetStrIdentifierStatic()
{
    return std::string("2D3NGAUSS12IPDETAIL");
}

#ifdef ENABLE_VISUALIZE
void NuTo::IntegrationType2D3NGauss12IpDetail::GetVisualizationCells(unsigned int& NumVisualizationPoints, std::vector<double>& VisualizationPointLocalCoordinates, unsigned int& NumVisualizationCells, std::vector<NuTo::eCellTypes>& VisualizationCellType,
        std::vector<unsigned int>& VisualizationCellsIncidence, std::vector<unsigned int>& VisualizationCellsIP) const
{

    NumVisualizationPoints = mVisualizationPointCoordinates.size() / 2;
    NumVisualizationCells = mVisualizationCellIPIndices.size();

    VisualizationPointLocalCoordinates = mVisualizationPointCoordinates;
    VisualizationCellsIncidence = mVisualizationCellIndices;
    VisualizationCellsIP = mVisualizationCellIPIndices;
    VisualizationCellType = mVisualizationCellTypes;

}
#endif // ENABLE_VISUALIZE

#ifdef ENABLE_SERIALIZATION
// serializes the class
template<class Archive>
void NuTo::IntegrationType2D3NGauss12IpDetail::load(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize IntegrationType2D3NGauss12Ip" << std::endl;
#endif
    int size = 0;
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(IntegrationType2D);

    ar & BOOST_SERIALIZATION_NVP(size);
    mIntegrationPointCoordinates.resize(size);
    ar & boost::serialization::make_nvp("mIntegrationPointCoordinates", boost::serialization::make_array(mIntegrationPointCoordinates.data(), size));

    ar & BOOST_SERIALIZATION_NVP(size);
    mIntegrationPointWeights.resize(size);
    ar & boost::serialization::make_array(mIntegrationPointWeights.data(), size);
#ifdef ENABLE_VISUALIZE
    ar & boost::serialization::make_array(mVisualizationPointCoordinates.data(), mVisualizationPointCoordinates.size());
    ar & boost::serialization::make_array(mVisualizationCellIndices.data(), mVisualizationCellIndices.size());
    ar & boost::serialization::make_array(mVisualizationCellIPIndices.data(), mVisualizationCellIPIndices.size());
    ar & boost::serialization::make_array(mVisualizationCellTypes.data(), mVisualizationCellTypes.size());
#endif // ENABLE_VISUALIZE
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize IntegrationType2D3NGauss12Ip" << std::endl;
#endif
}

template<class Archive>
void NuTo::IntegrationType2D3NGauss12IpDetail::save(Archive & ar, const unsigned int version)const
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize IntegrationType2D3NGauss12Ip" << std::endl;
#endif
    int size = 0;
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(IntegrationType2D);

    size = mIntegrationPointCoordinates.size();
    ar & BOOST_SERIALIZATION_NVP(size);
    ar & boost::serialization::make_array(mIntegrationPointCoordinates.data(), size);

    size = mIntegrationPointWeights.size();
    ar & BOOST_SERIALIZATION_NVP(size);
    ar & boost::serialization::make_nvp("mIntegrationPointCoordinates", boost::serialization::make_array(mIntegrationPointWeights.data(), size));

#ifdef ENABLE_VISUALIZE
    ar & boost::serialization::make_array(mVisualizationPointCoordinates.data(), mVisualizationPointCoordinates.size());
    ar & boost::serialization::make_array(mVisualizationCellIndices.data(), mVisualizationCellIndices.size());
    ar & boost::serialization::make_array(mVisualizationCellIPIndices.data(), mVisualizationCellIPIndices.size());
    ar & boost::serialization::make_array(mVisualizationCellTypes.data(), mVisualizationCellTypes.size());
#endif // ENABLE_VISUALIZE
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize IntegrationType2D3NGauss12Ip" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::IntegrationType2D3NGauss12IpDetail)
#endif // ENABLE_SERIALIZATION
