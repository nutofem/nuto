#pragma once

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/integrationtypes/IntegrationType2D.h"
#include <eigen3/Eigen/Dense>

namespace NuTo
{

//! @author Thomas Titscher
//! @date August 2015
//! @brief ... 6th order integration in 2D with 12 points and visualization for every integration point via voronoi cells
class IntegrationType2D3NGauss12IpDetail : public IntegrationType2D
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:
    //! @brief constructor
    IntegrationType2D3NGauss12IpDetail();

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION

    //! @brief returns the local coordinates of an integration point
    //! @param rIpNum integration point (counting from zero)
    //! @param rCoordinates (result)
    void GetLocalIntegrationPointCoordinates2D(int rIpNum, double rCoordinates[2])const;


    //! @brief returns the total number of integration points for this integration type
    //! @return number of integration points
    int GetNumIntegrationPoints()const;

    //! @brief returns the weight of an integration point
    //! @param rIpNum integration point (counting from zero)
    //! @return weight of integration points
    double GetIntegrationPointWeight(int rIpNum)const;

    //! @brief returns a string with the identifier of the integration type
    //! @return identifier
    std::string GetStrIdentifier()const;

    //! @brief returns a string with the identifier of the integration type
    //! @return identifier
    static std::string GetStrIdentifierStatic();

#ifdef ENABLE_VISUALIZE
    void GetVisualizationCells(
        unsigned int& NumVisualizationPoints,
        std::vector<double>& VisualizationPointLocalCoordinates,
        unsigned int& NumVisualizationCells,
        std::vector<NuTo::CellBase::eCellTypes>& VisualizationCellType,
        std::vector<unsigned int>& VisualizationCellsIncidence,
        std::vector<unsigned int>& VisualizationCellsIP) const;
#endif // ENABLE_VISUALIZE

protected:

    std::vector<Eigen::Vector2d> mIntegrationPointCoordinates;
    std::vector<double> mIntegrationPointWeights;
#ifdef ENABLE_VISUALIZE
    std::vector<double> mVisualizationPointCoordinates;
    std::vector<unsigned int> mVisualizationCellIndices;
    std::vector<unsigned int> mVisualizationCellIPIndices;
    std::vector<NuTo::CellBase::eCellTypes> mVisualizationCellTypes;
#endif // ENABLE_VISUALIZE
};
}
#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::IntegrationType2D3NGauss12IpDetail)
#endif // ENABLE_SERIALIZATION

