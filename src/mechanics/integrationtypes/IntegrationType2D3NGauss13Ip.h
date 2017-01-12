// $Id: IntegrationType2D3NGauss13Ip.h 309 2010-09-22 22:21:24Z unger3 $
#pragma once

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION

#include "mechanics/integrationtypes/IntegrationType2D.h"

namespace NuTo
{
//! @author Joerg F. Unger
//! @date June 2010
//! @brief ... integration types in 2D with three nodes Gauss integration and 13 integration points (degree 6)
class IntegrationType2D3NGauss13Ip : public IntegrationType2D
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:
    //! @brief constructor
    IntegrationType2D3NGauss13Ip();

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION

    //! @brief returns the local coordinates of an integration point
    //! @param rIpNum integration point (counting from zero)
    //! @return rCoordinates (result)
    Eigen::VectorXd GetLocalIntegrationPointCoordinates(int rIpNum) const override;

    //! @brief returns the total number of integration points for this integration type
    //! @return number of integration points
    int GetNumIntegrationPoints() const override;

    //! @brief returns the weight of an integration point
    //! @param rIpNum integration point (counting from zero)
    //! @return weight of integration points
    double GetIntegrationPointWeight(int rIpNum) const override;

    //! @brief returns an enum with the type of the integration type
    //! @return enum type
    eIntegrationType GetEnumType() const override
    {
        return eIntegrationType::IntegrationType2D3NGauss13Ip;
    }

#ifdef ENABLE_VISUALIZE
    void GetVisualizationCells(
        unsigned int& NumVisualizationPoints,
        std::vector<double>& VisualizationPointLocalCoordinates,
        unsigned int& NumVisualizationCells,
        std::vector<NuTo::eCellTypes>& VisualizationCellType,
        std::vector<unsigned int>& VisualizationCellsIncidence,
        std::vector<unsigned int>& VisualizationCellsIP) const override;
#endif // ENABLE_VISUALIZE

protected:


};
}
#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::IntegrationType2D3NGauss13Ip)
#endif // ENABLE_SERIALIZATION

