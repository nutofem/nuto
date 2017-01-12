// $Id$
#pragma once

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif //ENABLE_SERIALIZATION

#include "mechanics/integrationtypes/IntegrationType1D.h"

namespace NuTo
{
//! @author Jörg F. Unger, BAM
//! @date November 2013
//! @brief ... integration types in 1D with two nodes Lobatto integration and 3 integration points
class IntegrationType1D2NLobatto4Ip : public IntegrationType1D
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:
    //! @brief constructor
    IntegrationType1D2NLobatto4Ip();

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize IntegrationType1D2NLobatto4Ip" << std::endl;
#endif
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(IntegrationType1D);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize IntegrationType1D2NLobatto4Ip" << std::endl;
#endif
    }
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
        return eIntegrationType::IntegrationType1D2NLobatto4Ip;
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
private:
    //! @brief ... integration points coordinates
    double iPts[4];

    //! @brief ... weights for the integration
    double weights[4];
};
} // namespace

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::IntegrationType1D2NLobatto4Ip)
#endif

