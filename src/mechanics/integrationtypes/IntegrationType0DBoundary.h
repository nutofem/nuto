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

#include "mechanics/integrationtypes/IntegrationTypeBase.h"

namespace NuTo
{
//! @author Thomas Titscher
//! @date Jan 2015
//! @brief ... integration types in 0D, more like a dummy integration type
class IntegrationType0DBoundary : public IntegrationTypeBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:
    //! @brief constructor
    IntegrationType0DBoundary() = default;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION

    //! @brief returns the dimension of the integration type
    //! @return dimension
    int GetCoordinateDimension()const
    {
        return 0;
    }

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
    double GetIntegrationPointWeight(int rIpNum) const  override;

    //! @brief returns an enum with the type of the integration type
    //! @return enum type
    eIntegrationType GetEnumType() const override
    {
        return eIntegrationType::IntegrationType0DBoundary;
    }

    //! @brief ... check compatibility between element type and integration type
    //! @param rElementType ... element type (enum is defined in ElementBase, but forward declaration of enums not yet possible->int)
    //! @return ... <B>true</B> if the element is compatible with the constitutive relationship, <B>false</B> otherwise.
    bool CheckElementCompatibility(NuTo::Element::eElementType rElementType) const;

#ifdef ENABLE_VISUALIZE
    void GetVisualizationCells(
        unsigned int& NumVisualizationPoints,
        std::vector<double>& VisualizationPointLocalCoordinates,
        unsigned int& NumVisualizationCells,
        std::vector<NuTo::eCellTypes>& VisualizationCellType,
        std::vector<unsigned int>& VisualizationCellsIncidence,
        std::vector<unsigned int>& VisualizationCellsIP) const;
#endif // ENABLE_VISUALIZE
};
}
