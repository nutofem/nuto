#pragma once

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#include <boost/serialization/array.hpp>
#include "math/CustomBoostSerializationExtensions.h"
#endif // ENABLE_SERIALIZATION

#include "mechanics/integrationtypes/IntegrationType2D.h"

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
    //! @brief save (serializes) the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void save(Archive & ar, const unsigned int version) const;

    //! @brief load (deserializes) the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void load(Archive & ar, const unsigned int version);

    BOOST_SERIALIZATION_SPLIT_MEMBER()
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

#ifdef ENABLE_VISUALIZE
    IpCellInfo GetVisualizationCells() const override
    {
        return mIpCellInfo;
    }
#endif // ENABLE_VISUALIZE

protected:

    std::vector<Eigen::Vector2d> mIntegrationPointCoordinates;
    std::vector<double> mIntegrationPointWeights;
#ifdef ENABLE_VISUALIZE
    IntegrationTypeBase::IpCellInfo mIpCellInfo;
#endif // ENABLE_VISUALIZE
};
}
#ifdef ENABLE_SERIALIZATION
//namespace boost
//{
//    //! @brief tell boost how to serialize an Eigen::Matrix
//    template<class Archive, typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
//    inline void serialize(
//        Archive & ar,
//        Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> & t,
//        const unsigned int file_version
//    )
//    {
//        ar & boost::serialization::make_array(t.data(), t.size());
//    }
//}
BOOST_CLASS_EXPORT_KEY(NuTo::IntegrationType2D3NGauss12IpDetail)
#endif // ENABLE_SERIALIZATION

