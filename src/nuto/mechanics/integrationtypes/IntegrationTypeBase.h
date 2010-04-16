#ifndef INTEGRATIONTYPEBASE_H
#define INTEGRATIONTYPEBASE_H

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION

#include "nuto/mechanics/MechanicsException.h"
#ifdef ENABLE_VISUALIZE
#include "nuto/visualize/CellBase.h"
#endif // ENABLE_VISUALIZE
namespace NuTo
{
//! @author Jörg F. Unger, ISM
//! @date November 2009
//! @brief ... standard abstract class for all integration types
class IntegrationTypeBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:
#ifndef SWIG
    enum eIntegrationType
    {
        IntegrationType1D2NGauss1Ip=0,
        IntegrationType1D2NGauss2Ip,
        IntegrationType2D4NGauss1Ip,
        IntegrationType2D4NGauss4Ip,
        IntegrationType3D4NGauss1Ip,
        IntegrationType3D8NGauss1Ip,
        IntegrationType3D8NGauss2x2x2Ip,
        NumIntegrationTypes
    };
#endif

//! @brief constructor
    IntegrationTypeBase();

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {}
#endif // ENABLE_SERIALIZATION

    //! @brief returns the dimension of the integration type
    //! @return dimension = 1, 2 or 3
    virtual int GetCoordinateDimension()const=0;

    //! @brief returns the local coordinates of an integration point
    //! @param rIpNum integration point (counting from zero)
    //! @param rCoordinates (result)
    virtual void GetLocalIntegrationPointCoordinates1D(int rIpNum, double& rCoordinates)const;

    //! @brief returns the local coordinates of an integration point
    //! @param rIpNum integration point (counting from zero)
    //! @param rCoordinates (result)
    virtual void GetLocalIntegrationPointCoordinates2D(int rIpNum, double rCoordinates[2])const;

    //! @brief returns the local coordinates of an integration point
    //! @param rIpNum integration point (counting from zero)
    //! @param rCoordinates (result)
    virtual void GetLocalIntegrationPointCoordinates3D(int rIpNum, double rCoordinates[3])const;

    //! @brief returns the total number of integration points for this integration type
    //! @return number of integration points
    virtual int GetNumIntegrationPoints()const=0;

    //! @brief returns the weight of an integration point
    //! @param rIpNum integration point (counting from zero)
    //! @return weight of integration points
    virtual double GetIntegrationPointWeight(int rIpNum)const=0;

    //! @brief returns a string with the identifier of the integration type
    //! @return identifier
    virtual std::string GetStrIdentifier()const=0;

    //! @brief info about the integration type
    //! @param rVerboseLevel determines how detailed the information is
    void Info(int rVerboseLevel)const;

#ifdef ENABLE_VISUALIZE
    virtual void GetVisualizationCells(
        unsigned int& NumVisualizationPoints,
        std::vector<double>& VisualizationPointLocalCoordinates,
        unsigned int& NumVisualizationCells,
        std::vector<NuTo::CellBase::eCellTypes>& VisualizationCellType,
        std::vector<unsigned int>& VisualizationCellsIncidence,
        std::vector<unsigned int>& VisualizationCellsIP) const = 0;
#endif // ENABLE_VISUALIZE
protected:


};
}//namespace NuTo
#endif //INTEGRATIONTYPEBASE_H
