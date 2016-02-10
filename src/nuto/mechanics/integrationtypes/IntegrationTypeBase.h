#ifndef INTEGRATIONTYPEBASE_H
#define INTEGRATIONTYPEBASE_H

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/elements/ElementEnum.h"
#include "nuto/mechanics/MechanicsException.h"

#ifdef ENABLE_VISUALIZE
#include "nuto/visualize/CellBase.h"
#endif // ENABLE_VISUALIZE

#include  <vector>
#include <iostream>
namespace NuTo
{
class IntegrationPointBase;
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
#endif

    //! @brief constructor
    IntegrationTypeBase();

    //! @brief ... destructor
    virtual ~IntegrationTypeBase(){};

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version){}
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

    //! @brief ... check compatibility between element type and integration type
    //! @param rElementType ... element type (enum is defined in ElementBase, but forward declaration of enums not yet possible->int)
    //! @return ... <B>true</B> if the element is compatible with the constitutive relationship, <B>false</B> otherwise.
    virtual bool CheckElementCompatibility(NuTo::Element::eElementType rElementType) const = 0;

    //! @brief creates new integration-cells/order/area
    //! @param rArea (Input) polygonal surface of integration area
    //! @param rOrder (Input) integration order (or number of integration points)
    virtual void AddIntegrationPoints(std::vector< std::vector<double> > & rArea, const unsigned short rOrder);
    
    //! @brief adds a new integration point
    //! @param rIp (Input) integration point
    virtual void AddIntegrationPoint(const IntegrationPointBase & rIp);
    // @brief deletes an integration point
    // @param rIpNum (Input) integration point (counting from zero)
    virtual void DeleteIntegrationPoint(const int rIpNum);

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
#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::IntegrationTypeBase)
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::IntegrationTypeBase)
#endif // ENABLE_SERIALIZATION

#endif //INTEGRATIONTYPEBASE_H
