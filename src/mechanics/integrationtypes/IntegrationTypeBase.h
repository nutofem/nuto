#pragma once

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#include <iostream>
#endif // ENABLE_SERIALIZATION

#include "mechanics/MechanicsException.h"
#include "mechanics/integrationtypes/IntegrationTypeEnum.h"
#include <eigen3/Eigen/Core>
#include  <vector>
#include <memory>


namespace NuTo
{
#ifdef ENABLE_VISUALIZE
class CellBase;
enum class eCellTypes;
#endif //ENABLE_VISUALIZE


class IntegrationPointBase;

namespace Element
{
    enum class eElementType;
}// namespace Element


//! @author Jörg F. Unger, ISM
//! @date November 2009
//! @brief ... standard abstract class for all integration types
class IntegrationTypeBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:
    //! @brief constructor
    IntegrationTypeBase();

    //! @brief ... destructor
    virtual ~IntegrationTypeBase(){};

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize IntegrationTypeBase" << std::endl;
#endif
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize IntegrationTypeBase" << std::endl;
#endif

    }
#endif // ENABLE_SERIALIZATION

    virtual int GetDimension() const = 0;

    //! @brief returns the local coordinates of an integration point
    //! @param rIpNum integration point (counting from zero)
    //! @return rCoordinates (result)
    virtual Eigen::VectorXd GetLocalIntegrationPointCoordinates(int rIpNum) const = 0;


    virtual Eigen::MatrixXd GetNaturalIntegrationPointCoordinates() const
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "Not implemented in base class.");
    }

    //! @brief returns the total number of integration points for this integration type
    //! @return number of integration points
    virtual int GetNumIntegrationPoints() const = 0;

    //! @brief returns the weight of an integration point
    //! @param rIpNum integration point (counting from zero)
    //! @return weight of integration points
    virtual double GetIntegrationPointWeight(int rIpNum) const = 0;

    //! @brief returns a string with the identifier of the integration type
    //! @return identifier
    virtual eIntegrationType GetEnumType() const = 0;

    //! @brief info about the integration type
    //! @param rVerboseLevel determines how detailed the information is
    void Info(int rVerboseLevel)const;

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

    struct CellWithIpId
    {
        std::unique_ptr<CellBase> cell;
        int ipId;
    }

    struct IpCellInfo
    {
        std::vector<CellWithIpId> cells;
        std::vector<Eigen::VectorXd> cellVertices;
    }

    IpCellInfo GetVisualizationCells() const
    {
     unsigned int NumVisualizationPoints;
    std::vector<double> VisualizationPointLocalCoordinates;
    unsigned int NumVisualizationCells;
    std::vector<NuTo::eCellTypes> VisualizationCellType;
    std::vector<unsigned int> VisualizationCellsIncidence;
    std::vector<unsigned int> VisualizationCellsIP;

    GetVisualizationCells(NumVisualizationPoints, VisualizationPointLocalCoordinates, NumVisualizationCells, VisualizationCellType, VisualizationCellsIncidence, VisualizationCellsIP);
    

    }

    virtual void GetVisualizationCells(
        unsigned int& NumVisualizationPoints,
        std::vector<double>& VisualizationPointLocalCoordinates,
        unsigned int& NumVisualizationCells,
        std::vector<NuTo::eCellTypes>& VisualizationCellType,
        std::vector<unsigned int>& VisualizationCellsIncidence,
        std::vector<unsigned int>& VisualizationCellsIP) const = 0;
#endif // ENABLE_VISUALIZE
protected:
};
}//namespace NuTo

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::IntegrationTypeBase)
#endif // ENABLE_SERIALIZATION


