// $Id$
#pragma once

#include "mechanics/integrationtypes/IntegrationType3D.h"

namespace NuTo
{
//! @author , ISM
//! @date December 2009
//! @brief ... integration types in 3D with 4 nodes Gauss integration and 1 integration points
class IntegrationType3D4NGauss1Ip : public IntegrationType3D
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:
    //! @brief constructor
    IntegrationType3D4NGauss1Ip();




#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize IntegrationType3D4NGauss1Ip" << std::endl;
#endif
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(IntegrationType3D);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize IntegrationType3D4NGauss1Ip" << std::endl;
#endif
    }
#endif // ENABLE_SERIALIZATION



    //! @brief returns the local coordinates of an integration point
    //! @param rIpNum integration point (counting from zero)
    //! @param rCoordinates (result)
    void GetLocalIntegrationPointCoordinates3D(int rIpNum, double rCoordinates[3])const;


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
        std::vector<NuTo::eCellTypes>& VisualizationCellType,
        std::vector<unsigned int>& VisualizationCellsIncidence,
        std::vector<unsigned int>& VisualizationCellsIP) const;
#endif // ENABLE_VISUALIZE

protected:


};
}// namespace NuTo

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::IntegrationType3D4NGauss1Ip)
#endif
