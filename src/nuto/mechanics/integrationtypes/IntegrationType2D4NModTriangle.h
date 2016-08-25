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



#include "nuto/mechanics/integrationtypes/IntegrationType2DMod.h"

namespace NuTo
{
//! @author Daniel Arnold, ISM
//! @date February 2011
//! @brief ... integration types in 2D with four nodes and variable number of integration points
class IntegrationType2D4NModTriangle : public IntegrationType2DMod
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:
    //! @brief constructor
    IntegrationType2D4NModTriangle(const std::string rName);
    
	//! @brief creates new integration-cells/order/area
	//! @author Daniel Arnold, ISM
	//! @date May 2012
	//! @param rArea (Input) polygonal surface of integration area
	//! @param rOrder (Input) integration order (or number of integration points)
	void AddIntegrationPoints(std::vector< std::vector<double> > & rArea, const unsigned short rOrder);

	//! @brief calculates the shape functions for the triangular area
	//! @param rLocalCoordinates local coordinates of the integration point
	//! @param shape functions for all the nodes
    void CalculateShapeFunctions(const double rNaturalCoordinates[2], std::vector<double>& rShapeFunctions)const;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION

protected:

};
}

