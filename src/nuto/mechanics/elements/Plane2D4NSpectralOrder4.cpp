// $Id: Plane2D4NSpectralOrder2.cpp 647 2013-10-31 13:23:00Z unger3 $

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

#include <tuple>

#include <assert.h>

#include <boost/foreach.hpp>

#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/elements/Plane2D4NSpectral.h"
#include "nuto/mechanics/nodes/NodeBase.h"

//! @brief calculates the shape functions
//! @param rLocalCoordinates local coordinates of the integration point
//! @param shape functions for all the nodes
namespace NuTo
{
template<>
void NuTo::Plane2D4NSpectral<4>::CalculateShapeFunctionsField1D(double rNaturalCoordinate, std::vector<double>& rShapeFunctions)const
{
	assert(((int)rShapeFunctions.size())==GetNumNodesField1D());
	double s2 = rNaturalCoordinate*rNaturalCoordinate;
	double s3 = rNaturalCoordinate*s2;
	double s4 = rNaturalCoordinate*s3;
	rShapeFunctions[0] =    +0.375            *rNaturalCoordinate -0.375            *s2 -0.875*s3             + 0.875*s4;
	rShapeFunctions[1] =    -1.336584577695453*rNaturalCoordinate +2.041666666666666*s2 +1.336584577695453*s3 -2.041666666666666*s4;
	rShapeFunctions[2] =  1.                                      -3.333333333333333*s2                       +2.333333333333333*s4;
	rShapeFunctions[3] =    +1.336584577695453*rNaturalCoordinate +2.041666666666666*s2 -1.336584577695453*s3 -2.041666666666666*s4;
	rShapeFunctions[4] =    -0.375            *rNaturalCoordinate -0.375            *s2 + 0.875*s3            +0.875*s4;
}

//! @brief calculates the derivative of the shape functions
//! @param rLocalCoordinates local coordinates of the integration point
//! @param derivative of the shape functions for all the nodes,
//! first all the directions for a single node, and then for the next node
template<>
void NuTo::Plane2D4NSpectral<4>::CalculateDerivativeShapeFunctionsFieldNatural1D(double rNaturalCoordinate, std::vector<double>& rDerivativeShapeFunctions)const
{
	assert(((int)rDerivativeShapeFunctions.size())==GetNumNodesField1D());
	double s2 = rNaturalCoordinate*rNaturalCoordinate;
	double s3 = rNaturalCoordinate*s2;
	rDerivativeShapeFunctions[0] =  0.375             -0.75             *rNaturalCoordinate-2.625               *s2+3.5*s3;
	rDerivativeShapeFunctions[1] = -1.336584577695453 +4.083333333333333*rNaturalCoordinate+4.009753733086359517*s2-8.16666666666666*s3;
	rDerivativeShapeFunctions[2] =                    -6.666666666666666*rNaturalCoordinate                        +9.33333333333333*s3;
	rDerivativeShapeFunctions[3] =  1.336584577695453 +4.083333333333333*rNaturalCoordinate-4.009753733086359517*s2-8.16666666666666*s3;
	rDerivativeShapeFunctions[4] = -0.375             -0.75             *rNaturalCoordinate+2.625               *s2+3.5*s3;
}


//! @brief returns the enum of the standard integration type for this element
template<>
NuTo::IntegrationType::eIntegrationType NuTo::Plane2D4NSpectral<4>::GetStandardIntegrationType()
{
    return NuTo::IntegrationType::IntegrationType2D4NLobatto25Ip;
}

//! @brief returns the enum (type of the element)
//! @return enum
template<>
NuTo::Element::eElementType NuTo::Plane2D4NSpectral<4>::GetEnumType()const
{
    return NuTo::Element::PLANE2D4NSPECTRALORDER4;
}

} //end namespace nuto

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::Plane2D4NSpectral<4>)
#endif // ENABLE_SERIALIZATION
