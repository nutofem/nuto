// $Id: $

#ifndef CONSTITUTIVETANGENTLOCAL_H
#define CONSTITUTIVETANGENTLOCAL_H

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/constitutive/ConstitutiveTangentLocal_Def.h"
#include "nuto/mechanics/MechanicsException.h"

//! @brief ... constructor
template <int TNumRows, int TNumColumns>
NuTo::ConstitutiveTangentLocal<TNumRows,TNumColumns>::ConstitutiveTangentLocal() : ConstitutiveTangentBase::ConstitutiveTangentBase(), FullMatrix<double,TNumRows,TNumColumns>()
{}

//! @brief ... destructor
template <int TNumRows, int TNumColumns>
NuTo::ConstitutiveTangentLocal<TNumRows,TNumColumns>::~ConstitutiveTangentLocal()
{}

//! @brief ... get the number of rows of the tangent matrix
//! @return ... number of rows
template <int TNumRows, int TNumColumns>
unsigned int NuTo::ConstitutiveTangentLocal<TNumRows,TNumColumns>::GetNumberOfRows() const
{
	return TNumRows;
}

//! @brief ... get the number of columns of the tangent matrix
//! @return ... number of columns
template <int TNumRows, int TNumColumns>
unsigned int NuTo::ConstitutiveTangentLocal<TNumRows,TNumColumns>::GetNumberOfColumns() const
{
	return TNumColumns;
}


#endif// CONSTITUTIVETANGENTLOCAL_H

