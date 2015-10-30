

#include <nuto/mechanics/constitutive/ConstitutiveTangentLocal.h>
#include "nuto/mechanics/MechanicsException.h"

namespace NuTo
{

template class ConstitutiveTangentLocal<1,1>;
template class ConstitutiveTangentLocal<2,1>;
template class ConstitutiveTangentLocal<3,1>;
template class ConstitutiveTangentLocal<6,1>;
template class ConstitutiveTangentLocal<Eigen::Dynamic,1>;

template class ConstitutiveTangentLocal<2,2>;
template class ConstitutiveTangentLocal<3,3>;
template class ConstitutiveTangentLocal<6,6>;
template class ConstitutiveTangentLocal<Eigen::Dynamic,Eigen::Dynamic>;


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


//! @brief reinterpret as ConstitutiveTangentDynamic, otherwise throw an exception
template <int TNumRows, int TNumColumns>
NuTo::ConstitutiveTangentLocal<1,1>& NuTo::ConstitutiveTangentLocal<TNumRows,TNumColumns>::AsConstitutiveTangentLocal_1x1()
{
	throw MechanicsException(std::string(__PRETTY_FUNCTION__) + ":\t  wrong return type.");
}


//! @brief reinterpret as ConstitutiveTangentDynamic, otherwise throw an exception
template <int TNumRows, int TNumColumns>
NuTo::ConstitutiveTangentLocal<1,2>& NuTo::ConstitutiveTangentLocal<TNumRows,TNumColumns>::AsConstitutiveTangentLocal_1x2()
{
	throw MechanicsException(std::string(__PRETTY_FUNCTION__) + ":\t  wrong return type.");
}

//! @brief reinterpret as ConstitutiveTangentDynamic, otherwise throw an exception
template <int TNumRows, int TNumColumns>
NuTo::ConstitutiveTangentLocal<2,1>& NuTo::ConstitutiveTangentLocal<TNumRows,TNumColumns>::AsConstitutiveTangentLocal_2x1()
{
	throw MechanicsException(std::string(__PRETTY_FUNCTION__) + ":\t  wrong return type.");
}

//! @brief reinterpret as ConstitutiveTangentDynamic, otherwise throw an exception
template <int TNumRows, int TNumColumns>
NuTo::ConstitutiveTangentLocal<2,2>& NuTo::ConstitutiveTangentLocal<TNumRows,TNumColumns>::AsConstitutiveTangentLocal_2x2()
{
	throw MechanicsException(std::string(__PRETTY_FUNCTION__) + ":\t  wrong return type.");
}

//! @brief reinterpret as ConstitutiveTangentDynamic, otherwise throw an exception
template <int TNumRows, int TNumColumns>
NuTo::ConstitutiveTangentLocal<3,1>& NuTo::ConstitutiveTangentLocal<TNumRows,TNumColumns>::AsConstitutiveTangentLocal_3x1()
{
	throw MechanicsException(std::string(__PRETTY_FUNCTION__) + ":\t  wrong return type.");
}

//! @brief reinterpret as ConstitutiveTangentDynamic, otherwise throw an exception
template <int TNumRows, int TNumColumns>
NuTo::ConstitutiveTangentLocal<3,3>& NuTo::ConstitutiveTangentLocal<TNumRows,TNumColumns>::AsConstitutiveTangentLocal_3x3()
{
	throw MechanicsException(std::string(__PRETTY_FUNCTION__) + ":\t  wrong return type.");
}

//! @brief reinterpret as ConstitutiveTangentDynamic, otherwise throw an exception
template <int TNumRows, int TNumColumns>
NuTo::ConstitutiveTangentLocal<6,1>& NuTo::ConstitutiveTangentLocal<TNumRows,TNumColumns>::AsConstitutiveTangentLocal_6x1()
{
	throw MechanicsException(std::string(__PRETTY_FUNCTION__) + ":\t  wrong return type.");
}

//! @brief reinterpret as ConstitutiveTangentDynamic, otherwise throw an exception
template <int TNumRows, int TNumColumns>
NuTo::ConstitutiveTangentLocal<6,6>& NuTo::ConstitutiveTangentLocal<TNumRows,TNumColumns>::AsConstitutiveTangentLocal_6x6()
{
	throw MechanicsException(std::string(__PRETTY_FUNCTION__) + ":\t  wrong return type.");
}



////////////////////////////////////////////////////////////
// 			TEMPLATE SPECIALIZATION
////////////////////////////////////////////////////////////

//! @brief reinterpret as ConstitutiveTangentDynamic, otherwise throw an exception
template <>
NuTo::ConstitutiveTangentLocal<1,1>& NuTo::ConstitutiveTangentLocal<1,1>::AsConstitutiveTangentLocal_1x1()
{
	return *this;
}

//! @brief reinterpret as ConstitutiveTangentDynamic, otherwise throw an exception
template <>
NuTo::ConstitutiveTangentLocal<1,2>& NuTo::ConstitutiveTangentLocal<1,2>::AsConstitutiveTangentLocal_1x2()
{
	return *this;
}

//! @brief reinterpret as ConstitutiveTangentDynamic, otherwise throw an exception
template <>
NuTo::ConstitutiveTangentLocal<2,1>& NuTo::ConstitutiveTangentLocal<2,1>::AsConstitutiveTangentLocal_2x1()
{
	return *this;
}


//! @brief reinterpret as ConstitutiveTangentDynamic, otherwise throw an exception
template <>
NuTo::ConstitutiveTangentLocal<2,2>& NuTo::ConstitutiveTangentLocal<2,2>::AsConstitutiveTangentLocal_2x2()
{
	return *this;
}

//! @brief reinterpret as ConstitutiveTangentDynamic, otherwise throw an exception
template <>
NuTo::ConstitutiveTangentLocal<3,1>& NuTo::ConstitutiveTangentLocal<3,1>::AsConstitutiveTangentLocal_3x1()
{
    return *this;
}

//! @brief reinterpret as ConstitutiveTangentDynamic, otherwise throw an exception
template <>
NuTo::ConstitutiveTangentLocal<3,3>& NuTo::ConstitutiveTangentLocal<3,3>::AsConstitutiveTangentLocal_3x3()
{
	return *this;
}

//! @brief reinterpret as ConstitutiveTangentDynamic, otherwise throw an exception
template <>
NuTo::ConstitutiveTangentLocal<6,1>& NuTo::ConstitutiveTangentLocal<6,1>::AsConstitutiveTangentLocal_6x1()
{
	return *this;
}

//! @brief reinterpret as ConstitutiveTangentDynamic, otherwise throw an exception
template <>
NuTo::ConstitutiveTangentLocal<6,6>& NuTo::ConstitutiveTangentLocal<6,6>::AsConstitutiveTangentLocal_6x6()
{
	return *this;
}



}// namepsace NuTo
