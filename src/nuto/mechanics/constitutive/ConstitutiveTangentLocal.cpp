// $Id: $

#include "nuto/mechanics/constitutive/ConstitutiveTangentLocal.h"

//! @brief reinterpret as ConstitutiveTangentDynamic, otherwise throw an exception
namespace NuTo
{

//tangent 1x1
template <>
NuTo::ConstitutiveTangentLocal<1,1>& NuTo::ConstitutiveTangentLocal<1,1>::AsConstitutiveTangentLocal_1x1()
{
	return *this;
}

//! @brief reinterpret as ConstitutiveTangentDynamic, otherwise throw an exception
template <>
NuTo::ConstitutiveTangentLocal<2,1>& NuTo::ConstitutiveTangentLocal<1,1>::AsConstitutiveTangentLocal_2x1()
{
	throw MechanicsException("[NuTo::ConstitutiveTangentLocal<1,1>::AsConstitutiveTangentLocal_2x1] wrong return type.");
}

//! @brief reinterpret as ConstitutiveTangentDynamic, otherwise throw an exception
template <>
NuTo::ConstitutiveTangentLocal<2,2>& NuTo::ConstitutiveTangentLocal<1,1>::AsConstitutiveTangentLocal_2x2()
{
	throw MechanicsException("[NuTo::ConstitutiveTangentLocal<1,1>::AsConstitutiveTangentLocal_2x2] wrong return type.");
}

//! @brief reinterpret as ConstitutiveTangentDynamic, otherwise throw an exception
template <>
NuTo::ConstitutiveTangentLocal<3,3>& NuTo::ConstitutiveTangentLocal<1,1>::AsConstitutiveTangentLocal_3x3()
{
	throw MechanicsException("[NuTo::ConstitutiveTangentLocal<1,1>::AsConstitutiveTangentLocal_3x3] wrong return type.");
}

//! @brief reinterpret as ConstitutiveTangentDynamic, otherwise throw an exception
template <>
NuTo::ConstitutiveTangentLocal<6,1>& NuTo::ConstitutiveTangentLocal<1,1>::AsConstitutiveTangentLocal_6x1()
{
	throw MechanicsException("[NuTo::ConstitutiveTangentLocal<1,1>::AsConstitutiveTangentLocal_6x1] wrong return type.");
}

//! @brief reinterpret as ConstitutiveTangentDynamic, otherwise throw an exception
template <>
NuTo::ConstitutiveTangentLocal<6,6>& NuTo::ConstitutiveTangentLocal<1,1>::AsConstitutiveTangentLocal_6x6()
{
	throw MechanicsException("[NuTo::ConstitutiveTangentLocal<1,1>::AsConstitutiveTangentLocal_6x6] wrong return type.");
}

//tangent 2x1
template <>
NuTo::ConstitutiveTangentLocal<1,1>& NuTo::ConstitutiveTangentLocal<2,1>::AsConstitutiveTangentLocal_1x1()
{
	throw MechanicsException("[NuTo::ConstitutiveTangentLocal<2,1>::AsConstitutiveTangentLocal_2x1] wrong return type.");
}

//! @brief reinterpret as ConstitutiveTangentDynamic, otherwise throw an exception
template <>
NuTo::ConstitutiveTangentLocal<2,1>& NuTo::ConstitutiveTangentLocal<2,1>::AsConstitutiveTangentLocal_2x1()
{
	return *this;
}

//! @brief reinterpret as ConstitutiveTangentDynamic, otherwise throw an exception
template <>
NuTo::ConstitutiveTangentLocal<2,2>& NuTo::ConstitutiveTangentLocal<2,1>::AsConstitutiveTangentLocal_2x2()
{
	throw MechanicsException("[NuTo::ConstitutiveTangentLocal<2,1>::AsConstitutiveTangentLocal_2x2] wrong return type.");
}

//! @brief reinterpret as ConstitutiveTangentDynamic, otherwise throw an exception
template <>
NuTo::ConstitutiveTangentLocal<3,3>& NuTo::ConstitutiveTangentLocal<2,1>::AsConstitutiveTangentLocal_3x3()
{
	throw MechanicsException("[NuTo::ConstitutiveTangentLocal<2,1>::AsConstitutiveTangentLocal_3x3] wrong return type.");
}

//! @brief reinterpret as ConstitutiveTangentDynamic, otherwise throw an exception
template <>
NuTo::ConstitutiveTangentLocal<6,1>& NuTo::ConstitutiveTangentLocal<2,1>::AsConstitutiveTangentLocal_6x1()
{
	throw MechanicsException("[NuTo::ConstitutiveTangentLocal<2,1>::AsConstitutiveTangentLocal_6x1] wrong return type.");
}

//! @brief reinterpret as ConstitutiveTangentDynamic, otherwise throw an exception
template <>
NuTo::ConstitutiveTangentLocal<6,6>& NuTo::ConstitutiveTangentLocal<2,1>::AsConstitutiveTangentLocal_6x6()
{
	throw MechanicsException("[NuTo::ConstitutiveTangentLocal<2,1>::AsConstitutiveTangentLocal_6x6] wrong return type.");
}

//tangent 2x2
//! @brief reinterpret as ConstitutiveTangentDynamic, otherwise throw an exception
template <>
NuTo::ConstitutiveTangentLocal<1,1>& NuTo::ConstitutiveTangentLocal<2,2>::AsConstitutiveTangentLocal_1x1()
{
	throw MechanicsException("[NuTo::ConstitutiveTangentLocal<2,2>::AsConstitutiveTangentLocal_1x1] wrong return type.");
}

//! @brief reinterpret as ConstitutiveTangentDynamic, otherwise throw an exception
template <>
NuTo::ConstitutiveTangentLocal<2,1>& NuTo::ConstitutiveTangentLocal<2,2>::AsConstitutiveTangentLocal_2x1()
{
	throw MechanicsException("[NuTo::ConstitutiveTangentLocal<2,2>::AsConstitutiveTangentLocal_2x1] wrong return type.");
}

//! @brief reinterpret as ConstitutiveTangentDynamic, otherwise throw an exception
template <>
NuTo::ConstitutiveTangentLocal<2,2>& NuTo::ConstitutiveTangentLocal<2,2>::AsConstitutiveTangentLocal_2x2()
{
	return *this;
}

//! @brief reinterpret as ConstitutiveTangentDynamic, otherwise throw an exception
template <>
NuTo::ConstitutiveTangentLocal<3,3>& NuTo::ConstitutiveTangentLocal<2,2>::AsConstitutiveTangentLocal_3x3()
{
	throw MechanicsException("[NuTo::ConstitutiveTangentLocal<2,2>::AsConstitutiveTangentLocal_3x3] wrong return type.");
}

//! @brief reinterpret as ConstitutiveTangentDynamic, otherwise throw an exception
template <>
NuTo::ConstitutiveTangentLocal<6,1>& NuTo::ConstitutiveTangentLocal<2,2>::AsConstitutiveTangentLocal_6x1()
{
	throw MechanicsException("[NuTo::ConstitutiveTangentLocal<2,2>::AsConstitutiveTangentLocal_1x1] wrong return type.");
}

//! @brief reinterpret as ConstitutiveTangentDynamic, otherwise throw an exception
template <>
NuTo::ConstitutiveTangentLocal<6,6>& NuTo::ConstitutiveTangentLocal<2,2>::AsConstitutiveTangentLocal_6x6()
{
	throw MechanicsException("[NuTo::ConstitutiveTangentLocal<2,2>::AsConstitutiveTangentLocal_6x6] wrong return type.");
}

//tangent 3x3
//! @brief reinterpret as ConstitutiveTangentDynamic, otherwise throw an exception
template <>
NuTo::ConstitutiveTangentLocal<1,1>& NuTo::ConstitutiveTangentLocal<3,3>::AsConstitutiveTangentLocal_1x1()
{
	throw MechanicsException("[NuTo::ConstitutiveTangentLocal<3,3>::AsConstitutiveTangentLocal_1x1] wrong return type.");
}

//! @brief reinterpret as ConstitutiveTangentDynamic, otherwise throw an exception
template <>
NuTo::ConstitutiveTangentLocal<2,1>& NuTo::ConstitutiveTangentLocal<3,3>::AsConstitutiveTangentLocal_2x1()
{
	throw MechanicsException("[NuTo::ConstitutiveTangentLocal<3,3>::AsConstitutiveTangentLocal_2x2] wrong return type.");
}

//! @brief reinterpret as ConstitutiveTangentDynamic, otherwise throw an exception
template <>
NuTo::ConstitutiveTangentLocal<2,2>& NuTo::ConstitutiveTangentLocal<3,3>::AsConstitutiveTangentLocal_2x2()
{
	throw MechanicsException("[NuTo::ConstitutiveTangentLocal<3,3>::AsConstitutiveTangentLocal_2x2] wrong return type.");
}

//! @brief reinterpret as ConstitutiveTangentDynamic, otherwise throw an exception
template <>
NuTo::ConstitutiveTangentLocal<3,3>& NuTo::ConstitutiveTangentLocal<3,3>::AsConstitutiveTangentLocal_3x3()
{
	return *this;
}

//! @brief reinterpret as ConstitutiveTangentDynamic, otherwise throw an exception
template <>
NuTo::ConstitutiveTangentLocal<6,1>& NuTo::ConstitutiveTangentLocal<3,3>::AsConstitutiveTangentLocal_6x1()
{
	throw MechanicsException("[NuTo::ConstitutiveTangentLocal<2,2>::AsConstitutiveTangentLocal_6x1] wrong return type.");
}

//! @brief reinterpret as ConstitutiveTangentDynamic, otherwise throw an exception
template <>
NuTo::ConstitutiveTangentLocal<6,6>& NuTo::ConstitutiveTangentLocal<3,3>::AsConstitutiveTangentLocal_6x6()
{
	throw MechanicsException("[NuTo::ConstitutiveTangentLocal<2,2>::AsConstitutiveTangentLocal_6x6] wrong return type.");
}

//tangent 6x1
//! @brief reinterpret as ConstitutiveTangentDynamic, otherwise throw an exception
template <>
NuTo::ConstitutiveTangentLocal<1,1>& NuTo::ConstitutiveTangentLocal<6,1>::AsConstitutiveTangentLocal_1x1()
{
	throw MechanicsException("[NuTo::ConstitutiveTangentLocal<6,1>::AsConstitutiveTangentLocal_1x1] wrong return type.");
}

template <>
NuTo::ConstitutiveTangentLocal<2,1>& NuTo::ConstitutiveTangentLocal<6,1>::AsConstitutiveTangentLocal_2x1()
{
	throw MechanicsException("[NuTo::ConstitutiveTangentLocal<6,1>::AsConstitutiveTangentLocal_2x1] wrong return type.");
}

//! @brief reinterpret as ConstitutiveTangentDynamic, otherwise throw an exception
template <>
NuTo::ConstitutiveTangentLocal<2,2>& NuTo::ConstitutiveTangentLocal<6,1>::AsConstitutiveTangentLocal_2x2()
{
	throw MechanicsException("[NuTo::ConstitutiveTangentLocal<6,1>::AsConstitutiveTangentLocal_2x2] wrong return type.");
}

//! @brief reinterpret as ConstitutiveTangentDynamic, otherwise throw an exception
template <>
NuTo::ConstitutiveTangentLocal<3,3>& NuTo::ConstitutiveTangentLocal<6,1>::AsConstitutiveTangentLocal_3x3()
{
	throw MechanicsException("[NuTo::ConstitutiveTangentLocal<6,1>::AsConstitutiveTangentLocal_3x3] wrong return type.");
}

//! @brief reinterpret as ConstitutiveTangentDynamic, otherwise throw an exception
template <>
NuTo::ConstitutiveTangentLocal<6,1>& NuTo::ConstitutiveTangentLocal<6,1>::AsConstitutiveTangentLocal_6x1()
{
	return *this;
}

//! @brief reinterpret as ConstitutiveTangentDynamic, otherwise throw an exception
template <>
NuTo::ConstitutiveTangentLocal<6,6>& NuTo::ConstitutiveTangentLocal<6,1>::AsConstitutiveTangentLocal_6x6()
{
	throw MechanicsException("[NuTo::ConstitutiveTangentLocal<6,1>::AsConstitutiveTangentLocal_6x6] wrong return type.");
}

//tangent 6x6
//! @brief reinterpret as ConstitutiveTangentDynamic, otherwise throw an exception
template <>
NuTo::ConstitutiveTangentLocal<1,1>& NuTo::ConstitutiveTangentLocal<6,6>::AsConstitutiveTangentLocal_1x1()
{
	throw MechanicsException("[NuTo::ConstitutiveTangentLocal<6,6>::AsConstitutiveTangentLocal_1x1] wrong return type.");
}

template <>
NuTo::ConstitutiveTangentLocal<2,1>& NuTo::ConstitutiveTangentLocal<6,6>::AsConstitutiveTangentLocal_2x1()
{
	throw MechanicsException("[NuTo::ConstitutiveTangentLocal<6,6>::AsConstitutiveTangentLocal_2x1] wrong return type.");
}

//! @brief reinterpret as ConstitutiveTangentDynamic, otherwise throw an exception
template <>
NuTo::ConstitutiveTangentLocal<2,2>& NuTo::ConstitutiveTangentLocal<6,6>::AsConstitutiveTangentLocal_2x2()
{
	throw MechanicsException("[NuTo::ConstitutiveTangentLocal<6,6>::AsConstitutiveTangentLocal_2x2] wrong return type.");
}

//! @brief reinterpret as ConstitutiveTangentDynamic, otherwise throw an exception
template <>
NuTo::ConstitutiveTangentLocal<3,3>& NuTo::ConstitutiveTangentLocal<6,6>::AsConstitutiveTangentLocal_3x3()
{
	throw MechanicsException("[NuTo::ConstitutiveTangentLocal<6,6>::AsConstitutiveTangentLocal_3x3] wrong return type.");
}

//! @brief reinterpret as ConstitutiveTangentDynamic, otherwise throw an exception
template <>
NuTo::ConstitutiveTangentLocal<6,1>& NuTo::ConstitutiveTangentLocal<6,6>::AsConstitutiveTangentLocal_6x1()
{
	throw MechanicsException("[NuTo::ConstitutiveTangentLocal<6,6>::AsConstitutiveTangentLocal_6x1] wrong return type.");
}

//! @brief reinterpret as ConstitutiveTangentDynamic, otherwise throw an exception
template <>
NuTo::ConstitutiveTangentLocal<6,6>& NuTo::ConstitutiveTangentLocal<6,6>::AsConstitutiveTangentLocal_6x6()
{
	return *this;
}



}
