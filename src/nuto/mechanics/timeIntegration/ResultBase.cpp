
#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/timeIntegration/ResultBase.h"


//! @brief constructor
NuTo::ResultBase::ResultBase(const std::string& rIdent) : NuToObject()
{
	mIdent = rIdent;
	mCalculated = false;
}

//! @brief deconstructor
NuTo::ResultBase::~ResultBase()
{
}

void NuTo::ResultBase::SetIdent(const std::string& rIdent)
{
	mIdent = rIdent;
}

std::string NuTo::ResultBase::GetIdent()const
{
	return mIdent;
}

bool NuTo::ResultBase::IsCalculated() const
{
	return mCalculated;
}

void NuTo::ResultBase::SetCalculated(bool rCalculated)
{
	mCalculated = rCalculated;
}
