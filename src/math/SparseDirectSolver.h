// $Id$

#pragma once

#include "base/NuToObject.h"
#include "math/MathException.h"

namespace NuTo
{
class SparseDirectSolver : public NuToObject
{
public:
    //! @brief ... default constructor
    SparseDirectSolver() : NuToObject()
    {
    }

	virtual void Save (const std::string &filename, std::string rType )const
	{
		throw MathException("NuTo::SparseDirectSolver::Save] To be implemented.");
	}

	virtual void Restore (const std::string &filename, std::string rType )
	{
		throw MathException("NuTo::SparseDirectSolver::Restore] To be implemented.");
	}
	virtual void Info()const=0;
protected:
};
} // namespace NuTo
