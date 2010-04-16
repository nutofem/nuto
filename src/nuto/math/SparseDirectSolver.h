// $Id$

#ifndef SPARSE_DIRECT_SOLVER_H
#define SPARSE_DIRECT_SOLVER_H
namespace NuTo
{
class SparseDirectSolver
{
public:
    //! @brief ... default constructor
    SparseDirectSolver() : mVerboseLevel(0)
    {
    }

    //! @brief ... set the verbose level of the solver output
    //! @param rVerboseLevel ... verbose level
    inline void SetVerboseLevel(unsigned int rVerboseLevel)
    {
        this->mVerboseLevel = rVerboseLevel;
    }

    //! @brief ... print informations about the solver
    virtual void Info() = 0;
protected:
    //! @brief ... controls the verbose level of the solver output (0 ... no output)
    unsigned int mVerboseLevel;
};
} // namespace NuTo
#endif // SPARSE_DIRECT_SOLVER_H
