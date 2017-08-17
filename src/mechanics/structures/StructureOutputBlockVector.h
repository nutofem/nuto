#pragma once

#include "mechanics/dofSubMatrixStorage/BlockFullVector.h"
#include "mechanics/structures/StructureOutputBase.h"


namespace NuTo
{
// Foward declarations
class BlockSparseMatrix;


//! @author Volker Hirthammer
//! @date January 25, 2015
//! @brief ...
class StructureOutputBlockVector : public StructureOutputBase
{

public:
    //! @brief constructor
    //! @param rDofStatus: Reference to DofStatus needed for block vectors
    //! @param rAutomaticResize ... resizes the submatrices, if true
    StructureOutputBlockVector(const DofStatus& rDofStatus, bool rAutomaticResize = false);

    //! @brief copy constructor
    StructureOutputBlockVector(const StructureOutputBlockVector& rOther) = default;

#ifndef SWIG
    //! @brief move constructor
    //! @param rOther ... other StructureOutputBlockVector
    StructureOutputBlockVector(StructureOutputBlockVector&& rOther) = default;

    //! @brief destructor
    ~StructureOutputBlockVector();


    // Operator overloads
    // ------------------


    //! @brief copy assignment operator
    StructureOutputBlockVector& operator=(const StructureOutputBlockVector& rOther);


    //! @brief move assignment operator
    //! @param rOther ... other StructureOutputBlockVector
    StructureOutputBlockVector& operator=(StructureOutputBlockVector&& rOther) = default;

    StructureOutputBlockVector& operator+=(const StructureOutputBlockVector& rRhs);
    StructureOutputBlockVector& operator-=(const StructureOutputBlockVector& rRhs);
    StructureOutputBlockVector& operator*=(double rRhs);
    StructureOutputBlockVector& operator/=(double rRhs);


    friend StructureOutputBlockVector operator+(StructureOutputBlockVector rLhs, const StructureOutputBlockVector& rRhs)
    {
        return std::move(rLhs += rRhs);
    }
    friend StructureOutputBlockVector operator-(StructureOutputBlockVector rLhs, const StructureOutputBlockVector& rRhs)
    {
        return std::move(rLhs -= rRhs);
    }
    friend StructureOutputBlockVector operator*(StructureOutputBlockVector rLhs, double rRhs)
    {
        return std::move(rLhs *= rRhs);
    }
    friend StructureOutputBlockVector operator*(double rLhs, StructureOutputBlockVector rRhs)
    {
        return std::move(rRhs *= rLhs);
    }
    friend StructureOutputBlockVector operator/(StructureOutputBlockVector rLhs, double rRhs)
    {
        return std::move(rLhs /= rRhs);
    }
    friend StructureOutputBlockVector operator/(double rLhs, StructureOutputBlockVector rRhs)
    {
        return std::move(rRhs /= rLhs);
    }

    friend std::ostream& operator<<(std::ostream& rOut, const StructureOutputBlockVector& rStructureOutputBlockVector);

    // Member functions
    // ----------------

    void AddElementVector(const NuTo::BlockFullVector<double>& rElementVector,
                          const NuTo::BlockFullVector<int>& rGlobalRowDofNumbers);

    //! @brief resizes every member vector according to numActiveDofs/numDependentDofs
    void Resize(const std::map<Node::eDof, int>& rNumActiveDofsMap,
                const std::map<Node::eDof, int>& rNumDependentDofsMap);

    //! @brief sets each member vector to zero
    void SetZero() override;

    StructureOutputBlockVector& AsStructureOutputBlockVector() override
    {
        return *this;
    }

    Eigen::VectorXd ExportToEigenVector() const
    {
        Eigen::VectorXd j = J.Export();
        Eigen::VectorXd k = K.Export();

        Eigen::VectorXd vec = Eigen::VectorXd::Zero(j.rows() + k.rows());

        vec.head(j.rows()) = j;
        vec.tail(k.rows()) = k;

        return vec;
    }


#endif // SWIG

public:
    NuTo::BlockFullVector<double> J;
    NuTo::BlockFullVector<double> K;
};


} // namespace NuTo
