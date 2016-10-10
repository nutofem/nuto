#pragma once

namespace NuTo
{

class StructureOutputBlockMatrix;
class StructureOutputBlockVector;


//! @author Volker Hirthammer
//! @date June 25, 2015
//! @brief ...
class StructureOutputBase
{
public:

    StructureOutputBase();

#ifndef SWIG

    virtual ~StructureOutputBase();

    virtual StructureOutputBlockMatrix& AsStructureOutputBlockMatrix();

    virtual StructureOutputBlockVector& AsStructureOutputBlockVector();

    virtual void SetSymmetry(bool rSymmetric);

    virtual bool IsSymmetric()const;

    virtual void SetZero();

#endif
};

} // namespace NuTo
