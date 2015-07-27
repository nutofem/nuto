#ifndef STRUCTUREOUTPUTFULLVECTORDOUBLE_H
#define STRUCTUREOUTPUTFULLVECTORDOUBLE_H

#include <nuto/math/FullVector.h>
#include "nuto/mechanics/structures/StructureOutputBase.h"


namespace NuTo
{
//! @author Volker Hirthammer
//! @date June 29, 2015
//! @brief ...
class StructureOutputFullVectorDouble: public StructureOutputBase
{
    friend class Structure;
public:
    StructureOutputFullVectorDouble(FullVector<double,Eigen::Dynamic>& rVector)
        : StructureOutputBase(),
          mVector(rVector)
    {}

    virtual ~StructureOutputFullVectorDouble(){}

    virtual FullVector<double,Eigen::Dynamic>& GetFullVectorDouble() override
    {
        return mVector;
    }

    virtual int GetNumSubvectors() const override
    {
        return 1;
    }

    virtual bool GetConstant() const override
    {
        return mConstant;
    }

    virtual void SetConstant(bool rConstant) override
    {
        mConstant = rConstant;
    }

private:
    bool mConstant  = false;

    FullVector<double,Eigen::Dynamic>& mVector;
};


} // namespace NuTo

#endif // STRUCTUREOUTPUTFULLVECTORDOUBLE_H

