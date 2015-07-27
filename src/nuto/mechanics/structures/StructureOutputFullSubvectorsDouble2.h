#ifndef STRUCTUREOUTBUTFULLSUBVECTORSDOUBLE2
#define STRUCTUREOUTBUTFULLSUBVECTORSDOUBLE2

#include <nuto/math/FullVector.h>
#include "nuto/mechanics/structures/StructureOutputBase.h"


namespace NuTo
{
//! @author Volker Hirthammer
//! @date July 08, 2015
//! @brief ...
class StructureOutputFullSubvectorsDouble2: public StructureOutputBase
{
    friend class Structure;
public:
    StructureOutputFullSubvectorsDouble2(FullVector<double,Eigen::Dynamic>& rVectorJ,
                                         FullVector<double,Eigen::Dynamic>& rVectorK)
        : StructureOutputBase(),
          mVectorJ(rVectorJ),
          mVectorK(rVectorK)
    {}

    virtual ~StructureOutputFullSubvectorsDouble2(){}

    virtual FullVector<double,Eigen::Dynamic>& GetFullVectorDouble(StructureEnum::eSubVector rSubvectorEnum) override
    {
        if(rSubvectorEnum == StructureEnum::eSubVector::J)
        {
            return mVectorJ;
        }
        else
        {
            return mVectorK;
        }
    }

    virtual int GetNumSubvectors() const override
    {
        return 2;
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

    FullVector<double,Eigen::Dynamic>& mVectorJ;
    FullVector<double,Eigen::Dynamic>& mVectorK;
};


} // namespace NuTo



#endif // STRUCTUREOUTBUTFULLSUBVECTORSDOUBLE2

