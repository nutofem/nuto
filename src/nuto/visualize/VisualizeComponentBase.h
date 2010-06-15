// $ld: $ 
#ifndef VISUALIZECOMPONENTBASE_H_
#define VISUALIZECOMPONENTBASE_H_

#include "nuto/visualize/VisualizeBase.h"

namespace NuTo
{
class ElementBase;
//! @author Joerg F. Unger
//! @date Apr 27, 2010
//! @brief a base class to store additional information about the data to be plotted (e.g. the element and ip for nonlocal weights or the stress/strain/displacement component, if not all should be exported to the file
class VisualizeComponentBase
{
public:
    VisualizeComponentBase(){};

    virtual NuTo::VisualizeBase::eVisualizeWhat GetComponentEnum()const=0;

    virtual std::string GetComponentName()const=0;

    virtual int GetElementId()const;

    virtual const NuTo::ElementBase* GetElement()const;

    virtual int GetIp()const;

protected:

};
}
#endif /* VISUALIZECOMPONENTBASE_H_ */
