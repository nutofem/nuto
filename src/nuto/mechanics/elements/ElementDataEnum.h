// $ld: $ 
#ifndef ELEMENTDATAENUM_H_
#define ELEMENTDATAENUM_H_

namespace NuTo
{
namespace ElementData
{
enum eElementDataType
{
    NOELEMENTDATA=0,
    CONSTITUTIVELAWIP,                      //!< constitutive law and integration points
    CONSTITUTIVELAWIPNONLOCAL               //!< constitutive law, integration points and nonlocal data
};
}
}
#endif /* ELEMENTDATAENUM_H_ */
