// $Id$
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
    CONSTITUTIVELAWIPNONLOCAL,              //!< constitutive law, integration points and nonlocal data
    CONSTITUTIVELAWIPCRACK,                  //!< constitutive law, integration points and crack data
    VARIABLECONSTITUTIVELAWIP,              //!< variable constitutive laws and integration points
};
}
}
#endif /* ELEMENTDATAENUM_H_ */
