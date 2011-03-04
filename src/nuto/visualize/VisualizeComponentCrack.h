// $Id$

/*
 * VisualizeComponentCrack.h
 *
 *  Created on: October 2010
 *      Author: Daniel Arnold, Institut for Structural Mechanics, Bauhaus-Universität Weimar
 */

#ifndef VISUALIZECOMPONENTCRACK_H_
#define VISUALIZECOMPONENTCRACK_H_

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/visualize/VisualizeComponentBase.h"

namespace NuTo
{
//! @author Daniel Arnold
//! @date October 2010
//! @brief visualize the Crack-ID
class VisualizeComponentCrack : public VisualizeComponentBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:
	VisualizeComponentCrack();

    inline NuTo::VisualizeBase::eVisualizeWhat GetComponentEnum()const
    {
    	return NuTo::VisualizeBase::CRACK;
    }

    inline std::string GetComponentName()const
    {
    	return std::string("Crack");
    }

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION
protected:
};
}
#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::VisualizeComponentCrack)
#endif // ENABLE_SERIALIZATION

#endif /* VISUALIZECOMPONENTCRACK_H_ */
