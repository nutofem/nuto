// $Id$

#ifndef CrackBase_H
#define CrackBase_H

#include <string>

#ifdef ENABLE_VISUALIZE
#include "nuto/visualize/VisualizeUnstructuredGrid.h"
#endif // ENABLE_VISUALIZE

namespace NuTo
{
class ElementBase;
class NodeBase;

//! @author Daniel Arnold, ISM
//! @date October 2010
//! @brief ... standard abstract class for all cracks
class CrackBase
{

public:
	//! @brief constructor
	CrackBase();
	CrackBase(const std::string);

	//! @brief destructor
	~CrackBase();

	//! @brief gives the crack information
	virtual void Info(int rVerboseLevel)const=0;

	//! @brief initiate the crack into the structure
    //! @param rElements ... vector of element pointer
    //! @param rCrackedElements ... vector of cracked element pointers
	virtual void Initiate(std::vector<ElementBase*> & rElements, std::vector<ElementBase*> & rCrackedElements)=0;

    //! @brief appends a crack point to the end of the crack
    //! @param rMember new member
	virtual void PushBack(NuTo::NodeBase* rMember)=0;

    //! @brief appends a crack point to the begin of the crack
    //! @param rMember new member
	virtual void PushFront(NuTo::NodeBase* rMember)=0;

    //! @brief computes the intersection point of a crack and a line
    //! @return a bool if the line segment AB is cracked or not
    //! @param rNodeA (Input) ... const reference to first node of the line
    //! @param rNodeB (Input) ... const reference to second node of the line
    //! @param rNodeI (Output) ... reference to intersection node of the line with the crack
    //! @param rDist  (Output) ... relative coordinate of the intersection
    //! @param rSeg   (Output) ... id of the crack segment which intersects the line segment AB (max(size_t) if AB is not cracked)
    virtual bool Intersect(const NuTo::NodeBase* rNodeA, const NuTo::NodeBase* rNodeB, NuTo::NodeBase* rNodeI, double& rDist, size_t& rSeg)=0;

    //! @brief computes the intersection point between the end-ray of this crack and a line segment
    //! @return a bool if the line segment AB is intersected or not
    //! @param rNodeA (Input) ... const reference to first node of the line
    //! @param rNodeB (Input) ... const reference to second node of the line
    //! @param rNodeI (Output) ... reference to intersection node of the line with the crack
    //! @param rDist  (Output) ... relative coordinate of the segment intersection
    virtual bool ExtendEnd(const NuTo::NodeBase* rNodeA, const NuTo::NodeBase* rNodeB, NuTo::NodeBase* rNodeI, double& rDist)=0;

#ifdef ENABLE_VISUALIZE
    virtual void Visualize(VisualizeUnstructuredGrid& rVisualize) const;
#endif // ENABLE_VISUALIZE

protected:
    //! @brief ... name of the crack
	std::string mName;

}; //class definition
}
#endif //CrackBase_H
