#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include "nuto/mechanics/nodes/NodeBase.h"

#include "nuto/mechanics/nodes/NodeDof.h"
//#include "nuto/mechanics/constraints/ConstraintLinearDisplacementsPeriodic2D.h"
//#include "nuto/mechanics/constraints/ConstraintNode.h"

#include <fstream>

#include <string>


int main(int argc, char* argv[])
{
    std::string file = "NodeOut";
    std::cout << "\n\n\nWriting a node to " << file << "\n\n\n";

    // serialize it
    std::ofstream outFileStream(file);
    boost::archive::text_oarchive outArchive(outFileStream);

    NuTo::NodeDof<1, 0, 0, 0, 0, 0, 0, 0, 0, 0> myNodeDofDef;

    NuTo::NodeBase* myNodePtr = &myNodeDofDef;
    outArchive << myNodePtr;

//    NuTo::NodeBase& myNodeAdr = myNodeDofDef;
//    std::cout << "Reference to (reference): " << typeid(myNodeAdr).name() << std::endl;
//    outArchive << myNodeAdr;

    const NuTo::NodeBase* myNodePtrConst = &myNodeDofDef;
    outArchive << const_cast<NuTo::NodeBase*&>(myNodePtrConst);

    return 0;
}
