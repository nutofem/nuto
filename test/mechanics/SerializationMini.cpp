#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/nodes/NodeDof.h"

#include <fstream>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

int main(int argc, char* argv[])
{
    std::string file = "NodeOut";
    std::cout << "\n\n\nWriting a node to " << file << "\n\n\n";

    // serialize it
    std::ofstream outFileStream(file);
    boost::archive::text_oarchive outArchive(outFileStream);

    NuTo::NodeDof<1, 0, 0, 0, 0, 0, 0, 0, 0, 0> myNodeDofDef;

    const NuTo::NodeBase* myNodePtr2 = &myNodeDofDef;
    outArchive << const_cast<NuTo::NodeBase*&>(myNodePtr2);

    NuTo::NodeBase* myNodePtr1 = &myNodeDofDef;
    outArchive << myNodePtr1;


    return 0;
}
