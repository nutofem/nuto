#include "nuto/mechanics/structures/unstructured/StructureFETI.h"
#include "/usr/lib/openmpi/include/mpi.h"
#include <boost/mpi.hpp>

#include "nuto/mechanics/nodes/NodeEnum.h"
#include "nuto/mechanics/groups/GroupEnum.h"
#include "nuto/mechanics/sections/SectionEnum.h"
#include "nuto/mechanics/constitutive/ConstitutiveEnum.h"
#include "nuto/visualize/VisualizeEnum.h"
#include "nuto/mechanics/interpolationtypes/InterpolationTypeEnum.h"
#include "nuto/mechanics/elements/IpDataEnum.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeEnum.h"
#include "nuto/base/ErrorEnum.h"
#include "nuto/mechanics/sections/SectionEnum.h"



using std::cout;
using std::endl;
using NuTo::Constitutive::eConstitutiveType;
using NuTo::Constitutive::eConstitutiveParameter;
using NuTo::Node::eDof;
using NuTo::Interpolation::eTypeOrder;
using NuTo::Interpolation::eShapeType;
using NuTo::eGroupId;
using NuTo::eVisualizeWhat;
using NuTo::eSectionType;


constexpr int dim = 2;

int main(int argc, char* argv[])
{
    boost::mpi::environment  env(argc, argv);
    boost::mpi::communicator world;

    const int rank = world.rank();

    NuTo::StructureFETI structure(dim);

    std::string meshFile = "feti_rectangle_tri.msh" + std::to_string(rank);
    std::cout << meshFile << std::endl;


    const int interpolationTypeId = structure.InterpolationTypeCreate(eShapeType::TRIANGLE2D);
    structure.InterpolationTypeAdd(interpolationTypeId, eDof::COORDINATES,     eTypeOrder::EQUIDISTANT1);
    structure.InterpolationTypeAdd(interpolationTypeId, eDof::DISPLACEMENTS,   eTypeOrder::EQUIDISTANT1);

    structure.ImportMeshJson(meshFile,interpolationTypeId);



}
