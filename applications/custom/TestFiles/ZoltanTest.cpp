#include <Xpetra_MapFactory.hpp>
#include <Xpetra_CrsMatrixFactory.hpp>

#include <Amesos2.hpp>
#include <BelosTpetraAdapter.hpp>
#include <BelosSolverFactory.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosBlockGmresSolMgr.hpp>
#include <BelosPseudoBlockGmresSolMgr.hpp>
#include <BelosEpetraAdapter.hpp>
#include <Ifpack.h>
#include <Ifpack_AdditiveSchwarz.h>

#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_DefaultComm.hpp>

#include <zoltan.h>
#include <zoltan_types.h>

//#include "mechanics/mesh/MeshFem.h"
#include "mechanics/mesh/MeshFemDofConvert.h"
#include "mechanics/interpolation/InterpolationQuadLinear.h"
#include "mechanics/constitutive/LinearElastic.h"
#include "mechanics/cell/Cell.h"
#include "mechanics/cell/SimpleAssembler.h"
#include "mechanics/dofs/DofNumbering.h"
#include "mechanics/constraints/Constraints.h"
#include "mechanics/constraints/ConstraintCompanion.h"
#include "base/Group.h"

#include "../TestClasses/ZoltanMesh.h"


using Teuchos::RCP;
using Teuchos::rcp;

using namespace NuTo;


class MyMeshFem : public MeshFem
{
public:
    MyMeshFem() : MeshFem(){}
    MyMeshFem(int rDim) : dim(rDim){}
    int dim;
    ZOLTAN_ID_PTR myGlobalNodeIDs;
    ZOLTAN_ID_PTR myGlobalElementIDs;
    int numLocalNodes = 0;
    int numLocalElements = 0;
    int numGlobalNodes = 0;
    int numGlobalElements = 0;
    int numLocalPins = 0;
    int numGlobalPins = 0;
    ZOLTAN_ID_PTR myPinGIDs;
    int* neighborIndex;
};


MyMeshFem testMesh_1D(int rNodeCount)
{
    int dim = 1;
    MyMeshFem mesh(dim);
    int nodeCount = rNodeCount;
    int elementCount = (nodeCount > 1 ? nodeCount - 1 : 0);
    double h = 1;
    mesh.numLocalNodes = nodeCount;
    mesh.numLocalElements = elementCount;
    mesh.numGlobalNodes = nodeCount;
    mesh.numGlobalElements = elementCount;
    mesh.numLocalPins = 2*nodeCount - 3;
    mesh.numGlobalPins = 2*nodeCount - 3;

    const auto& interpolation = mesh.CreateInterpolation(InterpolationQuadLinear());
    for (int i = 0; i < elementCount; ++i)
    {
//        auto& n_i0 = mesh.Nodes.Add(Eigen::Vector2d(i*h, 0));
//        auto& n_i1 = mesh.Nodes.Add(Eigen::Vector2d((i+1)*h, 0));
        auto& n_i0 = mesh.Nodes.Add(i*h);
        auto& n_i1 = mesh.Nodes.Add((i+1)*h);

        mesh.Elements.Add({{{n_i0, n_i1}, interpolation}});
    }

    mesh.myGlobalNodeIDs = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * mesh.numGlobalNodes);
    for (int j = 0; j < nodeCount; ++j)
    {
        mesh.myGlobalNodeIDs[j] = j;
    }

    mesh.neighborIndex = new int[nodeCount];
    mesh.neighborIndex[0] = 0;
    for (int i = 1; i < nodeCount-1; ++i)
        mesh.neighborIndex[i] = 2*i-1;
    mesh.neighborIndex[nodeCount-1] = mesh.neighborIndex[nodeCount-2] + 1;

    mesh.myPinGIDs = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * mesh.numGlobalPins);
    int index = -1;
    for (int i = 0; i < mesh.numGlobalPins; ++i)
    {
        if (i % 2 == 0)
            ++index;
        mesh.myPinGIDs[i] = index;
    }

    mesh.myGlobalElementIDs = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * mesh.numGlobalElements);
    for (int j = 0; j < elementCount; ++j)
    {
        mesh.myGlobalElementIDs[j] = j;
    }

    return mesh;
}

MyMeshFem testMesh_1D()
{
    return testMesh_1D(0);
}

void generateMesh(int myRank, int numProcs, int globalNodeCount, MyMeshFem* mesh)
{
    double h = 1.0;
    int ackno_tag = 0;
    int localNum_tag = 1;
    int globalIDs_tag = 2;
    int coords_tag = 3;
    int nobj = 0;
    int ack = 0;
    ZOLTAN_ID_PTR gids;
    float* coords;
    MPI_Status status;


    if (myRank == 0)
    {
//        MyMeshFem mesh = std::move(testMesh_1D(globalNodeCount));
//        mesh = std::move(testMesh_1D(globalNodeCount));

        gids = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * (nobj + 1));
        coords = (float *)malloc(sizeof(float) * (nobj + 1));

        for (int k = 0; k < nobj; ++k)
        {
            gids[k] = k;
            coords[k] = k*h;
        }

        for (int i = 1; i < numProcs; ++i)
        {
            MPI_Send(&nobj, 1, MPI_INT, i, localNum_tag, MPI_COMM_WORLD);
            MPI_Recv(&ack, 1, MPI_INT, i, ackno_tag, MPI_COMM_WORLD, &status);
            if (nobj > 0)
            {
                MPI_Send(&gids, nobj, ZOLTAN_ID_MPI_TYPE, i, globalIDs_tag, MPI_COMM_WORLD);
                MPI_Send(&coords, nobj, MPI_FLOAT, i, coords_tag, MPI_COMM_WORLD);
            }
        }

        for (int i=1; i < numProcs; ++i){
            MPI_Send(&ack, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
        }

        free(gids);
        free(coords);
    }
    else
    {
//        MyMeshFem mesh = std::move(testMesh_1D());
//        mesh = std::move(testMesh_1D());
        MPI_Recv(&(mesh->numLocalNodes), 1, MPI_INT, 0, localNum_tag, MPI_COMM_WORLD, &status);

        ack = 0;
        if (mesh->numLocalNodes > 0)
        {
            MPI_Recv(&(mesh->myGlobalNodeIDs), mesh->numLocalNodes, ZOLTAN_ID_MPI_TYPE, 0, globalIDs_tag, MPI_COMM_WORLD, &status);
            MPI_Recv(&coords, mesh->numLocalNodes, MPI_FLOAT, 0, coords_tag, MPI_COMM_WORLD, &status);
        }
        else
            MPI_Send(&ack, 1, MPI_INT, 0, ackno_tag, MPI_COMM_WORLD);

        MPI_Recv(&ack, 1, MPI_INT, 0, ackno_tag, MPI_COMM_WORLD, &status);
    }
}


static int get_number_of_objects(void *data, int *ierr);
static void get_object_list(void *data, int sizeGID, int sizeLID,
            ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
                  int wgt_dim, float *obj_wgts, int *ierr);
static int get_num_geometry(void *data, int *ierr);
static void get_geometry_list(void *data, int sizeGID, int sizeLID,
             int num_obj, ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
             int num_dim, double *geom_vec, int *ierr);

void showSimpleMeshPartitions(int myProc, int numIDs, int numGlobalIDs, ZOLTAN_ID_PTR GIDs, int *parts);
void showSimpleMeshPartitions_1D(int numIDs, int numGlobalIDs, ZOLTAN_ID_PTR GIDs, int *parts);


void Zoltan_RCB_test(int argc, char **argv)
{
    Teuchos::GlobalMPISession mpiSession(&argc, &argv);

    RCP<const Teuchos::Comm<int>> commTeuchos = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
    int rc, i;
    int myRank = commTeuchos->getRank();
    int numProcs = commTeuchos->getSize();
    struct Zoltan_Struct *zz;
    int changes, numGidEntries, numLidEntries, numImport, numExport;
    ZOLTAN_ID_PTR importGlobalGids, importLocalGids, exportGlobalGids, exportLocalGids;
    int *importProcs, *importToPart, *exportProcs, *exportToPart;
    int *parts;
//    MESH_DATA myMesh;
    float ver;
    rc = Zoltan_Initialize(argc, argv, &ver);

    if (rc != ZOLTAN_OK){
        printf("sorry...\n");
        MPI_Finalize();
        exit(0);
    }

    int nodeCount = 5;
    MyMeshFem myMesh = testMesh_1D(nodeCount);
//    if (myRank == 0)
//        myMesh = std::move(testMesh_1D(nodeCount));
//    else
//        myMesh = std::move(testMesh_1D());

    generateMesh(myRank, numProcs, nodeCount, &myMesh);

    /******************************************************************
      ** Create a Zoltan library structure for this instance of load
      ** balancing.  Set the parameters and query functions that will
      ** govern the library's calculation.  See the Zoltan User's
      ** Guide for the definition of these and many other parameters.
      ******************************************************************/

    zz = Zoltan_Create(MPI_COMM_WORLD);

    /* General parameters */

    Zoltan_Set_Param(zz, "DEBUG_LEVEL", "1");
    Zoltan_Set_Param(zz, "LB_METHOD", "RCB");
    Zoltan_Set_Param(zz, "NUM_GID_ENTRIES", "1");
    Zoltan_Set_Param(zz, "NUM_LID_ENTRIES", "1");
    Zoltan_Set_Param(zz, "OBJ_WEIGHT_DIM", "0");
    Zoltan_Set_Param(zz, "RETURN_LISTS", "ALL");

    /* RCB parameters */

    Zoltan_Set_Param(zz, "RCB_OUTPUT_LEVEL", "0");
    Zoltan_Set_Param(zz, "RCB_RECTILINEAR_BLOCKS", "1");
    /*Zoltan_Set_Param(zz, "RCB_RECTILINEAR_BLOCKS", "0"); */
    Zoltan_Set_Param(zz, "CHECK_GEOM", "1");

    /* Query functions, to provide geometry to Zoltan */

    Zoltan_Set_Num_Obj_Fn(zz, get_number_of_objects, &myMesh);
    Zoltan_Set_Obj_List_Fn(zz, get_object_list, &myMesh);
    Zoltan_Set_Num_Geom_Fn(zz, get_num_geometry, &myMesh);
    Zoltan_Set_Geom_Multi_Fn(zz, get_geometry_list, &myMesh);

    /******************************************************************
      ** Zoltan can now partition the vertices in the simple mesh.
      ** In this simple example, we assume the number of partitions is
      ** equal to the number of processes.  Process rank 0 will own
      ** partition 0, process rank 1 will own partition 1, and so on.
      ******************************************************************/
    rc = Zoltan_LB_Partition(zz, /* input (all remaining fields are output) */
                             &changes,        /* 1 if partitioning was changed, 0 otherwise */
                             &numGidEntries,  /* Number of integers used for a global ID */
                             &numLidEntries,  /* Number of integers used for a local ID */
                             &numImport,      /* Number of vertices to be sent to me */
                             &importGlobalGids,  /* Global IDs of vertices to be sent to me */
                             &importLocalGids,   /* Local IDs of vertices to be sent to me */
                             &importProcs,    /* Process rank for source of each incoming vertex */
                             &importToPart,   /* New partition for each incoming vertex */
                             &numExport,      /* Number of vertices I must send to other processes*/
                             &exportGlobalGids,  /* Global IDs of the vertices I must send */
                             &exportLocalGids,   /* Local IDs of the vertices I must send */
                             &exportProcs,    /* Process to which I send each of the vertices */
                             &exportToPart);  /* Partition to which each vertex will belong */


    if (myRank == 0)
    {
        std::cout << "Partition summary\n-----------------" << std::endl;
        std::cout << "changes: " << changes << std::endl;
        std::cout << "numGidEntries: " << numGidEntries << std::endl;
        std::cout << "numLidEntries: " << numLidEntries << std::endl;
        std::cout << "numImport: " << numImport << std::endl;
        std::cout << "importGlobalGids: " << std::endl;
        for (int i = 0; i < numImport; ++i)
        {
            std::cout << importGlobalGids[i] << ", ";
        }
        std::cout << std::endl;
        std::cout << "importLocalGids: " << std::endl;
        for (int i = 0; i < numImport; ++i)
        {
            std::cout << importLocalGids[i] << ", ";
        }
        std::cout << std::endl;
        std::cout << "importProcs: " << std::endl;
        for (int i = 0; i < numImport; ++i)
        {
            std::cout << importProcs[i] << ", ";
        }
        std::cout << std::endl;
        std::cout << "importToPart: " << std::endl;
        for (int i = 0; i < numImport; ++i)
        {
            std::cout << importToPart[i] << ", ";
        }
        std::cout << "numExport: " << numExport << std::endl;
        std::cout << "exportGlobalGids: " << std::endl;
        for (int i = 0; i < numExport; ++i)
        {
            std::cout << exportGlobalGids[i] << ", ";
        }
        std::cout << std::endl;
        std::cout << "exportLocalGids: " << std::endl;
        for (int i = 0; i < numExport; ++i)
        {
            std::cout << exportLocalGids[i] << ", ";
        }
        std::cout << std::endl;
        std::cout << "exportProcs: " << std::endl;
        for (int i = 0; i < numExport; ++i)
        {
            std::cout << exportProcs[i] << ", ";
        }
        std::cout << std::endl;
        std::cout << "exportToPart: " << std::endl;
        for (int i = 0; i < numExport; ++i)
        {
            std::cout << exportToPart[i] << ", ";
        }
        std::cout << std::endl;
    }

    if (rc != ZOLTAN_OK){
        printf("sorry...\n");
        MPI_Finalize();
        Zoltan_Destroy(&zz);
        exit(0);
    }

    /******************************************************************
      ** Visualize the mesh partitioning before and after calling Zoltan.
      ******************************************************************/

    parts = (int *)malloc(sizeof(int) * myMesh.numLocalNodes);

    for (i=0; i < myMesh.numLocalNodes; i++){
        parts[i] = myRank;
    }

    if (myRank== 0){
        printf("\nMesh partition assignments before calling Zoltan\n");
    }

//    showSimpleMeshPartitions_1D(myMesh.numLocalNodes, myMesh.numGlobalNodes, myMesh.myGlobalNodeIDs, parts);

    for (i=0; i < numExport; i++){
        parts[exportLocalGids[i]] = exportToPart[i];
    }

    if (myRank == 0){
        printf("Mesh partition assignments after calling Zoltan\n");
    }

//    showSimpleMeshPartitions_1D(myMesh.numLocalNodes, myMesh.numGlobalNodes, myMesh.myGlobalNodeIDs, parts);

    free(parts);

    /******************************************************************
      ** Free the arrays allocated by Zoltan_LB_Partition, and free
      ** the storage allocated for the Zoltan structure.
      ******************************************************************/

    Zoltan_LB_Free_Part(&importGlobalGids, &importLocalGids,
                        &importProcs, &importToPart);
    Zoltan_LB_Free_Part(&exportGlobalGids, &exportLocalGids,
                        &exportProcs, &exportToPart);

    Zoltan_Destroy(&zz);

    /**********************
      ** all done ***********
      **********************/

    if (myMesh.numLocalNodes > 0){
//        free(myMesh.myGlobalNodeIDs);
    }
}


static int get_number_of_objects(void *data, int *ierr)
{
    MyMeshFem* mesh = (MyMeshFem*)data;
    *ierr = ZOLTAN_OK;
    return mesh->numLocalNodes;
}

static void get_object_list(void *data, int sizeGID, int sizeLID,
                            ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
                            int wgt_dim, float *obj_wgts, int *ierr)
{
    MyMeshFem *mesh= (MyMeshFem *)data;
    *ierr = ZOLTAN_OK;

    for (int i = 0; i < mesh->numLocalNodes; ++i){
        globalID[i] = mesh->myGlobalNodeIDs[i];
        localID[i] = i;
    }
}

static int get_num_geometry(void *data, int *ierr)
{
    *ierr = ZOLTAN_OK;
    MyMeshFem* mesh = (MyMeshFem*)data;
    return mesh->dim;
}

static void get_geometry_list(void *data, int sizeGID, int sizeLID,
                              int num_obj,
                              ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
                              int num_dim, double *geom_vec, int *ierr)
{
    MyMeshFem* mesh = (MyMeshFem*)data;

    if ( (sizeGID != 1) || (sizeLID != 1)){
        *ierr = ZOLTAN_FATAL;
        return;
    }

    *ierr = ZOLTAN_OK;

    for (int i = 0; i < num_obj; ++i){
        geom_vec[i] = (double)mesh->Nodes[i].mValues[0];
    }

    return;
}


void showSimpleMeshPartitions(int myProc, int numIDs, int numGlobalIDs, ZOLTAN_ID_PTR GIDs, int *parts)
{
int partAssign[25], allPartAssign[25];
int i, j, part;

  memset(partAssign, 0, sizeof(int) * 25);

  for (i=0; i < numIDs; i++){
    partAssign[GIDs[i]-1] = parts[i];
  }

  MPI_Reduce(partAssign, allPartAssign, 25, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);

  if (myProc == 0){

    for (i=20; i >= 0; i-=5){
      for (j=0; j < 5; j++){
        part = allPartAssign[i + j];
        if (j < 4)
          printf("%d-----",part);
        else
          printf("%d\n",part);
      }
      if (i > 0)
        printf("|     |     |     |     |\n");
    }
    printf("\n");
  }
}

void showSimpleMeshPartitions_1D(int numIDs, int numGlobalIDs, ZOLTAN_ID_PTR GIDs, int *parts)
{
    int partAssign[numGlobalIDs], allPartAssign[numGlobalIDs];
    int i, j, part;

    memset(partAssign, 0, sizeof(int) * numGlobalIDs);

    for (i=0; i < numIDs; i++){
        partAssign[GIDs[i]] = parts[i];
    }

    MPI_Reduce(partAssign, allPartAssign, numGlobalIDs, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);

    for (j=0; j < numGlobalIDs; ++j){
        part = allPartAssign[j];
        if (j < numGlobalIDs-1)
            printf("%d-----",part);
        else
            printf("%d\n",part);
    }
    printf("\n");
}


void Zoltan_HyperGraphPartitioning_test(int argc, char** argv)
{
    Teuchos::GlobalMPISession mpiSession(&argc, &argv);

    RCP<const Teuchos::Comm<int>> commTeuchos = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();

    int i, rc;
    float ver;
    struct Zoltan_Struct *zz;
    struct Zoltan_DD_Struct *dd;
    int changes, numGidEntries, numLidEntries, numImport, numExport;
    int myRank = commTeuchos->getRank();
    int numProcs = commTeuchos->getSize();
    ZOLTAN_ID_PTR importGlobalGids, importLocalGids, exportGlobalGids, exportLocalGids;
    int *importProcs, *importToPart, *exportProcs, *exportToPart;
    int *parts;
    ZOLTAN_ID_PTR lids;

    int globalNodeCount = 5;
    if (argc > 1)
        globalNodeCount = atoi(argv[1]);

//    MyMeshFem myMesh = testMesh_1D(globalNodeCount);
//    generateMesh_HyperPartition(myRank, numProcs, &myMesh);

//    ZoltanMesh myMesh_z = ZoltanMesh::create_1D_mesh(globalNodeCount);
//    ZoltanMesh<zoltanElementFEM_2D_QUAD> myMesh_z = ZoltanMesh<zoltanElementFEM_2D_QUAD>::create_2D_mesh();
    ZoltanMesh<zoltanElementFEM_1D> myMesh_z = ZoltanMesh<zoltanElementFEM_1D>::create_1D_mesh(globalNodeCount);
    myMesh_z.convert2NuToEntities_1D();
//    NuTo::NodeSimple n1(4);
//    NuTo::NodeSimple n2(14);
//    std::vector<NuTo::NodeSimple> nodes;
//    nodes.push_back(n1);
//    nodes.push_back(n2);
//    auto& interp = myMesh_z.CreateInterpolation(NuTo::InterpolationQuadLinear());
//    zoltanElementFEM elem(1, nodes, interp);
//    zoltanElementFEM elem2(elem.GetID(), elem.GetNodes(), elem.Interpolation());
//    std::cout << elem2.GetNodes().size() << std::endl;
//    myMesh_z.zoltanFEMElements.push_back(elem);

//    ZoltanMesh myMesh_z = ZoltanMesh::create_2D_mesh();
//    generateMesh_HyperPartition_z(myRank, numProcs, &myMesh_z);
//    ZoltanMesh::distribute(myRank, numProcs, &myMesh_z);

//    ZoltanMesh<zoltanElementFEM_2D_QUAD>::distributeToRoot(myRank, numProcs, &myMesh_z);
    ZoltanMesh<zoltanElementFEM_1D>::distributeToRoot(myRank, numProcs, &myMesh_z);
    rc = Zoltan_Initialize(argc, argv, &ver);

    if (rc != ZOLTAN_OK){
        printf("sorry...\n");
        MPI_Finalize();
        exit(0);
    }

    parts = (int*)malloc(sizeof(int)* myMesh_z.numLocalElements);
    lids = (ZOLTAN_ID_PTR)malloc(sizeof(ZOLTAN_ID_TYPE) * myMesh_z.numLocalElements);

//    rc = Zoltan_DD_Create(&dd, MPI_COMM_WORLD, 1, 1, 0, myMesh_z.numLocalElements, 0);

    for (int i = 0; i < myMesh_z.numLocalElements; ++i)
    {
        parts[i] = myRank;
        lids[i] = (ZOLTAN_ID_TYPE)i;
    }

//    rc = Zoltan_DD_Update(dd, myMesh_z.myGlobalElementIDs, lids, NULL, parts, myMesh_z.numLocalElements);
//    myMesh_z.zoltan_dd = dd;

//    read_input_file(myRank, numProcs, fname, &hg);

    zz = Zoltan_Create(MPI_COMM_WORLD);

    /* General parameters */

    Zoltan_Set_Param(zz, "DEBUG_LEVEL", "0");
    Zoltan_Set_Param(zz, "LB_METHOD", "HYPERGRAPH");   /* partitioning method */
    Zoltan_Set_Param(zz, "HYPERGRAPH_PACKAGE", "PHG"); /* version of method */
    Zoltan_Set_Param(zz, "NUM_GID_ENTRIES", "1");/* global IDs are integers */
    Zoltan_Set_Param(zz, "NUM_LID_ENTRIES", "1");/* local IDs are integers */
    Zoltan_Set_Param(zz, "RETURN_LISTS", "ALL"); /* export AND import lists */
    Zoltan_Set_Param(zz, "OBJ_WEIGHT_DIM", "0"); /* use Zoltan default vertex weights */
    Zoltan_Set_Param(zz, "EDGE_WEIGHT_DIM", "0");/* use Zoltan default hyperedge weights */

    /* PHG parameters  - see the Zoltan User's Guide for many more
    *   (The "REPARTITION" approach asks Zoltan to create a partitioning that is
    *    better but is not too far from the current partitioning, rather than partitioning
    *    from scratch.  It may be faster but of lower quality that LB_APPROACH=PARTITION.)
    */

//    Zoltan_Set_Param(zz, "LB_APPROACH", "REPARTITION");
    Zoltan_Set_Param(zz, "LB_APPROACH", "PARTITION");
    Zoltan_Set_Param(zz, "AUTO_MIGRATE", "TRUE");
//    Zoltan_Set_Param(zz, "ORDER_METHOD", "PARMETIS");

    /* Application defined query functions */
//    Zoltan_Set_Num_Obj_Fn(zz, get_number_of_vertices, &myMesh);
//    Zoltan_Set_Obj_List_Fn(zz, get_vertex_list, &myMesh);
//    Zoltan_Set_HG_Size_CS_Fn(zz, get_hypergraph_size, &myMesh);
//    Zoltan_Set_HG_CS_Fn(zz, get_hypergraph, &myMesh);

    Zoltan_Set_Num_Obj_Fn(zz, myMesh_z.get_number_of_localElements, &myMesh_z);
    Zoltan_Set_Obj_List_Fn(zz, myMesh_z.get_localElement_list, &myMesh_z);
    Zoltan_Set_HG_Size_CS_Fn(zz, myMesh_z.get_number_of_localNodes_localPins, &myMesh_z);
    Zoltan_Set_HG_CS_Fn(zz, myMesh_z.get_hypergraph, &myMesh_z);

    //Query functions for migration
//    Zoltan_Set_Obj_Size_Multi_Fn(zz, myMesh_z.get_element_sizes, &myMesh_z);
//    Zoltan_Set_Pack_Obj_Multi_Fn(zz, myMesh_z.pack_elements, &myMesh_z);
//    Zoltan_Set_Unpack_Obj_Multi_Fn(zz, myMesh_z.unpack_elements, &myMesh_z);

//    Zoltan_Set_Obj_Size_Fn(zz, myMesh_z.get_element_size, &myMesh_z);
//    Zoltan_Set_Pack_Obj_Fn(zz, myMesh_z.pack_element, &myMesh_z);
//    Zoltan_Set_Unpack_Obj_Fn(zz, myMesh_z.unpack_element, &myMesh_z);
//    Zoltan_Set_Mid_Migrate_PP_Fn(zz, myMesh_z.mid_migrate, &myMesh_z);

    Zoltan_Set_Obj_Size_Fn(zz, myMesh_z.get_element_size_elementFEM, &myMesh_z);
    Zoltan_Set_Pack_Obj_Fn(zz, myMesh_z.pack_element_elementFEM, &myMesh_z);
    Zoltan_Set_Unpack_Obj_Fn(zz, myMesh_z.unpack_element_elementFEM, &myMesh_z);
    Zoltan_Set_Mid_Migrate_PP_Fn(zz, myMesh_z.mid_migrate_elementFEM, &myMesh_z);

//    ZoltanMesh<zoltanElementFEM_2D_QUAD>::printMesh(myRank, &myMesh_z);
    ZoltanMesh<zoltanElementFEM_1D>::printMesh(myRank, &myMesh_z);
    rc = Zoltan_LB_Partition(zz, /* input (all remaining fields are output) */
                             &changes,        /* 1 if partitioning was changed, 0 otherwise */
                             &numGidEntries,  /* Number of integers used for a global ID */
                             &numLidEntries,  /* Number of integers used for a local ID */
                             &numImport,      /* Number of vertices to be sent to me */
                             &importGlobalGids,  /* Global IDs of vertices to be sent to me */
                             &importLocalGids,   /* Local IDs of vertices to be sent to me */
                             &importProcs,    /* Process rank for source of each incoming vertex */
                             &importToPart,   /* New partition for each incoming vertex */
                             &numExport,      /* Number of vertices I must send to other processes*/
                             &exportGlobalGids,  /* Global IDs of the vertices I must send */
                             &exportLocalGids,   /* Local IDs of the vertices I must send */
                             &exportProcs,    /* Process to which I send each of the vertices */
                             &exportToPart);  /* Partition to which each vertex will belong */

    if (rc != ZOLTAN_OK){
        printf("sorry...\n");
        MPI_Finalize();
        Zoltan_Destroy(&zz);
        exit(0);
    }

    myMesh_z.convert2NuToEntities_1D();

    if (myRank == 0)
    {
        std::cout << "Partition summary\n-----------------" << std::endl;
        std::cout << "changes: " << changes << std::endl;
        std::cout << "numGidEntries: " << numGidEntries << std::endl;
        std::cout << "numLidEntries: " << numLidEntries << std::endl;
        std::cout << "numImport: " << numImport << std::endl;
        std::cout << "importGlobalGids: " << std::endl;
        for (int i = 0; i < numImport; ++i)
        {
            std::cout << importGlobalGids[i] << ", ";
        }
        std::cout << std::endl;
        std::cout << "importLocalGids: " << std::endl;
        for (int i = 0; i < numImport; ++i)
        {
            std::cout << importLocalGids[i] << ", ";
        }
        std::cout << std::endl;
        std::cout << "importProcs: " << std::endl;
        for (int i = 0; i < numImport; ++i)
        {
            std::cout << importProcs[i] << ", ";
        }
        std::cout << std::endl;
        std::cout << "importToPart: " << std::endl;
        for (int i = 0; i < numImport; ++i)
        {
            std::cout << importToPart[i] << ", ";
        }
        std::cout << "numExport: " << numExport << std::endl;
        std::cout << "exportGlobalGids: " << std::endl;
        for (int i = 0; i < numExport; ++i)
        {
            std::cout << exportGlobalGids[i] << ", ";
        }
        std::cout << std::endl;
        std::cout << "exportLocalGids: " << std::endl;
        for (int i = 0; i < numExport; ++i)
        {
            std::cout << exportLocalGids[i] << ", ";
        }
        std::cout << std::endl;
        std::cout << "exportProcs: " << std::endl;
        for (int i = 0; i < numExport; ++i)
        {
            std::cout << exportProcs[i] << ", ";
        }
        std::cout << std::endl;
        std::cout << "exportToPart: " << std::endl;
        for (int i = 0; i < numExport; ++i)
        {
            std::cout << exportToPart[i] << ", ";
        }
        std::cout << std::endl;
    }

    /******************************************************************
    ** Visualize the hypergraph partitioning before and after calling Zoltan.
    ******************************************************************/

//    parts = (int *)malloc(sizeof(int) * hg.numMyVertices);
//    parts = (int*)malloc(sizeof(int)* myMesh.numLocalElements);


//    for (i=0; i < hg.numMyVertices; i++){
//    for (i=0; i < myMesh.numLocalElements; i++){
//        parts[i] = myRank;
//    }
    for (i=0; i < myMesh_z.numLocalElements; i++){
        parts[i] = myRank;
    }

    if (myRank== 0){
        printf("\nHypergraph partition before calling Zoltan\n");
    }
//    showHypergraph_1D(myRank, numProcs, myMesh.numLocalElements, myMesh.numGlobalElements, myMesh.myGlobalElementIDs, parts);
//    showHypergraph_1D(myRank, numProcs, myMesh_z.numLocalElements, myMesh_z.numGlobalElements, myMesh_z.myGlobalElementIDs, parts);
//    ZoltanMesh::showHypergraph_1D(myRank, numProcs, myMesh_z.numLocalElements, myMesh_z.numGlobalElements, myMesh_z.myGlobalElementIDs, parts);

    for (i=0; i < numExport; i++){
        parts[exportLocalGids[i]] = exportToPart[i];
    }


//    rc = Zoltan_DD_Update(dd, myMesh_z.myGlobalElementIDs, lids, NULL, parts, myMesh_z.numLocalElements);

//    rc = Zoltan_Invert_Lists(zz, numImport, importGlobalGids, importLocalGids, importProcs, importToPart, &numExport, &exportGlobalGids, &exportLocalGids, &exportProcs, &exportToPart);
//    rc = Zoltan_Migrate(zz, numImport, importGlobalGids, importLocalGids, importProcs, importToPart, numExport,
//                        exportGlobalGids, exportLocalGids, exportProcs, exportToPart);

    if (myRank == 0){
        printf("Hypergraph partition after calling Zoltan\n");
    }
//    showHypergraph_1D(myRank, numProcs, myMesh.numLocalElements, myMesh.numGlobalElements, myMesh.myGlobalElementIDs, parts);
//    showHypergraph_1D(myRank, numProcs, myMesh_z.numLocalElements, myMesh_z.numGlobalElements, myMesh_z.myGlobalElementIDs, parts);
//    ZoltanMesh::showHypergraph_1D(myRank, numProcs, myMesh_z.numLocalElements, myMesh_z.numGlobalElements, myMesh_z.myGlobalElementIDs, parts);
//    ZoltanMesh<zoltanElementFEM_2D_QUAD>::printMesh(myRank, &myMesh_z);
    ZoltanMesh<zoltanElementFEM_1D>::printMesh(myRank, &myMesh_z);

//    std::cout << "NumElements = " << myMesh_z.Elements.Size() << std::endl;
//    std::cout << "NumNodes = " << myMesh_z.Elements[0].CoordinateElement().GetNumNodes() << std::endl;
//    std::cout << "Node(0) = " << myMesh_z.Elements[0].CoordinateElement().GetNode(0).GetNumValues() << std::endl;

    Zoltan_LB_Free_Part(&importGlobalGids, &importLocalGids,
                        &importProcs, &importToPart);
    Zoltan_LB_Free_Part(&exportGlobalGids, &exportLocalGids,
                        &exportProcs, &exportToPart);

    Zoltan_Destroy(&zz);
}


static int get_number_of_vertices(void *data, int *ierr)
{
    MyMeshFem* mesh = (MyMeshFem*)data;
    *ierr = ZOLTAN_OK;
    return mesh->numLocalElements;
}


static void get_vertex_list(void *data, int sizeGID, int sizeLID, ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID, int wgt_dim, float *obj_wgts, int *ierr)
{
    MyMeshFem* mesh = (MyMeshFem*)data;
    *ierr = ZOLTAN_OK;

    /* In this example, return the IDs of our vertices, but no weights.
       * Zoltan will assume equally weighted vertices.
       */

    for (int i = 0; i < mesh->numLocalElements; ++i){
        globalID[i] = mesh->myGlobalElementIDs[i];
        localID[i] = i;
    }
}


static void get_hypergraph_size(void *data, int *num_lists, int *num_nonzeroes, int *format, int *ierr)
{
    MyMeshFem* mesh = (MyMeshFem*)data;
    *ierr = ZOLTAN_OK;

//    *num_lists = hg->numMyHEdges;
//    *num_nonzeroes = hg->numAllNbors;
    *num_lists = mesh->numLocalNodes;
    *num_nonzeroes = mesh->numLocalPins;

    /* We will provide compressed hyperedge (row) format.  The alternative is
       * is compressed vertex (column) format: ZOLTAN_COMPRESSED_VERTEX.
       */

    *format = ZOLTAN_COMPRESSED_EDGE;
}


static void get_hypergraph(void *data, int sizeGID, int num_edges, int num_nonzeroes, int format, ZOLTAN_ID_PTR edgeGID, int *vtxPtr, ZOLTAN_ID_PTR vtxGID, int *ierr)
{
    MyMeshFem* mesh = (MyMeshFem*)data;
    *ierr = ZOLTAN_OK;

    if ( (num_edges != mesh->numLocalNodes) || (num_nonzeroes != mesh->numLocalPins) ||
         (format != ZOLTAN_COMPRESSED_EDGE)) {
        *ierr = ZOLTAN_FATAL;
        return;
    }

    for (int i = 0; i < num_edges; ++i){
        edgeGID[i] = mesh->myGlobalNodeIDs[i];
//        vtxPtr[i] = hg->nborIndex[i];
        vtxPtr[i] = mesh->neighborIndex[i];
    }

    for (int i = 0; i < num_nonzeroes; ++i){
        vtxGID[i] = mesh->myPinGIDs[i];

    }
}


static void showHypergraph_1D(int myProc, int numProcs, int numLocalIDs, int numGlobalIDs, ZOLTAN_ID_TYPE *GIDs, int *parts)
{
    int partAssign[numGlobalIDs], allPartAssign[numGlobalIDs];
    int i, j, part;

    memset(partAssign, 0, sizeof(int) * numGlobalIDs);

    for (i=0; i < numLocalIDs; i++){
        partAssign[GIDs[i]] = parts[i];
    }
    MPI_Reduce(partAssign, allPartAssign, numGlobalIDs, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);

    printf("[");
    for (j=0; j < numGlobalIDs; ++j){
        part = allPartAssign[j];
        if (j < numGlobalIDs-1)
            printf("]---%d---[",part);
        else
            printf("]---%d---[]\n",part);
    }
    printf("\n");

//    std::cout << "Proc: " << myProc << "\nProc's: " << numProcs << "\nLocalElements: " << numIDs << "\nGlobalElements: " << num_vert << "\nGlobalNodes: " << num_ed << std::endl;
//    std::cout << "Local elements on Proc [" << myProc << "]" << std::endl;
//    for (int i = 0; i < numIDs; ++i)
//    {
////        if (parts[i] == myProc)
////        {
//            std::cout << GIDs[i] << " ";
////        }
//    }
//    std::cout << std::endl;
}


std::map<int, int> local2GlobalDofNumbering(std::map<int, std::vector<int>> globalDofs, std::map<int, std::vector<int>> localDofs)
{
    std::map<int,int> local2Global;
    for (std::pair<int, std::vector<int>> nodeDofs : localDofs)
    {
        for (int i = 0; i < nodeDofs.second.size(); ++i)
        {
            local2Global[nodeDofs.second[i]] = globalDofs[nodeDofs.first][i];
        }
    }
    return local2Global;
}


using Teuchos::RCP;
using Teuchos::rcp;

void AssemblerZoltanTest_1D(int argc, char** argv)
{
    Teuchos::GlobalMPISession mpiSession(&argc, &argv);
    RCP<const Teuchos::Comm<int>> commTeuchos = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();

    int myRank = commTeuchos->getRank();
    int nodeCount = 5;

    if (argc > 1)
        nodeCount = atoi(argv[1]);

    std::cout << "Before partition" << std::endl;
//    ZoltanMesh<zoltanElementFEM_1D> myMesh_z = ZoltanMesh<zoltanElementFEM_1D>::create_1D_mesh(nodeCount);
//    myMesh_z.convert2NuToEntities_1D();
    ZoltanMesh<zoltanElementFEM_2D_QUAD> myMesh_z = ZoltanMesh<zoltanElementFEM_2D_QUAD>::create_2D_mesh();
    myMesh_z.convert2NuToEntities_2D();

    std::cout << myMesh_z.Elements[0].CoordinateElement().GetNode(0).GetValues() << std::endl;
//    std::map<int, std::vector<int>> globalDofs = myMesh_z.getDofNumbering_1D();
    std::map<int, std::vector<int>> globalDofs = myMesh_z.getDofNumbering_2D();

    DofType displ("displacements", 2);
    const auto& interpolation = myMesh_z.CreateInterpolation(InterpolationQuadLinear());

    AddDofInterpolation(&myMesh_z, displ, interpolation);

    Constraint::Constraints constraints;
//    Group<NodeSimple> nodesConstrainedInX = mesh.NodesAtAxis(eDirection::X, dof);
    Group<NodeSimple> nodesConstrainedInY = Group<NodeSimple>(myMesh_z.NodeAtCoordinate(Eigen::Vector2d(0, 0), displ));

//    constraints.Add(dof, Constraint::Component(nodesConstrainedInX, {eDirection::X}));
    constraints.Add(displ, Constraint::Component(nodesConstrainedInY, {eDirection::Y}));

    DofNumbering::DofInfo dofInfo = DofNumbering::Build(myMesh_z.NodesTotal(displ), displ, constraints);

    //-------------- mesh partitioning

    int i, rc;
    float ver;
    struct Zoltan_Struct *zz;
    struct Zoltan_DD_Struct *dd;
    int changes, numGidEntries, numLidEntries, numImport, numExport;
    int numProcs = commTeuchos->getSize();
    ZOLTAN_ID_PTR importGlobalGids, importLocalGids, exportGlobalGids, exportLocalGids;
    int *importProcs, *importToPart, *exportProcs, *exportToPart;
    int *parts;
    ZOLTAN_ID_PTR lids;

//    ZoltanMesh<zoltanElementFEM_1D>::distributeToRoot(myRank, numProcs, &myMesh_z);
    ZoltanMesh<zoltanElementFEM_2D_QUAD>::distributeToRoot(myRank, numProcs, &myMesh_z);
    rc = Zoltan_Initialize(argc, argv, &ver);

    if (rc != ZOLTAN_OK){
        printf("sorry...\n");
        MPI_Finalize();
        exit(0);
    }

    zz = Zoltan_Create(MPI_COMM_WORLD);

    /* General parameters */
    Zoltan_Set_Param(zz, "DEBUG_LEVEL", "0");
    Zoltan_Set_Param(zz, "LB_METHOD", "HYPERGRAPH");   /* partitioning method */
    Zoltan_Set_Param(zz, "HYPERGRAPH_PACKAGE", "PHG"); /* version of method */
    Zoltan_Set_Param(zz, "NUM_GID_ENTRIES", "1");/* global IDs are integers */
    Zoltan_Set_Param(zz, "NUM_LID_ENTRIES", "1");/* local IDs are integers */
    Zoltan_Set_Param(zz, "RETURN_LISTS", "ALL"); /* export AND import lists */
    Zoltan_Set_Param(zz, "OBJ_WEIGHT_DIM", "0"); /* use Zoltan default vertex weights */
    Zoltan_Set_Param(zz, "EDGE_WEIGHT_DIM", "0");/* use Zoltan default hyperedge weights */

//    Zoltan_Set_Param(zz, "LB_APPROACH", "REPARTITION");
    Zoltan_Set_Param(zz, "LB_APPROACH", "PARTITION");
    Zoltan_Set_Param(zz, "AUTO_MIGRATE", "TRUE");
//    Zoltan_Set_Param(zz, "ORDER_METHOD", "PARMETIS");

    /* Application defined query functions */
    Zoltan_Set_Num_Obj_Fn(zz, myMesh_z.get_number_of_localElements, &myMesh_z);
    Zoltan_Set_Obj_List_Fn(zz, myMesh_z.get_localElement_list, &myMesh_z);
    Zoltan_Set_HG_Size_CS_Fn(zz, myMesh_z.get_number_of_localNodes_localPins, &myMesh_z);
    Zoltan_Set_HG_CS_Fn(zz, myMesh_z.get_hypergraph, &myMesh_z);

    //Query functions for migration
    Zoltan_Set_Obj_Size_Fn(zz, myMesh_z.get_element_size_elementFEM, &myMesh_z);
    Zoltan_Set_Pack_Obj_Fn(zz, myMesh_z.pack_element_elementFEM, &myMesh_z);
    Zoltan_Set_Unpack_Obj_Fn(zz, myMesh_z.unpack_element_elementFEM, &myMesh_z);
    Zoltan_Set_Mid_Migrate_PP_Fn(zz, myMesh_z.mid_migrate_elementFEM, &myMesh_z);

    rc = Zoltan_LB_Partition(zz, /* input (all remaining fields are output) */
                             &changes,        /* 1 if partitioning was changed, 0 otherwise */
                             &numGidEntries,  /* Number of integers used for a global ID */
                             &numLidEntries,  /* Number of integers used for a local ID */
                             &numImport,      /* Number of vertices to be sent to me */
                             &importGlobalGids,  /* Global IDs of vertices to be sent to me */
                             &importLocalGids,   /* Local IDs of vertices to be sent to me */
                             &importProcs,    /* Process rank for source of each incoming vertex */
                             &importToPart,   /* New partition for each incoming vertex */
                             &numExport,      /* Number of vertices I must send to other processes*/
                             &exportGlobalGids,  /* Global IDs of the vertices I must send */
                             &exportLocalGids,   /* Local IDs of the vertices I must send */
                             &exportProcs,    /* Process to which I send each of the vertices */
                             &exportToPart);  /* Partition to which each vertex will belong */

    if (rc != ZOLTAN_OK){
        printf("sorry...\n");
        MPI_Finalize();
        Zoltan_Destroy(&zz);
        exit(0);
    }

    std::cout << "After partition" << std::endl;
//    myMesh_z.convert2NuToEntities_1D();
    myMesh_z.convert2NuToEntities_2D();

    const auto& interpolation_sub = myMesh_z.CreateInterpolation(InterpolationQuadLinear());
    AddDofInterpolation(&myMesh_z, displ, interpolation_sub);

    DofNumbering::DofInfo currDofInfo = DofNumbering::Build(myMesh_z.NodesTotal(displ), displ, constraints);
//    std::map<int, std::vector<int>> localDofs = myMesh_z.getDofNumbering_1D();
    std::map<int, std::vector<int>> localDofs = myMesh_z.getDofNumbering_2D();
//    std::vector<std::map<int, int>> local2GlobalDofs = local2GlobalNumbering(&mesh, &subMesh0, &subMesh1, displ);
    std::map<int, int> local2GlobalDofs = local2GlobalDofNumbering(globalDofs, localDofs);

//    //--------- do the following steps for every part of the mesh

//    // ************************************************************************
//    //                 add continuum cells
//    // ************************************************************************
//    constexpr double E = 20000;
//    constexpr double nu = 0.2;
//    Laws::LinearElastic<2> linearElasticLaw(E, nu, ePlaneState::PLANE_STRESS);
//    Integrands::TimeDependent::MomentumBalance<2> momentumBalance(displ, linearElasticLaw);

//    boost::ptr_vector<CellInterface> cellContainer;
//    IntegrationTypeTensorProduct<2> integrationType(2, eIntegrationMethod::GAUSS);

//    NuTo::Groups::Group<CellInterface> cellGroup;
//    for (auto& element : subMeshes[rank].Elements)
//    {
//        cellContainer.push_back(new Cell(element, integrationType, momentumBalance));
//        cellGroup.Add(cellContainer.back());
//    }

//    // ************************************************************************
//    //                 add boundary cells
//    // ************************************************************************

//    // manually add the boundary element
//    const auto& interpolationBc = subMeshes[rank].CreateInterpolation(InterpolationTrussLinear(2));

//    IntegrationTypeTensorProduct<1> integrationTypeBc(1, eIntegrationMethod::GAUSS);
//    Eigen::Vector2d pressureBC(1, 0);
//    Integrands::TimeDependent::NeumannBc<2> neumannBc(displ, pressureBC);

//    // extract existing nodes
//    NuTo::Groups::Group<NodeSimple> boundaryCoordNodes = subMeshes[rank].NodesAtAxis(eDirection::X, 10);
//    if (boundaryCoordNodes.Size() == 2)
//    {
//        NodeSimple& nc1 = *boundaryCoordNodes.begin();
//        NodeSimple& nc2 = *(boundaryCoordNodes.begin() + 1);

//        auto boundaryDisplNodes = subMeshes[rank].NodesAtAxis(eDirection::X, displ, 10);
//        NodeSimple& nd1 = *boundaryDisplNodes.begin();
//        NodeSimple& nd2 = *(boundaryDisplNodes.begin() + 1);

//        // add the boundary element
//        auto& boundaryElement = subMeshes[rank].Elements.Add({{{nc1, nc2}, interpolationBc}});
//        boundaryElement.AddDofElement(displ, {{nd1, nd2}, interpolationBc});

//        cellContainer.push_back(new Cell(boundaryElement, integrationTypeBc, neumannBc));
//        cellGroup.Add(cellContainer.back());
//    }

//    // ************************************************************************
//    //                  assemble
//    // ************************************************************************
//    SimpleAssembler assembler(currDofInfo.numIndependentDofs, currDofInfo.numDependentDofs);

//    auto gradient = assembler.BuildVector(cellGroup, {&displ}, Integrands::TimeDependent::Gradient());
//    auto hessian = assembler.BuildMatrix(cellGroup, {&displ}, Integrands::TimeDependent::Hessian0());
//    Eigen::SparseMatrix<double> A_JJ = hessian.JJ(displ,displ);
//    Eigen::VectorXd r_J = gradient.J[displ];

//    //-----------

//    //******************************************
//    //*     create overlapping index map       *
//    //******************************************
//    std::vector<int> myGlobalDofIDs = getAllDofNumbers(local2GlobalDofs, rank);
//    int* myGlobalDofIDs_arr = &myGlobalDofIDs[0];
//    RCP<Tpetra::Map<int, int>> overlappingMap_tpetra = rcp(new Tpetra::Map<int, int>(Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(), myGlobalDofIDs_arr, myGlobalDofIDs.size(), 0, commTeuchos));


//    //******************************************
//    //*       create owning index map          *
//    //******************************************
//    std::vector<int> myOwningGlobalActiveDofIDs = getGlobalOwningActiveDofNumbers(local2GlobalDofs, rank, currDofInfo.activeDofs);
//    int* myOwningGlobalActiveDofIDs_arr = &myOwningGlobalActiveDofIDs[0];
//    RCP<Tpetra::Map<int, int>> owningMap_tpetra = rcp(new Tpetra::Map<int, int>(Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(), myOwningGlobalActiveDofIDs_arr, myOwningGlobalActiveDofIDs.size(), 0, commTeuchos));

////    std::ostream &out = std::cout;
////    RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));
////    *fos << "OwningMap :" << std::endl;
////    owningMap_tpetra->describe(*fos,Teuchos::VERB_EXTREME);
////    *fos << std::endl;
////    *fos << "OverlappingMap :" << std::endl;
////    overlappingMap_tpetra->describe(*fos,Teuchos::VERB_EXTREME);
////    *fos << std::endl;

//    //******************************************
//    //*         create index graphs            *
//    //******************************************
//    int maxNonZeros = 18;   //got by interpolation type (order), e.g. QUAD2 + EQUIDISTANT1 => 18
//    RCP<Tpetra::CrsGraph<int, int>> owningGraph_tpetra = rcp(new Tpetra::CrsGraph<int,int>(owningMap_tpetra, maxNonZeros, Tpetra::ProfileType::DynamicProfile));
//    RCP<Tpetra::CrsGraph<int, int>> overlappingGraph_tpetra = rcp(new Tpetra::CrsGraph<int, int>(overlappingMap_tpetra, maxNonZeros, Tpetra::ProfileType::DynamicProfile));
//    std::vector<int> columnIndices;
//    for (int k=0; k<A_JJ.outerSize(); ++k)
//    {
//        columnIndices.clear();
//        for (Eigen::SparseMatrix<double>::InnerIterator it(A_JJ,k); it; ++it)
//        {
//            // describe position of entries
//            overlappingGraph_tpetra->insertGlobalIndices(myGlobalDofIDs_arr[it.row()], 1, &myGlobalDofIDs_arr[it.col()]);
//        }
//    }

//    //******************************************
//    //*   define inter-process communication   *
//    //*      for local-to-global indices       *
//    //******************************************
//    RCP<const Tpetra::Export<int, int>> exporter_tpetra = rcp(new Tpetra::Export<int, int>(overlappingMap_tpetra, owningMap_tpetra));
//    overlappingGraph_tpetra->fillComplete();
//    owningGraph_tpetra->doExport(*overlappingGraph_tpetra.get(), *exporter_tpetra.get(), Tpetra::CombineMode::INSERT);
//    owningGraph_tpetra->fillComplete();


//    //******************************************
//    //*  initialize Trilinos matrix and vector *
//    //******************************************
//    RCP<Tpetra::CrsMatrix<double, int, int>> globalA_JJ_tpetra = rcp(new Tpetra::CrsMatrix<double, int, int>(owningGraph_tpetra));
//    RCP<Tpetra::Vector<double, int, int>> globalRhsVector_tpetra = rcp(new Tpetra::Vector<double, int, int>(owningMap_tpetra));
//    globalRhsVector_tpetra->putScalar(0.0);

//    Eigen::SparseMatrix<double, Eigen::RowMajor> A_JJ_rowMajor(A_JJ);

//    //******************************************
//    //*    conversion from NuTo to Trilinos    *
//    //******************************************
//    ConversionTools converter2;
//    Teuchos::RCP<Tpetra::CrsMatrix<double, int, int>> localA_JJ_tpetra = converter2.convertEigen2TpetraCrsMatrix(A_JJ_rowMajor, overlappingGraph_tpetra, true);
//    Teuchos::RCP<Tpetra::Vector<double, int, int>> localRhsVector_tpetra = converter2.convertEigen2TpetraVector(r_J, overlappingMap_tpetra);
//    globalA_JJ_tpetra->doExport(*localA_JJ_tpetra.get(), *exporter_tpetra.get(), Tpetra::CombineMode::ADD);
//    globalRhsVector_tpetra->doExport(*localRhsVector_tpetra.get(), *exporter_tpetra.get(), Tpetra::CombineMode::INSERT);
//    globalRhsVector_tpetra->scale(-1.);

//    //******************************************
//    //*        solve complete problem          *
//    //******************************************
//    Teuchos::RCP<Tpetra::MultiVector<double, int, int>> sol_tpetra = solveSystem_tpetra(globalA_JJ_tpetra, globalRhsVector_tpetra, false);
}


void AssemblerZoltanTest_2D()
{

}


int main(int argc, char **argv)
{
//    Zoltan_RCB_test(argc, argv);

//    Zoltan_HyperGraphPartitioning_test(argc, argv);

    AssemblerZoltanTest_1D(argc, argv);
    return 0;
}

















