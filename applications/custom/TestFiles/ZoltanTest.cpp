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

#include "mechanics/mesh/MeshFem.h"
#include "mechanics/mesh/MeshFemDofConvert.h"
#include "mechanics/interpolation/InterpolationQuadLinear.h"

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

    const auto& interpolation = mesh.CreateInterpolation(InterpolationQuadLinear(dim));
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


static int get_number_of_vertices(void *data, int *ierr);
static void get_vertex_list(void *data, int sizeGID, int sizeLID, ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID, int wgt_dim, float *obj_wgts, int *ierr);
static void get_hypergraph_size(void *data, int *num_lists, int *num_nonzeroes, int *format, int *ierr);
static void get_hypergraph(void *data, int sizeGID, int num_edges, int num_nonzeroes, int format, ZOLTAN_ID_PTR edgeGID, int *vtxPtr, ZOLTAN_ID_PTR vtxGID, int *ierr);
static void showHypergraph_1D(int myProc, int numProcs, int numLocalIDs, int numGlobalIDs, ZOLTAN_ID_TYPE *GIDs, int *parts);

void generateMesh_HyperPartition(int myRank, int numProcs, MyMeshFem* mesh)
{
    int numGlobalVertices, numGlobalEdges, numGlobalNZ;
    int num, count, nnbors, ack=0;
    int to=-1, from, remaining;
    int vGID;
    int i, j;
    int vals[128], send_count[3];
    ZOLTAN_ID_TYPE *idx;
    unsigned int id;
    MPI_Status status;
    int ack_tag = 5, count_tag = 10, id_tag = 15;
    MyMeshFem* send_hg;
    MyMeshFem* send_mesh;
    int numGlobalElements = 0;
    int numGlobalNodes = 0;
    numGlobalVertices = mesh->numGlobalElements;
    numGlobalEdges = mesh->numGlobalNodes;
    numGlobalNZ = mesh->numGlobalPins;
//    numGlobalElements = mesh->numGlobalElements;
//    numGlobalNodes = mesh->numGlobalNodes;

    if (myRank == 0)
    {
        /* Create a sub graph for each process */

        send_hg = (MyMeshFem *)calloc(sizeof(MyMeshFem) , numProcs);

        /*
        * Divide the vertices across the processes
        *            (elements)*/

        remaining = numGlobalVertices;
//        remaining = numGlobalEdges;
        count = (numGlobalVertices / numProcs) + 1;
        idx = mesh->myGlobalElementIDs;

        for (i=0; i < numProcs; i++){

            if (remaining == 0) count = 0;
            if (count > remaining) count = remaining;

            send_hg[i].numLocalElements = count;

            if (count){

                send_hg[i].myGlobalElementIDs = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * count);

                for (j=0; j < count; j++){
                    send_hg[i].myGlobalElementIDs[j] = *idx++;
                }
            }

            remaining -= count;
        }

        /*
* Assign hyperedges to processes, and create a sub-hypergraph for each process.
*/

        remaining = numGlobalEdges;
        count = (numGlobalEdges / numProcs) + 1;
        from = 0;

        for (i=0; i < numProcs; i++){

            if (remaining == 0) count = 0;
            if (count > remaining) count = remaining;

            send_hg[i].numLocalNodes = count;
            send_hg[i].numLocalPins = 0;
//            send_hg[i].numLocalElements = 0;

            if (count > 0){

                to = from + count;
//                to = from + count - 1;

//                nnbors = global_hg.nborIndex[to] - global_hg.nborIndex[from];
                nnbors = mesh->neighborIndex[to] - mesh->neighborIndex[from];

                send_hg[i].numLocalPins = nnbors;

                send_hg[i].myGlobalNodeIDs = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * count);
                memcpy(send_hg[i].myGlobalNodeIDs, mesh->myGlobalNodeIDs + from, sizeof(ZOLTAN_ID_TYPE) * count);

                send_hg[i].neighborIndex = (int *)malloc(sizeof(int) * (count + 1));
                send_hg[i].neighborIndex[0] = 0;

                if (nnbors > 0){

//                    num = global_hg.nborIndex[from];
                    num = mesh->neighborIndex[from];

                    for (j=1; j <= count; j++){
                        send_hg[i].neighborIndex[j] = mesh->neighborIndex[from+j] - num;
                    }

//                    send_hg[i].nborGID = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * nnbors);
                    send_hg[i].myPinGIDs = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * nnbors);
//                    memcpy(send_hg[i].nborGID,
//                           global_hg.nborGID + global_hg.nborIndex[from],
//                           sizeof(ZOLTAN_ID_TYPE) * nnbors);
                    memcpy(send_hg[i].myPinGIDs,
                           mesh->myPinGIDs + mesh->neighborIndex[from],
                           sizeof(ZOLTAN_ID_TYPE) * nnbors);
                }
            }

            remaining -= count;
            from = to;
        }

        /* Send each process its hyperedges and the vertices in its partition */

//        *hg = send_hg[0];
        *mesh = std::move(send_hg[0]);
        mesh->numGlobalElements = numGlobalVertices;
        mesh->numGlobalNodes = numGlobalEdges;
        mesh->numGlobalPins = numGlobalNZ;

        for (i=1; i < numProcs; i++){
//            send_count[0] = send_hg[i].numMyVertices;
//            send_count[1] = send_hg[i].numMyHEdges;
//            send_count[2] = send_hg[i].numAllNbors;
            send_count[0] = send_hg[i].numLocalElements;
            send_count[1] = send_hg[i].numLocalNodes;
            send_count[2] = send_hg[i].numLocalPins;

            MPI_Send(send_count, 3, MPI_INT, i, count_tag, MPI_COMM_WORLD);
            MPI_Recv(&ack, 1, MPI_INT, i, ack_tag, MPI_COMM_WORLD, &status);

            if (send_count[0] > 0){
//                MPI_Send(send_hg[i].vtxGID, send_count[0], ZOLTAN_ID_MPI_TYPE, i, id_tag, MPI_COMM_WORLD);
//                free(send_hg[i].vtxGID);
                MPI_Send(send_hg[i].myGlobalElementIDs, send_count[0], ZOLTAN_ID_MPI_TYPE, i, id_tag, MPI_COMM_WORLD);
                free(send_hg[i].myGlobalElementIDs);
            }

            if (send_count[1] > 0){
//                MPI_Send(send_hg[i].edgeGID, send_count[1], ZOLTAN_ID_MPI_TYPE, i, id_tag + 1, MPI_COMM_WORLD);
//                free(send_hg[i].edgeGID);
                MPI_Send(send_hg[i].myGlobalNodeIDs, send_count[1], ZOLTAN_ID_MPI_TYPE, i, id_tag + 1, MPI_COMM_WORLD);
                free(send_hg[i].myGlobalNodeIDs);

//                MPI_Send(send_hg[i].nborIndex, send_count[1] + 1, MPI_INT, i, id_tag + 2, MPI_COMM_WORLD);
//                free(send_hg[i].nborIndex);
                MPI_Send(send_hg[i].neighborIndex, send_count[1] + 1, MPI_INT, i, id_tag + 2, MPI_COMM_WORLD);
                free(send_hg[i].neighborIndex);

                if (send_count[2] > 0){
//                    MPI_Send(send_hg[i].nborGID, send_count[2], ZOLTAN_ID_MPI_TYPE, i, id_tag + 3, MPI_COMM_WORLD);
//                    free(send_hg[i].nborGID);
                    MPI_Send(send_hg[i].myPinGIDs, send_count[2], ZOLTAN_ID_MPI_TYPE, i, id_tag + 3, MPI_COMM_WORLD);
                    free(send_hg[i].myPinGIDs);
                }
            }
        }

        free(send_hg);

        /* signal all procs it is OK to go on */
        ack = 0;
        for (i=1; i < numProcs; i++){
            MPI_Send(&ack, 1, MPI_INT, i, ack_tag, MPI_COMM_WORLD);
        }
    }
    else{

        MPI_Recv(send_count, 3, MPI_INT, 0, count_tag, MPI_COMM_WORLD, &status);

        if (send_count[0] < 0){
            MPI_Finalize();
            exit(1);
        }

        ack = 0;

        memset(mesh, 0, sizeof(MyMeshFem));
        mesh->numGlobalElements = numGlobalVertices;
        mesh->numGlobalNodes = numGlobalEdges;
        mesh->numGlobalPins = numGlobalNZ;

//        hg->numMyVertices = send_count[0];
//        hg->numMyHEdges = send_count[1];
//        hg->numAllNbors = send_count[2];
        mesh->numLocalElements = send_count[0];
        mesh->numLocalNodes = send_count[1];
        mesh->numLocalPins = send_count[2];

        if (send_count[0] > 0){
//            hg->vtxGID = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * send_count[0]);
            mesh->myGlobalElementIDs = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * send_count[0]);
        }

        if (send_count[1] > 0){
//            hg->edgeGID = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * send_count[1]);
//            hg->nborIndex = (int *)malloc(sizeof(int) * (send_count[1] + 1));
            mesh->myGlobalNodeIDs = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * send_count[1]);
            mesh->neighborIndex = (int *)malloc(sizeof(int) * (send_count[1] + 1));

            if (send_count[2] > 0){
//                hg->nborGID = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * send_count[2]);
                mesh->myPinGIDs = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * send_count[2]);
            }
        }

        MPI_Send(&ack, 1, MPI_INT, 0, ack_tag, MPI_COMM_WORLD);

        if (send_count[0] > 0){
//            MPI_Recv(hg->vtxGID,send_count[0], ZOLTAN_ID_MPI_TYPE, 0, id_tag, MPI_COMM_WORLD, &status);
            MPI_Recv(mesh->myGlobalElementIDs,send_count[0], ZOLTAN_ID_MPI_TYPE, 0, id_tag, MPI_COMM_WORLD, &status);

            if (send_count[1] > 0){
                MPI_Recv(mesh->myGlobalNodeIDs,send_count[1], ZOLTAN_ID_MPI_TYPE, 0, id_tag + 1, MPI_COMM_WORLD, &status);
                MPI_Recv(mesh->neighborIndex,send_count[1] + 1, MPI_INT, 0, id_tag + 2, MPI_COMM_WORLD, &status);

                if (send_count[2] > 0){
//                    MPI_Recv(hg->nborGID,send_count[2], ZOLTAN_ID_MPI_TYPE, 0, id_tag + 3, MPI_COMM_WORLD, &status);
                    MPI_Recv(mesh->myPinGIDs,send_count[2], ZOLTAN_ID_MPI_TYPE, 0, id_tag + 3, MPI_COMM_WORLD, &status);
                }
            }
        }

        /* ok to go on? */

        MPI_Recv(&ack, 1, MPI_INT, 0, ack_tag, MPI_COMM_WORLD, &status);
        if (ack < 0){
            MPI_Finalize();
            exit(1);
        }
    }
}

void generateMesh_HyperPartition_z(int myRank, int numProcs, ZoltanMesh* mesh)
{
    int numGlobalVertices, numGlobalEdges, numGlobalNZ;
    int num, count, nnbors, ack=0;
    int to=-1, from, remaining;
    int vGID;
    int i, j;
    int vals[128], send_count[3];
    ZOLTAN_ID_TYPE *idx;
    unsigned int id;
    MPI_Status status;
    int ack_tag = 5, count_tag = 10, id_tag = 15;
    ZoltanMesh* send_hg;
    ZoltanMesh* send_mesh;
    int numGlobalElements = 0;
    int numGlobalNodes = 0;
    numGlobalVertices = mesh->numGlobalElements;
    numGlobalEdges = mesh->numGlobalNodes;
    numGlobalNZ = mesh->numGlobalPins;
//    numGlobalElements = mesh->numGlobalElements;
//    numGlobalNodes = mesh->numGlobalNodes;

    if (myRank == 0)
    {
        /* Create a sub graph for each process */

        send_hg = (ZoltanMesh *)calloc(sizeof(ZoltanMesh) , numProcs);

        /*
        * Divide the vertices across the processes
        *            (elements)*/

        remaining = numGlobalVertices;
//        remaining = numGlobalEdges;
        count = (numGlobalVertices / numProcs) + 1;
        idx = mesh->myGlobalElementIDs;

        for (i=0; i < numProcs; i++){

            if (remaining == 0) count = 0;
            if (count > remaining) count = remaining;

            send_hg[i].numLocalElements = count;

            if (count){

                send_hg[i].myGlobalElementIDs = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * count);

                for (j=0; j < count; j++){
                    send_hg[i].myGlobalElementIDs[j] = *idx++;
                }
            }

            remaining -= count;
        }

        /*
* Assign hyperedges to processes, and create a sub-hypergraph for each process.
*/

        remaining = numGlobalEdges;
        count = (numGlobalEdges / numProcs) + 1;
        from = 0;

        for (i=0; i < numProcs; i++){

            if (remaining == 0) count = 0;
            if (count > remaining) count = remaining;

            send_hg[i].numLocalNodes = count;
            send_hg[i].numLocalPins = 0;
//            send_hg[i].numLocalElements = 0;

            if (count > 0){

                to = from + count;
//                to = from + count - 1;

//                nnbors = global_hg.nborIndex[to] - global_hg.nborIndex[from];
                nnbors = mesh->neighborIndex[to] - mesh->neighborIndex[from];

                send_hg[i].numLocalPins = nnbors;

                send_hg[i].myGlobalNodeIDs = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * count);
                memcpy(send_hg[i].myGlobalNodeIDs, mesh->myGlobalNodeIDs + from, sizeof(ZOLTAN_ID_TYPE) * count);

                send_hg[i].neighborIndex = (int *)malloc(sizeof(int) * (count + 1));
                send_hg[i].neighborIndex[0] = 0;

                if (nnbors > 0){

//                    num = global_hg.nborIndex[from];
                    num = mesh->neighborIndex[from];

                    for (j=1; j <= count; j++){
                        send_hg[i].neighborIndex[j] = mesh->neighborIndex[from+j] - num;
                    }

//                    send_hg[i].nborGID = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * nnbors);
                    send_hg[i].myPinGIDs = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * nnbors);
//                    memcpy(send_hg[i].nborGID,
//                           global_hg.nborGID + global_hg.nborIndex[from],
//                           sizeof(ZOLTAN_ID_TYPE) * nnbors);
                    memcpy(send_hg[i].myPinGIDs,
                           mesh->myPinGIDs + mesh->neighborIndex[from],
                           sizeof(ZOLTAN_ID_TYPE) * nnbors);
                }
            }

            remaining -= count;
            from = to;
        }

        /* Send each process its hyperedges and the vertices in its partition */

//        *hg = send_hg[0];
        *mesh = std::move(send_hg[0]);
        mesh->numGlobalElements = numGlobalVertices;
        mesh->numGlobalNodes = numGlobalEdges;
        mesh->numGlobalPins = numGlobalNZ;

        for (i=1; i < numProcs; i++){
//            send_count[0] = send_hg[i].numMyVertices;
//            send_count[1] = send_hg[i].numMyHEdges;
//            send_count[2] = send_hg[i].numAllNbors;
            send_count[0] = send_hg[i].numLocalElements;
            send_count[1] = send_hg[i].numLocalNodes;
            send_count[2] = send_hg[i].numLocalPins;

            MPI_Send(send_count, 3, MPI_INT, i, count_tag, MPI_COMM_WORLD);
            MPI_Recv(&ack, 1, MPI_INT, i, ack_tag, MPI_COMM_WORLD, &status);

            if (send_count[0] > 0){
//                MPI_Send(send_hg[i].vtxGID, send_count[0], ZOLTAN_ID_MPI_TYPE, i, id_tag, MPI_COMM_WORLD);
//                free(send_hg[i].vtxGID);
                MPI_Send(send_hg[i].myGlobalElementIDs, send_count[0], ZOLTAN_ID_MPI_TYPE, i, id_tag, MPI_COMM_WORLD);
                free(send_hg[i].myGlobalElementIDs);
            }

            if (send_count[1] > 0){
//                MPI_Send(send_hg[i].edgeGID, send_count[1], ZOLTAN_ID_MPI_TYPE, i, id_tag + 1, MPI_COMM_WORLD);
//                free(send_hg[i].edgeGID);
                MPI_Send(send_hg[i].myGlobalNodeIDs, send_count[1], ZOLTAN_ID_MPI_TYPE, i, id_tag + 1, MPI_COMM_WORLD);
                free(send_hg[i].myGlobalNodeIDs);

//                MPI_Send(send_hg[i].nborIndex, send_count[1] + 1, MPI_INT, i, id_tag + 2, MPI_COMM_WORLD);
//                free(send_hg[i].nborIndex);
                MPI_Send(send_hg[i].neighborIndex, send_count[1] + 1, MPI_INT, i, id_tag + 2, MPI_COMM_WORLD);
                free(send_hg[i].neighborIndex);

                if (send_count[2] > 0){
//                    MPI_Send(send_hg[i].nborGID, send_count[2], ZOLTAN_ID_MPI_TYPE, i, id_tag + 3, MPI_COMM_WORLD);
//                    free(send_hg[i].nborGID);
                    MPI_Send(send_hg[i].myPinGIDs, send_count[2], ZOLTAN_ID_MPI_TYPE, i, id_tag + 3, MPI_COMM_WORLD);
                    free(send_hg[i].myPinGIDs);
                }
            }
        }

        free(send_hg);

        /* signal all procs it is OK to go on */
        ack = 0;
        for (i=1; i < numProcs; i++){
            MPI_Send(&ack, 1, MPI_INT, i, ack_tag, MPI_COMM_WORLD);
        }
    }
    else{

        MPI_Recv(send_count, 3, MPI_INT, 0, count_tag, MPI_COMM_WORLD, &status);

        if (send_count[0] < 0){
            MPI_Finalize();
            exit(1);
        }

        ack = 0;

        memset(mesh, 0, sizeof(ZoltanMesh));
        mesh->numGlobalElements = numGlobalVertices;
        mesh->numGlobalNodes = numGlobalEdges;
        mesh->numGlobalPins = numGlobalNZ;

//        hg->numMyVertices = send_count[0];
//        hg->numMyHEdges = send_count[1];
//        hg->numAllNbors = send_count[2];
        mesh->numLocalElements = send_count[0];
        mesh->numLocalNodes = send_count[1];
        mesh->numLocalPins = send_count[2];

        if (send_count[0] > 0){
//            hg->vtxGID = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * send_count[0]);
            mesh->myGlobalElementIDs = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * send_count[0]);
        }

        if (send_count[1] > 0){
//            hg->edgeGID = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * send_count[1]);
//            hg->nborIndex = (int *)malloc(sizeof(int) * (send_count[1] + 1));
            mesh->myGlobalNodeIDs = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * send_count[1]);
            mesh->neighborIndex = (int *)malloc(sizeof(int) * (send_count[1] + 1));

            if (send_count[2] > 0){
//                hg->nborGID = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * send_count[2]);
                mesh->myPinGIDs = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * send_count[2]);
            }
        }

        MPI_Send(&ack, 1, MPI_INT, 0, ack_tag, MPI_COMM_WORLD);

        if (send_count[0] > 0){
//            MPI_Recv(hg->vtxGID,send_count[0], ZOLTAN_ID_MPI_TYPE, 0, id_tag, MPI_COMM_WORLD, &status);
            MPI_Recv(mesh->myGlobalElementIDs,send_count[0], ZOLTAN_ID_MPI_TYPE, 0, id_tag, MPI_COMM_WORLD, &status);

            if (send_count[1] > 0){
                MPI_Recv(mesh->myGlobalNodeIDs,send_count[1], ZOLTAN_ID_MPI_TYPE, 0, id_tag + 1, MPI_COMM_WORLD, &status);
                MPI_Recv(mesh->neighborIndex,send_count[1] + 1, MPI_INT, 0, id_tag + 2, MPI_COMM_WORLD, &status);

                if (send_count[2] > 0){
//                    MPI_Recv(hg->nborGID,send_count[2], ZOLTAN_ID_MPI_TYPE, 0, id_tag + 3, MPI_COMM_WORLD, &status);
                    MPI_Recv(mesh->myPinGIDs,send_count[2], ZOLTAN_ID_MPI_TYPE, 0, id_tag + 3, MPI_COMM_WORLD, &status);
                }
            }
        }

        /* ok to go on? */

        MPI_Recv(&ack, 1, MPI_INT, 0, ack_tag, MPI_COMM_WORLD, &status);
        if (ack < 0){
            MPI_Finalize();
            exit(1);
        }
    }
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

    ZoltanMesh myMesh_z = ZoltanMesh::create_1D_mesh(globalNodeCount);
//    ZoltanMesh myMesh_z = ZoltanMesh::create_2D_mesh();
//    generateMesh_HyperPartition_z(myRank, numProcs, &myMesh_z);
//    ZoltanMesh::distribute(myRank, numProcs, &myMesh_z);
    ZoltanMesh::distributeToRoot(myRank, numProcs, &myMesh_z);
    myMesh_z.copyElements();
    std::cout << "Elements = " << myMesh_z.femElements.size() << std::endl;
    rc = Zoltan_Initialize(argc, argv, &ver);

    if (rc != ZOLTAN_OK){
        printf("sorry...\n");
        MPI_Finalize();
        exit(0);
    }

    parts = (int*)malloc(sizeof(int)* myMesh_z.numLocalElements);
    lids = (ZOLTAN_ID_PTR)malloc(sizeof(ZOLTAN_ID_TYPE) * myMesh_z.numLocalElements);

    rc = Zoltan_DD_Create(&dd, MPI_COMM_WORLD, 1, 1, 0, myMesh_z.numLocalElements, 0);

    for (int i = 0; i < myMesh_z.numLocalElements; ++i)
    {
        parts[i] = myRank;
        lids[i] = (ZOLTAN_ID_TYPE)i;
    }

    rc = Zoltan_DD_Update(dd, myMesh_z.myGlobalElementIDs, lids, NULL, parts, myMesh_z.numLocalElements);
    myMesh_z.zoltan_dd = dd;

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
    Zoltan_Set_Obj_Size_Fn(zz, myMesh_z.get_element_size, &myMesh_z);
//    Zoltan_Set_Pack_Obj_Multi_Fn(zz, myMesh_z.pack_elements, &myMesh_z);
    Zoltan_Set_Pack_Obj_Fn(zz, myMesh_z.pack_element, &myMesh_z);
//    Zoltan_Set_Unpack_Obj_Multi_Fn(zz, myMesh_z.unpack_elements, &myMesh_z);
    Zoltan_Set_Unpack_Obj_Fn(zz, myMesh_z.unpack_element, &myMesh_z);
    Zoltan_Set_Mid_Migrate_PP_Fn(zz, myMesh_z.mid_migrate, &myMesh_z);

    ZoltanMesh::printMesh(myRank, &myMesh_z);

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


    rc = Zoltan_DD_Update(dd, myMesh_z.myGlobalElementIDs, lids, NULL, parts, myMesh_z.numLocalElements);
//    rc = Zoltan_Invert_Lists(zz, numImport, importGlobalGids, importLocalGids, importProcs, importToPart, &numExport, &exportGlobalGids, &exportLocalGids, &exportProcs, &exportToPart);
//    rc = Zoltan_Migrate(zz, numImport, importGlobalGids, importLocalGids, importProcs, importToPart, numExport,
//                        exportGlobalGids, exportLocalGids, exportProcs, exportToPart);

    if (myRank == 0){
        printf("Graph partition after calling Zoltan\n");
    }
//    showHypergraph_1D(myRank, numProcs, myMesh.numLocalElements, myMesh.numGlobalElements, myMesh.myGlobalElementIDs, parts);
//    showHypergraph_1D(myRank, numProcs, myMesh_z.numLocalElements, myMesh_z.numGlobalElements, myMesh_z.myGlobalElementIDs, parts);
//    ZoltanMesh::showHypergraph_1D(myRank, numProcs, myMesh_z.numLocalElements, myMesh_z.numGlobalElements, myMesh_z.myGlobalElementIDs, parts);
    ZoltanMesh::printMesh(myRank, &myMesh_z);

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


int main(int argc, char **argv)
{
//    Zoltan_RCB_test(argc, argv);

    Zoltan_HyperGraphPartitioning_test(argc, argv);
    return 0;
}

















