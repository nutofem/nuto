#include <zoltan.h>

#include "mechanics/mesh/MeshFem.h"
#include "mechanics/elements/ElementCollection.h"


struct zoltanElement
{
    int id;
    int numNodes;

    zoltanElement() : id(0), numNodes(0) {}
    zoltanElement(int rID, int rNumNodes) : id(rID), numNodes(rNumNodes) {}
};


class ZoltanMesh : public NuTo::MeshFem
{
public:
    ZoltanMesh() : NuTo::MeshFem(){}
    ZoltanMesh(int rDim) : dim(rDim){}
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
    struct Zoltan_DD_Struct* zoltan_dd;
    std::vector<zoltanElement> zoltanElements;
    std::vector<NuTo::ElementCollectionFem> femElements;

    void copyElements();

    //Query functions for hyperpartitioning
    static int get_number_of_localElements(void *data, int *ierr);
    static void get_localElement_list(void *data, int sizeGID, int sizeLID, ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID, int wgt_dim, float *obj_wgts, int *ierr);
    static void get_number_of_localNodes_localPins(void *data, int *num_lists, int *num_nonzeroes, int *format, int *ierr);
    static void get_hypergraph(void *data, int sizeGID, int num_edges, int num_nonzeroes, int format, ZOLTAN_ID_PTR edgeGID, int *vtxPtr, ZOLTAN_ID_PTR vtxGID, int *ierr);

    //Visualization
    static void showHypergraph_1D(int myProc, int numProcs, int numLocalIDs, int numGlobalIDs, ZOLTAN_ID_TYPE *GIDs, int *parts);
    static void printMesh(int myProc, ZoltanMesh* mesh);

    //Query functions for migration
    static int get_element_size(void *data, int gidSize, int lidSize, ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID, int* ierr);
    static void get_element_sizes(void *data, int gidSize, int lidSize, int num_ids, ZOLTAN_ID_PTR globalIDs, ZOLTAN_ID_PTR localIDs, int* sizes, int* ierr);
    static void pack_element(void *data, int gidSize, int lidSize, ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID, int dest, int size, char *buf, int *ierr);
    static void pack_elements(void *data, int gidSize, int lidSize, int num_ids, ZOLTAN_ID_PTR globalIDs, ZOLTAN_ID_PTR localIDs, int *dests, int *sizes, int *idx, char *buf, int *ierr);
    static void unpack_element(void *data, int gidSize, ZOLTAN_ID_PTR globalID, int size, char *buf, int *ierr);
    static void unpack_elements(void *data, int gidSize, int num_ids, ZOLTAN_ID_PTR globalIDs, int *size, int *idx, char *buf, int *ierr);

    static void mid_migrate(void *data, int gidSize, int lidSize, int numImport, ZOLTAN_ID_PTR importGlobalID, ZOLTAN_ID_PTR importLocalID, int *importProc, int *importPart, int numExport, ZOLTAN_ID_PTR exportGlobalID, ZOLTAN_ID_PTR exportLocalID, int *exportProc, int *exportPart, int *ierr);

    //Default creation
    static ZoltanMesh create_1D_mesh(int rNodeCount);
    static ZoltanMesh create_2D_mesh();

    //Default distribution
    static void distributeToRoot(int rRank, int rNumProc, ZoltanMesh* rMesh);
    static void distribute(int rRank, int rNumProc, ZoltanMesh* rMesh);
};
