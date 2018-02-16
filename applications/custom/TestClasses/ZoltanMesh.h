#include <zoltan.h>

#include <map>

#include "mechanics/mesh/MeshFem.h"
#include "mechanics/elements/ElementCollection.h"
#include "mechanics/interpolation/InterpolationQuadLinear.h"

#include "ZoltanMeshEntities.h"


template <typename zoltanElementType>
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
//    std::vector<zoltanElement> zoltanElements;
//    std::vector<zoltanElementFEM<2>> zoltanFEMElements;
//    std::vector<zoltanElement> zoltanFEMElements;
//    std::vector<zoltanElementFEM_1D> zoltanFEMElements;
    std::vector<zoltanElementType> zoltanFEMElements;

    void convert2NuToEntities_1D()
    {
        const auto& interpolation = this->CreateInterpolation(NuTo::InterpolationQuadLinear());

        this->Nodes.Erase(this->Nodes.begin(), this->Nodes.end());
        this->Elements.Erase(this->Elements.begin(), this->Elements.end());

        for (zoltanElementFEM_1D elem : this->zoltanFEMElements)
        {
//            NuTo::NodeSimple n_i0(elem.getNode(0).getValues());
//            NuTo::NodeSimple n_i1(elem.getNode(1).getValues());

//            this->Nodes.Add(n_i0);
//            this->Nodes.Add(n_i1);
            NuTo::NodeSimple& n_i0 = this->Nodes.Add(elem.getNode(0).getValues());
            NuTo::NodeSimple& n_i1 = this->Nodes.Add(elem.getNode(1).getValues());

            this->Elements.Add({{{n_i0, n_i1}, interpolation}});
        }
    }

    void convert2NuToEntities_2D()
    {
        const auto& interpolation = this->CreateInterpolation(NuTo::InterpolationQuadLinear());

        this->Nodes.Erase(this->Nodes.begin(), this->Nodes.end());
        this->Elements.Erase(this->Elements.begin(), this->Elements.end());

        std::map<int, bool> nodeCreatedMap;
        std::map<int, NuTo::NodeSimple*> nodeMap;


        NuTo::NodeSimple* n_i0;
        NuTo::NodeSimple* n_i1;
        NuTo::NodeSimple* n_i2;
        NuTo::NodeSimple* n_i3;

        for (zoltanElementFEM_2D_QUAD elem : this->zoltanFEMElements)
        {
            if (nodeCreatedMap[elem.getNode(0).getID()])
                n_i0 = nodeMap[elem.getNode(0).getID()];
            else
            {
                n_i0 = &(this->Nodes.Add(elem.getNode(0).getValues()));
                nodeMap[elem.getNode(0).getID()] = n_i0;
                nodeCreatedMap[elem.getNode(0).getID()] = true;
            }
            if (nodeCreatedMap[elem.getNode(1).getID()])
                n_i1 = nodeMap[elem.getNode(1).getID()];
            else
            {
                n_i1 = &(this->Nodes.Add(elem.getNode(1).getValues()));
                nodeMap[elem.getNode(1).getID()] = n_i1;
                nodeCreatedMap[elem.getNode(1).getID()] = true;
            }
            if (nodeCreatedMap[elem.getNode(2).getID()])
                n_i2 = nodeMap[elem.getNode(2).getID()];
            else
            {
                n_i2 = &(this->Nodes.Add(elem.getNode(2).getValues()));
                nodeMap[elem.getNode(2).getID()] = n_i2;
                nodeCreatedMap[elem.getNode(2).getID()] = true;
            }
            if (nodeCreatedMap[elem.getNode(3).getID()])
                n_i3 = nodeMap[elem.getNode(3).getID()];
            else
            {
                n_i3 = &(this->Nodes.Add(elem.getNode(3).getValues()));
                nodeMap[elem.getNode(3).getID()] = n_i3;
                nodeCreatedMap[elem.getNode(3).getID()] = true;
            }

            this->Elements.Add({{{n_i0, n_i1, n_i2, n_i3}, interpolation}});
        }

//        for (zoltanElementFEM_2D_QUAD elem : this->zoltanFEMElements)
//        {
//            NuTo::NodeSimple& n_i0 = this->Nodes.Add(elem.getNode(0).getValues());
//            NuTo::NodeSimple& n_i1 = this->Nodes.Add(elem.getNode(1).getValues());
//            NuTo::NodeSimple& n_i2 = this->Nodes.Add(elem.getNode(2).getValues());
//            NuTo::NodeSimple& n_i3 = this->Nodes.Add(elem.getNode(3).getValues());
//            this->Elements.Add({{{n_i0, n_i1, n_i2, n_i3}, interpolation}});
//        }
    }


    void copyDofNumbering(NuTo::DofType rDof)
    {
        for (int iElem = 0; iElem < this->Elements.Size(); ++iElem)
        {
            NuTo::ElementCollectionFem& elem = this->Elements[iElem];
            for (int i = 0; i < elem.DofElement(rDof).GetNumNodes(); ++i)
            {
                NuTo::NodeSimple currNode = elem.DofElement(rDof).GetNode(i);
                Eigen::VectorXi dofs = currNode.mDofNumbers;
                for (int iDofs = 0; iDofs < dofs.rows(); ++iDofs)
                {
                    this->zoltanFEMElements[iElem].getNode(i).setDofNumber(iDofs, dofs(iDofs));
                }
            }
        }
    }


    std::map<int, std::vector<int>> getDofNumbering_1D()
    {
        std::map<int, std::vector<int>> dofMapping;
        for (zoltanElementType elem : this->zoltanFEMElements)
        {
            for (zoltanNode<1, 1> node : elem.getNodes())
            {
                int id = node.getID();
                auto dofs = node.getDofNumbers();
                std::vector<int> currDofs;

                for (int i = 0; i < dofs.rows(); ++i)
                {
                    currDofs.push_back(dofs[i]);
                }
                dofMapping[id] = currDofs;
            }
        }

        return dofMapping;
    }

    std::map<int, std::vector<int>> getDofNumbering_2D()
    {
        std::map<int, std::vector<int>> dofMapping;
        for (zoltanElementType elem : this->zoltanFEMElements)
        {
            for (zoltanNode<2, 2> node : elem.getNodes())
            {
                int id = node.getID();
                auto dofs = node.getDofNumbers();
                std::vector<int> currDofs;

                for (int i = 0; i < dofs.rows(); ++i)
                {
                    currDofs.push_back(dofs[i]);
                }
                dofMapping[id] = currDofs;
            }
        }

        return dofMapping;
    }

    //Query functions for hyperpartitioning
    static int get_number_of_localElements(void *data, int *ierr)
    {
        ZoltanMesh* mesh = (ZoltanMesh*)data;
        *ierr = ZOLTAN_OK;
        return mesh->numLocalElements;
    }

    static void get_localElement_list(void *data, int sizeGID, int sizeLID, ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID, int wgt_dim, float *obj_wgts, int *ierr)
    {
        ZoltanMesh* mesh = (ZoltanMesh*)data;
        *ierr = ZOLTAN_OK;

        /* In this example, return the IDs of our vertices, but no weights.
           * Zoltan will assume equally weighted vertices.
           */

        for (int i = 0; i < mesh->numLocalElements; ++i){
            globalID[i] = mesh->myGlobalElementIDs[i];
            localID[i] = i;
        }
    }

    static void get_number_of_localNodes_localPins(void *data, int *num_lists, int *num_nonzeroes, int *format, int *ierr)
    {
        ZoltanMesh* mesh = (ZoltanMesh*)data;
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
        ZoltanMesh* mesh = (ZoltanMesh*)data;
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

    //Visualization
    static void showHypergraph_1D(int myProc, int numProcs, int numLocalIDs, int numGlobalIDs, ZOLTAN_ID_TYPE *GIDs, int *parts)
    {
        int partAssign[numGlobalIDs], allPartAssign[numGlobalIDs];
        int part;

        memset(partAssign, 0, sizeof(int) * numGlobalIDs);

        std::cout << "numLocalIDs (on " << myProc << ") " << numLocalIDs << std::endl;
        for (int i = 0; i < numLocalIDs; ++i){
            partAssign[GIDs[i]] = parts[i];
        }
        MPI_Reduce(partAssign, allPartAssign, numGlobalIDs, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);

        printf("[");
        for (int j = 0; j < numGlobalIDs; ++j){
            part = allPartAssign[j];
            if (j < numGlobalIDs-1)
                printf("]---%d---[",part);
            else
                printf("]---%d---[]\n",part);
        }
        printf("\n");
    }

    static void printMesh(int myProc, ZoltanMesh *mesh)
    {
        printf("Proc(%d) NumElements = %d\n", myProc, mesh->numLocalElements);
        int numNodes = 0;

        std::cout << "Elements = " << mesh->Elements.Size() << std::endl;
        std::cout << "Nodes = " << mesh->Nodes.Size() << std::endl;

        for (int i = 0; i < mesh->numLocalElements; ++i) {
    //        printf("Proc(%d) LocalElement [%d] with GID [%d]: %d nodes", myProc, i, mesh->myGlobalElementIDs[i], mesh->zoltanFEMElements[i].getNumNodes());
    //        printf("Proc(%d) LocalElement [%d] with GID [%d]: %d nodes", myProc, i, mesh->zoltanFEMElements[i].getID(), mesh->zoltanFEMElements[i].getNumNodes());
            numNodes = mesh->zoltanFEMElements[i].getNumNodes();
            printf("Proc(%d) LocalElement [%d] with GID [%d]: %d nodes ", myProc, i, mesh->zoltanFEMElements[i].getID(), numNodes);
            for (int j = 0; j < numNodes; ++j)
            {
                printf("{");
                printf(mesh->zoltanFEMElements[i].getNode(j).toString().c_str());
                if (j < numNodes-1)
                    printf("} ");
               else
                    printf("}");
    //            printf("%d", mesh->zoltanFEMElements[i].getNode(j).getID());
    //            if (j < numNodes-1)
    //            {
    //                printf(", ");
    //            }
            }
            printf("\n");
            fflush(stdout);
        }
    }

    //Query functions for migration
    static int get_element_size(void *data, int gidSize, int lidSize, ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID, int *ierr)
    {
        return sizeof(zoltanElement);
        //    return sizeof(std::string);
    }

    static int get_element_size_elementFEM(void *data, int gidSize, int lidSize, ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID, int *ierr)
    {
    //    return sizeof(NuTo::ElementCollectionFem);
        ZoltanMesh* mesh = (ZoltanMesh*)data;
        if (mesh->dim == 1)
            return sizeof(zoltanElementFEM_1D);
        else
            return sizeof(zoltanElementFEM_2D_QUAD);
    }

    static void get_element_sizes(void *data, int gidSize, int lidSize, int num_ids, ZOLTAN_ID_PTR globalIDs, ZOLTAN_ID_PTR localIDs, int *sizes, int *ierr)
    {
        ZoltanMesh *mesh;
        int len;

        mesh = (ZoltanMesh *)data;
        *ierr = ZOLTAN_OK;

        //save complete elements
        for (int i = 0; i < num_ids; ++i)
        {
            sizes[i] = sizeof(NuTo::ElementCollectionFem);
    //        sizes[i] = sizeof(int);
    //        sizes[i] = sizeof(std::string);
        }
    }

    static void pack_element(void *data, int gidSize, int lidSize, ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID, int dest, int size, char *buf, int *ierr)
    {
    //    NuTo::ElementCollectionFem* elem;
        ZoltanMesh* mesh = (ZoltanMesh*)data;
        *ierr = ZOLTAN_OK;

    //    elem = (NuTo::ElementCollectionFem*)buf;
    //    memcpy(elem, &(mesh->Elements[*localID]), sizeof(NuTo::ElementCollectionFem));

        zoltanElement* elem = (zoltanElement*)buf;
//        elem[0] = mesh->zoltanElements[*localID];

    //    NuTo::ElementCollectionFem* elem = (NuTo::ElementCollectionFem*)buf;
    ////    elem[0] = mesh->femElements[*localID];
    //    elem = new (buf) NuTo::ElementCollectionFem(mesh->femElements[*localID]);

    //    std::cout << &elem << " " << &elem[0] << std::endl;

    //    auto it = mesh->Elements.begin();
    //    for (int i = 0; i < mesh->numLocalElements; ++i)
    //    {
    //        std::cout << "myGlobalElementIDs[i] = " << mesh->myGlobalElementIDs[i] << std::endl;
    //        if (mesh->myGlobalElementIDs[i] == *globalID)
    //        {
    //            mesh->Elements.Erase(it);
    //            for (int j = i; j < mesh->numLocalElements-1; ++j)
    //            {
    //                mesh->myGlobalElementIDs[j] = mesh->myGlobalElementIDs[j+1];
    //            }
    //            break;
    //        }
    //        else
    //        {
    //            it++;
    //        }
    //    }

    //    mesh->numLocalElements -= 1;
    //    mesh->myGlobalElementIDs = (ZOLTAN_ID_PTR)realloc(mesh->myGlobalElementIDs, sizeof(ZOLTAN_ID_TYPE)*mesh->numLocalElements);

    //    std::string* str = new (buf) std::string("hallo");
        //    std::cout << &str << " " << &str [0] << std::endl;
    }

    static void pack_element_elementFEM(void *data, int gidSize, int lidSize, ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID, int dest, int size, char *buf, int *ierr)
    {
        ZoltanMesh* mesh = (ZoltanMesh*)data;
    //    zoltanElementFEM* elem = (zoltanElementFEM*)buf;
    //    elem[0] = mesh->zoltanFEMElements[*localID];

        zoltanElementType* elem = new (buf) zoltanElementType(mesh->zoltanFEMElements[*localID]);
//        if (mesh->dim == 1)
//            zoltanElementFEM_1D* elem = new (buf) zoltanElementFEM_1D(mesh->zoltanFEMElements[*localID]);
//        else
//            zoltanElementFEM_2D_QUAD* elem = new (buf) zoltanElementFEM_2D_QUAD(mesh->zoltanFEMElements[*localID]);
    }

    static void pack_elements(void *data, int gidSize, int lidSize, int num_ids, ZOLTAN_ID_PTR globalIDs, ZOLTAN_ID_PTR localIDs, int *dests, int *sizes, int *idx, char *buf, int *ierr)
    {
        int num_nbors;
        ZOLTAN_ID_TYPE *nbors=NULL, *ibuf=NULL;
        NuTo::ElementCollectionFem* elems = NULL;
    //    elems = (NuTo::ElementCollectionFem*)malloc(num_ids * sizeof(NuTo::ElementCollectionFem));
        Zoltan_Memory_Debug(1);
        elems = (NuTo::ElementCollectionFem*) Zoltan_Malloc(num_ids * sizeof(NuTo::ElementCollectionFem), __FILE__, __LINE__);
        int* someValues;
        std::string* someText;
        std::vector<NuTo::ElementCollectionFem> elemVec;
    //    NuTo::ValueVector<NuTo::ElementCollectionFem> els;
        ZOLTAN_ID_TYPE lid;
        ZoltanMesh *mesh;
        *ierr = ZOLTAN_OK;
        mesh = (ZoltanMesh *)data;

        //For each element write complete elements (nodes, interpolation) to the buffer
    //    elems = (NuTo::ElementCollectionFem *)calloc(num_ids, sizeof(NuTo::ElementCollectionFem));

        for (int i = 0; i < num_ids; ++i){
            lid = localIDs[i];
    //        nbors = graph->nborGID + graph->nborIndex[lid];
    //        num_nbors = graph->nborIndex[lid+1] - graph->nborIndex[lid];

    //        ibuf = (ZOLTAN_ID_TYPE *)(buf + idx[i]);
    //        ibuf[0] = num_nbors;

    //        for (j=1; j <= num_nbors; j++){
    //            ibuf[j] = *nbors++;
    //        }

    //        someText = (std::string*)(buf + idx[i]);
    //        std::cout << "someText = " << &someText[0] << std::endl;
    //        someText = new (buf + idx[i]) std::string("aha");

    //        std::string* st = (std::string*)(buf + idx[i]);
    //        std::cout << "st = " << *st << std::endl;

    //        someValues = (int*)(buf + idx[i]);
    //        someValues[0] = i;
    //        std::cout << "pack someValues address " << &someValues << std::endl;
    //        std::cout << "pack someValues = " << someValues[0] << std::endl;

    //        elemVec.push_back(mesh->Elements[lid]);
    //        elems = new (buf + idx[i]) NuTo::ElementCollectionFem(mesh->Elements[lid]);
            elems = (NuTo::ElementCollectionFem*)(buf);
            std::cout << "pack elems address " << &elems << std::endl;
            std::cout << "pack elems address " << &elems[0] << std::endl;
            elems = (NuTo::ElementCollectionFem*)(buf + idx[i]);
            std::cout << "pack buf = " << atoi(buf) << std::endl;
            std::cout << "pack idx[i] = " << idx[i] << std::endl;
            std::cout << "pack elems address " << &elems << std::endl;
            std::cout << "pack elems address " << &elems[0] << std::endl;
            std::cout << "pack elements[lid] = " << mesh->Elements[lid].CoordinateElement().GetNumNodes() << std::endl;
    //        memcpy(elems, &(mesh->Elements[lid]), sizeof(NuTo::ElementCollectionFem));
            elems[0] = mesh->Elements[lid];
    //        std::cout << "pack elements[lid] = " << elems[0].CoordinateElement().GetNumNodes() << std::endl;
            std::cout << "pack elems address " << &elems << std::endl;
            std::cout << "pack elems address " << &elems[0] << std::endl;
    //        elems[0] = mesh->Elements[lid];
    //        elems = &elemVec[0];

        }

        //    Zoltan_Memory_Stats();
    }

    static void unpack_element(void *data, int gidSize, ZOLTAN_ID_PTR globalID, int size, char *buf, int *ierr)
    {
        ZoltanMesh<zoltanElementType>* mesh = (ZoltanMesh<zoltanElementType>*)data;
    //    NuTo::ElementCollectionFem* newElem;
        int numElems = mesh->numLocalElements + 1;

        mesh->myGlobalElementIDs = (ZOLTAN_ID_TYPE*) realloc(mesh->myGlobalElementIDs, sizeof(ZOLTAN_ID_TYPE) * numElems);
        mesh->myGlobalElementIDs[mesh->numLocalElements] = *globalID;

        zoltanElement* newElem = (zoltanElement*)buf;
//        mesh->zoltanElements.push_back(*newElem);

        mesh->numLocalElements += 1;
    }

    static void unpack_element_elementFEM(void *data, int gidSize, ZOLTAN_ID_PTR globalID, int size, char *buf, int *ierr)
    {
        ZoltanMesh* mesh = (ZoltanMesh*)data;
        int numElems = mesh->numLocalElements + 1;

        mesh->myGlobalElementIDs = (ZOLTAN_ID_TYPE*) realloc(mesh->myGlobalElementIDs, sizeof(ZOLTAN_ID_TYPE) * numElems);
        mesh->myGlobalElementIDs[mesh->numLocalElements] = *globalID;

        zoltanElementType* newElem = (zoltanElementType*)buf;
        mesh->zoltanFEMElements.push_back(newElem[0]);
//        if (mesh->dim == 1)
//        {
//            zoltanElementFEM_1D* newElem = (zoltanElementFEM_1D*)buf;
//            mesh->zoltanFEMElements.push_back(newElem[0]);
//        }
//        else
//        {
//            zoltanElementFEM_2D_QUAD* newElem = (zoltanElementFEM_2D_QUAD*)buf;
//            mesh->zoltanFEMElements.push_back(newElem[0]);
//        }
    //    zoltanElementFEM elem(newElem->GetID(), newElem->GetNodes(), newElem->Interpolation());
    //    mesh->zoltanFEMElements[mesh->numLocalElements] = newElem[0];

        mesh->numLocalElements += 1;
    }

    static void unpack_elements(void *data, int gidSize, int num_ids, ZOLTAN_ID_PTR globalIDs, int *size, int *idx, char *buf, int *ierr)
    {
        int len, next_vertex, next_nbor;
        ZOLTAN_ID_TYPE *ibuf=NULL;
        NuTo::ElementCollectionFem* elems;
        Zoltan_Memory_Debug(1);
    //    elems = (NuTo::ElementCollectionFem*)Zoltan_Malloc(num_ids * sizeof(NuTo::ElementCollectionFem), __FILE__, __LINE__);
        NuTo::ElementCollectionFem* newElem;
    //    newElem = (NuTo::ElementCollectionFem*)malloc(sizeof(NuTo::ElementCollectionFem));
        int* someValues;
        std::string* someText;
        ZoltanMesh *mesh;
        *ierr = ZOLTAN_OK;
        mesh = (ZoltanMesh *)data;

        //Add incoming elements
        int num_myElems = mesh->numLocalElements;

        int num_elems = 0;
        num_elems = num_myElems + num_ids;

        mesh->myGlobalElementIDs = (ZOLTAN_ID_TYPE*) realloc(mesh->myGlobalElementIDs, sizeof(ZOLTAN_ID_TYPE) * num_elems);
    //    mesh->Elements = (NuTo::ElementCollectionFem*) realloc(mesh->Elements, sizeof(NuTo::ElementCollectionFem) * num_elems);

        for (int i = 0; i < num_ids; ++i){
            mesh->myGlobalElementIDs[num_myElems + i] = globalIDs[i];
    //        someText = (std::string*)(buf);
    //        std::cout << "&someText[0] = " << &someText[0] << std::endl;
    //        someText = (std::string*)(buf + idx[i]);
    //        std::cout << "received values" << std::endl;
    //        std::cout << "&someText[0] = " << &someText[0] << std::endl;
    //        std::cout << "someText[0] = " << someText[0] << std::endl;

    //        someValues = (int*)(buf + idx[i]);
    //        std::cout << "unpack elems address " << &someValues << std::endl;
    //        std::cout << "unpack someValues = " << someValues[0] << std::endl;
    //        std::cout << "unpack buf = " << atoi(buf) << std::endl;
    //        std::cout << "unpack idx[i] = " << idx[i] << std::endl;


            elems = (NuTo::ElementCollectionFem *)(buf + idx[i]);
            std::cout << "unpack elems address " << &elems << std::endl;
            std::cout << "unpack elems address " << &elems[0] << std::endl;
            std::cout << "num nodes unpack" << elems[0].CoordinateElement().GetNumNodes() << std::endl;
    //        mesh->Elements.Add(elems[0]);
    //        mesh->Elements.Add(*elems);
        }

    //    graph->numMyVertices += num_ids;
        Zoltan_Memory_Stats();
        mesh->numLocalElements += num_ids;
    }


    static void mid_migrate(void *data, int gidSize, int lidSize, int numImport, ZOLTAN_ID_PTR importGlobalID, ZOLTAN_ID_PTR importLocalID, int *importProc, int *importPart,
                            int numExport, ZOLTAN_ID_PTR exportGlobalID, ZOLTAN_ID_PTR exportLocalID, int *exportProc, int *exportPart, int *ierr)
    {
        ZoltanMesh *mesh;
        int *exports;
    //    int* imports;

        *ierr = ZOLTAN_OK;
        mesh = (ZoltanMesh *)data;

        ZOLTAN_ID_PTR elem_GIDs = mesh->myGlobalElementIDs;
        int newElemCount = mesh->numLocalElements - numExport;

        exports = (int*)calloc(sizeof(int), mesh->numLocalElements);
        for (int i = 0; i < numExport; ++i)
        {
            exports[exportLocalID[i]] = 1;
        }

    //    auto it = mesh->Elements.begin();
//        auto it = mesh->zoltanElements.begin();

        //  mesh->myGlobalElementIDs = (ZOLTAN_ID_TYPE*) realloc(mesh->myGlobalElementIDs, sizeof(ZOLTAN_ID_TYPE) * (mesh->numLocalElements - numExport));
        mesh->myGlobalElementIDs = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * newElemCount);

        int indexCounter = 0;
        for (int i = 0; i < mesh->numLocalElements; ++i)
        {
            if (exports[i] == 1)
            {
    //            mesh->Elements.Erase(it);
//                mesh->zoltanElements.erase(it);
//                it--;
            }
            else
            {
                mesh->myGlobalElementIDs[indexCounter] = elem_GIDs[i];
                indexCounter++;
            }
//            it++;
        }

        free(exports);
        //  free(imports);
        //  graph->numMyVertices = next_vertex;
        mesh->numLocalElements = newElemCount;
    }

    static void mid_migrate_elementFEM(void *data, int gidSize, int lidSize, int numImport, ZOLTAN_ID_PTR importGlobalID, ZOLTAN_ID_PTR importLocalID, int *importProc, int *importPart,
                                       int numExport, ZOLTAN_ID_PTR exportGlobalID, ZOLTAN_ID_PTR exportLocalID, int *exportProc, int *exportPart, int *ierr)
    {
        ZoltanMesh *mesh;
        int *exports;

        *ierr = ZOLTAN_OK;
        mesh = (ZoltanMesh *)data;

        ZOLTAN_ID_PTR elem_GIDs = mesh->myGlobalElementIDs;
        int newElemCount = mesh->numLocalElements - numExport;

        exports = (int*)calloc(sizeof(int), mesh->numLocalElements);
        for (int i = 0; i < numExport; ++i)
        {
            exports[exportLocalID[i]] = 1;
        }

        auto it = mesh->Elements.begin();
        auto it2 = mesh->zoltanFEMElements.begin();


        mesh->myGlobalElementIDs = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * newElemCount);
    //    std::cout << "numExport = " << numExport << ": " << "zoltanFEMElements.size() = " << mesh->zoltanFEMElements.size() << std::endl;

        int j = 0;
        int indexCounter = 0;
        for (int i = 0; i < mesh->numLocalElements; ++i)
        {
            if (exports[i] == 1)
            {
    //            std::cout << "exports[" << i << "] = " << exports[i] << std::endl;
                mesh->Elements.Erase(it);
                mesh->zoltanFEMElements.erase(mesh->zoltanFEMElements.begin() + i-(j++));
    //            mesh->zoltanFEMElements.erase(it2);
                it--;
    //            it2--;
            }
            else
            {
                mesh->myGlobalElementIDs[indexCounter] = elem_GIDs[i];
                indexCounter++;
            }
            it++;
    //        it2++;
        }
    //    std::cout << "numExport = " << numExport << ": " << "zoltanFEMElements.size() = " << mesh->zoltanFEMElements.size() << std::endl;


        free(exports);
        mesh->numLocalElements = newElemCount;
    }

    //Default creation
    static ZoltanMesh<zoltanElementFEM_1D> create_1D_mesh(int rNodeCount)
    {
        ZoltanMesh<zoltanElementFEM_1D> mesh(1);
        int nodeCount = rNodeCount;
        int elementCount = (nodeCount > 1 ? nodeCount - 1 : 0);
        double h = 0.1;
        mesh.numLocalNodes = nodeCount;
        mesh.numLocalElements = elementCount;
        mesh.numGlobalNodes = nodeCount;
        mesh.numGlobalElements = elementCount;
        mesh.numLocalPins = 2*nodeCount - 2;
        mesh.numGlobalPins = 2*nodeCount - 2;

//        mesh.zoltanElements.resize(mesh.numLocalElements);
//        for (int i = 0; i < mesh.numLocalElements; ++i)
//        {
//            mesh.zoltanElements[i].mId = i;
//        }

        //Geometry, not necessary for hyperpartitioning
//        const auto& interpolation = mesh.CreateInterpolation(NuTo::InterpolationQuadLinear());
//        for (int i = 0; i < elementCount; ++i)
//        {
//            auto& n_i0 = mesh.Nodes.Add(i*h);
//            auto& n_i1 = mesh.Nodes.Add((i+1)*h);

//            mesh.Elements.Add({{{n_i0, n_i1}, interpolation}});
//        }

        for (int i = 0; i < elementCount; ++i)
        {
            zoltanNode<1,1> n0(i, i*h);
            n0.setDofNumber(0, i);
            zoltanNode<1,1> n1(i+1, (i+1)*h);
            n1.setDofNumber(0, i+1);
            mesh.zoltanFEMElements.push_back({i, {n0, n1}});
        }


        mesh.myGlobalNodeIDs = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * mesh.numGlobalNodes);
        for (int j = 0; j < nodeCount; ++j)
        {
            mesh.myGlobalNodeIDs[j] = j;
        }

        mesh.neighborIndex = new int[nodeCount];
        mesh.neighborIndex[0] = 0;
        for (int i = 1; i < nodeCount; ++i)
            mesh.neighborIndex[i] = 2*i-1;

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

    static ZoltanMesh<zoltanElementFEM_2D_QUAD> create_2D_mesh()
    {
        /* Something like this:
         *
         *    3-----------------------2
         * /| | - _        e2       / | -->
         * /| |     -7------------6   | -->
         * /| | e4  /     e3      |   | -->
         * /| |    /              |e1 | -->
         * /| |   /    _____------5   | -->
         * /| |  4-----            \  | --> p
         * /| | /      e0           \ | -->
         * /| |/                     \| -->
         *    0-----------------------1
         *   /_\
         *   ///
         *              (c) ttitsche :)
         */


        int dim = 2;
        ZoltanMesh<zoltanElementFEM_2D_QUAD> mesh(dim);
    //    int nodeCount = rNodesX*rNodesY;
    //    int elementsX = (rNodesX > 1 ? rNodesX - 1 : 0);
    //    int elementsY = (rNodesY > 1 ? rNodesY - 1 : 0);
    //    int elementCount = elementsX*elementsY;
    //    double h = 1;

        int nodeCount = 8;
        int elementCount = 5;
        int pinCount = 20;

        mesh.numLocalNodes = nodeCount;
        mesh.numLocalElements = elementCount;
        mesh.numLocalPins = pinCount;
        mesh.numGlobalNodes = nodeCount;
        mesh.numGlobalElements = elementCount;
        mesh.numGlobalPins = pinCount;

        auto& n0 = mesh.Nodes.Add(Eigen::Vector2d(0, 0));
        auto& n1 = mesh.Nodes.Add(Eigen::Vector2d(10, 0));
        auto& n2 = mesh.Nodes.Add(Eigen::Vector2d(10, 10));
        auto& n3 = mesh.Nodes.Add(Eigen::Vector2d(0, 10));

        auto& n4 = mesh.Nodes.Add(Eigen::Vector2d(2, 2));
        auto& n5 = mesh.Nodes.Add(Eigen::Vector2d(8, 3));
        auto& n6 = mesh.Nodes.Add(Eigen::Vector2d(8, 7));
        auto& n7 = mesh.Nodes.Add(Eigen::Vector2d(4, 7));

        const auto& interpolation = mesh.CreateInterpolation(NuTo::InterpolationQuadLinear());

        mesh.Elements.Add({{{n0, n1, n5, n4}, interpolation}});
        mesh.Elements.Add({{{n1, n2, n6, n5}, interpolation}});
        mesh.Elements.Add({{{n7, n6, n2, n3}, interpolation}});
        mesh.Elements.Add({{{n4, n5, n6, n7}, interpolation}});
        mesh.Elements.Add({{{n0, n4, n7, n3}, interpolation}});

        zoltanNode<2,2> n_0(0, Eigen::Vector2d(0, 0));
        zoltanNode<2,2> n_1(1, Eigen::Vector2d(10, 0));
        zoltanNode<2,2> n_2(2, Eigen::Vector2d(10, 10));
        zoltanNode<2,2> n_3(3, Eigen::Vector2d(0, 10));
        zoltanNode<2,2> n_4(4, Eigen::Vector2d(2, 2));
        zoltanNode<2,2> n_5(5, Eigen::Vector2d(8, 3));
        zoltanNode<2,2> n_6(6, Eigen::Vector2d(8, 7));
        zoltanNode<2,2> n_7(7, Eigen::Vector2d(4, 7));
        n_0.setDofNumber(0,0);
        n_1.setDofNumber(0,1);
        n_2.setDofNumber(0,2);
        n_3.setDofNumber(0,3);
        n_4.setDofNumber(0,4);
        n_5.setDofNumber(0,5);
        n_6.setDofNumber(0,6);
        n_7.setDofNumber(0,7);

        zoltanElementFEM_2D_QUAD e_0(0, {n_0, n_1, n_5, n_4});
        zoltanElementFEM_2D_QUAD e_1(1, {n_1, n_2, n_6, n_5});
        zoltanElementFEM_2D_QUAD e_2(2, {n_7, n_6, n_2, n_3});
        zoltanElementFEM_2D_QUAD e_3(3, {n_4, n_5, n_6, n_7});
        zoltanElementFEM_2D_QUAD e_4(4, {n_0, n_4, n_7, n_3});

        mesh.zoltanFEMElements.push_back(e_0);
        mesh.zoltanFEMElements.push_back(e_1);
        mesh.zoltanFEMElements.push_back(e_2);
        mesh.zoltanFEMElements.push_back(e_3);
        mesh.zoltanFEMElements.push_back(e_4);

    //    mesh.zoltanFEMElements.push_back({0, {n_0, n_1, n_5, n_4}, interpolation});
    //    mesh.zoltanFEMElements.push_back({1, {n_1, n_2, n_6, n_5}, interpolation});
    //    mesh.zoltanFEMElements.push_back({2, {n_7, n_6, n_2, n_3}, interpolation});
    //    mesh.zoltanFEMElements.push_back({3, {n_4, n_5, n_6, n_7}, interpolation});
    //    mesh.zoltanFEMElements.push_back({4, {n_0, n_4, n_7, n_3}, interpolation});


        mesh.myGlobalNodeIDs = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * nodeCount);
        for (int j = 0; j < nodeCount; ++j)
        {
            mesh.myGlobalNodeIDs[j] = j;
        }

        mesh.neighborIndex = new int[nodeCount];
        int nborIndex[] = {0, 2, 4, 6, 8, 11, 14, 17};
        std::copy(std::begin(nborIndex), std::end(nborIndex), mesh.neighborIndex);

        mesh.myPinGIDs = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * pinCount);
        unsigned int pinGids[] = {0, 4, 0, 1, 1, 2, 2, 4, 0, 3, 4, 0, 1, 3, 1, 2, 3, 2, 3, 4};
        std::copy(std::begin(pinGids), std::end(pinGids), mesh.myPinGIDs);

        mesh.myGlobalElementIDs = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * elementCount);
        for (int j = 0; j < elementCount; ++j)
        {
            mesh.myGlobalElementIDs[j] = j;
        }

        return mesh;
    }

    static ZoltanMesh<zoltanElementFEM_2D_QUAD> create_2D_mesh(int mx, int my, double lx)
    {
        ZoltanMesh<zoltanElementFEM_2D_QUAD> mesh(2);
        int nodesX = mx + 1;
        int nodesY = my + 1;
        int nodeCount = nodesX * nodesY;
        int elementCount = mx * my;
        double h = lx/mx;

        mesh.numLocalNodes = nodeCount;
        mesh.numLocalElements = elementCount;
        mesh.numGlobalNodes = nodeCount;
        mesh.numGlobalElements = elementCount;

        int cornerNodes = 4;
        int boundaryNodes = 2*(nodesX + nodesY) - 4 - cornerNodes;
        int innerNodes = nodeCount - cornerNodes - boundaryNodes;
        int pinCount = cornerNodes + boundaryNodes*2 + innerNodes*4;

        mesh.numLocalPins = pinCount;
        mesh.numGlobalPins = pinCount;

        std::map<std::pair<double, double>, zoltanNode<2,1>> coord2NodeIDs;
        int indexCounter = 0;
        //number boundary at first
        for (int i = 0; i < nodesX; ++i)
        {
            zoltanNode<2,1> node(indexCounter++, Eigen::Vector2d(i*h, 0));
            std::pair<double, double> coords(i*h, 0);
            coord2NodeIDs[coords] = node;
        }

        for (int i = 1; i < nodesY; ++i)
        {
            zoltanNode<2,1> node(indexCounter++, Eigen::Vector2d(lx, i*h));
            std::pair<double, double> coords(lx, i*h);
            coord2NodeIDs[coords] = node;
        }

        for (int i = 1; i < nodesX; ++i)
        {
            zoltanNode<2,1> node(indexCounter++, Eigen::Vector2d(lx-i*h, lx));
            std::pair<double, double> coords(lx-i*h, lx);
            coord2NodeIDs[coords] = node;
        }

        for (int i = 1; i < nodesY-1; ++i)
        {
            zoltanNode<2,1> node(indexCounter++, Eigen::Vector2d(0, lx-i*h));
            std::pair<double, double> coords(0, lx-i*h);
            coord2NodeIDs[coords] = node;
        }

        //number inner nodes
        for (int i = 1; i < nodesX-1; ++i)
        {
            for (int j = 1; j < nodesY-1; ++j)
            {
                zoltanNode<2,1> node(indexCounter++, Eigen::Vector2d(i*h, i*h));
                std::pair<double, double> coords(i*h, i*h);
                coord2NodeIDs[coords] = node;
            }
        }

        int elementIndexCounter = 0;
        for (int iy = 0; iy < my; ++iy)
        {
            for (int ix = 0; ix < mx; ++ix)
            {
                std::pair<double, double> coords_0(ix*h, iy*h);
                std::pair<double, double> coords_1((ix+1)*h, iy*h);
                std::pair<double, double> coords_2((ix+1)*h, (iy+1)*h);
                std::pair<double, double> coords_3(ix*h, (iy+1)*h);


                zoltanNode<2,1> n_i0 = coord2NodeIDs[coords_0];
                zoltanNode<2,1> n_i1 = coord2NodeIDs[coords_1];
                zoltanNode<2,1> n_i2 = coord2NodeIDs[coords_2];
                zoltanNode<2,1> n_i3 = coord2NodeIDs[coords_3];

                zoltanElementFEM_2D_QUAD e_i(elementIndexCounter++, {n_i0, n_i1, n_i2, n_i3});
                mesh.zoltanFEMElements.push_back(e_i);
            }
        }

        mesh.myGlobalNodeIDs = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * mesh.numGlobalNodes);
        for (int j = 0; j < nodeCount; ++j)
        {
            mesh.myGlobalNodeIDs[j] = j;
        }

        mesh.neighborIndex = new int[nodeCount];
        mesh.neighborIndex[0] = 0;
        for (int i = 1; i < nodeCount; ++i)
            mesh.neighborIndex[i] = 2*i-1;

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
//            mesh.myGlobalElementIDs[j] = j;
            mesh.myGlobalElementIDs[j] = mesh.zoltanFEMElements[j].getID();
        }

        return mesh;
    }


    //Default distribution
    static void distributeToRoot(int rRank, int rNumProc, ZoltanMesh *rMesh)
    {
        int num, count, nnbors, ack=0;
        int to=-1, from, remaining;
        int send_count[3];
        ZOLTAN_ID_PTR myElementGIDs;
        MPI_Status status;
        int ack_tag = 5, count_tag = 10, id_tag = 15;
        ZoltanMesh* send_mesh;
        int numGlobalHyperVertices = rMesh->numGlobalElements;
        int numGlobalHyperEdges = rMesh->numGlobalNodes;
        int numGlobalHyperPins = rMesh->numGlobalPins;
    //    NuTo::ValueVector<NuTo::ElementCollectionFem> myElements = rMesh->Elements;
        NuTo::ElementCollectionFem* myEls;
        myEls = (NuTo::ElementCollectionFem*)calloc(sizeof(NuTo::ElementCollectionFem), rMesh->numLocalElements);
        for (int i = 0; i < rMesh->numLocalElements; ++i)
            myEls[i] = rMesh->Elements[i];

        if (rRank == 0)
        {
            /* Create a sub graph for each process */
            send_mesh = (ZoltanMesh *)calloc(sizeof(ZoltanMesh) , rNumProc);

            /*
            * Divide the vertices across the processes
            *            (elements)*/

            remaining = numGlobalHyperVertices;
            count = numGlobalHyperVertices;
            myElementGIDs = rMesh->myGlobalElementIDs;

            for (int i = 0; i < rNumProc; ++i)
            {
                if (i == 0)
                    count = remaining;
                else
                    count = 0;

                send_mesh[i].numLocalElements = count;

                if (count)
                {
                    send_mesh[i].myGlobalElementIDs = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * count);

                    for (int j = 0; j < count; ++j)
                    {
    //                    send_mesh[i].Elements.Add(myElements[*myElementGIDs]);
                        send_mesh[i].Elements.Add(myEls[*myElementGIDs]);
                        send_mesh[i].myGlobalElementIDs[j] = *myElementGIDs++;
    //                    send_mesh[i].zoltanElements.push_back(*(new zoltanElement(j, j)));
//                        send_mesh[i].zoltanElements.push_back(rMesh->zoltanElements[j]);
                        send_mesh[i].zoltanFEMElements.push_back(rMesh->zoltanFEMElements[j]);
                    }
                }
            }

            /*
            * Assign hyperedges to processes, and create a sub-hypergraph for each process.
            */

            remaining = numGlobalHyperEdges;
            count = numGlobalHyperEdges;
            from = 0;

            for (int i=0; i < rNumProc; ++i)
            {
                if (i == 0)
                    count = remaining;
                else
                    count = 0;

                send_mesh[i].numLocalNodes = count;
                send_mesh[i].numLocalPins = 0;

                if (count > 0)
                {
    //                to = from + count;
                    to = from + count - 1;
                    nnbors = rMesh->neighborIndex[to] - rMesh->neighborIndex[from];
                    std::cout << "rMesh->neighborIndex[to=" << to << "] = " << rMesh->neighborIndex[to] << std::endl;
                    std::cout << "rMesh->neighborIndex[from=" << from << "] = " << rMesh->neighborIndex[from] << std::endl;
                    std::cout << "nnbors = " << nnbors << std::endl;
                    std::cout << "count = " << count << std::endl;
                    send_mesh[i].numLocalPins = nnbors;

                    send_mesh[i].myGlobalNodeIDs = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * count);
                    memcpy(send_mesh[i].myGlobalNodeIDs, rMesh->myGlobalNodeIDs + from, sizeof(ZOLTAN_ID_TYPE) * count);

                    send_mesh[i].neighborIndex = (int *)malloc(sizeof(int) * (count + 1));
                    send_mesh[i].neighborIndex[0] = 0;

                    if (nnbors > 0)
                    {
                        num = rMesh->neighborIndex[from];

                        for (int j=1; j <= count; ++j){
                            send_mesh[i].neighborIndex[j] = rMesh->neighborIndex[from+j] - num;
                        }

                        send_mesh[i].myPinGIDs = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * nnbors);
                        memcpy(send_mesh[i].myPinGIDs, rMesh->myPinGIDs + rMesh->neighborIndex[from], sizeof(ZOLTAN_ID_TYPE) * nnbors);
                    }
                }
                from = to;
            }

            /* Send each process its hyperedges and the vertices in its partition */


            *rMesh = std::move(send_mesh[0]);
            rMesh->numGlobalElements = numGlobalHyperVertices;
            rMesh->numGlobalNodes = numGlobalHyperEdges;
            rMesh->numGlobalPins = numGlobalHyperPins;

    //        rMesh->numLocalElements = numGlobalHyperVertices;
    //        rMesh->numLocalNodes = numGlobalHyperEdges;
    //        rMesh->numLocalPins = numGlobalHyperPins;

            for (int i=1; i < rNumProc; ++i)
            {
                send_count[0] = send_mesh[i].numLocalElements;
                send_count[1] = send_mesh[i].numLocalNodes;
                send_count[2] = send_mesh[i].numLocalPins;

                MPI_Send(send_count, 3, MPI_INT, i, count_tag, MPI_COMM_WORLD);
                MPI_Recv(&ack, 1, MPI_INT, i, ack_tag, MPI_COMM_WORLD, &status);

                if (send_count[0] > 0){
                    MPI_Send(send_mesh[i].myGlobalElementIDs, send_count[0], ZOLTAN_ID_MPI_TYPE, i, id_tag, MPI_COMM_WORLD);
                    free(send_mesh[i].myGlobalElementIDs);
                }

                if (send_count[1] > 0){
                    MPI_Send(send_mesh[i].myGlobalNodeIDs, send_count[1], ZOLTAN_ID_MPI_TYPE, i, id_tag + 1, MPI_COMM_WORLD);
                    free(send_mesh[i].myGlobalNodeIDs);

                    MPI_Send(send_mesh[i].neighborIndex, send_count[1] + 1, MPI_INT, i, id_tag + 2, MPI_COMM_WORLD);
                    free(send_mesh[i].neighborIndex);

                    if (send_count[2] > 0){
                        MPI_Send(send_mesh[i].myPinGIDs, send_count[2], ZOLTAN_ID_MPI_TYPE, i, id_tag + 3, MPI_COMM_WORLD);
                        free(send_mesh[i].myPinGIDs);
                    }
                }
            }

            free(send_mesh);

            /* signal all procs it is OK to go on */
            ack = 0;
            for (int i=1; i < rNumProc; ++i){
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

            memset(rMesh, 0, sizeof(ZoltanMesh));
            rMesh->numGlobalElements = numGlobalHyperVertices;
            rMesh->numGlobalNodes = numGlobalHyperEdges;
            rMesh->numGlobalPins = numGlobalHyperPins;

            rMesh->numLocalElements = send_count[0];
            rMesh->numLocalNodes = send_count[1];
            rMesh->numLocalPins = send_count[2];

            if (send_count[0] > 0){
                rMesh->myGlobalElementIDs = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * send_count[0]);
            }

            if (send_count[1] > 0){
                rMesh->myGlobalNodeIDs = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * send_count[1]);
                rMesh->neighborIndex = (int *)malloc(sizeof(int) * (send_count[1] + 1));

                if (send_count[2] > 0){
                    rMesh->myPinGIDs = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * send_count[2]);
                }
            }

            MPI_Send(&ack, 1, MPI_INT, 0, ack_tag, MPI_COMM_WORLD);

            if (send_count[0] > 0){
                MPI_Recv(rMesh->myGlobalElementIDs,send_count[0], ZOLTAN_ID_MPI_TYPE, 0, id_tag, MPI_COMM_WORLD, &status);

                if (send_count[1] > 0){
                    MPI_Recv(rMesh->myGlobalNodeIDs,send_count[1], ZOLTAN_ID_MPI_TYPE, 0, id_tag + 1, MPI_COMM_WORLD, &status);
                    MPI_Recv(rMesh->neighborIndex,send_count[1] + 1, MPI_INT, 0, id_tag + 2, MPI_COMM_WORLD, &status);

                    if (send_count[2] > 0){
                        MPI_Recv(rMesh->myPinGIDs,send_count[2], ZOLTAN_ID_MPI_TYPE, 0, id_tag + 3, MPI_COMM_WORLD, &status);
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

    static void distribute(int rRank, int rNumProc, ZoltanMesh *rMesh)
    {
        int num, count, nnbors, ack=0;
        int to=-1, from, remaining;
        int send_count[3];
        ZOLTAN_ID_PTR myElementGIDs;
        MPI_Status status;
        int ack_tag = 5, count_tag = 10, id_tag = 15;
        ZoltanMesh* send_mesh;
        int numGlobalHyperVertices = rMesh->numGlobalElements;
        int numGlobalHyperEdges = rMesh->numGlobalNodes;
        int numGlobalHyperPins = rMesh->numGlobalPins;
    //    NuTo::ValueVector<NuTo::ElementCollectionFem> myElements = rMesh->Elements;

        if (rRank == 0)
        {
            /* Create a sub graph for each process */
            send_mesh = (ZoltanMesh *)calloc(sizeof(ZoltanMesh) , rNumProc);

            /*
            * Divide the vertices across the processes
            *            (elements)*/

            remaining = numGlobalHyperVertices;
            count = (numGlobalHyperVertices / rNumProc) + 1;
            myElementGIDs = rMesh->myGlobalElementIDs;

            for (int i = 0; i < rNumProc; ++i)
            {
                if (remaining == 0)
                    count = 0;

                if (count > remaining)
                    count = remaining;

                send_mesh[i].numLocalElements = count;

                if (count)
                {
                    send_mesh[i].myGlobalElementIDs = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * count);

                    for (int j = 0; j < count; ++j)
                    {
    //                    send_mesh[i].Elements.Add(myElements[*myElementGIDs]);
                        send_mesh[i].myGlobalElementIDs[j] = *myElementGIDs++;
                    }
                }
                remaining -= count;
            }

            /*
            * Assign hyperedges to processes, and create a sub-hypergraph for each process.
            */

            remaining = numGlobalHyperEdges;
            count = (numGlobalHyperEdges / rNumProc) + 1;
            from = 0;

            for (int i=0; i < rNumProc; ++i)
            {
                if (remaining == 0)
                    count = 0;

                if (count > remaining)
                    count = remaining;

                send_mesh[i].numLocalNodes = count;
                send_mesh[i].numLocalPins = 0;

                if (count > 0)
                {
                    to = from + count;
    //                to = from + count - 1;
                    nnbors = rMesh->neighborIndex[to] - rMesh->neighborIndex[from];
                    send_mesh[i].numLocalPins = nnbors;

                    send_mesh[i].myGlobalNodeIDs = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * count);
                    memcpy(send_mesh[i].myGlobalNodeIDs, rMesh->myGlobalNodeIDs + from, sizeof(ZOLTAN_ID_TYPE) * count);

                    send_mesh[i].neighborIndex = (int *)malloc(sizeof(int) * (count + 1));
                    send_mesh[i].neighborIndex[0] = 0;

                    if (nnbors > 0)
                    {
                        num = rMesh->neighborIndex[from];

                        for (int j=1; j <= count; ++j){
                            send_mesh[i].neighborIndex[j] = rMesh->neighborIndex[from+j] - num;
                        }

                        send_mesh[i].myPinGIDs = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * nnbors);
                        memcpy(send_mesh[i].myPinGIDs, rMesh->myPinGIDs + rMesh->neighborIndex[from], sizeof(ZOLTAN_ID_TYPE) * nnbors);
                    }
                }

                remaining -= count;
                from = to;
            }

            /* Send each process its hyperedges and the vertices in its partition */

            *rMesh = std::move(send_mesh[0]);
            rMesh->numGlobalElements = numGlobalHyperVertices;
            rMesh->numGlobalNodes = numGlobalHyperEdges;
            rMesh->numGlobalPins = numGlobalHyperPins;

            for (int i=1; i < rNumProc; ++i)
            {
                send_count[0] = send_mesh[i].numLocalElements;
                send_count[1] = send_mesh[i].numLocalNodes;
                send_count[2] = send_mesh[i].numLocalPins;

                MPI_Send(send_count, 3, MPI_INT, i, count_tag, MPI_COMM_WORLD);
                MPI_Recv(&ack, 1, MPI_INT, i, ack_tag, MPI_COMM_WORLD, &status);

                if (send_count[0] > 0){
                    MPI_Send(send_mesh[i].myGlobalElementIDs, send_count[0], ZOLTAN_ID_MPI_TYPE, i, id_tag, MPI_COMM_WORLD);
                    free(send_mesh[i].myGlobalElementIDs);
                }

                if (send_count[1] > 0){
                    MPI_Send(send_mesh[i].myGlobalNodeIDs, send_count[1], ZOLTAN_ID_MPI_TYPE, i, id_tag + 1, MPI_COMM_WORLD);
                    free(send_mesh[i].myGlobalNodeIDs);

                    MPI_Send(send_mesh[i].neighborIndex, send_count[1] + 1, MPI_INT, i, id_tag + 2, MPI_COMM_WORLD);
                    free(send_mesh[i].neighborIndex);

                    if (send_count[2] > 0){
                        MPI_Send(send_mesh[i].myPinGIDs, send_count[2], ZOLTAN_ID_MPI_TYPE, i, id_tag + 3, MPI_COMM_WORLD);
                        free(send_mesh[i].myPinGIDs);
                    }
                }
            }

            free(send_mesh);

            /* signal all procs it is OK to go on */
            ack = 0;
            for (int i=1; i < rNumProc; ++i){
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

            memset(rMesh, 0, sizeof(ZoltanMesh));
            rMesh->numGlobalElements = numGlobalHyperVertices;
            rMesh->numGlobalNodes = numGlobalHyperEdges;
            rMesh->numGlobalPins = numGlobalHyperPins;

            rMesh->numLocalElements = send_count[0];
            rMesh->numLocalNodes = send_count[1];
            rMesh->numLocalPins = send_count[2];

            if (send_count[0] > 0){
                rMesh->myGlobalElementIDs = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * send_count[0]);
            }

            if (send_count[1] > 0){
                rMesh->myGlobalNodeIDs = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * send_count[1]);
                rMesh->neighborIndex = (int *)malloc(sizeof(int) * (send_count[1] + 1));

                if (send_count[2] > 0){
                    rMesh->myPinGIDs = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * send_count[2]);
                }
            }

            MPI_Send(&ack, 1, MPI_INT, 0, ack_tag, MPI_COMM_WORLD);

            if (send_count[0] > 0){
                MPI_Recv(rMesh->myGlobalElementIDs,send_count[0], ZOLTAN_ID_MPI_TYPE, 0, id_tag, MPI_COMM_WORLD, &status);

                if (send_count[1] > 0){
                    MPI_Recv(rMesh->myGlobalNodeIDs,send_count[1], ZOLTAN_ID_MPI_TYPE, 0, id_tag + 1, MPI_COMM_WORLD, &status);
                    MPI_Recv(rMesh->neighborIndex,send_count[1] + 1, MPI_INT, 0, id_tag + 2, MPI_COMM_WORLD, &status);

                    if (send_count[2] > 0){
                        MPI_Recv(rMesh->myPinGIDs,send_count[2], ZOLTAN_ID_MPI_TYPE, 0, id_tag + 3, MPI_COMM_WORLD, &status);
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

};
