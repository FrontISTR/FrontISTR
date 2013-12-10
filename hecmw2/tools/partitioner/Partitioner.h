/*
 ----------------------------------------------------------
|
| Software Name :part Ver 0.1 beta
|
|
|                     Written by T.Takeda,    2012/06/01
|                                Y.Sato,      2012/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include "HEC_MW3.h"
extern "C" {
#include "func.h"
#include <metis.h>
}

using namespace pmw;
using namespace std;

class CPartitioner
{
private:
    CMW *pMW;

    // metis data structure
    int ne, nn;
    idxtype* elmnts;
    int etype, numflag, nparts;
    int edgecut;
    idxtype* epart;
    idxtype* npart;
    uiint numGlobalCommMesh;
    vector<vector<uiint> > npartTable;
    vector<vector<uiint> > epartTable;
    vector<uiint> numNodesTable;
    vector<uiint> numElementsTable;
    vector<uiint> globalComm;
    vector<vector<uiint> > commFaceNumTable;
    vector<vector<vector<uiint> > > commNodeTable;
    vector<vector<vector<uiint> > > commNodeIndexTable;
    vector<vector<uiint> > commFaceElementTable;
    vector<vector<uiint> > commFaceIndexTable;
    vector<vector<uiint> > commFaceElementVertTable;
    vector<vector<uiint> > commFaceIndexVertTable;
    vector<vector<uiint> > commFaceRegionTable;
    vector<vector<uiint> > commFaceTypeTable; //point, beam, beam2, triangle, triangle2, quad, quad2
    uiint iMeshMax;
public:
    CPartitioner(CMW *pMW, uiint n);
    ~CPartitioner();

    void prepare();
    void partition();
    void post();
    void write(string prefix);
    void writeUCD(string filename);
    void printReport();
    uiint getCommMeshIX(uiint myrank, uiint rank);
};
