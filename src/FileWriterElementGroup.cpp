//
// FileWriterElementGroup.cpp
//
//              2010.10.26
//              k.Takeda
#include "FileWriterElementGroup.h"
using namespace FileIO;


CFileWriterElementGroup::CFileWriterElementGroup()
{
    ;
}
CFileWriterElementGroup::~CFileWriterElementGroup()
{
    ;
}


void CFileWriterElementGroup::Write(ofstream& ofs, const uint& mgLevel)
{
    uint iMesh, iGrp;
    uint nNumOfMesh, nNumOfGrp;

    pmw::CAssyModel *pAssyModel = mpGMGModel->getAssyModel(mgLevel);

    nNumOfMesh = pAssyModel->getNumOfMesh();

    for(iMesh=0; iMesh < nNumOfMesh; iMesh++){
        pmw::CMesh *pMesh = pAssyModel->getMesh(iMesh);

        nNumOfGrp = pMesh->getNumOfElemGrp();
        pmw::CElementGroup *pElemGrp;
        for(iGrp=0; iGrp < nNumOfGrp; iGrp++){
            pElemGrp = pMesh->getElemGrpIX(iGrp);

            ofs << "Element Group ID : " << pElemGrp->getID() << "  Name : " << pElemGrp->getName() << endl;

            uint nNumOfElemID = pElemGrp->getNumOfElementID();
            uint i;
            for(i=0; i < nNumOfElemID; i++){
                uint ElemID = pElemGrp->getElementID(i);

                ofs << ElemID << " ";

                if( (i+1)%5==0 ) ofs << endl;
            };
            ofs << endl;
        };
    };
}










