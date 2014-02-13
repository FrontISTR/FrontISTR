/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/FileReaderBoundaryNodeMesh.cpp
|
|                     Written by T.Takeda,    2013/03/26
|                                Y.Sato,      2013/03/26
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include "FileReaderBoundaryNodeMesh.h"
using namespace FileIO;

CFileReaderBoundaryNodeMesh::CFileReaderBoundaryNodeMesh()
{
    ;
}
CFileReaderBoundaryNodeMesh::~CFileReaderBoundaryNodeMesh()
{
    ;
}
string CFileReaderBoundaryNodeMesh::Name()
{
    return  "FileReaderBoundaryNodeMesh";
}

bool CFileReaderBoundaryNodeMesh::Read(ifstream& ifs, string& sLine)
{
    istringstream iss;
    uiint mgLevel(0), bnd_id, bnd_type, mesh_id, numOfBoundary;
    string s_bnd_type, s_bnd_name;

    if( TagCheck(sLine, FileBlockName::StartBoundaryNodeMesh()) ) {
        while(!ifs.eof()) {
            sLine = getLine(ifs);
            if(TagCheck(sLine, FileBlockName::EndBoundaryNodeMesh()) ) break;

            iss.clear();
            iss.str(sLine.c_str());
            iss >> mesh_id >> numOfBoundary;

            mpFactory->reserveBoundaryNodeMesh(mgLevel, mesh_id, numOfBoundary);

            uiint ibound;
            for(ibound=0; ibound < numOfBoundary; ibound++) {
                sLine= getLine(ifs);
                iss.clear();
                iss.str(sLine.c_str());
                uiint nCount(0);
                while(iss) {
                    if(nCount==0) iss >> bnd_id;
                    if(nCount==1) iss >> s_bnd_type;
                    if(nCount==2) iss >> s_bnd_name;
                    nCount++;
                    if(nCount > 2) break;
                }
                bnd_type= IntBndType(s_bnd_type);

                map<uiint,string> mStrForm;//---- 数式<dof,数式>
                // 数式入力
                if( mpFileManage->isAvailable_NumFrom() ) {
                    uiint nNumOfForm;
                    iss >> nNumOfForm;
                    for(uiint i=0; i < nNumOfForm; i++) {
                        uiint nFormDOF;
                        string sForm;
                        iss >> nFormDOF >> sForm;

                        mStrForm[nFormDOF]= sForm;
                    };
                }
                mpFactory->GeneBoundaryNodeMesh(mgLevel, mesh_id, bnd_id, bnd_type, s_bnd_name, mStrForm);
            };
        };
        return true;
    } else {
        return false;
    }
}
bool CFileReaderBoundaryNodeMesh::Read_bin(ifstream& ifs)
{
    CFileReaderBinCheck *pBinCheck= CFileReaderBinCheck::Instance();
    bool bOrder= pBinCheck->isByteOrderSwap();
    bool b32, bCheck;
    string sClassName("FileReaderBoundaryNodeMesh");
    if( !Check_IntSize(b32, bCheck, sClassName) ) return false;
    char cHead='B';
    if( !TagCheck_Bin(ifs, bCheck, cHead, FileBlockName::StartBoundaryNodeMesh(), FileBlockName::BoundaryNodeMesh_Len())) return false;
    uiint mgLevel(0), bnd_id, bnd_type, mesh_id, nNumOfBoundary;
    string s_bnd_type, s_bnd_name;
    map<uiint,string> mStrForm;//数式<dof,数式>

    while(!ifs.eof()) {
        if( Check_End(ifs) ) break;
        ifs.read((char*)&mesh_id, sizeof(uiint));
        if(bOrder) pBinCheck->ByteOrderSwap(mesh_id);
        ifs.read((char*)&nNumOfBoundary, sizeof(uiint));
        if(bOrder) pBinCheck->ByteOrderSwap(nNumOfBoundary);
        mpFactory->reserveBoundaryNodeMesh(mgLevel, mesh_id, nNumOfBoundary);
        uiint ibound;
        for(ibound=0; ibound < nNumOfBoundary; ibound++) {
            ifs.read((char*)&bnd_id, sizeof(uiint));
            if(bOrder) pBinCheck->ByteOrderSwap(bnd_id);
            Read_BndType(ifs, s_bnd_type);
            Read_AnyName(ifs, s_bnd_name);
            bnd_type= IntBndType(s_bnd_type);
            //
            //TODO:数式 読み込み
            //
            mpFactory->GeneBoundaryNodeMesh(mgLevel, mesh_id, bnd_id, bnd_type, s_bnd_name, mStrForm);
        };
    };
    return true;
}
