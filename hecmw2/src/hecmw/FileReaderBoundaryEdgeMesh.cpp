/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/FileReaderBoundaryEdgeMesh.cpp
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
#include "FileReaderBoundaryEdgeMesh.h"
using namespace FileIO;

CFileReaderBoundaryEdgeMesh::CFileReaderBoundaryEdgeMesh()
{
    ;
}
CFileReaderBoundaryEdgeMesh::~CFileReaderBoundaryEdgeMesh()
{
    ;
}
string CFileReaderBoundaryEdgeMesh::Name()
{
    return  "FileReaderBoundaryEdgeMesh";
}

bool CFileReaderBoundaryEdgeMesh::Read(ifstream& ifs, string& sLine)
{
    uiint mgLevel(0);
    uiint bnd_id, bnd_type, mesh_id, numOfBoundary, numOfDOF;
    vuint vDOF;
    string s_bnd_type, s_bnd_name;

    istringstream iss;
    if( TagCheck(sLine, FileBlockName::StartBoundaryEdgeMesh()) ) {
        while(!ifs.eof()) {
            sLine = getLine(ifs);
            if(TagCheck(sLine, FileBlockName::EndBoundaryEdgeMesh()) ) break;
            iss.clear();
            iss.str(sLine);
            iss >> mesh_id >> numOfBoundary;

            mpFactory->reserveBoundaryEdgeMesh(mgLevel, mesh_id, numOfBoundary);

            uiint ibound;
            for(ibound=0; ibound < numOfBoundary; ibound++) {
                sLine= getLine(ifs);
                iss.clear();
                iss.str(sLine);
                iss >> bnd_id >> s_bnd_type >> s_bnd_name >> numOfDOF;

                bnd_type = IntBndType(s_bnd_type);
                vDOF.clear();
                vDOF.resize(numOfDOF);
                for(uiint i=0; i < numOfDOF; i++) {
                    iss >> vDOF[i];
                };

                map<uiint,string> mStrForm;//---- 数式<dof,数式>
                // 数式入力
                if( mpFileManage->isAvailable_NumFrom() ) {
                    uiint nNumOfForm;
                    iss >> nNumOfForm;
                    for(uiint i=0; i < nNumOfForm; i++) {
                        uiint nFormDOF;
                        string sForm;
                        iss >> nFormDOF >> sForm;

                        ////debug
                        //cout << "FileReaderBoundaryEdgeMesh::Read,  sForm:" << sForm << endl;

                        mStrForm[nFormDOF]= sForm;
                    };
                }
                mpFactory->GeneBoundaryEdgeMesh(mgLevel, mesh_id, bnd_id, bnd_type, s_bnd_name, numOfDOF, vDOF, mStrForm);
            };
        };
        return true;
    } else {
        return false;
    }
}
bool CFileReaderBoundaryEdgeMesh::Read_bin(ifstream& ifs)
{
    CFileReaderBinCheck *pBinCheck= CFileReaderBinCheck::Instance();
    bool bOrder= pBinCheck->isByteOrderSwap();
    bool b32, bCheck;
    string sClassName("FileReaderBoundaryEdgeMesh");
    if( !Check_IntSize(b32, bCheck, sClassName) ) return false;
    char cHead='B';
    if( !TagCheck_Bin(ifs, bCheck, cHead, FileBlockName::StartBoundaryEdgeMesh(), FileBlockName::BoundaryEdgeMesh_Len())) return false;
    uiint mgLevel(0);
    uiint bnd_id, bnd_type, mesh_id, nNumOfBoundary, nNumOfDOF;
    vuint vDOF;
    string s_bnd_type, s_bnd_name;
    map<uiint,string> mStrForm;//数式<dof,数式>

    while(!ifs.eof()) {
        ifs.read((char*)&mesh_id, sizeof(uiint));
        if(bOrder) pBinCheck->ByteOrderSwap(mesh_id);
        ifs.read((char*)&nNumOfBoundary, sizeof(uiint));
        if(bOrder) pBinCheck->ByteOrderSwap(nNumOfBoundary);

        mpFactory->reserveBoundaryEdgeMesh(mgLevel, mesh_id, nNumOfBoundary);

        for(uiint ibound=0; ibound < nNumOfBoundary; ibound++) {
            ifs.read((char*)&bnd_id, sizeof(uiint));
            if(bOrder) pBinCheck->ByteOrderSwap(bnd_id);
            Read_BndType(ifs, s_bnd_type);
            Read_AnyName(ifs, s_bnd_name);
            ifs.read((char*)&nNumOfDOF, sizeof(uiint));
            if(bOrder) pBinCheck->ByteOrderSwap(nNumOfDOF);
            bnd_type = IntBndType(s_bnd_type);
            vDOF.clear();
            vDOF.resize(nNumOfDOF);
            for(uiint i=0; i < nNumOfDOF; i++) {
                ifs.read((char*)&vDOF[i], sizeof(uiint));
                if(bOrder) pBinCheck->ByteOrderSwap(vDOF[i]);
            };
            uiint nNumOfForm,nFormDOF;
            ifs.read((char*)&nNumOfForm, sizeof(uiint));
            if(bOrder) pBinCheck->ByteOrderSwap(nNumOfForm);
            for(uiint i=0; i < nNumOfForm; i++) {
                ifs.read((char*)&nFormDOF,sizeof(uiint));
                if(bOrder) pBinCheck->ByteOrderSwap(nFormDOF);
                //バイナリの場合文字列数
                //
                // TODO: 数式
                //
            };
            mpFactory->GeneBoundaryEdgeMesh(mgLevel, mesh_id, bnd_id, bnd_type, s_bnd_name, nNumOfDOF, vDOF, mStrForm);
        };
    };
    return true;
}
