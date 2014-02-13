/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/FileReaderBoundaryFaceMesh.cpp
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
#include "FileReaderBoundaryFaceMesh.h"
using namespace FileIO;

CFileReaderBoundaryFaceMesh::CFileReaderBoundaryFaceMesh()
{
    ;
}
CFileReaderBoundaryFaceMesh::~CFileReaderBoundaryFaceMesh()
{
    ;
}
string CFileReaderBoundaryFaceMesh::Name()
{
    return  "FileReaderBoundaryFaceMesh";
}

bool CFileReaderBoundaryFaceMesh::Read(ifstream& ifs, string& sLine)
{
    uiint mgLevel(0);
    uiint bnd_id, bnd_type, mesh_id, numOfBoundary, numOfDOF;
    vuint vDOF;
    string s_bnd_type, s_bnd_name;

    istringstream iss;
    if( TagCheck(sLine, FileBlockName::StartBoundaryFaceMesh()) ) {
        while(!ifs.eof()) {
            sLine = getLine(ifs);
            if(TagCheck(sLine, FileBlockName::EndBoundaryFaceMesh()) ) break;
            iss.clear();
            iss.str(sLine);
            iss >> mesh_id >> numOfBoundary;
            mpFactory->reserveBoundaryFaceMesh(mgLevel, mesh_id, numOfBoundary);

            uiint ibound;
            for(ibound=0; ibound < numOfBoundary; ibound++) {
                sLine= getLine(ifs);
                iss.clear();
                iss.str(sLine);

                iss >> bnd_id >> s_bnd_type >> s_bnd_name >> numOfDOF;

                bnd_type= IntBndType(s_bnd_type);
                vDOF.clear();
                vDOF.resize(numOfDOF);
                for(uiint i=0; i < numOfDOF; i++) {
                    iss >> vDOF[i];
                }

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
                mpFactory->GeneBoundaryFaceMesh(mgLevel, mesh_id, bnd_id, bnd_type, s_bnd_name, numOfDOF, vDOF, mStrForm);
            };
        };
        return true;
    } else {
        return false;
    }
}
bool CFileReaderBoundaryFaceMesh::Read_bin(ifstream& ifs)
{
    CFileReaderBinCheck *pBinCheck= CFileReaderBinCheck::Instance();
    bool bOrder= pBinCheck->isByteOrderSwap();
    bool b32, bCheck;
    string sClassName("FileReaderBoundaryFaceMesh");
    if( !Check_IntSize(b32, bCheck, sClassName) ) return false;
    char cHead='B';
    if( !TagCheck_Bin(ifs, bCheck, cHead, FileBlockName::StartBoundaryFaceMesh(), FileBlockName::BoundaryFaceMesh_Len())) return false;
    uiint mgLevel(0);
    uiint bnd_id, bnd_type, mesh_id, nNumOfBoundary, nNumOfDOF;
    vuint vDOF;
    string s_bnd_type, s_bnd_name;
    map<uiint,string> mStrForm;//数式<dof,数式>

    while(!ifs.eof()) {
        if( CFileReader::Check_End(ifs) )  break;
        ifs.read((char*)&mesh_id, sizeof(uiint));
        if(bOrder) pBinCheck->ByteOrderSwap(mesh_id);
        ifs.read((char*)&nNumOfBoundary, sizeof(uiint));
        if(bOrder) pBinCheck->ByteOrderSwap(nNumOfBoundary);
        mpFactory->reserveBoundaryFaceMesh(mgLevel, mesh_id, nNumOfBoundary);
        uiint ibound;
        for(ibound=0; ibound < nNumOfBoundary; ibound++) {
            ifs.read((char*)&bnd_id, sizeof(uiint));
            if(bOrder) pBinCheck->ByteOrderSwap(bnd_id);
            Read_BndType(ifs, s_bnd_type);
            Read_AnyName(ifs, s_bnd_name);
            ifs.read((char*)&nNumOfDOF, sizeof(uiint));
            if(bOrder) pBinCheck->ByteOrderSwap(nNumOfDOF);
            bnd_type= IntBndType(s_bnd_type);
            vDOF.clear();
            vDOF.resize(nNumOfDOF);
            for(uiint i=0; i < nNumOfDOF; i++) {
                ifs.read((char*)&vDOF[i], sizeof(uiint));
                if(bOrder) pBinCheck->ByteOrderSwap(vDOF[i]);
            }
            //
            //TODO:数式 読み込み
            //
            mpFactory->GeneBoundaryFaceMesh(mgLevel, mesh_id, bnd_id, bnd_type, s_bnd_name, nNumOfDOF, vDOF, mStrForm);
        };
    };
    return true;
}
