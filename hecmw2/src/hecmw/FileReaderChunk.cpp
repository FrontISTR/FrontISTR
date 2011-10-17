/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.0beta
|
|   ../src/FileReaderChunk.cpp
|
|                     Written by T.Takeda,    2011/06/01
|                                Y.Sato       2011/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include "HEC_MPI.h"
#include "FileReaderChunk.h"
#include "FileReaderContactMesh.h"
using namespace FileIO;
CFileReaderChunk::CFileReaderChunk()
{
    mpLogger = Utility::CLogger::Instance();
    mvReader.reserve(21);
    mvReader.push_back(new CFileReaderNode());
    mvReader.push_back(new CFileReaderElement());
    mvReader.push_back(new CFileReaderAssyModel);
    mvReader.push_back(new CFileReaderBoundaryNode);
    mvReader.push_back(new CFileReaderBoundaryFace);
    mvReader.push_back(new CFileReaderBoundaryVolume);
    mvReader.push_back(new CFileReaderBoundaryEdge);
    mvReader.push_back(new CFileReaderBoundaryNodeMesh);
    mvReader.push_back(new CFileReaderBoundaryFaceMesh);
    mvReader.push_back(new CFileReaderBoundaryVolumeMesh);
    mvReader.push_back(new CFileReaderBoundaryEdgeMesh);
    mvReader.push_back(new CFileReaderMaterial);
    mvReader.push_back(new CFileReaderCommMesh);
    mvReader.push_back(new CFileReaderCommNode);
    mvReader.push_back(new CFileReaderCommElement);
    mvReader.push_back(new CFileReaderContactMesh);
    mvReader.push_back(new CFileReaderCommMesh2);
    mvReader.push_back(new CFileReaderCommFace);
    mvReader.push_back(new CFileReaderCommNodeCM2);
    mvReader.push_back(new CFileReaderElementGroup);
    mvReader.push_back(new CFileReaderElementGroupEntity);
    mpAlgebraReader = new CFileReaderAlgebra;
    mb_fstr= false;
}
CFileReaderChunk::~CFileReaderChunk()
{
    for_each(mvReader.begin(), mvReader.end(), pmw::DeleteObject());
}
void CFileReaderChunk::setFactory(pmw::CMeshFactory *pFactory)
{
    uiint i;
    for(i=0; i<mvReader.size(); i++){
            mvReader[i]->setFactory(pFactory);
    };
}
void CFileReaderChunk::Read(string filename, bool bBinary)
{
    char c_Line[BUFFERLENGTH];
    string s_Line;
    mpLogger->Info(Utility::LoggerMode::MWDebug,"Input Mesh Filename => ",filename);
    ifstream ifs;
    if(bBinary){
        ifs.open(filename.c_str(), ios::in|ios::binary);
    }else{
        ifs.open(filename.c_str(), ios::in);
    }
    CFileReaderBinCheck *pBinCheck = CFileReaderBinCheck::Instance();
    uiint iBlock;
    uiint nNumOfBlock= mvReader.size();
    if(ifs && !bBinary){
        while(!ifs.eof()){
            ifs.getline(c_Line,sizeof(c_Line),'\n');
            s_Line = c_Line;
            for(iBlock=0; iBlock < nNumOfBlock; iBlock++){
                    mvReader[iBlock]->Read(ifs, s_Line);
            };
        };
    }else if(ifs && bBinary){
        bool bEndian(false);
        while(!ifs.eof()){
            if(!bEndian) bEndian= pBinCheck->Read_bin(ifs);
            if(!bEndian) break;
            if(bEndian){
                for(iBlock=0; iBlock < nNumOfBlock; iBlock++)
                    mvReader[iBlock]->Read_bin(ifs);
            }
            break;
        };
    }else{
        mpLogger->Info(Utility::LoggerMode::Error, "Mesh_file not found, filename => ", filename);
    }
    ifs.close();
}
void CFileReaderChunk::markingFstrStyle()
{
    mb_fstr=true;
}
bool CFileReaderChunk::ReadAlgebra(const uiint& nStep, string filename, bool bBinary)
{
    bool bState(false);
    stringstream ss;
    ss << nStep;
    string sFileName;
    if(mb_fstr){
        sFileName= filename + "." + ss.str();
    }else{
        sFileName= filename + "." + ss.str() + ".res";
    }
    ifstream ifs;
    if(bBinary){
        ifs.open(sFileName.c_str(), ios::in|ios::binary);
    }else{
        ifs.open(sFileName.c_str(), ios::in);
    }
    string sLine;
    if(ifs){
        if(!bBinary){
            while(!ifs.eof()){
                getline(ifs, sLine);
                mpAlgebraReader->Read(ifs, sLine);
            };
        }else{
            CFileReaderBinCheck *pBinCheck= CFileReaderBinCheck::Instance();
            pBinCheck->Read_bin(ifs);
            mpAlgebraReader->Read_bin(ifs);
        }
        bState=true;
    }else{
        mpLogger->Info(Utility::LoggerMode::Warn, "Res_file not found(Algebra), filename => ", sFileName);
        bState=false;
    }
    ifs.close();
    return bState;
}
uiint CFileReaderChunk::getNumOfEquation()
{
    return mpAlgebraReader->getNumOfEquation();
}
uiint& CFileReaderChunk::getEquationDOF(const uiint& ieq)
{
    return mpAlgebraReader->getEquationDOF(ieq);
}
bool CFileReaderChunk::ReadRes(const uiint& nStep, string filename, bool bBinary)
{
    bool bState(false);
    stringstream ss;
    ss << nStep;
    string sFileName;
    if(mb_fstr){
        sFileName= filename + "." + ss.str();
    }else{
        sFileName= filename + "." + ss.str() + ".res";
    }
    ifstream ifs;
    if(bBinary){
        ifs.open(sFileName.c_str(), ios::in|ios::binary);
    }else{
        ifs.open(sFileName.c_str(), ios::in);
    }
    CFileReaderRes oReaderRes;
    string s_Line;
    if(ifs){
        if(!bBinary){
            while(!ifs.eof()){
                getline(ifs, s_Line);
                oReaderRes.Read(ifs, s_Line);
            };
        }else{
            oReaderRes.Read_bin(ifs);
        }
        bState=true;
    }else{
        mpLogger->Info(Utility::LoggerMode::Warn, "Res_file not found(Res), filename => ", sFileName);
        bState=false;
    }
    ifs.close();
    return bState;
}
bool CFileReaderChunk::ReadCnt()
{
    bool bSuccess(false);
    ifstream ifs;
    string s_Line;
    ifs.open("mw3.cnt", ios::in);
    if(ifs){
        while(!ifs.eof()){
            getline(ifs, s_Line);
            bSuccess = mpCntReader->Read(ifs, s_Line);
        };
    }else{
        mpLogger->Info(Utility::LoggerMode::Info, "MW3 manage filename");
    }
    ifs.close();
    return bSuccess;
}
bool CFileReaderChunk::Read_fstr(string& ctrlname)
{
    if( 0==ctrlname.length() ) ctrlname = "hecmw_ctrl.dat";
    ifstream ifs;
    ifs.open(ctrlname.c_str(), ios::in);
    if(ifs){
        bool bCheck = mpCntReader->Read_fstr_ctrl_file(ifs);
        return bCheck;
    }else{
        mpLogger->Info(Utility::LoggerMode::Error, "hecmw_ctrl not found");
        return false;
    }
    ifs.close();
}
