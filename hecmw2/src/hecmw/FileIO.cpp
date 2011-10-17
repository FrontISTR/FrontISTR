/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.0beta
|
|   ../src/FileIO.cpp
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
#include "FileIO.h"
using namespace FileIO;
CFileIO::CFileIO()
{
    moReader.setCntReader(&moCntReader);
}
CFileIO::~CFileIO()
{
    ;
}
void CFileIO::setBaseName(char base[], const uiint& nLength)
{
    moCntReader.setBaseName(base, nLength);
}
void CFileIO::setBaseName(const string& base)
{
    moCntReader.setBaseName(base);
}
void CFileIO::setFstr_MeshName(char name[], const uiint& nLength)
{
    moCntReader.setFstr_MeshName(name, nLength);
}
void CFileIO::setFstr_ControlName(char name[], const uiint& nLength)
{
    moCntReader.setFstr_ControlName(name, nLength);
}
void CFileIO::setFstr_ResultName(char name[], const uiint& nLength)
{
    moCntReader.setFstr_ResultName(name, nLength);
}
void CFileIO::setFstr_RestartName(char name[], const uiint& nLength)
{
    moCntReader.setFstr_RestartName(name, nLength);
}
void CFileIO::setFstr_VisName_Mesh(char name[], const uiint& nLength)
{
    moCntReader.setFstr_VisName_Mesh(name, nLength);
}
void CFileIO::setFstr_VisName_IN(char name[], const uiint& nLength)
{
    moCntReader.setFstr_VisName_IN(name, nLength);
}
void CFileIO::setFstr_VisName_OUT(char name[], const uiint& nLength)
{
    moCntReader.setFstr_VisName_OUT(name, nLength);
}
void CFileIO::setFstr_PartName_IN(char name[], const uiint& nLength)
{
    moCntReader.setFstr_PartName_IN(name, nLength);
}
void CFileIO::setFstr_PartName_OUT(char name[], const uiint& nLength)
{
    moCntReader.setFstr_PartName_OUT(name, nLength);
}
void CFileIO::setFstr_FileName(char name[], const uiint& nLength, int nType)
{
    switch(nType){
        case( FileTypeMW2::Mesh ):
            moCntReader.setFstr_MeshName(name, nLength);
            break;
        case( FileTypeMW2::Control ):
            moCntReader.setFstr_ControlName(name, nLength);
            break;
        case( FileTypeMW2::Restart ):
            moCntReader.setFstr_RestartName(name, nLength);
            break;
        case( FileTypeMW2::Result ):
            moCntReader.setFstr_ResultName(name, nLength);
            break;
        case( FileTypeMW2::PartIN ):
            moCntReader.setFstr_PartName_IN(name, nLength);
            break;
        case( FileTypeMW2::PartOUT ):
            moCntReader.setFstr_PartName_OUT(name, nLength);
            break;
        case( FileTypeMW2::VisMesh ):
            moCntReader.setFstr_VisName_Mesh(name, nLength);
            break;
        case( FileTypeMW2::VisIN ):
            moCntReader.setFstr_VisName_IN(name, nLength);
            break;
        case( FileTypeMW2::VisOUT ):
            moCntReader.setFstr_VisName_OUT(name, nLength);
            break;
        default:
            break;
    }
}
void CFileIO::setFactory(pmw::CMeshFactory* pFactory)
{
    moReader.setFactory(pFactory);
}
void CFileIO::setLogger(Utility::CLogger* pLogger)
{
    moReader.setLogger(pLogger);
}
void CFileIO::setSolutionType(const uiint& nSolutionType)
{
    mnSolutionType = nSolutionType;
    moWriter.setSolutionType(nSolutionType);
}
bool CFileIO::ReadCntFile()
{
    bool bSuccess(false);
    bSuccess = moReader.ReadCnt();
    return bSuccess;
}
bool CFileIO::Read_fstr_CntFile(string& ctrlname)
{
    bool bSuccess(false);
    bSuccess = moReader.Read_fstr(ctrlname);
    return bSuccess;
}
string& CFileIO::getFstr_MeshFileName(){ return  moCntReader.getFstr_MeshFileName();}
string& CFileIO::getFstr_ControlFileName(){ return  moCntReader.getFstr_ControlFileName();}
string& CFileIO::getFstr_ResultFileName(){ return  moCntReader.getFstr_ResultFileName();}
string& CFileIO::getFstr_RestartFileName(){ return  moCntReader.getFstr_RestartFileName();}
string& CFileIO::getFstr_VisFileName_Mesh(){ return  moCntReader.getFstr_VisFileName_Mesh();}
string& CFileIO::getFstr_VisFileName_IN(){ return  moCntReader.getFstr_VisFileName_IN();}
string& CFileIO::getFstr_VisFileName_OUT(){ return  moCntReader.getFstr_VisFileName_OUT();}
string& CFileIO::getFstr_PartFileName_IN(){ return  moCntReader.getFstr_PartFileName_IN();}
string& CFileIO::getFstr_PartFileName_OUT(){ return  moCntReader.getFstr_PartFileName_OUT();}
string& CFileIO::getFstr_FileName(int nType)
{
    switch(nType){
        case( FileTypeMW2::Mesh ):
            return moCntReader.getFstr_MeshFileName();
        case( FileTypeMW2::Control ):
            return moCntReader.getFstr_ControlFileName();
        case( FileTypeMW2::Restart ):
            return moCntReader.getFstr_RestartFileName();
        case( FileTypeMW2::Result ):
            return moCntReader.getFstr_ResultFileName();
        case( FileTypeMW2::PartIN ):
            return moCntReader.getFstr_PartFileName_IN();
        case( FileTypeMW2::PartOUT ):
            return moCntReader.getFstr_PartFileName_OUT();
        case( FileTypeMW2::VisMesh):
            return moCntReader.getFstr_VisFileName_Mesh();
        case( FileTypeMW2::VisIN ):
            return moCntReader.getFstr_VisFileName_IN();
        case( FileTypeMW2::VisOUT ):
            return moCntReader.getFstr_VisFileName_OUT();
        default:
            break;
    }
}
void CFileIO::ReadFile(string filename, bool bBinary)
{
    moReader.Read(filename, bBinary);
}
void CFileIO::WriteFile_Debug(string filename, const uiint& nNumOfLevel)
{
    moWriter.WriteDebug(filename, nNumOfLevel);
}
void CFileIO::markingFstrStyle()
{
    moReader.markingFstrStyle();
    moWriter.markingFstrStyle();
}
bool CFileIO::ReadAlgebraBlock(const uiint& nStep, string filename, bool bBinary)
{
    return moReader.ReadAlgebra(nStep, filename, bBinary);
}
uiint CFileIO::getNumOfEquation()
{
    return moReader.getNumOfEquation();
}
uiint& CFileIO::getEquationDOF(const uiint& ieq)
{
    return moReader.getEquationDOF(ieq);
}
bool CFileIO::ReadResBlock(const uiint& nStep, string filename, bool bBinary)
{
    return moReader.ReadRes(nStep, filename, bBinary);
}
void CFileIO::WriteResFile(const uiint& nStep, string filename, bool bBinary)
{
    moWriter.WriteRes(nStep, filename, bBinary);
}
void CFileIO::PrintResult_Start(const uiint& nStep, string filename, bool bBinary)
{
    moWriter.PrintResult_Start(nStep, filename, bBinary);
}
void CFileIO::PrintResult(const uiint& width, char* format, vector<void*>& param)
{
    moWriter.PrintResult(width, format, param);
}
void CFileIO::PrintResult_End()
{
    moWriter.PrintResult_End();
}
void CFileIO::WriteAVS_Basis(string filename, const uiint& iMesh, const uiint& nLevel)
{
    CFileWriterAVS* pAVS = CFileWriterAVS::Instance();
    ofstream ofs;
    ofs.open(filename.c_str());
    pAVS->WriteBasis(ofs, iMesh, nLevel);
    ofs.close();
}
void CFileIO::recAVS_Label(const uiint& iMesh, char* cLabel, char* cUnit, const uiint& nNumOfDOF)
{
    CFileWriterAVS *pAVS= CFileWriterAVS::Instance();
    pAVS->recAVS_Label(iMesh, cLabel, cUnit, nNumOfDOF);
}
void CFileIO::recAVS_Variable(const uiint& iMesh, const uiint& nNumOfNode, char* cLabel, double* pvValue)
{
    CFileWriterAVS *pAVS= CFileWriterAVS::Instance();
    pAVS->recAVS_Variable(iMesh, nNumOfNode, cLabel, pvValue);
}
void CFileIO::WriteAVS_FEM(string& filename, const uiint& iMesh, const uiint& nLevel)
{
    CFileWriterAVS *pAVS= CFileWriterAVS::Instance();
    ofstream ofs;
    ofs.open(filename.c_str());
    pAVS->WriteFEM(ofs, iMesh, nLevel);
    ofs.close();
}
