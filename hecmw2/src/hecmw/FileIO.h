/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/FileIO.h
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
#include "FileWriterChunk.h"
#include "CommonStd.h"
#include "TypeDef.h"
#include "FileReaderChunk.h"
#include "FileReaderCnt.h"
#include "FileWriterAVS.h"
#include "FileWriterVTK.h"
#include "FileWriterUNS.h"

namespace FileIO
{
#ifndef FILE_IO_HH_4D6193F1_360F_4f01_A219_725098D1C2FC
#define FILE_IO_HH_4D6193F1_360F_4f01_A219_725098D1C2FC
class CFileIO
{
private:
    CFileIO();
public:
    static CFileIO* Instance() {
        static CFileIO  file_io;
        return &file_io;
    }
    virtual ~CFileIO();
private:
    CFileReaderChunk moReader;
    CFileWriterChunk moWriter;
    CFileReaderCnt   moCntReader;
    ofstream moRltStream;
    uiint mnSolutionType;
public:
    void setFactory(pmw::CMeshFactory *pFactory);
    void setLogger(Utility::CLogger *pLogger);
    void setSolutionType(const uiint& nSolutionType);
    void markingFstrStyle();
    bool ReadCntFile();
    string& getMeshFileBaseName() {
        return moCntReader.getMeshFileBaseName();
    }
    bool Read_fstr_CntFile(string& ctrlname);

    string& getFstr_MeshFileName();
    string& getFstr_ControlFileName();
    string& getFstr_ResultFileName();
    string& getFstr_RestartFileName();
    string& getFstr_VisFileName_Mesh();
    string& getFstr_VisFileName_IN();
    string& getFstr_VisFileName_OUT();
    string& getFstr_PartFileName_IN();
    string& getFstr_PartFileName_OUT();
    string& getFstr_RefineCADFitName();
    string& getFstr_FileName(int nType);

    uiint&  getFstr_RefineNum();
    uiint&  getFstr_RefineType();


    void setBaseName(char base[], const uiint& nLength);
    void setBaseName(const string& base);
    void setFstr_MeshName(char name[], const uiint& nLength);
    void setFstr_ControlName(char name[], const uiint& nLength);
    void setFstr_ResultName(char name[], const uiint& nLength);
    void setFstr_RestartName(char name[], const uiint& nLength);
    void setFstr_VisName_Mesh(char name[], const uiint& nLength);
    void setFstr_VisName_IN(char name[], const uiint& nLength);
    void setFstr_VisName_OUT(char name[], const uiint& nLength);
    void setFstr_PartName_IN(char name[], const uiint& nLength);
    void setFstr_PartName_OUT(char name[], const uiint& nLength);
    void setFstr_CADFitFileName(char name[], const uiint& nLength);
    void setFstr_FileName(char name[], const uiint& nLength, int nType);

    void setFstr_RefineNum(const uiint& nRefineNum);
    void setFstr_RefineType(const uiint& nRefineType);

    void ReadFile(string filename, bool bBinary);
    void WriteFile_Debug(string filename, const uiint& nNumOfLevel);

    bool ReadAlgebraBlock(const uiint& nStep, string filename, bool bBinary);
    uiint  getNumOfLevel();
    uiint  getNumOfEquation();
    uiint  getNumOfParts();
    uiint& getEquationDOF(const uiint& ieq, const uiint& ipart);
    bool ReadResBlock(const uiint& nStep, string filename, bool bBinary);
    void WriteResFile(const uiint& nStep, string filename, bool bBinary);

    void PrintResult_Start(const uiint& nStep, string filename, bool bBinary);
    void PrintResult(const uiint& width, string format, vector<void*>& param);
    void PrintResult_End();
    void WriteAVS_Basis(string filename, const uiint& ieq, const uiint& iMesh, const uiint& nLevel);
    void recAVS_Label(const uiint& iMesh, char* cLabel, char* cUnit, const uiint& nNumOfDOF);
    void recAVS_Variable(const uiint& iMesh, const uiint& nNumOfNode, char* cLabel, double* pvValue);
    void WriteAVS_FEM(string& filename, const uiint& iMesh, const uiint& nLevel);

    void recVTK_Label(const uiint& iMesh, char* cLabel, char* cUnit, const uiint& nNumOfDOF);
    void recVTK_Variable(const uiint& iMesh, const uiint& nNumOfNode, char* cLabel, double* pvValue);
    void WriteVTK_FEM_Rank0(string& filename, const uiint& iMesh, const uiint& nLevel, string& basename, const uiint& nNumOfProcs);
    void WriteVTK_FEM(string& filename, const uiint& iMesh, const uiint& nLevel);

    void recUNS_Label(const uiint& iMesh, char* cLabel, char* cUnit, const uiint& nNumOfDOF);
    void recUNS_Variable(const uiint& iMesh, const uiint& nNumOfNode, char* cLabel, double* pvValue);
    void WriteUNS_FEM(string& filename, const uiint& iMesh, const uiint& nLevel);
};
#endif  //FILE_IO_HH_4D6193F1_360F_4f01_A219_725098D1C2FC
}

