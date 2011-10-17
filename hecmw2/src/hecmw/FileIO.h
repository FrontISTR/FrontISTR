/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.0beta
|
|   ../src/FileIO.h
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
#ifndef FILE_IO_HH_4D6193F1_360F_4f01_A219_725098D1C2FC
#define FILE_IO_HH_4D6193F1_360F_4f01_A219_725098D1C2FC
#include "CommonStd.h"
#include "TypeDef.h"
#include "FileReaderChunk.h"
#include "FileWriterChunk.h"
#include "FileReaderCnt.h"
#include "FileWriterAVS.h"
namespace FileIO{
class CFileIO
{
private:
    CFileIO();
public:
    static CFileIO* Instance(){
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
    string& getMeshFileBaseName(){return moCntReader.getMeshFileBaseName();}
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
    string& getFstr_FileName(int nType);
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
    void setFstr_FileName(char name[], const uiint& nLength, int nType);
    void ReadFile(string filename, bool bBinary);
    void WriteFile_Debug(string filename, const uiint& nNumOfLevel);
    bool ReadAlgebraBlock(const uiint& nStep, string filename, bool bBinary);
    uiint  getNumOfEquation();              
    uiint& getEquationDOF(const uiint& ieq); 
    bool ReadResBlock(const uiint& nStep, string filename, bool bBinary);    
    void WriteResFile(const uiint& nStep, string filename, bool bBinary);
    void PrintResult_Start(const uiint& nStep, string filename, bool bBinary);
    void PrintResult(const uiint& width, char* format, vector<void*>& param);
    void PrintResult_End();
    void WriteAVS_Basis(string filename, const uiint& iMesh, const uiint& nLevel);
    void recAVS_Label(const uiint& iMesh, char* cLabel, char* cUnit, const uiint& nNumOfDOF);
    void recAVS_Variable(const uiint& iMesh, const uiint& nNumOfNode, char* cLabel, double* pvValue);
    void WriteAVS_FEM(string& filename, const uiint& iMesh, const uiint& nLevel);
};
}
#endif
