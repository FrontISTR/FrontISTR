/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.0beta
|
|   ../src/FileReaderCnt.h
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
#include "CommonFile.h"
#include "FileReader.h"
#include "FileReaderBinCheck.h" 
#include "FileBlockNameMW2.h"
namespace FileIO{
#ifndef _FILEREADERCNT_H
#define	_FILEREADERCNT_H
class CFileReaderCnt:public CFileReader{
public:
    CFileReaderCnt();
    virtual ~CFileReaderCnt();
private:
    string msMeshFileBaseName;
    string ms_fstrMsh;    
    bool mb_fstrDistType; 
    string ms_fstrCnt;    
    string ms_fstrResult; 
    string ms_fstrRestart;
    string ms_fstrVisMesh;
    string ms_fstrVisIn;  
    string ms_fstrVisOut; 
    string ms_fstrPartIn; 
    string ms_fstrPartOut;
    short mnMSH, mnCNT, mnResult, mnRestart;
    short mnPartIN, mnPartOUT;
    short mnVisMesh, mnVisIN, mnVisOUT;
    bool fstr_line_check(string& sLine);
    string fstr_tag_split(string& sTag);
    string Read_fstr_mesh(ifstream& ifs, string& sLine);
    string Read_fstr_control(ifstream& ifs, string& sLine);
    string Read_fstr_result(ifstream& ifs, string& sLine);
    string Read_fstr_restart(ifstream& ifs, string& sLine);
public:
    virtual bool Read(ifstream& ifs, string& sLine);
    virtual bool Read_bin(ifstream& ifs);
    string& getMeshFileBaseName(){ return msMeshFileBaseName;}
    bool Read_fstr_ctrl_file(ifstream& ifs);
    string& getFstr_MeshFileName(){ return ms_fstrMsh;}
    string& getFstr_ControlFileName() { return ms_fstrCnt;}
    string& getFstr_ResultFileName() { return ms_fstrResult;}
    string& getFstr_RestartFileName() { return ms_fstrRestart;}
    string& getFstr_VisFileName_Mesh(){ return ms_fstrVisMesh;}
    string& getFstr_VisFileName_IN() { return ms_fstrVisIn;}
    string& getFstr_VisFileName_OUT() { return ms_fstrVisOut;}
    string& getFstr_PartFileName_IN() { return ms_fstrPartIn;}
    string& getFstr_PartFileName_OUT(){ return ms_fstrPartOut;}
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
};
#endif	/* _FILEREADERCNT_H */
}
