/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/FileReaderCnt.h
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
#include "CommonFile.h"
#include "FileReader.h"
#include "FileReaderBinCheck.h"
#include "FileBlockNameMW2.h"

namespace FileIO
{
#ifndef _FILEREADERCNT_H
#define	_FILEREADERCNT_H
//--
// FrontISTR専用 全体制御ファイル 読み込み
//--
class CFileReaderCnt:public CFileReader
{
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
    string ms_fstrCADFitFile;//リファイン用CADFittingFile名称
    uiint  mnRefineNum; //リファイン数(fstr全体制御ファイル)
    uiint  mnRefineType;//リファイン・タイプを表す整数

    // ブロックの出現の有無(整数)
    short mnMSH, mnCNT, mnResult, mnRestart;
    short mnPartIN, mnPartOUT;
    short mnVisMesh, mnVisIN, mnVisOUT;
    short mnRefineBlock;

    bool fstr_line_check(string& sLine);
    string fstr_tag_split(string& sTag);
    string fstr_comma_remove(string& sToken);

    string Read_fstr_mesh(ifstream& ifs, string& sLine);
    string Read_fstr_control(ifstream& ifs, string& sLine);
    string Read_fstr_result(ifstream& ifs, string& sLine);
    string Read_fstr_restart(ifstream& ifs, string& sLine);
    string Read_fstr_refine(ifstream& ifs, string& sLine);
public:
    virtual bool Read(ifstream& ifs, string& sLine);
    virtual bool Read_bin(ifstream& ifs);

    virtual string Name();

    string& getMeshFileBaseName() {
        return msMeshFileBaseName;
    }

    bool Read_fstr_ctrl_file(ifstream& ifs);//------- fstr全体制御ファイル 読み込み

    string& getFstr_MeshFileName() {
        return ms_fstrMsh;
    }
    string& getFstr_ControlFileName() {
        return ms_fstrCnt;
    }
    string& getFstr_ResultFileName() {
        return ms_fstrResult;
    }
    string& getFstr_RestartFileName() {
        return ms_fstrRestart;
    }
    string& getFstr_VisFileName_Mesh() {
        return ms_fstrVisMesh;
    }
    string& getFstr_VisFileName_IN() {
        return ms_fstrVisIn;
    }
    string& getFstr_VisFileName_OUT() {
        return ms_fstrVisOut;
    }
    string& getFstr_PartFileName_IN() {
        return ms_fstrPartIn;
    }
    string& getFstr_PartFileName_OUT() {
        return ms_fstrPartOut;
    }
    string& getFstr_RefineCADFitFileName() {
        return ms_fstrCADFitFile;
    }

    uiint&  getFstr_RefineNum() {
        return mnRefineNum;
    }
    uiint&  getFstr_RefineType() {
        return mnRefineType;
    }


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

    void setFstr_RefineNum(const uiint& nRefineNum);
    void setFstr_RefineType(const uiint& nRefineType);
};
#endif	/* _FILEREADERCNT_H */
}
