//
//	FileIO.h
//
//			2008.12.08
//			2008.12.08
//			k.Takeda
#ifndef FILE_IO_HH_4D6193F1_360F_4f01_A219_725098D1C2FC
#define FILE_IO_HH_4D6193F1_360F_4f01_A219_725098D1C2FC

#include "CommonStd.h"
#include "TypeDef.h"

#include "FileReaderChunk.h"
#include "FileWriterChunk.h"

#include "FileReaderCnt.h"//全体制御ファイル hecmw_ctrl

#include "FileWriterAVS.h"//MicroAVS *.inp 出力

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
    //string msCntPathName;// ctrlファイルのパス
    //string msDatPathName;// データファイルパス

    CFileReaderChunk moReader;
    CFileWriterChunk moWriter;

    CFileReaderCnt   moCntReader;
    ofstream moRltStream;

    uiint mnSolutionType;// FEM .or. FVM

public:
    void setFactory(pmw::CMeshFactory *pFactory);
    void setLogger(Utility::CLogger *pLogger);

    void setSolutionType(const uiint& nSolutionType);

    void markingFstrStyle();//リスタートの拡張子の付け方 => fstr仕様

    // mw3.cnt(テスト)からメッシュのベースネームを取得
    bool ReadCntFile();
    string& getMeshFileBaseName(){return moCntReader.getMeshFileBaseName();}

    // hecmw_ctrl 読み込み
    bool Read_fstr_CntFile(string& ctrlname);
    
    // hecmw_ctrl (FrontISTR全体制御ファイル)記述のファイル名を取得
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


    //// 全体ctrlファイルへのパス
    //void setCntPathName(const char* path);
    //string& getCntPathName(){ return msCntPathName;}

    //// データファイルへのパス
    //void setDatPathName(const char* path);
    //string& getDatPathName(){ return msDatPathName;}


    // MW3固有 : MPI経由でのファイルベース名セット(rank!=0)
    void setBaseName(char base[], const uiint& nLength);//mw3.cntファイルのBaseNameをセット(テスト)
    void setBaseName(const string& base);               //main関数からファイルベース名セット(一般)
    
    // fstr関連: MPI経由でファイル名をセット
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
    
    
    // メッシュ
    void ReadFile(string filename, bool bBinary);//メッシュ(*.msh)
    void WriteFile_Debug(string filename, const uiint& nNumOfLevel);//Data_Check(*.out)

    // リスタート
    bool ReadAlgebraBlock(const uiint& nStep, string filename, bool bBinary);//Algebraブロック(リスタートファイル)
    uiint  getNumOfEquation();              //方程式の個数  :Algebraブロック
    uiint& getEquationDOF(const uiint& ieq); //各方程式のDOF :Algebraブロック
    bool ReadResBlock(const uiint& nStep, string filename, bool bBinary);    //Resブロック(リスタートファイル)
    void WriteResFile(const uiint& nStep, string filename, bool bBinary);

    // リザルト
    void PrintResult_Start(const uiint& nStep, string filename, bool bBinary);
    void PrintResult(const uiint& width, char* format, vector<void*>& param);
    void PrintResult_End();

    // MicroAVS *.inp 出力 : 基礎変数   nLevel:出力する階層
    void WriteAVS_Basis(string filename, const uiint& iMesh, const uiint& nLevel);
    void recAVS_Label(const uiint& iMesh, char* cLabel, char* cUnit, const uiint& nNumOfDOF);//出力ラベルの登録
    void recAVS_Variable(const uiint& iMesh, const uiint& nNumOfNode, char* cLabel, double* pvValue);//出力変数の登録
    void WriteAVS_FEM(string& filename, const uiint& iMesh, const uiint& nLevel);//登録変数の出力
};
}
#endif
