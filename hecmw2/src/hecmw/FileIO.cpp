//
//  FileIO.cpp
//              2008.12.08
//              2008.12.08
//              k.Takeda

#include "FileIO.h"
using namespace FileIO;

// construct && destruct
// --
CFileIO::CFileIO()
{
    // CntReaderのセット
    moReader.setCntReader(&moCntReader);
}
CFileIO::~CFileIO()
{
    ;
}


//// pathの設定(ctrlファイルパス)
////
//void CFileIO::setCntPathName(const char* cpath)
//{
//    msCntPathName= cpath;
//}
//void CFileIO::setDatPathName(const char* path)
//{
//    msDatPathName= path;
//}


// MPI経由でのファイルベース名セット(rank!=0):mw3.cntのデータ
void CFileIO::setBaseName(char base[], const uiint& nLength)
{
    moCntReader.setBaseName(base, nLength);
}
// main関数からファイルベース名セット(一般)
void CFileIO::setBaseName(const string& base)
{
    moCntReader.setBaseName(base);
}



// MPI経由でfstr関連ファイル名をセット
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
//
// MW Initialize の BCast で使用
//
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

// Factory を 各Reader に設置
//
void CFileIO::setFactory(pmw::CMeshFactory* pFactory)
{
    moReader.setFactory(pFactory);
}


// Logger を 各Reader に設置
//
void CFileIO::setLogger(Utility::CLogger* pLogger)
{
    moReader.setLogger(pLogger);
}

// Solution Type
void CFileIO::setSolutionType(const uiint& nSolutionType)
{
    mnSolutionType = nSolutionType;

    moWriter.setSolutionType(nSolutionType);
}

// --
// mw3.cnt(テスト)を読み込む -> ベースネームを取得
// --
bool CFileIO::ReadCntFile()//テスト
{
    bool bSuccess(false);

    //moReader.setCntPath(msCntPathName);
    bSuccess = moReader.ReadCnt();

    return bSuccess;
}
// --
// hecmw_ctrl.dat(FrontISTR全体制御ファイル)の読み込み
// --
bool CFileIO::Read_fstr_CntFile(string& ctrlname)
{
    bool bSuccess(false);

    //moReader.setCntPath(msCntPathName);
    bSuccess = moReader.Read_fstr(ctrlname);
    
    return bSuccess;
}
// --
// hecmw_ctrl (FrontISTR全体制御ファイル)記述のファイル名を取得
// --
string& CFileIO::getFstr_MeshFileName(){ return  moCntReader.getFstr_MeshFileName();}
string& CFileIO::getFstr_ControlFileName(){ return  moCntReader.getFstr_ControlFileName();}
string& CFileIO::getFstr_ResultFileName(){ return  moCntReader.getFstr_ResultFileName();}
string& CFileIO::getFstr_RestartFileName(){ return  moCntReader.getFstr_RestartFileName();}
string& CFileIO::getFstr_VisFileName_Mesh(){ return  moCntReader.getFstr_VisFileName_Mesh();}
string& CFileIO::getFstr_VisFileName_IN(){ return  moCntReader.getFstr_VisFileName_IN();}
string& CFileIO::getFstr_VisFileName_OUT(){ return  moCntReader.getFstr_VisFileName_OUT();}
string& CFileIO::getFstr_PartFileName_IN(){ return  moCntReader.getFstr_PartFileName_IN();}
string& CFileIO::getFstr_PartFileName_OUT(){ return  moCntReader.getFstr_PartFileName_OUT();}

//
// MW Initialize の BCastで使用
//
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


// メッシュ
//
void CFileIO::ReadFile(string filename, bool bBinary)
{
    moReader.Read(filename, bBinary);
}
// Data_Check(Debug)ファイル
//
void CFileIO::WriteFile_Debug(string filename, const uiint& nNumOfLevel)
{
    moWriter.WriteDebug(filename, nNumOfLevel);
}

// ----
// Resファイル入力
// ----
//
// リスタートファイル名の拡張子の付け方 => fstr仕様
//
void CFileIO::markingFstrStyle()
{
    moReader.markingFstrStyle();
    moWriter.markingFstrStyle();
}
//
// Algebraブロック(リスタートファイル)
//
bool CFileIO::ReadAlgebraBlock(const uiint& nStep, string filename, bool bBinary)
{
    return moReader.ReadAlgebra(nStep, filename, bBinary);
}
//
// 各方程式のDOF :Algebraブロック
//
uiint CFileIO::getNumOfEquation()
{
    return moReader.getNumOfEquation();
}
//
// 各方程式のDOF :Algebraブロック
//
uiint& CFileIO::getEquationDOF(const uiint& ieq)
{
    return moReader.getEquationDOF(ieq);
}
//
// Resブロック(リスタートファイル)
//
bool CFileIO::ReadResBlock(const uiint& nStep, string filename, bool bBinary)
{
    return moReader.ReadRes(nStep, filename, bBinary);
}
// ----
// Resファイル出力
// ----
void CFileIO::WriteResFile(const uiint& nStep, string filename, bool bBinary)
{
    moWriter.WriteRes(nStep, filename, bBinary);
}

//
// Resultブロック
//
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

// ----
// MicroAVS *.inp 出力 
// ----
//
// 基礎変数   nLevel:出力する階層
//
void CFileIO::WriteAVS_Basis(string filename, const uiint& iMesh, const uiint& nLevel)
{
    CFileWriterAVS* pAVS = CFileWriterAVS::Instance();

    ofstream ofs;
    ofs.open(filename.c_str());

    pAVS->WriteBasis(ofs, iMesh, nLevel);

    ofs.close();
}
//
//出力ラベルの登録
//
void CFileIO::recAVS_Label(const uiint& iMesh, char* cLabel, char* cUnit, const uiint& nNumOfDOF)
{
    CFileWriterAVS *pAVS= CFileWriterAVS::Instance();
    pAVS->recAVS_Label(iMesh, cLabel, cUnit, nNumOfDOF);
}
//
// 変数登録(FEM)
//
void CFileIO::recAVS_Variable(const uiint& iMesh, const uiint& nNumOfNode, char* cLabel, double* pvValue)
{
    CFileWriterAVS *pAVS= CFileWriterAVS::Instance();
    pAVS->recAVS_Variable(iMesh, nNumOfNode, cLabel, pvValue);
}
//
// 登録変数の出力（FEM)
//
void CFileIO::WriteAVS_FEM(string& filename, const uiint& iMesh, const uiint& nLevel)
{
    CFileWriterAVS *pAVS= CFileWriterAVS::Instance();

    ofstream ofs;
    ofs.open(filename.c_str());

    pAVS->WriteFEM(ofs, iMesh, nLevel);

    ofs.close();
}


