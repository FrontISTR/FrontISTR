//
//  FileReaderChank.cpp
//
//			2009.05.01
//			2008.12.09
//			k.Takeda
#include "FileReaderChunk.h"
#include "FileReaderContactMesh.h"
using namespace FileIO;

// Constructor
//
CFileReaderChunk::CFileReaderChunk()
{
    // --
    // テスト用(プロセス０専用)
    //  cntファイル名:"mw3.cnt" 2009.09.22
    // --
    // msCntFileName= "mw3.cnt";

    // --
    // FrontISTRの全体制御ファイル
    // --
    // msCntFileName= "hecmw_ctrl.dat";

    
    mpLogger = Utility::CLogger::Instance();

    mvReader.reserve(21);

    mvReader.push_back(new CFileReaderNode());
    mvReader.push_back(new CFileReaderElement());
    mvReader.push_back(new CFileReaderAssyModel);
    //// mvReader.push_back(new CFileReaderRefine);

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

    
    mpAlgebraReader = new CFileReaderAlgebra;// resファイルの線形方程式ブロック(各方程式のDOF)

    mb_fstr= false;//標準スタイル(リスタート拡張子)

}
// Destructor
//
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

// MW3メッシュ・ファイル
// --
void CFileReaderChunk::Read(string filename, bool bBinary)
{
    char c_Line[BUFFERLENGTH];
    string s_Line;

    mpLogger->Info(Utility::LoggerMode::MWDebug,"Input Mesh Filename => ",filename);//debug

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
        // ASCII
        while(!ifs.eof()){
            ifs.getline(c_Line,sizeof(c_Line),'\n');
            s_Line = c_Line;

            for(iBlock=0; iBlock < nNumOfBlock; iBlock++){
                    mvReader[iBlock]->Read(ifs, s_Line);
            };
        };
    }else if(ifs && bBinary){
        bool bEndian(false);
        // Binary
        while(!ifs.eof()){
            //エンディアン・ブロック 有無(先頭に存在しなければ中断)
            if(!bEndian) bEndian= pBinCheck->Read_bin(ifs);
            if(!bEndian) break;

            if(bEndian){
                for(iBlock=0; iBlock < nNumOfBlock; iBlock++)
                    mvReader[iBlock]->Read_bin(ifs);//其々,ifsを全てチェック
            }
            break;
        };
    }else{
        mpLogger->Info(Utility::LoggerMode::Error, "Mesh_file not found, filename => ", filename);
    }
    

    ifs.close();
}

// リスタートの拡張子の付け方管理のマーキング
void CFileReaderChunk::markingFstrStyle()
{
    mb_fstr=true;
}
// MW3リスタート・ファイル(Algebraブロック)
//
bool CFileReaderChunk::ReadAlgebra(const uiint& nStep, string filename, bool bBinary)
{
    bool bState(false);
    
    stringstream ss;
    ss << nStep;
    
    string sFileName;
    if(mb_fstr){
        sFileName= filename + "." + ss.str();//ステップ番号を追加
    }else{
        sFileName= filename + "." + ss.str() + ".res";//"."+ステップ番号+拡張子
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
//
// Algebraブロック：方程式の個数
// 
uiint CFileReaderChunk::getNumOfEquation()
{
    return mpAlgebraReader->getNumOfEquation();
}
//
// Algebraブロック：各方程式のDOF
//
uiint& CFileReaderChunk::getEquationDOF(const uiint& ieq)
{
    return mpAlgebraReader->getEquationDOF(ieq);
}

// MW3リスタート・ファイル(Resブロック)
// --
bool CFileReaderChunk::ReadRes(const uiint& nStep, string filename, bool bBinary)
{
    bool bState(false);
    
    stringstream ss;
    ss << nStep;
    
    string sFileName;
    if(mb_fstr){
        sFileName= filename + "." + ss.str();//ステップ番号を追加
    }else{
        sFileName= filename + "." + ss.str() + ".res";//"."+ステップ番号+拡張子
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

// --
// テスト (mw3.cntファイル)
// --
bool CFileReaderChunk::ReadCnt()
{
    bool bSuccess(false);
    ifstream ifs;
    string s_Line;

    ifs.open("mw3.cnt", ios::in);
    if(ifs){
        while(!ifs.eof()){
            //ifs.getline(c_Line,sizeof(c_Line),'\n');
            getline(ifs, s_Line);

            bSuccess = mpCntReader->Read(ifs, s_Line);
        };
    }else{
        mpLogger->Info(Utility::LoggerMode::Info, "MW3 manage filename");//main()からのベース名セットを利用して、ファイル名を決定
    }
    ifs.close();

    return bSuccess;
}

// --
// FrontISTR全体制御ファイル
// --
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
//// --
//// cntファイル名をフルパス名に変更
//// --
//void CFileReaderChunk::setCntPath(string& cntpath)
//{
//    string sFullPathName;
//
//    sFullPathName = cntpath + msCntFileName;
//    msCntFileName = sFullPathName;
//}









