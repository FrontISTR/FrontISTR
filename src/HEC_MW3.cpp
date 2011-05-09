//
//  HEC_MW3.cpp
//
//			2009.04.20
//			2008.11.19
//			k.Takeda
#include "HEC_MW3.h"
using namespace pmw;


// Constructor
//
CMWMain::CMWMain(void)
{
    // GMGModel
    mpGMGModel= CGMGModel::Instance();

    // Factory
    mpFactory = CMeshFactory::Instance();
    mpFactory->setGMGModel(mpGMGModel);

    // FileIO
    mpFileIO = FileIO::CFileIO::Instance();
    mb_file = false;// ファイル読み込み後？

    // Factory => FileIO
    mpFileIO->setFactory(mpFactory);
    
    // MPI
    mpMPI = pmw::CHecMPI::Instance();

    mpLogger = Utility::CLogger::Instance();
}
// Destructor
//
CMWMain::~CMWMain(void)
{
    // ↓これらは,シングルトンなのでExeの破棄のときに破棄される
    //
    //  mpGMGModel;
    //  mpFactory;
    //  mpFileIO;
    //  mpMPI
    //  mpLogger;
}

// MPI, Logger の初期化
//
int CMWMain::Initialize(int argc, char** argv)
{
    // MPI 初期化
    mpMPI->Initialize(argc, argv);
    // rank 取得
    uint rank= mpMPI->getRank();

    //-- このLogger設定はテスト用 --
    //
    mpLogger->setMode(Utility::LoggerMode::MWDebug);

    mpLogger->setProperty(Utility::LoggerMode::MWDebug, Utility::LoggerDevice::Disk);
    mpLogger->setProperty(Utility::LoggerMode::Debug, Utility::LoggerDevice::Display);
    mpLogger->setProperty(Utility::LoggerMode::Error, Utility::LoggerDevice::Display);
    mpLogger->setProperty(Utility::LoggerMode::Warn, Utility::LoggerDevice::Display);
    mpLogger->setProperty(Utility::LoggerMode::Info, Utility::LoggerDevice::Display);
    //-- テスト用途 end --

    // Logger初期化
    mpLogger->initializeLogFile(rank);// ログファイル-オープン
    mpLogger->InfoDisplay();
    mpLogger->Info(Utility::LoggerMode::Info," HEC_MW3 Initialized");
    

    //cntファイルから入力ファイルのベースネームを取得
    mpFileIO->ReadCntFile();

    // ファイルパス
    msInputFileName= mpFileIO->getPathName();
    msOutputFileName=mpFileIO->getPathName();

    // ファイルパスにファイルベース名を追加, アンダースコアを追加
    msInputFileName += mpFileIO->getMeshFileBaseName() + "_";
    msOutputFileName+= mpFileIO->getMeshFileBaseName() + "_";
    
    // 自分のrankを"ファイルベース名"に追加
    msInputFileName  += boost::lexical_cast<string>(rank);
    msOutputFileName += boost::lexical_cast<string>(rank);

    // 拡張子を追加
    msInputFileName  += ".msh";
    msOutputFileName += ".out";

    return 1;
}

//  MPI, Logger 後処理
//
int CMWMain::Finalize()
{
    // MPI
    mpMPI->Finalize();

    // Logger
    mpLogger->Info(Utility::LoggerMode::Info," HEC_MW3 Finalized");
    mpLogger->finalizeLogFile();//ログファイル-クローズ
    return 1;
}


// FileIO内からFactoryコール、
// 1.Factory内でAssyModel(マルチグリッド階層数分の生成)
// 2.Mesh(Level=0)を生成
//
int CMWMain::FileRead()
{
    mpLogger->Info(Utility::LoggerMode::Info," HEC_MW3 FileRead");

    mpFileIO->ReadFile(msInputFileName);
    mb_file = true;

    return 1;
}

int CMWMain::FileWrite()// const uint& nmgLevel
{
    mpLogger->Info(Utility::LoggerMode::Info," HEC_MW3 FileWrite");

    uint nLevel= mpFactory->getMGLevel();
    mpFileIO->WriteFile(msOutputFileName, nLevel+1);//階層数+1 <= 出力数

    return 1;
}

// マルチグリッドのための,Refine(prolongation)メッシュ生成
//
int CMWMain::Refine()
{
    if(mb_file){
        // ファイル読み込み後
        mpFactory->refineMesh();// MultiGridデータの生成(通信Meshも含む)
        return 1;
    }else{
        mpLogger->Info(Utility::LoggerMode::Error," not read");
        return 0;
    }
    return 1;
}





















