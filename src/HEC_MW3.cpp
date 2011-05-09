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
    //mpGMGModel = new CGMGModel();
    mpGMGModel= CGMGModel::Instance();

    //mpFactory = new CMeshFactory();
    mpFactory = CMeshFactory::Instance();
    mpFactory->setGMGModel(mpGMGModel);

    //mpFileIO = new FileIO::CFileIO();
    mpFileIO = FileIO::CFileIO::Instance();
    mb_file = false;// ファイル読み込み後？

    // Factory => FileIO
    mpFileIO->setFactory(mpFactory);

    mpLogger = Utility::CLogger::Instance();
}
// Destructor
//
CMWMain::~CMWMain(void)
{
    // ↓これらは,シングルトンなのでExeの破棄のときに破棄される
    //
    //    delete mpGMGModel;
    //    delete mpFactory;
    //    delete mpFileIO;
    //    delete mpLogger;
}

// MPI, Logger
//
int CMWMain::Initialize()
{
    //-- このLogger設定はテスト用 --
    //
    mpLogger->setMode(Utility::LoggerMode::MWDebug);
    mpLogger->setProperty(Utility::LoggerMode::MWDebug, Utility::LoggerDevice::Disk);
    mpLogger->setProperty(Utility::LoggerMode::Debug, Utility::LoggerDevice::Display);
    mpLogger->setProperty(Utility::LoggerMode::Error, Utility::LoggerDevice::Display);
    mpLogger->setProperty(Utility::LoggerMode::Warn, Utility::LoggerDevice::Display);
    mpLogger->setProperty(Utility::LoggerMode::Info, Utility::LoggerDevice::Display);
    //-- テスト用途 end --


    mpLogger->initializeLogFile();
    mpLogger->InfoDisplay();
    mpLogger->Info(Utility::LoggerMode::Info," HEC_MW3 Initialized");
    
    return 1;
}

//  
//
int CMWMain::Finalize()
{
    mpLogger->Info(Utility::LoggerMode::Info," HEC_MW3 Finalized");
    mpLogger->finalizeLogFile();//ログファイルをclose
    return 1;
}


// FileIO内からFactoryコール、
// 1.Factory内でAssyModel(マルチグリッド階層数分の生成)
// 2.Mesh(Level=0)を生成
//
int CMWMain::FileRead(string& file_name)
{
    mpLogger->Info(Utility::LoggerMode::Info," HEC_MW3 FileRead");

    mpFileIO->ReadFile(file_name);
    mb_file = true;

    return 1;
}

int CMWMain::FileWrite(string& file_name)// ,const uint& nmgLevel)
{
    mpLogger->Info(Utility::LoggerMode::Info," HEC_MW3 FileWrite");

    uint nLevel= mpFactory->getMGLevel();
    mpFileIO->WriteFile(file_name, nLevel+1);//階層数+1 <= 出力数

    return 1;
}

// マルチグリッドのための,Refine(prolongation)メッシュ生成
//
int CMWMain::Refine()
{
    if(mb_file){
        mpFactory->refineMesh();// MultiGridデータの生成 (ファイル読み込み後)
        return 1;
    }else{
        mpLogger->Info(Utility::LoggerMode::Error," not read");
        return 0;
    }
    return 1;
}





















