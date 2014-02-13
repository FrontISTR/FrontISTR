/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/FileReaderChunk.cpp
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
#include "FileReaderChunk.h"
#include "FileReaderContactMesh.h"
using namespace FileIO;
CFileReaderChunk::CFileReaderChunk()
{
    mpLogger = Utility::CLogger::Instance();
    mvReader.reserve(21);
    mvReader.push_back(new CFileReaderNode());
    mvReader.push_back(new CFileReaderElement());
    mvReader.push_back(new CFileReaderAssyModel);
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
    mvReader.push_back(new CFileReaderVersion);

    mpAlgebraReader = new CFileReaderAlgebra;

    mb_fstr= false;

    mpFileManage = CFileIOManage::Instance();
}
CFileReaderChunk::~CFileReaderChunk()
{
    for_each(mvReader.begin(), mvReader.end(), pmw::DeleteObject());
}
void CFileReaderChunk::setFactory(pmw::CMeshFactory *pFactory)
{
    uiint i;
    for(i=0; i<mvReader.size(); i++) {
        mvReader[i]->setFactory(pFactory);
    };
}

void CFileReaderChunk::setFileVersion(string sVer)
{
    // ファイルから取得したバージョン番号
    double dVersion;
    stringstream ss( sVer );
    ss >> dVersion;

    // MW3の現行ファイル・バージョン番号
    string sCrrVer = FileBlockName::Version_Num();
    double dCrrVer;
    ss.clear();
    ss.str( sCrrVer );
    ss >> dCrrVer;

    if(dVersion > dCrrVer) {
        mpLogger->Info(Utility::LoggerMode::Error, "invalid mesh version");
        exit(0);//------------------------------------------------------ exit(0)
    }

    mpFileManage->setFileVersion(dVersion);
}
void CFileReaderChunk::version_Check(CFileReader* pReader, bool& bVerFlg, bool& bBlockFlg, bool& bPrevOtherBlock)
{
    // RTTI(もどき) & ブロック通過判定
    if( pReader->Name()=="FileReaderVersion" && bBlockFlg && !bPrevOtherBlock) {
        bVerFlg=true;
        CFileReaderVersion *pFileVer= dynamic_cast<CFileReaderVersion*>(pReader);
        if(pFileVer) {
            string sVerNum = FileBlockName::Version_Num();
            //ファイル・バージョン チェック
            if( sVerNum == pFileVer->getVersion() ) {
                mpLogger->Info( Utility::LoggerMode::Info, "mesh ver "+pFileVer->getVersion() );
                setFileVersion( pFileVer->getVersion() );
            } else {
                mpLogger->Info( Utility::LoggerMode::Error, "no current mesh version "+pFileVer->getVersion() );
                setFileVersion( pFileVer->getVersion() );
            }
        }
    }

    // Versionブロックの前に他のブロックが入ってる
    if(bBlockFlg && !bVerFlg && !bPrevOtherBlock) {
        bPrevOtherBlock = true;
        mpLogger->Info(Utility::LoggerMode::Warn, "could not find version block.");
        setFileVersion( FileBlockName::Default_Version_Num() );
    }

    // 他のブロックが来た後でVersion登場 => exit(0)
    if(bVerFlg && bPrevOtherBlock) {
        mpLogger->Info(Utility::LoggerMode::Error, "The first block, need \"version\" block");
        exit(0);//--------------------------------------------------------------------exit(0)
    }
}
void CFileReaderChunk::Read(string filename, bool bBinary)
{
    pmw::CHecMPI* pMPI= pmw::CHecMPI::Instance();

    char c_Line[BUFFERLENGTH];
    string s_Line;
    mpLogger->Info(Utility::LoggerMode::MWDebug,"Input Mesh Filename => ",filename);

    ifstream ifs;
    if(bBinary) {
        ifs.open(filename.c_str(), ios::in|ios::binary);//---------------- open
    } else {
        ifs.open(filename.c_str(), ios::in);//---------------------------- open
    }
    CFileReaderBinCheck *pBinCheck = CFileReaderBinCheck::Instance();

    bool  bVerFlg(false),bBlockFlg(false),bPrevOtherBlock(false);
    uiint iBlock;
    uiint nNumOfBlock= mvReader.size();
    if(ifs && !bBinary) {
        //--
        // ascii
        //--
        while(!ifs.eof()) {
            ifs.getline(c_Line,sizeof(c_Line),'\n');
            s_Line = c_Line;
            for(iBlock=0; iBlock < nNumOfBlock; iBlock++) {

                bBlockFlg = mvReader[iBlock]->Read(ifs, s_Line);//-------------------------- each reader(ascii)

                //////debug : 読み込みBlock表示
                ////if(bBlockFlg){
                ////    uiint rank= pMPI->getRank();
                ////    string sMssg= "FileBlock Name: " + mvReader[iBlock]->Name() + " rank:";
                ////    vector<void*> vParam;
                ////    vParam.push_back((void*)sMssg.c_str()); vParam.push_back(&rank);
                ////    mpLogger->Info( Utility::LoggerMode::Info, "%s%u", vParam );
                ////}

                version_Check(mvReader[iBlock], bVerFlg, bBlockFlg, bPrevOtherBlock);
                bBlockFlg=false;
            };
        };// while end

    } else if(ifs && bBinary) {
        //--
        // binary
        //--
        bool bEndian(false);
        while(!ifs.eof()) {
            if(!bEndian) bEndian= pBinCheck->Read_bin(ifs);
            if(!bEndian) break;
            if(bEndian) {
                for(iBlock=0; iBlock < nNumOfBlock; iBlock++) {

                    bBlockFlg = mvReader[iBlock]->Read_bin(ifs);//----------------------------- each reader(bin)

                    version_Check(mvReader[iBlock], bVerFlg, bBlockFlg, bPrevOtherBlock);
                    bBlockFlg=false;
                };
            }
        };// while end

    } else {
        mpLogger->Info(Utility::LoggerMode::Error, "mesh file not found, filename => ", filename);
    }

    ifs.close();//--------------------------------------------------------close
}

void CFileReaderChunk::markingFstrStyle()
{
    mb_fstr=true;
}
//--
// リスタートデータ
//--
bool CFileReaderChunk::ReadAlgebra(const uiint& nStep, string filename, bool bBinary)
{
    bool bState(false);
    stringstream ss;
    ss << nStep;
    string sFileName;
    if(mb_fstr) {
        sFileName= filename + "." + ss.str();
    } else {
        sFileName= filename + "." + ss.str() + ".res";
    }
    ifstream ifs;
    if(bBinary) {
        ifs.open(sFileName.c_str(), ios::in|ios::binary);
    } else {
        ifs.open(sFileName.c_str(), ios::in);
    }
    string sLine;
    if(ifs) {
        if(!bBinary) {
            while(!ifs.eof()) {
                getline(ifs, sLine);
                mpAlgebraReader->Read(ifs, sLine);
            };
        } else {
            CFileReaderBinCheck *pBinCheck= CFileReaderBinCheck::Instance();
            pBinCheck->Read_bin(ifs);
            mpAlgebraReader->Read_bin(ifs);
        }
        bState=true;
    } else {
        mpLogger->Info(Utility::LoggerMode::Warn, "Res_file not found(Algebra), filename => ", sFileName);
        bState=false;
    }
    ifs.close();
    return bState;
}
uiint CFileReaderChunk::getNumOfLevel()
{
    return mpAlgebraReader->getNumOfLevel();
}
uiint CFileReaderChunk::getNumOfEquation()
{
    return mpAlgebraReader->getNumOfEquation();
}
uiint CFileReaderChunk::getNumOfParts()
{
    return mpAlgebraReader->getNumOfParts();
}
uiint& CFileReaderChunk::getEquationDOF(const uiint& ieq, const uiint& ipart)
{
    return mpAlgebraReader->getEquationDOF(ieq, ipart);
}
bool CFileReaderChunk::ReadRes(const uiint& nStep, string filename, bool bBinary)
{
    bool bState(false);
    stringstream ss;
    ss << nStep;
    string sFileName;
    if(mb_fstr) {
        sFileName= filename + "." + ss.str();
    } else {
        sFileName= filename + "." + ss.str() + ".res";
    }
    ifstream ifs;
    if(bBinary) {
        ifs.open(sFileName.c_str(), ios::in|ios::binary);
    } else {
        ifs.open(sFileName.c_str(), ios::in);
    }
    CFileReaderRes oReaderRes;
    string s_Line;
    if(ifs) {
        if(!bBinary) {
            while(!ifs.eof()) {
                getline(ifs, s_Line);
                oReaderRes.Read(ifs, s_Line);
            };
        } else {
            oReaderRes.Read_bin(ifs);
        }
        bState=true;
    } else {
        mpLogger->Info(Utility::LoggerMode::Warn, "Res_file not found(Res), filename => ", sFileName);
        bState=false;
    }
    ifs.close();
    return bState;
}
bool CFileReaderChunk::ReadCnt()//FrontISTR全体制御ファイル(FrontISTR専用)
{
    bool bSuccess(false);
    ifstream ifs;
    string s_Line;
    ifs.open("mw3.cnt", ios::in);
    if(ifs) {
        while(!ifs.eof()) {
            getline(ifs, s_Line);
            bSuccess = mpCntReader->Read(ifs, s_Line);
        };
    } else {
        mpLogger->Info(Utility::LoggerMode::Info, "MW3 manage filename");
    }
    ifs.close();
    return bSuccess;
}
bool CFileReaderChunk::Read_fstr(string& ctrlname)
{
    if( 0==ctrlname.length() ) ctrlname = "hecmw_ctrl.dat";
    ifstream ifs;
    ifs.open(ctrlname.c_str(), ios::in);
    if(ifs) {
        bool bCheck = mpCntReader->Read_fstr_ctrl_file(ifs);
        return bCheck;
    } else {
        mpLogger->Info(Utility::LoggerMode::Error, "hecmw_ctrl not found");
        return false;
    }
    ifs.close();
}
