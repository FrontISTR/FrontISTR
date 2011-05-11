//
//  FileReaderCnt.cpp
//
//
//
//                  2009.09.22
//                  2009.09.22
//                  k.Takeda
#include "FileReaderCnt.h"
#include "FileBlockName.h"
using namespace FileIO;

// construct & destruct
// --
CFileReaderCnt::CFileReaderCnt()
{
    mb_fstrDistType=false;
}
CFileReaderCnt::~CFileReaderCnt()
{
    ;
}
// MPI経由でのファイルベース名セット(rank!=0):mw3.cnt ファイルのデータ
void CFileReaderCnt::setBaseName(char base[], const uiint& nLength)
{
    msMeshFileBaseName.resize(nLength);

    for(uiint i=0; i < nLength; i++) msMeshFileBaseName[i] = base[i];
}
// main関数からファイルベース名セット(一般)
void CFileReaderCnt::setBaseName(const string& base)
{
    msMeshFileBaseName = base;
}

// MPI経由でfstr関連ファイル名をセット
void CFileReaderCnt::setFstr_MeshName(char name[], const uiint& nLength)
{
    ms_fstrMsh.resize(nLength);
    for(uiint i=0; i < nLength; i++) ms_fstrMsh[i] = name[i];
}
void CFileReaderCnt::setFstr_ControlName(char name[], const uiint& nLength)
{
    ms_fstrCnt.resize(nLength);
    for(uiint i=0; i < nLength; i++) ms_fstrCnt[i] = name[i];
}
void CFileReaderCnt::setFstr_ResultName(char name[], const uiint& nLength)
{
    ms_fstrResult.resize(nLength);
    for(uiint i=0; i < nLength; i++) ms_fstrResult[i] = name[i];
}
void CFileReaderCnt::setFstr_RestartName(char name[], const uiint& nLength)
{
    ms_fstrRestart.resize(nLength);
    for(uiint i=0; i < nLength; i++) ms_fstrRestart[i] = name[i];
}
void CFileReaderCnt::setFstr_VisName_Mesh(char name[], const uiint& nLength)
{
    ms_fstrVisMesh.resize(nLength);
    for(uiint i=0; i < nLength; i++) ms_fstrVisMesh[i] = name[i];
}
void CFileReaderCnt::setFstr_VisName_IN(char name[], const uiint& nLength)
{
    ms_fstrVisIn.resize(nLength);
    for(uiint i=0; i < nLength; i++) ms_fstrVisIn[i] = name[i];
}
void CFileReaderCnt::setFstr_VisName_OUT(char name[], const uiint& nLength)
{
    ms_fstrVisOut.resize(nLength);
    for(uiint i=0; i < nLength; i++) ms_fstrVisOut[i] = name[i];
}
void CFileReaderCnt::setFstr_PartName_IN(char name[], const uiint& nLength)
{
    ms_fstrPartIn.resize(nLength);
    for(uiint i=0; i < nLength; i++) ms_fstrPartIn[i] = name[i];
}
void CFileReaderCnt::setFstr_PartName_OUT(char name[], const uiint& nLength)
{
    ms_fstrPartOut.resize(nLength);
    for(uiint i=0; i < nLength; i++) ms_fstrPartOut[i] = name[i];
}

// --
// mw3.cnt(テストファイル) ＊他のReadメソッドとは,ReaderChunkでの使い方が異なっているので注意.
// --
bool CFileReaderCnt::Read(ifstream& ifs, string& sLine)
{
    // --
    // mw3.cnt (テストファイル)
    // --
    // MeshFileName
    if(TagCheck(sLine, FileBlockName::StartMeshFileName()) ){

        while(true){
            sLine = getLineSt(ifs);
            if(sLine==FileBlockName::EndMeshFileName()) break;

            istringstream iss(sLine.c_str());
            iss >> msMeshFileBaseName;

            mpLogger->Info(Utility::LoggerMode::Debug,"Basefile name => ",msMeshFileBaseName);
        };
        return true;
    }else{
        return false;
    }
}
bool CFileReaderCnt::Read_bin(ifstream& ifs)
{
    mpLogger->Info(Utility::LoggerMode::Error, "invalid method, FileReaderCnt::Read_bin(ifstream& ifs)");
    return false;
}

// --
// hecmw_ctrl.dat ( FrontISTR全体制御ファイル )
// --
bool CFileReaderCnt::Read_fstr_ctrl_file(ifstream& ifs)
{
    // fstr NAME check
    mnMSH=0, mnCNT=0, mnResult=0, mnRestart=0;
    mnPartIN=0, mnPartOUT=0;
    mnVisMesh=0, mnVisIN=0, mnVisOUT=0;

    string sLine;
    
    while(!ifs.eof()){
        sLine = getLine(ifs);

        sLine= Read_fstr_mesh( ifs, sLine );
        sLine= Read_fstr_control( ifs, sLine );
        sLine= Read_fstr_restart( ifs, sLine );
        sLine= Read_fstr_result( ifs, sLine );
    }
    //    cout << "mnMSH " << mnMSH << " mnCNT " << mnCNT << flush;
    //    cout << " mnResutl " << mnResult  << " mnRestart " << mnRestart << flush;
    //    cout << " mnPartIN " << mnPartIN  << " mnPartOUT " << mnPartOUT << flush;
    //    cout << " mnVisMesh " << mnVisMesh << " mnVisIN " << mnVisIN << " mnVisOUT " << mnVisOUT << endl;

    short nTotal= mnMSH + mnCNT + mnResult + mnRestart + mnPartIN + mnPartOUT + mnVisMesh + mnVisIN + mnVisOUT;
    
    if( nTotal > 0 ){
        return true;
    }else{
        return false;
    }
}
// --
// 1行の内容 check
// --
bool CFileReaderCnt::fstr_line_check(string& sLine)
{
    // 何か別のブロックが始まっている
    if( sLine[0]==FileBlockNameMW2::Exclamation() && sLine[1]!=FileBlockNameMW2::Exclamation() ){// !M**,!C**,!R**
        mpLogger->Info(Utility::LoggerMode::Error, "hecmw_ctrl, no filename");
        return false;
    }
    // 空行
    if( sLine.length() == 0){
        return false;
    }
    //  # 記号 コメント
    if( sLine[0]==FileBlockNameMW2::HashMark() ){
        return false;
    }
    // !! 記号 コメント
    if( sLine[0]==FileBlockNameMW2::Exclamation() && sLine[1]==FileBlockNameMW2::Exclamation() ){
        return false;
    }
    return true;
}
// --
// タグを"="で分割
// --
string CFileReaderCnt::fstr_tag_split(string& sTag)
{
    uiint i, nLength= sTag.length();
    for(i=0; i < nLength; i++){
        if(sTag[i]=='=') sTag[i]=' ';
    };
    string sToken;
    istringstream iss(sTag);

    iss >> sToken;//から読み
    iss >> sToken;//パラメータ

    return sToken;
}
// --
// 1.fstr メッシュ・ブロック
// --
string CFileReaderCnt::Read_fstr_mesh(ifstream& ifs, string& sLine)
{
    istringstream iss;
    string sBlock, sName, sType;
    
    iss.clear(); iss.str(sLine);

    short nCount;
    for(nCount=0; nCount < 3; nCount++){
        if(iss){
            if(nCount==0) iss >> sBlock;
            if(nCount==1) iss >> sName;
            if(nCount==2) iss >> sType;
        }
    };
    sName = fstr_tag_split(sName);
    sType = fstr_tag_split(sType);
    
    // !MESH ブロック
    if( TagCheck(sBlock, FileBlockNameMW2::Mesh()) ){

        // Name : fstrMSH
        if( TagCheck(sName, FileBlockNameMW2::Name_Mesh())){

            if( TagCheck(sType, FileBlockNameMW2::Type_Dist()) ){//HECMW-DIST 分散データ
                mb_fstrDistType=true;
            }
            if( TagCheck(sType, FileBlockNameMW2::Type_Entire()) ){//HECMW-ENTIRE 単一領域データ
                mb_fstrDistType=false;
            }
            bool bCheck(false);
            
            while(!bCheck){
                sLine= getLine(ifs);
                bCheck= fstr_line_check(sLine);
                if(bCheck){
                    iss.clear();
                    iss.str(sLine);
                    iss >> ms_fstrMsh;
                    ++mnMSH;
                    mpLogger->Info(Utility::LoggerMode::MWDebug, "hecmw_ctrl fstrMesh ", ms_fstrMsh);
                }
            };
        }
        // Name : part_in (パーティショナー) 入力 : HECMW-ENTIRE 単一領域データ
        if( TagCheck(sName, FileBlockNameMW2::Name_Part_IN()) ){
            if( TagCheck(sType, FileBlockNameMW2::Type_Dist()) ){
                mpLogger->Info(Utility::LoggerMode::Warn, "fstr Part_IN Type is HECMW-ENTIRE");
            }
            bool bCheck(false);
            
            while(!bCheck){
                sLine= getLine(ifs);
                bCheck= fstr_line_check(sLine);
                if(bCheck){
                    iss.clear();
                    iss.str(sLine);
                    iss >> ms_fstrPartIn;
                    ++mnPartIN;
                    mpLogger->Info(Utility::LoggerMode::MWDebug, "hecmw_ctrl part_in ", ms_fstrPartIn);
                }
            };
        }
        // Name : part_in (パーティショナー) 出力 : HECMW-DIST   分散データ
        if( TagCheck(sName, FileBlockNameMW2::Name_Part_OUT()) ){
            if( TagCheck(sType, FileBlockNameMW2::Type_Entire()) ){
                mpLogger->Info(Utility::LoggerMode::Warn, "fstr Part_OUT Type is HECMW-DIST");
            }
            bool bCheck(false);
            
            while(!bCheck){
                sLine= getLine(ifs);
                bCheck= fstr_line_check(sLine);
                if(bCheck){
                    iss.clear();
                    iss.str(sLine);
                    iss >> ms_fstrPartOut;
                    ++mnPartOUT;
                    mpLogger->Info(Utility::LoggerMode::MWDebug, "hecmw_ctrl part_out ", ms_fstrPartOut);
                }
            };
        }
        // Name : mesh
        if( TagCheck(sName, FileBlockNameMW2::Name_VisMesh()) ){
            bool bCheck(false);
            while(!bCheck){
                sLine= getLine(ifs);
                bCheck= fstr_line_check(sLine);
                if(bCheck){
                    iss.clear();
                    iss.str(sLine);
                    iss >> ms_fstrVisMesh;
                    ++mnVisMesh;
                    mpLogger->Info(Utility::LoggerMode::MWDebug, "hecmw_ctrl  visualizer_mesh ", ms_fstrVisMesh);
                }
            };
        }
    }
    return sLine;
}
// --
// 2.fstr コントロール・ブロック
// --
string CFileReaderCnt::Read_fstr_control(ifstream& ifs, string& sLine)
{
    istringstream iss;
    string sBlock, sName, sType="";

    iss.clear(); iss.str(sLine);

    short nCount;
    for(nCount=0; nCount < 3; nCount++){
        if(iss){
            if(nCount==0) iss >> sBlock;
            if(nCount==1) iss >> sName;
        }
    };
    sName = fstr_tag_split(sName);

    // !CONTROLブロック
    if( TagCheck(sBlock, FileBlockNameMW2::Control()) ){

        // Name : fstrCNT
        if( TagCheck( sName, FileBlockNameMW2::Name_Control()) ){
            bool bCheck(false);
            
            while(!bCheck){
                sLine= getLine(ifs);
                bCheck= fstr_line_check(sLine);
                if(bCheck){
                    iss.clear();
                    iss.str(sLine);
                    iss >> ms_fstrCnt;
                    ++mnCNT;
                    mpLogger->Info(Utility::LoggerMode::MWDebug, "hecmw_ctrl fstrCnt ", ms_fstrCnt);
                }
            };

        }else{
            mpLogger->Info(Utility::LoggerMode::Error, "!CONTROL,  mismatch name", sName);
        }
    }
    return sLine;
}
// --
// 3.fstr リザルト・ブロック
// --
string CFileReaderCnt::Read_fstr_result(ifstream& ifs, string& sLine)
{
    istringstream iss;
    string sBlock, sName, sType="", sIO;

    iss.clear(); iss.str(sLine);

    short nCount;
    for(nCount=0; nCount < 3; nCount++){
        if(iss){
            if(nCount==0) iss >> sBlock;
            if(nCount==1) iss >> sName;
            if(nCount==2) iss >> sIO;
        }
    };
    sName = fstr_tag_split(sName);

    // !RESULTブロック
    if( TagCheck(sBlock, FileBlockNameMW2::Result()) ){

        // Name : fstrRES
        if( TagCheck( sName, FileBlockNameMW2::Name_Result()) ){
            bool bCheck(false);

            while(!bCheck){
                sLine= getLine(ifs);
                bCheck= fstr_line_check(sLine);
                if(bCheck){
                    iss.clear();
                    iss.str(sLine);
                    iss >> ms_fstrResult;
                    ++mnResult;
                    mpLogger->Info(Utility::LoggerMode::MWDebug, "hecmw_ctrl fstrRES ", ms_fstrResult);
                }
            };
        }
        // Name : result (visualizer in)
        if( TagCheck( sName, FileBlockNameMW2::Name_VisIn()) ){
            bool bCheck(false);

            while(!bCheck){
                sLine= getLine(ifs);
                bCheck= fstr_line_check(sLine);
                if(bCheck){
                    iss.clear();
                    iss.str(sLine);
                    iss >> ms_fstrVisIn;
                    ++mnVisIN;
                    mpLogger->Info(Utility::LoggerMode::MWDebug, "hecmw_ctrl visualizer_in ", ms_fstrVisIn);
                }
            };
        }
        // Name : vis_out (visualizer out)
        if( TagCheck( sName, FileBlockNameMW2::Name_VisOut()) ){
            bool bCheck(false);

            while(!bCheck){
                sLine= getLine(ifs);
                bCheck= fstr_line_check(sLine);
                if(bCheck){
                    iss.clear();
                    iss.str(sLine);
                    iss >> ms_fstrVisOut;
                    ++mnVisOUT;
                    mpLogger->Info(Utility::LoggerMode::MWDebug, "hecmw_ctrl visualizer_out ", ms_fstrVisOut);
                }
            };
        }
    }
    return sLine;
}
// --
// 4.fstr リスタート・ブロック
// --
string CFileReaderCnt::Read_fstr_restart(ifstream& ifs, string& sLine)
{
    istringstream iss;
    string sBlock, sName, sType="", sIO;

    iss.clear(); iss.str(sLine);

    short nCount;
    for(nCount=0; nCount < 3; nCount++){
        if(iss){
            if(nCount==0) iss >> sBlock;
            if(nCount==1) iss >> sName;
            if(nCount==2) iss >> sIO;
        }
    };
    sName = fstr_tag_split(sName);

    // !RESTARTブロック
    if( TagCheck(sBlock, FileBlockNameMW2::Restart()) ){

        // Name : restart
        if( TagCheck( sName, FileBlockNameMW2::Name_Restart()) ){
            bool bCheck(false);

            while(!bCheck){
                sLine= getLine(ifs);
                bCheck= fstr_line_check(sLine);
                if(bCheck){
                    iss.clear();
                    iss.str(sLine);
                    iss >> ms_fstrRestart;
                    ++mnRestart;
                    mpLogger->Info(Utility::LoggerMode::MWDebug, "hecmw_ctrl restart ", ms_fstrRestart);
                }
            };
        }
    }
    return sLine;
}









