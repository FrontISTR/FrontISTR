/* 
 * File:   FileReaderCnt.h
 * Author: ktakeda
 *
 * Created on 2009/09/22, 16:39
 */
#include "CommonFile.h"
#include "FileReader.h"
#include "FileReaderBinCheck.h" //入力ファイルとシステムのエンディアン相違判定

#include "FileBlockNameMW2.h"// hecmw_ctrl.dat タグ

namespace FileIO{
#ifndef _FILEREADERCNT_H
#define	_FILEREADERCNT_H
class CFileReaderCnt:public CFileReader{
public:
    CFileReaderCnt();
    virtual ~CFileReaderCnt();

private:
    // MW3
    string msMeshFileBaseName;// メッシュデータのベースネーム(一般)

    // fstr
    string ms_fstrMsh;    // fstr メッシュファイルベース名 + PE番号(内部付与)
    bool mb_fstrDistType; // fstr メッシュファイル(解析用) true=分散 false=単一 : 分散ならばPE#はつける.

    string ms_fstrCnt;    // fstr コントロールファイル名(そのまま)
    string ms_fstrResult; // fstr リザルトファイル名 + PE番号(内部付与) + ステップ番号(内部付与)
    string ms_fstrRestart;// fstr リスタートファイル名 + PE番号(内部付与)

    string ms_fstrVisMesh;// fstr ビジュアライズ(入力メッシュ名) + PE番号(内部付与)
    string ms_fstrVisIn;  // fstr ビジュアライズ(入力リザルト名) + PE番号(内部付与) + ステップ番号(内部付与)
    string ms_fstrVisOut; // fstr ビジュアライズ(出力ファイル名) + ステップ番号(内部付与) + 拡張子(内部付与)
    
    string ms_fstrPartIn; // fstr メッシュファイル(パーティショナー) 入力
    string ms_fstrPartOut;// fstr メッシュファイル(パーティショナー) 出力 : PE#をつける

    // fstr NAME check
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
    //--
    //mw3.cnt(テストファイル)
    //--
    virtual bool Read(ifstream& ifs, string& sLine);
    virtual bool Read_bin(ifstream& ifs);
    string& getMeshFileBaseName(){ return msMeshFileBaseName;}

    //--
    //hecmw_ctrl.dat(FrontISTR全体制御ファイル)
    //--
    bool Read_fstr_ctrl_file(ifstream& ifs);

    //--
    // FrontISTR固有 : fstr関連ファイル名 提供
    //--
    string& getFstr_MeshFileName(){ return ms_fstrMsh;}
    string& getFstr_ControlFileName() { return ms_fstrCnt;}
    string& getFstr_ResultFileName() { return ms_fstrResult;}
    string& getFstr_RestartFileName() { return ms_fstrRestart;}

    string& getFstr_VisFileName_Mesh(){ return ms_fstrVisMesh;}
    string& getFstr_VisFileName_IN() { return ms_fstrVisIn;}
    string& getFstr_VisFileName_OUT() { return ms_fstrVisOut;}

    string& getFstr_PartFileName_IN() { return ms_fstrPartIn;}
    string& getFstr_PartFileName_OUT(){ return ms_fstrPartOut;}

    
    //--
    // MW3固有 (ファイルベース名)
    //--
    void setBaseName(char base[], const uiint& nLength);//MPI経由でのファイルベース名セット(rank!=0):mw3.cntファイルのデータ
    void setBaseName(const string& base);               //main関数からファイルベース名セット(一般)

    //--
    // FrontISTR固有 : MPI経由でfstr関連ファイル名をセット(rank!=0)
    //--
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



