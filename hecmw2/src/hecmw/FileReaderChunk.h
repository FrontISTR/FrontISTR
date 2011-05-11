//
//  FileReaderChank.h
//  HEC_MW3 -> each block Reader
//
//			2009.09.22
//			2008.12.08
//			k.Takeda

/*
        << old HEC_MW >>
	<<  FileTagName.h  >>

	<< cnt file >>
	!CONTROL	 
	!MESH	     
	!MESH GROUP	 
	!RESTART
	!RESULT	     

	<< mesh file >>
	!HEADER	       
	!ZERO	       
	!NODE	       
	!ELEMENT	   
	!NGROUP	       
	!EGROUP	       
	!SGROUP        
	!EQUATION	   
	!AMPLITUDE	   
	!SECTION	   
	!MATERIAL	   
	!INITIAL  CONDITION
	!INCLUDE	   
	!CONNECTIVITY  
	!END	       
*/
#ifndef FILE_READER_HH_BBAE6774_525B_49ec_8315_9FD9F5052070
#define FILE_READER_HH_BBAE6774_525B_49ec_8315_9FD9F5052070

//pmwデータ型
#include "TypeDef.h"

//ロガー
#include "Logger.h"

//階層型メッシュデータ構造
#include "FileReaderNode.h"
#include "FileReaderElement.h"
#include "FileReaderAssyModel.h"
////#include "FileReaderRefine.h"

//材料データ
#include "FileReaderMaterial.h"

//要素共有型 通信
#include "FileReaderCommMesh.h"
#include "FileReaderCommNode.h"
#include "FileReaderCommElement.h"

//境界条件オブジェクト
#include "FileReaderBoundaryNode.h"
#include "FileReaderBoundaryFace.h"
#include "FileReaderBoundaryVolume.h"
#include "FileReaderBoundaryEdge.h"
//境界条件メッシュ
#include "FileReaderBoundaryNodeMesh.h"
#include "FileReaderBoundaryFaceMesh.h"
#include "FileReaderBoundaryVolumeMesh.h"
#include "FileReaderBoundaryEdgeMesh.h"

//mw3.cntファイル
#include "FileReaderCnt.h"

//MPCデータ(ContactMesh)
#include "FileReaderContactMesh.h"

//節点共有型 通信
#include "FileReaderCommMesh2.h"
#include "FileReaderCommFace.h"
#include "FileReaderCommNode_CM2.h"


//グループ
#include "FileReaderElementGroup.h"
#include "FileReaderElementGroupEntity.h"

//リスタート
#include "FileReaderAlgebra.h"//Algebraブロック(各線形方程式のDOF)
#include "FileReaderRes.h"    //Resブロック

//ヘッダー(入力ファイルのエンディアン判定)
#include "FileReaderBinCheck.h"

namespace FileIO{
class CFileReaderChunk
{
public:
    CFileReaderChunk();
    CFileReaderChunk(pmw::CMeshFactory *pFactory);
    virtual ~CFileReaderChunk();

private:
    vector<CFileReader*> mvReader;// Mesh本体: Node,Element,..etc
    
    CFileReaderCnt     *mpCntReader;    // cntファイル(hecmw_ctrl.dat) Reader
    CFileReaderAlgebra *mpAlgebraReader;// resファイルの線形方程式ブロック(各方程式のDOF)

    // CFileReaderHeader *mpHeadReader;// エンディアン判定 入力ファイル

    // string msCntFileName; //cntファイル名:テスト= mw3.cnt, FrontISTR = hecmw_ctrl.dat

    Utility::CLogger *mpLogger;

    bool mb_fstr;//拡張子の付け方管理(*resのステップ番号付け方)

public:
    void setCntReader(CFileReaderCnt* pReader){ mpCntReader= pReader;}

    // テスト (mw3.cntファイル)
    bool ReadCnt();
    
    // FrontISTR全体制御ファイル
    bool Read_fstr(string& ctrlname);
    
    // メッシュ
    void Read(string filename, bool bBinary);// メッシュ(*.msh)ファイル

    // リスタートの拡張子の付け方管理のマーキング
    void markingFstrStyle();
    
    // リスタート
    bool ReadAlgebra(const uiint& nStep, string filename, bool bBinary);// Algebraブロック:MW3書式のリスタート(*.res)ファイル
    uiint  getNumOfEquation();              // Algebraブロック:方程式の個数
    uiint& getEquationDOF(const uiint& ieq);// Algebraブロック:各方程式のDOF
    bool ReadRes(const uiint& nStep, string filename, bool bBinary);   // Resブロック    :MW3書式のリスタート(*.res)ファイル

    //// エンディアン判定 入力ファイル
    // bool isLittleEndian_File(){ return mpHeadReader->isLittleEndian_File();}
    // bool isBigEndian_File(){ return mpHeadReader->isBigEndian_File();}

    //// パス
    //void setCntPath(string& cntpath);

    
    void setFactory(pmw::CMeshFactory* pFactory);
    void setLogger(Utility::CLogger *pLogger){ mpLogger = pLogger;}
};
}
#endif
