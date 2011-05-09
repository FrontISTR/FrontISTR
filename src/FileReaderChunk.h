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

#include "CommonStd.h"
#include "TypeDef.h"

#include "Logger.h"

#include "FileReaderNode.h"
#include "FileReaderElement.h"

#include "FileReaderAssyModel.h"
#include "FileReaderRefine.h"

#include "FileReaderMaterial.h"

#include "FileReaderCommMesh.h"
#include "FileReaderCommNode.h"
#include "FileReaderCommElement.h"

#include "FileReaderBoundaryNode.h"
#include "FileReaderBoundaryFace.h"
#include "FileReaderBoundaryVolume.h"

#include "FileReaderCnt.h"//mw3.cntファイル

namespace FileIO{
class CFileReaderChunk
{
public:
    CFileReaderChunk();
    CFileReaderChunk(pmw::CMeshFactory *pFactory);
    virtual ~CFileReaderChunk();

private:
    vector<CFileReader*> mvReader;// Node,Element,etc...Reader
    CFileReaderCnt    *mpCntReader;// cntファイルReader

    string  msCntFileName;//cntファイル名:名称は,固定名

    Utility::CLogger *mpLogger;

public:
    void setCntReader(CFileReaderCnt* pReader){ mpCntReader= pReader;}

    void ReadCnt();// mw3.cnt読み込み
    void Read(string filename);// pMW書式のメッシュファイル読み込み
    
    void setPath(string& filepath);

    void setFactory(pmw::CMeshFactory* pFactory);
    void setLogger(Utility::CLogger *pLogger){ mpLogger = pLogger;}
};
}
#endif
