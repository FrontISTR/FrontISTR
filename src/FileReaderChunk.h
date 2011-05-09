//
//  FileReaderChank.h
//  HEC_MW3 -> each block Reader
//
//			2009.03.23
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

namespace FileIO{
class CFileReaderChunk
{
public:
	CFileReaderChunk();
	CFileReaderChunk(pmw::CMeshFactory *pFactory);
	virtual ~CFileReaderChunk();

protected:
	vector<CFileReader*> mvReader;// Node,Element,etc...Reader
        Utility::CLogger *mpLogger;

public:
	void Read(string filename);// pMW

        void setFactory(pmw::CMeshFactory* pFactory);
        void setLogger(Utility::CLogger *pLogger){ mpLogger = pLogger;}
};
}
#endif
