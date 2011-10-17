/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.0beta
|
|   ../src/FileReaderChunk.h
|
|                     Written by T.Takeda,    2011/06/01
|                                Y.Sato       2011/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#ifndef FILE_READER_HH_BBAE6774_525B_49ec_8315_9FD9F5052070
#define FILE_READER_HH_BBAE6774_525B_49ec_8315_9FD9F5052070
#include "TypeDef.h"
#include "Logger.h"
#include "FileReaderNode.h"
#include "FileReaderElement.h"
#include "FileReaderAssyModel.h"
#include "FileReaderMaterial.h"
#include "FileReaderCommMesh.h"
#include "FileReaderCommNode.h"
#include "FileReaderCommElement.h"
#include "FileReaderBoundaryNode.h"
#include "FileReaderBoundaryFace.h"
#include "FileReaderBoundaryVolume.h"
#include "FileReaderBoundaryEdge.h"
#include "FileReaderBoundaryNodeMesh.h"
#include "FileReaderBoundaryFaceMesh.h"
#include "FileReaderBoundaryVolumeMesh.h"
#include "FileReaderBoundaryEdgeMesh.h"
#include "FileReaderCnt.h"
#include "FileReaderContactMesh.h"
#include "FileReaderCommMesh2.h"
#include "FileReaderCommFace.h"
#include "FileReaderCommNode_CM2.h"
#include "FileReaderElementGroup.h"
#include "FileReaderElementGroupEntity.h"
#include "FileReaderAlgebra.h"
#include "FileReaderRes.h"    
#include "FileReaderBinCheck.h"
namespace FileIO{
class CFileReaderChunk
{
public:
    CFileReaderChunk();
    CFileReaderChunk(pmw::CMeshFactory *pFactory);
    virtual ~CFileReaderChunk();
private:
    vector<CFileReader*> mvReader;
    CFileReaderCnt     *mpCntReader;    
    CFileReaderAlgebra *mpAlgebraReader;
    Utility::CLogger *mpLogger;
    bool mb_fstr;
public:
    void setCntReader(CFileReaderCnt* pReader){ mpCntReader= pReader;}
    bool ReadCnt();
    bool Read_fstr(string& ctrlname);
    void Read(string filename, bool bBinary);
    void markingFstrStyle();
    bool ReadAlgebra(const uiint& nStep, string filename, bool bBinary);
    uiint  getNumOfEquation();              
    uiint& getEquationDOF(const uiint& ieq);
    bool ReadRes(const uiint& nStep, string filename, bool bBinary);   
    void setFactory(pmw::CMeshFactory* pFactory);
    void setLogger(Utility::CLogger *pLogger){ mpLogger = pLogger;}
};
}
#endif
