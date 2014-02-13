/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/FileReaderChunk.h
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
#include "FileReaderVersion.h"
#include "FileIOManage.h"

//#include <boost/lexical_cast.hpp>

#include "HEC_MPI.h"

namespace FileIO
{
#ifndef FILE_READER_HH_BBAE6774_525B_49ec_8315_9FD9F5052070
#define FILE_READER_HH_BBAE6774_525B_49ec_8315_9FD9F5052070
class CFileReaderChunk
{
public:
    CFileReaderChunk();
    CFileReaderChunk(pmw::CMeshFactory *pFactory);
    virtual ~CFileReaderChunk();
private:
    vector<CFileReader*> mvReader;
    CFileReaderCnt     *mpCntReader;    //FrontISTR全体制御ファイル(FrontISTRタイプ専用)
    CFileReaderAlgebra *mpAlgebraReader;//リスタート時の線形方程式情報 読み込み

    Utility::CLogger *mpLogger;
    bool mb_fstr;

    CFileIOManage *mpFileManage;
    void setFileVersion(string sVer);
    void version_Check(CFileReader* pReader, bool& bVerFlg, bool& bBlockFlg, bool& bPrevOtherBlock);
public:
    void setCntReader(CFileReaderCnt* pReader) {
        mpCntReader= pReader;
    }
    bool ReadCnt();
    bool Read_fstr(string& ctrlname);
    void Read(string filename, bool bBinary);
    void markingFstrStyle();

    //--
    // リスタート
    //--
    bool ReadAlgebra(const uiint& nStep, string filename, bool bBinary);
    uiint  getNumOfLevel();
    uiint  getNumOfEquation();
    uiint  getNumOfParts();
    uiint& getEquationDOF(const uiint& ieq, const uiint& ipart);
    bool ReadRes(const uiint& nStep, string filename, bool bBinary);

    void setFactory(pmw::CMeshFactory* pFactory);
    void setLogger(Utility::CLogger *pLogger) {
        mpLogger = pLogger;
    }
};
#endif //FILE_READER_HH_BBAE6774_525B_49ec_8315_9FD9F5052070
}

