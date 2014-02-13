/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/FileWriterChunk.h
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
#ifndef _FILEWRITERCHUNK_H_607f58ad_5a55_440c_9eab_a348cf76863e
#define	_FILEWRITERCHUNK_H_607f58ad_5a55_440c_9eab_a348cf76863e
#include <vector>
#include "FileWriterNode.h"
#include "FileWriterElement.h"
#include "FileWriterBoundaryNode.h"
#include "FileWriterBoundaryFace.h"
#include "FileWriterBoundaryEdge.h"
#include "FileWriterBoundaryVolume.h"
#include "FileWriterContactMesh.h"
#include "FileWriterCommMesh2.h"
#include "FileWriterElementGroup.h"
#include "FileWriterRes.h"
#include <cstring>
namespace FileIO
{
class CFileWriterChunk
{
public:
    CFileWriterChunk();
    virtual ~CFileWriterChunk();
protected:
    vector<CFileWriter*> mvWriter;
    ofstream m_ofs;
    bool mb_fstr;
public:
    void setSolutionType(const uiint& nSolutionType);
    void markingFstrStyle();
    void WriteDebug(string& filename, const uiint& nNumOfLevel);
    void WriteRes(const uiint& nStep, string& filename, bool bBinary);
    void PrintResult_Start(const uiint& nStep, string filename, bool bBinary);
    void PrintResult(const uiint& width, string format, vector<void*>& param);
    void PrintResult_End();
};
}
#endif	/* _FILEWRITERCHUNK_H */
