/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   FileWriterChunk.h
|
|                     Written by T.Takeda,    2010/06/01
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
#include "FileWriterContactMesh.h"
#include "FileWriterCommMesh2.h"
namespace FileIO{
class CFileWriterChunk{
public:
    CFileWriterChunk();
    virtual ~CFileWriterChunk();
protected:
    vector<CFileWriter*> mvWriter;
public:
    void Write(string& filename, const uint& numOfLevel);
};
}
#endif	/* _FILEWRITERCHUNK_H */
