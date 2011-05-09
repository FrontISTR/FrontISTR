/* 
 * File:   FileWriterChunk.h
 * Author: ktakeda
 *
 * Created on 2009/07/23, 19:40
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

namespace FileIO{
class CFileWriterChunk{
public:
    CFileWriterChunk();
    virtual ~CFileWriterChunk();

protected:
    vector<CFileWriter*> mvWriter;

public:
    void setSolutionType(const uint& nSolutionType);
    void Write(string& filename, const uint& numOfLevel);
};
}

#endif	/* _FILEWRITERCHUNK_H */

