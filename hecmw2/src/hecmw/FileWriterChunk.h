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

#include "FileWriterElementGroup.h"

#include "FileWriterRes.h"

#include <cstring> //strlen

namespace FileIO{
class CFileWriterChunk{
public:
    CFileWriterChunk();
    virtual ~CFileWriterChunk();

protected:
    vector<CFileWriter*> mvWriter;

    ofstream m_ofs;//リザルト

    bool mb_fstr;//拡張子の付け方管理(*resのステップ番号付け方)

public:
    void setSolutionType(const uiint& nSolutionType);

    // リスタートの拡張子の付け方管理のマーキング
    void markingFstrStyle();

    // Data Check(*.out)
    void WriteDebug(string& filename, const uiint& nNumOfLevel);//Debug Write (Data Check)

    // リスタート
    void WriteRes(const uiint& nStep, string& filename, bool bBinary);//Res(Restart) File output

    // リザルト
    void PrintResult_Start(const uiint& nStep, string filename, bool bBinary);
    void PrintResult(const uiint& width, char* format, vector<void*>& param);
    void PrintResult_End();
};
}

#endif	/* _FILEWRITERCHUNK_H */

