/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/FileWriterNode.h
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
#ifndef _FILEWRITERNODE_H_f576c116_dcdf_4d04_bc13_e214549933b4
#define	_FILEWRITERNODE_H_f576c116_dcdf_4d04_bc13_e214549933b4
#include "AssyModel.h"
#include "Mesh.h"
#include "Node.h"
#include "FileWriter.h"
#include "SolutionType.h"
namespace FileIO
{
class CFileWriterNode:public CFileWriter
{
public:
    CFileWriterNode();
    virtual ~CFileWriterNode();
public:
    virtual void WriteDebug(ofstream& ofs, const uiint& mgLevel);
};
}
#endif	/* _FILEWRITERNODE_H */
