/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/FileWriter.h
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
#ifndef _FILEWRITER_H_6e7d31b0_0854_4c88_af18_efc72922f25e
#define	_FILEWRITER_H_6e7d31b0_0854_4c88_af18_efc72922f25e
#include "HEC_MPI.h"
#include "CommonFile.h"
#include "GMGModel.h"
#include "SolutionType.h"
#include "ElementType.h"
namespace FileIO
{
class CFileWriter
{
public:
    CFileWriter();
    virtual ~CFileWriter();
protected:
    pmw::CGMGModel *mpGMGModel;
    uiint mnSolutionType;
public:
    void setSolutionType(const uiint& nSolutionType);
    virtual void WriteDebug(ofstream& ofs, const uiint& mgLevel)=0;
};
}
#endif	/* _FILEWRITER_H */
