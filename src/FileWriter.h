/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   FileWriter.h
|
|                     Written by T.Takeda,    2010/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#ifndef _FILEWRITER_H_6e7d31b0_0854_4c88_af18_efc72922f25e
#define	_FILEWRITER_H_6e7d31b0_0854_4c88_af18_efc72922f25e
#include "CommonFile.h"
#include "GMGModel.h"
#include "HEC_MPI.h"
namespace FileIO{
class CFileWriter{
public:
    CFileWriter();
    virtual ~CFileWriter();
protected:
    pmw::CGMGModel *mpGMGModel;
public:
    virtual void Write(ofstream& ofs, const uint& mgLevel)=0;
};
}
#endif	/* _FILEWRITER_H */
