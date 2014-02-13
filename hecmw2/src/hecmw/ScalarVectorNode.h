/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/ScalarVectorNode.h
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
#ifndef _SCALARVECTORNODE_H_7ba23d5d_34f7_425c_be80_3ba32504ab2f
#define	_SCALARVECTORNODE_H_7ba23d5d_34f7_425c_be80_3ba32504ab2f
#include "Node.h"
namespace pmw
{
class CScalarVectorNode:public CNode
{
public:
    CScalarVectorNode();
    virtual ~CScalarVectorNode();
private:
    static uiint mnType;
protected:
    uiint mnNumOfSDOF;
    uiint mnNumOfVDOF;
public:
    virtual uiint& getType() {
        return mnType;
    }
    virtual void setScalarDOF(const uiint& nNDOF);
    virtual void setVectorDOF(const uiint& nNDOF);
    virtual uiint& getScalarDOF();
    virtual uiint& getVectorDOF();
    virtual uiint getTotalDOF();
};
}
#endif	/* _SCALARVECTORNODE_H */
