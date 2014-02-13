/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.3 beta
|
|   ./src/ScalarNode.h
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
#ifndef _SCALARNODE_H_ffd42ed1_dfca_4d02_891c_908db139b697
#define	_SCALARNODE_H_ffd42ed1_dfca_4d02_891c_908db139b697
#include "Node.h"
namespace pmw
{
class CScalarNode:public CNode
{
public:
    CScalarNode();
    virtual ~CScalarNode();
private:
    static uiint mnType;
protected:
    uiint mnNumOfDOF;
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
#endif	/* _SCALARNODE_H */
