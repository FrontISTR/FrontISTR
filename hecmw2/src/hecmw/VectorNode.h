/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.0beta
|
|   ../src/VectorNode.h
|
|                     Written by T.Takeda,    2011/06/01
|                                Y.Sato       2011/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#ifndef _VECTORNODE_H_dc47e55f_2af8_4606_a36b_820ed5128260
#define	_VECTORNODE_H_dc47e55f_2af8_4606_a36b_820ed5128260
#include "Node.h"
namespace pmw{
class CVectorNode:public CNode{
public:
    CVectorNode();
    virtual ~CVectorNode();
private:
    static uiint mnType;
protected:
    uiint mnNumOfDOF;
public:
    virtual uiint& getType(){ return mnType;}
    virtual void setScalarDOF(const uiint& nNDOF);
    virtual void setVectorDOF(const uiint& nNDOF);
    virtual uiint& getScalarDOF();
    virtual uiint& getVectorDOF();
    virtual uiint getTotalDOF();
};
}
#endif	/* _VECTORNODE_H */
