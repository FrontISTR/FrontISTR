/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   VectorNode.h
|
|                     Written by T.Takeda,    2010/06/01
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
    static uint mnType;
protected:
    vdouble mvParam;
public:
    virtual uint& getType(){ return mnType;}
    virtual void resizeScalar(const uint& res_size);
    virtual void resizeVector(const uint& res_size);
    virtual void setScalar(const double& val, const uint& index);
    virtual void setVector(const double& val, const uint& index);
    virtual double& getScalar(const uint& i);
    virtual vdouble& getVector(){ return mvParam;}
    virtual uint numOfScalarParam();
    virtual uint numOfVectorParam();
    virtual uint numOfTotalParam();
};
}
#endif	/* _VECTORNODE_H */
