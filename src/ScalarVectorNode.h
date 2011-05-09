/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   ScalarVectorNode.h
|
|                     Written by T.Takeda,    2010/06/01
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
namespace pmw{
class CScalarVectorNode:public CNode{
public:
    CScalarVectorNode();
    virtual ~CScalarVectorNode();
private:
    void Message_getScalar();
    static uint mnType;
protected:
    vdouble mvScalarParam;
    vdouble mvVectorParam;
public:
    virtual uint& getType(){ return mnType;}
    virtual void resizeScalar(const uint& res_size);
    virtual void resizeVector(const uint& res_size);
    virtual void setScalar(const double& val, const uint& index);
    virtual void setVector(const double& val, const uint& index);
    virtual double& getScalar(const uint& i){
        if(i>= mvScalarParam.size()){
            Message_getScalar();
            return mvScalarParam[0];
        }
        return mvScalarParam[i];
    }
    virtual vdouble& getVector(){ return mvVectorParam;}
    virtual uint numOfScalarParam();
    virtual uint numOfVectorParam();
    virtual uint numOfTotalParam();
};
}
#endif	/* _SCALARVECTORNODE_H */
