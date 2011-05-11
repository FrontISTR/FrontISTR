/* 
 * File:   ScalarVectorNode.h
 * Author: ktakeda
 *
 * Created on 2009/05/26, 16:02
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
    static uiint mnType;

protected:
    ////    vdouble mvScalarParam;
    ////    vdouble mvVectorParam;
    uiint mnNumOfSDOF;
    uiint mnNumOfVDOF;

public:
    virtual uiint& getType(){ return mnType;}

    // parameter accessor
    //
    virtual void setScalarDOF(const uiint& nNDOF);
    virtual void setVectorDOF(const uiint& nNDOF);
    virtual uiint& getScalarDOF();
    virtual uiint& getVectorDOF();
    virtual uiint getTotalDOF();
    
////    virtual void resizeScalar(const uiint& res_size);
////    virtual void resizeVector(const uiint& res_size);
////
////    virtual void setScalar(const double& val, const uiint& index);
////    virtual void setVector(const double& val, const uiint& index);
////
////    virtual double& getScalar(const uiint& i){
////        if(i>= mvScalarParam.size()){
////            Message_getScalar();
////            return mvScalarParam[0];
////        }
////        return mvScalarParam[i];
////    }
////    virtual double& getVector(const uiint& i){ return mvVectorParam[i];}
////
////    virtual uiint numOfScalarParam();
////    virtual uiint numOfVectorParam();
////    virtual uiint numOfTotalParam();
};
}
#endif	/* _SCALARVECTORNODE_H */

