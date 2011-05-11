/* 
 * File:   ScalarNode.h
 * Author: ktakeda
 *
 * Created on 2009/05/26, 13:04
 */

#ifndef _SCALARNODE_H_ffd42ed1_dfca_4d02_891c_908db139b697
#define	_SCALARNODE_H_ffd42ed1_dfca_4d02_891c_908db139b697

#include "Node.h"

namespace pmw{
class CScalarNode:public CNode{
public:
    CScalarNode();
    virtual ~CScalarNode();
private:
    static uiint mnType;

protected:
    ////    vdouble mvParam;
    uiint mnNumOfDOF;
    
public:
    virtual uiint& getType(){return mnType;}

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
////        if(i>= mvParam.size()){
////            Message_getScalar();
////            return mvParam[0];
////        }
////        return mvParam[i];
////    }
////    virtual double& getVector(const uiint& i);
////
////    virtual uiint numOfScalarParam();
////    virtual uiint numOfVectorParam();
////    virtual uiint numOfTotalParam();
};
}

#endif	/* _SCALARNODE_H */

