/* 
 * File:   VectorNode.h
 * Author: ktakeda
 *
 * Created on 2009/05/26, 15:17
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
    ////    vdouble mvParam;
    uiint mnNumOfDOF;

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
////    virtual double& getScalar(const uiint& i);
////    virtual double& getVector(const uiint& i){ return mvParam[i];}
////
////    virtual uiint numOfScalarParam();
////    virtual uiint numOfVectorParam();
////    virtual uiint numOfTotalParam();
};
}
#endif	/* _VECTORNODE_H */

