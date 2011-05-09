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
    static uint mnType;

protected:
    vdouble mvParam;

public:
    virtual uint& getType(){ return mnType;}

    // parameter accessor
    //
    virtual void reserveScalar(const uint& res_size);
    virtual void reserveVector(const uint& res_size);

    virtual void setScalar(const double& val);
    virtual void setScalar(const double& val, const uint& index);
    virtual void setVector(const vdouble& vVal);
    virtual void setVector(const double& val);
    virtual void setVector(const double& val, const uint& index);

    virtual double& getScalar(const uint& i);
    virtual vdouble& getVector(){ return mvParam;}

    virtual uint numOfScalarParam();
    virtual uint numOfVectorParam();

};
}
#endif	/* _VECTORNODE_H */

