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
    void Message_getScalar();
    static uint mnType;
protected:
    vdouble mvParam;

public:
    virtual uint& getType(){return mnType;}

    // parameter accessor
    //
    //virtual void reserveScalar(const uint& res_size);
    virtual void resizeScalar(const uint& res_size);
    //virtual void reserveVector(const uint& res_size);
    virtual void resizeVector(const uint& res_size);

    //virtual void setScalar(const double& val);
    virtual void setScalar(const double& val, const uint& index);
    //virtual void setVector(const vdouble& vVal);
    //virtual void setVector(const double& val);
    virtual void setVector(const double& val, const uint& index);

    virtual double& getScalar(const uint& i){
        if(i>= mvParam.size()){
            Message_getScalar();
            return mvParam[0];
        }
        return mvParam[i];
    }
    virtual vdouble& getVector();

    virtual uint numOfScalarParam();
    virtual uint numOfVectorParam();
    virtual uint numOfTotalParam();
};
}

#endif	/* _SCALARNODE_H */

