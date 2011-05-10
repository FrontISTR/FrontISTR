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
    void Message_getScalar();
    static uint mnType;
protected:
    vdouble mvScalarParam;
    vdouble mvVectorParam;

public:
    virtual uint& getType(){ return mnType;}

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
        if(i>= mvScalarParam.size()){
            Message_getScalar();
            return mvScalarParam[0];
        }
        return mvScalarParam[i];
    }
    virtual double& getVector(const uint& i){ return mvVectorParam[i];}

    virtual uint numOfScalarParam();
    virtual uint numOfVectorParam();
    virtual uint numOfTotalParam();
};
}
#endif	/* _SCALARVECTORNODE_H */

