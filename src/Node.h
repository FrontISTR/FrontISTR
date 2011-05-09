/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   Node.h
|
|                     Written by T.Takeda,    2010/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#ifndef NODE_HH_5DBF3A6F_8A65_4ddd_B389_2CB193B87807
#define NODE_HH_5DBF3A6F_8A65_4ddd_B389_2CB193B87807
#include "CommonStd.h"
#include "TypeDef.h"
#include "NodeType.h"
#include "Logger.h"
#include "Vertex.h"
namespace pmw{
class CNode:public CVertex{
public:
    CNode();
    virtual ~CNode();
protected:
    uint mMGLevel;
    vector<CNode*> mvParentNode;
    vector<CNode*> mvChildNode;
public:
    void  setMGLevel(const uint& level){ mMGLevel=level;}
    uint& getMGLevel(){return mMGLevel;}
    virtual uint& getType()=0;
    virtual void resizeScalar(const uint& res_size)=0;
    virtual void resizeVector(const uint& res_size)=0;
    virtual void setScalar(const double& val, const uint& index)=0;
    virtual void setVector(const double& val, const uint& index)=0;
    virtual double& getScalar(const uint& i)=0;
    virtual vdouble& getVector()=0;
    virtual uint numOfScalarParam()=0;
    virtual uint numOfVectorParam()=0;
    virtual uint numOfTotalParam()=0;
    vector<CNode*>& getParentNode(){ return mvParentNode;}
    CNode* getParentNode(const uint& index){ return mvParentNode[index];}
    void reserveParentNode(const uint& res_size){ mvParentNode.reserve(res_size);}
    void addParentNode(CNode* pNode){ mvParentNode.push_back(pNode);}
    uint getNumOfParentNode(){ return mvParentNode.size();}
    vector<CNode*>& getChildNode(){ return mvChildNode;}
    CNode* getChildNode(const uint& index){ return mvChildNode[index];}
    void reserveChildNode(const uint& res_size){ mvChildNode.reserve(res_size);}
    void addChildNode(CNode* pNode){ mvChildNode.push_back(pNode);}
    uint getNumOfChildNode(){ return mvChildNode.size();}
};
}
#endif
