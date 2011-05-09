//
// Node.h <- Vertex
//
//
//			2008.05.25
//			2008.12.18
//			k.Takeda
#include "TypeDef.h"

#include "NodeType.h"
#include "Logger.h"

#include "Vertex.h"
namespace pmw{
#ifndef NODE_HH_
#define NODE_HH_
class CNode:public CVertex{
public:
    // NodeID==VertexID
    CNode();
    virtual ~CNode();
    
protected:

    //// Arbitrary DOF (任意自由度) 変数コンテナー
    // vvdouble mvvArbitVar;
    //
    // Solid :tx,ty,tz
    // Shell :tx,ty,tz, rx,ry,rz
    //
    // Solid: x, y, z, xy, yz, zx
    // Shell: Mx,My,Mxy,Qx,Qy
    //
    // ---> VectorNode,ScalarNode,VectorScalarNode に移設.

    
    // prolongater プロロンゲーター用
    // 
    vector<CNode*> mvParentNode;//refine時にNodeの親となったNode

public:
    // Nodeタイプ
    virtual uint& getType()=0;


    // parameter accessor
    //
    virtual void resizeScalar(const uint& res_size)=0;
    virtual void resizeVector(const uint& res_size)=0;

    virtual void setScalar(const double& val, const uint& index)=0;
    virtual void setVector(const double& val, const uint& index)=0;

    virtual double& getScalar(const uint& i)=0;
    virtual double& getVector(const uint& i)=0;

    
    virtual uint numOfScalarParam()=0;
    virtual uint numOfVectorParam()=0;
    virtual uint numOfTotalParam()=0;

    
    // for prolongate
    vector<CNode*>& getParentNode(){ return mvParentNode;}
    CNode* getParentNode(const uint& index){ return mvParentNode[index];}
    void reserveParentNode(const uint& res_size){ mvParentNode.reserve(res_size);}
    void addParentNode(CNode* pNode){ mvParentNode.push_back(pNode);}
    uint getNumOfParentNode(){ return mvParentNode.size();}
};
#endif
}

