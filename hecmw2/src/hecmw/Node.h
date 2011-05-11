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

    bool mbSComm;//通信相手が小RankのCommNode : Refineの新NodeのID決定時に利用

public:
    // Nodeタイプ
    virtual uiint& getType()=0;

    // 通信相手が小RankのCommNode に所属しているか.否か.
    void markingSCommNode();
    bool isSCommNode();


    // parameter accessor
    //
    virtual void setScalarDOF(const uiint& nNDOF)=0;
    virtual void setVectorDOF(const uiint& nNDOF)=0;
    virtual uiint& getScalarDOF()=0;
    virtual uiint& getVectorDOF()=0;
    virtual uiint getTotalDOF()=0;
////    virtual void resizeScalar(const uiint& res_size)=0;
////    virtual void resizeVector(const uiint& res_size)=0;
////
////    virtual void setScalar(const double& val, const uiint& index)=0;
////    virtual void setVector(const double& val, const uiint& index)=0;
////
////    virtual double& getScalar(const uiint& i)=0;
////    virtual double& getVector(const uiint& i)=0;
////
////    virtual uiint numOfScalarParam()=0;
////    virtual uiint numOfVectorParam()=0;
////    virtual uiint numOfTotalParam()=0;

    
    // for prolongate
    vector<CNode*>& getParentNode(){ return mvParentNode;}
    CNode* getParentNode(const uiint& index){ return mvParentNode[index];}
    void reserveParentNode(const uiint& res_size){ mvParentNode.reserve(res_size);}
    void addParentNode(CNode* pNode){ mvParentNode.push_back(pNode);}
    uiint getNumOfParentNode(){ return mvParentNode.size();}
};
#endif
}

